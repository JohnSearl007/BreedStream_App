gc()

library(shiny)
library(data.table)
library(bigmemory)

predictUI <- function(id) {
  ns <- NS(id)
  fluidPage(
    titlePanel("Make Predictions"),
    sidebarLayout(
      sidebarPanel(
        fileInput(ns("geno_csv"), "Transposed Genotype CSV File", accept = ".csv"),
        fileInput(ns("add_effects_csv"), "Additive Marker Effects CSV File", accept = ".csv"),
        checkboxInput(ns("dominance"), "Include Dominance", value = FALSE),
        conditionalPanel(
          condition = sprintf("input['%s'] == true", ns("dominance")),
          fileInput(ns("dom_effects_csv"), "Dominant Marker Effects CSV File", accept = ".csv")
        ),
        numericInput(ns("ploidy"), "Ploidy", value = 2, min = 2, step = 2),
        numericInput(ns("batch_size"), "Batch Size", value = 100, min = 1, step = 1),
        numericInput(ns("min_minor_allele"), "Minor Allele Threshold", value = 1, min = 0, step = 1),
        actionButton(ns("process_btn"), "Process Predictions")
      ),
      mainPanel(
        h4("Using transposed genotype data and marker effects (additive and optionally dominant), RR-BLUP predictions are made."),
        tableOutput(ns("prediction_table")),
        uiOutput(ns("file_status_ui")),
        verbatimTextOutput(ns("error_log"))
      )
    )
  )
}

predictServer <- function(input, output, session, workdir) {
  ns <- session$ns
  processed_data <- reactiveVal(NULL)
  output_path <- reactiveVal(NULL)
  error_message <- reactiveVal(NULL)
  
  observeEvent(input$process_btn, {
    # Clear previous errors and outputs
    error_message(NULL)
    processed_data(NULL)
    output_path(NULL)
    
    # Validate inputs
    req(input$geno_csv, input$add_effects_csv)
    if (input$dominance) req(input$dom_effects_csv)
    work_dir <- workdir()
    if (!is.character(work_dir) || !dir.exists(work_dir) || !file.access(work_dir, 2) == 0) {
      error_message("Invalid or non-writable working directory")
      return()
    }
    
    if (!is.numeric(input$ploidy) || input$ploidy < 2 || input$ploidy %% 2 != 0) {
      error_message("Ploidy must be an even positive integer")
      return()
    }
    if (!is.numeric(input$batch_size) || input$batch_size < 1 || input$batch_size %% 1 != 0) {
      error_message("Batch size must be a positive integer")
      return()
    }
    if (!is.numeric(input$min_minor_allele) || input$min_minor_allele < 0 || input$min_minor_allele %% 1 != 0) {
      error_message("Min minor allele must be a non-negative integer")
      return()
    }
    
    # Set up progress
    withProgress(message = "Processing predictions", value = 0, {
      
        # Step 1: Process genotype matrix (no parallel setup)
        incProgress(0.2, detail = "Processing genotype data...")
        test <- try({
          process_geno_matrix_transposed_columnwise(
            filename = input$geno_csv$datapath,
            ploidy = input$ploidy,
            dominance = input$dominance,
            min_minor_allele = input$min_minor_allele,
            batch_size = input$batch_size
          )
        }, silent = FALSE)
        
        if (inherits(test, "try-error")) {
          stop("Failed to process genotype matrix: ", conditionMessage(attr(test, "condition")))
        }
        
        # Step 2: Attach bigmemory matrices
        incProgress(0.3, detail = "Attaching coefficient matrices...")
        coeff <- bigmemory::attach.big.matrix(test$coeff)
        coeff.D <- if (input$dominance) bigmemory::attach.big.matrix(test$coeff.D) else NULL
        id <- test$id
        M <- ncol(coeff)
        N <- nrow(coeff)
        
        # Step 3: Load and validate effects files
        incProgress(0.4, detail = "Loading marker effects...")
        add_effects <- data.table::fread(input$add_effects_csv$datapath, select = c("marker", "effect"))
        if (!all(c("marker", "effect") %in% names(add_effects)) || nrow(add_effects) == 0) {
          stop("Additive effects CSV must contain non-empty 'marker' and 'effect' columns")
        }
        
        add_effects_vec <- add_effects$effect[1:min(nrow(add_effects), M)]
        if (nrow(add_effects) != M) {
          showNotification(sprintf("Additive effects CSV has %d markers, expected %d", nrow(add_effects), M), type = "warning")
        }
        if (length(add_effects_vec) < M) add_effects_vec <- c(add_effects_vec, rep(0, M - length(add_effects_vec)))
        
        if (input$dominance) {
          dom_effects <- data.table::fread(input$dom_effects_csv$datapath, select = c("marker", "effect"))
          if (!all(c("marker", "effect") %in% names(dom_effects)) || nrow(dom_effects) == 0) {
            stop("Dominant effects CSV must contain non-empty 'marker' and 'effect' columns")
          }
          
          dom_effects_vec <- dom_effects$effect[1:min(nrow(dom_effects), M)]
          if (nrow(dom_effects) != M) {
            showNotification(sprintf("Dominant effects CSV has %d markers, expected %d", nrow(dom_effects), M), type = "warning")
          }
          if (length(dom_effects_vec) < M) dom_effects_vec <- c(dom_effects_vec, rep(0, M - length(dom_effects_vec)))
        }
        
        # Step 4: Compute predictions in chunks and write directly to output file
        incProgress(0.5, detail = "Computing and writing predictions...")
        output_file <- file.path(work_dir, "HybridPredictions.csv")
        # Write header only (no data row)
        data.table::fwrite(data.table::data.table(id = character(), value = numeric()), output_file, col.names = TRUE, eol = ifelse(.Platform$OS.type == "windows", "\r\n", "\n"))
        
        num_chunks <- ceiling(N / input$batch_size)
        for (i in 1:num_chunks) {
          start_row <- (i - 1) * input$batch_size + 1
          end_row <- min(i * input$batch_size, N)
          id_chunk <- id[start_row:end_row]
          coeff_chunk <- coeff[start_row:end_row, , drop = FALSE]
          
          pred_add_chunk <- as.numeric(coeff_chunk %*% matrix(add_effects_vec, ncol = 1))
          
          if (input$dominance) {
            coeff.D_chunk <- coeff.D[start_row:end_row, , drop = FALSE]
            pred_dom_chunk <- as.numeric(coeff.D_chunk %*% matrix(dom_effects_vec, ncol = 1))
            value_chunk <- pred_add_chunk + pred_dom_chunk
          } else {
            value_chunk <- pred_add_chunk
          }
          
          output_chunk <- data.frame(id = id_chunk, value = value_chunk)
          # Append to output file without header
          data.table::fwrite(output_chunk, output_file, append = TRUE, col.names = FALSE, eol = ifelse(.Platform$OS.type == "windows", "\r\n", "\n"))
          incProgress(0.5 / num_chunks, detail = sprintf("Processed chunk %d/%d", i, num_chunks))
        }
        
        # Step 5: Finalize
        incProgress(1.0, detail = "Finalizing...")
        prediction_df <- data.table::fread(output_file, nrows = 10)[order(-value)]
        processed_data(prediction_df)
        showNotification(paste("Predictions saved to:", output_file), type = "message")
        output_path(output_file)
        })
    
    output$file_status_ui <- renderUI({
      req(output_path())
      tagList(
        h4("File Status"),
        p(paste("Saved:", output_path())),
        actionButton(ns("open_dir"), "Open Working Directory")
      )
    })
    
    observeEvent(input$open_dir, {
      req(workdir())
      if (!dir.exists(workdir())) {
        showNotification("Working directory does not exist.", type = "error")
        return()
      }
      if (.Platform$OS.type == "windows") {
        system(paste("explorer", shQuote(gsub("/", "\\", workdir(), fixed = TRUE))))
      } else {
        system(paste("xdg-open", shQuote(workdir())))
      }
    })
    
    output$prediction_table <- renderTable({
      req(processed_data())
      processed_data()
    })
    
    output$error_log <- renderText({
      req(error_message())
      error_message()
    })
    
    # Clean up on session end
    onStop(function() {
      gc()
      if (exists("stopImplicitCluster")) doParallel::stopImplicitCluster()
    })
  })
}

shinyApp(predictUI, predictServer)