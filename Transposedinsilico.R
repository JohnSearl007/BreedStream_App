gc()

library(shiny)
library(data.table)
library(dplyr)
library(BreedStream)

insilico_transposedUI <- function(id) {
  ns <- NS(id)
  fluidPage(
    titlePanel("in silico Hybrid Genotypes"),
    sidebarLayout(
      sidebarPanel(
        textInput(ns("set_names"), "Enter Hybrid Set Names (comma-separated)", placeholder = "e.g., SS,NSS"),
        uiOutput(ns("file_inputs")),
        actionButton(ns("run_analysis"), "Generate Hybrids"),
        uiOutput(ns("save_ui"))
      ),
      mainPanel(
        h4("Using parental genotypes, in silico genotypes are generated. The output format is transposed to that of the standard 'in silico Hybrids' panel in order to minimize the downstream computational burden for 'RR-BLUP -> Prediction'."),
        verbatimTextOutput(ns("status_message")),
        uiOutput(ns("file_status_ui"))
      )
    )
  )
}

insilico_transposedServer <- function(input, output, session, workdir) {
  ns <- session$ns
  output_path <- reactiveVal(NULL)
  
  set_names <- reactive({
    req(input$set_names)
    strsplit(input$set_names, ",")[[1]] %>% trimws()
  })
  
  output$file_inputs <- renderUI({
    sets <- set_names()
    do.call(tagList, lapply(sets, function(set) {
      list(
        h4(set),
        fileInput(ns(paste0("geno_file_female_", set)), paste("Upload Female Parent Genotype Data for", set, "(CSV)"), accept = ".csv"),
        fileInput(ns(paste0("geno_file_male_", set)), paste("Upload Male Parent Genotype Data for", set, "(CSV)"), accept = ".csv")
      )
    }))
  })
  
  observeEvent(input$run_analysis, {
    req(set_names())
    sets <- set_names()
    temp_files <- list()
    # metadata_files <- list()
    
    withProgress(message = "Generating hybrids", value = 0, {
      tryCatch({
        # Generate hybrid data for each set
        for (set in sets) {
          female_file <- input[[paste0("geno_file_female_", set)]]
          male_file <- input[[paste0("geno_file_male_", set)]]
          req(female_file, male_file)
          
          temp_output <- tempfile(fileext = ".csv")
          insilico_transposed(
            geno_file_female = female_file$datapath,
            geno_file_male = male_file$datapath,
            output_file = temp_output
          )
          temp_files[[set]] <- temp_output
          incProgress(1 / length(sets))
        }
        
        if (length(temp_files) == 0) {
          stop("No hybrid files generated.")
        }
        
        # Combine files
        output_file <- file.path(workdir(), "insilico_transposedhybrids.csv")
        if (length(temp_files) == 1) {
          system(paste("cp", shQuote(temp_files[[1]]), shQuote(output_file)))
        } else {
          cmd <- paste(
            "cat", shQuote(temp_files[[1]]), ">", shQuote(output_file),
            ";",
            paste(
              sapply(temp_files[-1], function(f) paste("tail -n +2", shQuote(f), ">>", shQuote(output_file))),
              collapse = "; "
            )
          )
          system(paste("bash -c", shQuote(cmd)))
        }
        
        output$save_ui <- renderUI({
          actionButton(ns("save_hybrids"), "Save Hybrids")
        })
        
        output$status_message <- renderText({
          "Hybrids generated and combined successfully."
        })
        
        showNotification("Hybrids generated successfully.", type = "message")
      }, error = function(e) {
        showNotification(paste("Error generating hybrids:", e$message), type = "error")
        output$status_message <- renderText({
          paste("Error:", e$message)
        })
      }, finally = {
        # Clean up temporary files
        unlink(unlist(temp_files))
      })
    })
  })
  
  output$file_status_ui <- renderUI({
    req(output_path())
    tagList(
      h4("File Status"),
      p(paste("Saved:", output_path())),
      actionButton(ns("open_dir"), "Open Working Directory")
    )
  })
  
  observeEvent(input$save_hybrids, {
    req(workdir())
    if (!dir.exists(workdir())) {
      showNotification("Selected working directory does not exist.", type = "error")
      return()
    }
    output_path_val <- file.path(workdir(), "insilico_transposedhybrids.csv")
    tryCatch({
      # File already saved during processing; confirm existence
      if (!file.exists(output_path_val)) {
        stop("Output file not found.")
      }
      output_path(output_path_val)
      showNotification(paste("File saved to:", output_path_val), type = "message")
    }, error = function(e) {
      showNotification(paste("Error confirming file save:", e$message), type = "error")
    })
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
  
  # Ensure garbage collection when server function exits
  on.exit(gc())
}

shinyApp(insilico_transposedUI, insilico_transposedServer)