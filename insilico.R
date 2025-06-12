gc()

library(shiny)
library(data.table)
library(dplyr)
library(BreedStream)

insilicoUI <- function(id) {
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
        h4("Using parental genotypes, in silico genotypes are generated."),
        verbatimTextOutput(ns("status_message")),
        uiOutput(ns("file_status_ui"))
      )
    )
  )
}

insilicoServer <- function(input, output, session, workdir) {
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
    metadata_files <- list()
    
    withProgress(message = "Generating hybrids", value = 0, {
      tryCatch({
        # Generate hybrid data for each set
        for (set in sets) {
          female_file <- input[[paste0("geno_file_female_", set)]]
          male_file <- input[[paste0("geno_file_male_", set)]]
          req(female_file, male_file)
          
          temp_output <- tempfile(fileext = ".csv")
          insilico(
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
        
        # Extract metadata for consistency check
        for (set in sets) {
          temp_file <- temp_files[[set]]
          metadata_file <- tempfile(fileext = ".csv")
          system(paste("cut -d, -f1-3", shQuote(temp_file), ">", shQuote(metadata_file)))
          metadata_files[[set]] <- metadata_file
        }
        
        # Validate metadata files
        for (set in sets) {
          if (!file.exists(metadata_files[[set]]) || file.info(metadata_files[[set]])$size == 0) {
            stop(paste("Metadata file for set", set, "is missing or empty"))
          }
        }
        
        # Check metadata consistency
        all_identical <- TRUE
        if (length(metadata_files) > 1) {
          base_metadata <- metadata_files[[1]]
          for (i in 2:length(metadata_files)) {
            exit_status <- system(
              paste("diff", shQuote(base_metadata), shQuote(metadata_files[[i]])),
              ignore.stdout = TRUE,
              ignore.stderr = TRUE
            )
            if (exit_status != 0) {
              all_identical <- FALSE
              break
            }
          }
        } else {
          # Single set: no comparison needed, but ensure metadata is valid
          base_metadata <- metadata_files[[1]]
          metadata_check <- read.csv(base_metadata, nrows = 1)
          if (!all(c("marker", "chrom", "position") %in% colnames(metadata_check))) {
            stop("Metadata file for single set does not contain required columns: marker, chrom, position")
          }
        }
        
        if (!all_identical) {
          stop("Metadata columns (marker, chrom, position) are not identical across all sets.")
        }
        
        # Combine files
        output_file <- file.path(workdir(), "insilicohybrids.csv")
        if (length(temp_files) == 1) {
          # Single set: copy directly
          system(paste("cp", shQuote(temp_files[[1]]), shQuote(output_file)))
        } else {
          # Multiple sets: combine with paste
          paste_cmd <- paste(
            "paste -d,",
            paste("<(cut -d, -f1-3", shQuote(temp_files[[1]]), ")"),
            paste(
              sapply(temp_files, function(f) paste("<(cut -d, -f4-", shQuote(f), ")")),
              collapse = " "
            ),
            ">", shQuote(output_file)
          )
          system(paste("bash -c", shQuote(paste_cmd)))
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
        unlink(unlist(metadata_files))
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
    output_path_val <- file.path(workdir(), "insilicohybrids.csv")
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

shinyApp(insilicoUI, insilicoServer)