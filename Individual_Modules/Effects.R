gc()

library(shiny)
library(data.table)
library(StageWise)

effectsUI <- function(id) {
  ns <- NS(id)
  fluidPage(
    titlePanel("Marker Effects"),
    sidebarLayout(
      sidebarPanel(
        fileInput(ns("geno_file"), "Genotype Data File (.csv)", accept = ".csv"),
        fileInput(ns("weight_file"), "Optimal Weight File (.txt)", accept = ".txt"),
        fileInput(ns("ped_file"), "Pedigree File", accept = ".csv"),
        fileInput(ns("tuned_file"), "Upload Two Stage Results .rds file", accept = ".rds"),
        fileInput(ns("index_file"), "Upload Multi-trait Index Coefficients .rds file (skip for single trait)", accept = ".rds"),
        numericInput(ns("ploidy"), "Ploidy Level:", value = 2, min = 2, step = 2),
        radioButtons(ns("map"), "Map Included:", choices = list("Yes" = "TRUE", "No" = "FALSE"), selected = "TRUE"),
        numericInput(ns("min.minor.allele"), "Minor Allele Threshold:", value = 1, min = 1, step = 1),
        radioButtons(ns("dominance"), "Include Dominance:", choices = list("Yes" = "TRUE", "No" = "FALSE"), selected = "TRUE"),
        actionButton(ns("processBtn"), "Calculate Marker Effects"),
        uiOutput(ns("save_ui"))
      ),
      mainPanel(
        h4("Calucation of marker effects (additive and optionally dominant)."),
        tableOutput(ns("add_effects_table")),
        uiOutput(ns("file_status_ui"))
      )
    )
  )
}

effectsServer <- function(input, output, session, workdir) {
  ns <- session$ns
  results <- reactiveVal(NULL)
  output_paths <- reactiveVal(NULL)
  
  observeEvent(input$processBtn, {
    req(input$geno_file, input$weight_file, input$ped_file, input$tuned_file)
    tryCatch({
      output$status <- renderText("Processing...")
      optimal_weight <- as.numeric(fread(input$weight_file$datapath))
      ped_data <- as.data.frame(fread(input$ped_file$datapath))
      tuned_data <- readRDS(input$tuned_file$datapath)
      index_coeff <- if (is.null(input$index_file)) NULL else readRDS(input$index_file$datapath)
      map_logical <- as.logical(input$map)
      dominance_logical <- as.logical(input$dominance)
      genoH <- StageWise::read_geno(
        filename = input$geno_file$datapath,
        ploidy = input$ploidy,
        map = map_logical,
        min.minor.allele = input$min.minor.allele,
        ped = ped_data,
        w = optimal_weight,
        dominance = dominance_logical
      )
      temp_results <- list(genoH = genoH)
      temp_results$add_effects <- blup(data = tuned_data, geno = genoH, what = "AM", index.coeff = index_coeff)
      if (dominance_logical) {
        temp_results$dom_effects <- blup(data = tuned_data, geno = genoH, what = "DM", index.coeff = index_coeff)
        message("Dominance effects calculated")
      }
      results(temp_results)
      output$save_ui <- renderUI({
        actionButton(ns("save_all"), "Save All Results")
      })
      showNotification("Marker effects calculated successfully.", type = "message")
    }, error = function(e) {
      showNotification(paste("Error calculating marker effects:", e$message), type = "error")
    })
  })
  
  output$add_effects_table <- renderTable({
    req(results())
    head(results()$add_effects)
  })
  
  output$file_status_ui <- renderUI({
    req(output_paths())
    paths <- output_paths()
    tagList(
      h4("File Status"),
      lapply(names(paths), function(file) {
        p(paste("Saved:", paths[[file]]))
      }),
      actionButton(ns("open_dir"), "Open Working Directory")
    )
  })
  
  observeEvent(input$save_all, {
    req(workdir(), results())
    if (!dir.exists(workdir())) {
      showNotification("Selected working directory does not exist.", type = "error")
      return()
    }
    paths <- list()
    tryCatch({
      # Save marker_effects.rds
      rds_path <- file.path(workdir(), "marker_effects.rds")
      saveRDS(results(), rds_path)
      paths[["marker_effects.rds"]] <- rds_path
      
      # Save add_effects.csv
      add_path <- file.path(workdir(), "add_effects.csv")
      fwrite(results()$add_effects, add_path)
      paths[["add_effects.csv"]] <- add_path
      
      # Save dom_effects.csv if applicable
      if ("dom_effects" %in% names(results())) {
        dom_path <- file.path(workdir(), "dom_effects.csv")
        fwrite(results()$dom_effects, dom_path)
        paths[["dom_effects.csv"]] <- dom_path
      }
      
      output_paths(paths)
      showNotification("All files saved successfully.", type = "message")
    }, error = function(e) {
      showNotification(paste("Error saving files:", e$message), type = "error")
      output_paths(NULL)
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

shinyApp(effectsUI, effectsServer)
