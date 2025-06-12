gc()

library(shiny)
library(data.table)
library(dplyr)
library(BreedStream)

twoStageUI <- function(id) {
  ns <- NS(id)
  fluidPage(
    titlePanel("Two-Stage Genomic Analysis with H.Matrix_Optimized"),
    sidebarLayout(
      sidebarPanel(
        h4("Data Inputs"),
        fileInput(ns("pheno_file"), "Upload Phenotypic Data (CSV)", accept = ".csv"),
        fileInput(ns("geno_file"), "Upload Genotypic Data (CSV)", accept = ".csv"),
        fileInput(ns("ped_file"), "Upload Pedigree Data (CSV)", accept = ".csv"),
        h4("Model"),
        selectInput(ns("traits"), "Select Traits:", multiple = TRUE, choices = NULL),
        selectInput(ns("fixed_effects"), "Fixed Effects:", multiple = TRUE, choices = NULL),
        selectInput(ns("random_effects"), "Random Effects:", multiple = TRUE, choices = NULL),
        h4("H-Matrix Optimization"),
        numericInput(ns("ploidy"), "Ploidy Level:", value = 2, min = 2, step = 2),
        numericInput(ns("max.iter"), "Max Iterations:", value = 100, min = 1),
        numericInput(ns("blend_lower"), "Blending Lower Bound:", value = 1e-5, min = 1e-5, max = 1-1e-5),
        numericInput(ns("blend_upper"), "Blending Upper Bound:", value = 1e-5, min = 1e-5, max = 1-1e-5),
        radioButtons(ns("dominance"), "Include Dominance:", choices = list("TRUE" = "TRUE", "FALSE" = "FALSE"), selected = "TRUE"),
        radioButtons(ns("map"), "Map Included:", choices = list("TRUE" = "TRUE", "FALSE" = "FALSE"), selected = "TRUE"),
        textInput(ns("workspace"), "Workspace Size:", value = "12Gb"),
        textInput(ns("pworkspace"), "PWorkspace Size:", value = "12Gb"),
        numericInput(ns("min.minor.allele"), "Minor Allele Threshold:", value = 1, min = 1, step = 1),
        textInput(ns("fix.eff.marker"), "Fixed Effect Markers (comma-separated, optional):", value = ""),
        textInput(ns("covariates"), "Covariates (comma-separated, optional):", value = ""),
        fileInput(ns("mask_file"), "Data to Mask (CSV, optional):", accept = ".csv"),
        selectInput(ns("method"), "Method (optional):", choices = c("None" = "", "MME" = "MME", "Vinv" = "Vinv"), selected = ""),
        textInput(ns("non.add"), "Non-Additive Effects:", value = "dom"),
        numericInput(ns("tol"), "Tolerance:", value = 1e-4),
        h4("Run Analysis"),
        actionButton(ns("run_analysis"), "Run Two-Stage Analysis"),
        uiOutput(ns("save_ui"))
      ),
      mainPanel(
        h4("A Two-Stage Analysis approach is implemented with user specification of parameters such as inclusion of a dominance effect. The user optionally has the ability to utilize a blended H Matrix by setting upper and lower blending limits. If an H Matrix is requested, the optimal blending will be identified based on the AIC values of the tested models."),
        h4("Analysis Results"),
        verbatimTextOutput(ns("results_summary")),
        uiOutput(ns("file_status_ui"))
      )
    )
  )
}

twoStageServer <- function(input, output, session, workdir) {
  ns <- session$ns
  analysis_results <- reactiveVal(NULL)
  output_path <- reactiveVal(NULL)
  
  observe({
    req(input$pheno_file)
    tryCatch({
      df <- data.table::fread(input$pheno_file$datapath)
      updateSelectInput(session, "traits", choices = names(df))
      updateSelectInput(session, "fixed_effects", choices = names(df))
      updateSelectInput(session, "random_effects", choices = names(df))
    }, error = function(e) {
      showNotification(paste("Error loading phenotypic data:", e$message), type = "error")
    })
  })
  
  mask_data <- reactive({
    if (is.null(input$mask_file)) return(NULL)
    tryCatch({
      df <- data.table::fread(input$mask_file$datapath)
      if (!all(colnames(df) %in% c("id", "env", "trait"))) stop("Mask CSV must contain only 'id', 'env', and/or 'trait' columns.")
      return(df)
    }, error = function(e) {
      showNotification(paste("Error loading mask file:", e$message), type = "error")
      stop(e$message)
    })
  })
  
  observeEvent(input$run_analysis, {
    req(input$pheno_file, input$geno_file, input$ped_file, input$traits)
    tryCatch({
      ped_data <- as.data.frame(data.table::fread(input$ped_file$datapath))
      fix.eff.marker <- if (nchar(trimws(input$fix.eff.marker)) > 0) strsplit(trimws(input$fix.eff.marker), ",\\s*")[[1]] else NULL
      covariates <- if (nchar(trimws(input$covariates)) > 0) strsplit(trimws(input$covariates), ",\\s*")[[1]] else NULL
      method <- if (input$method == "" || input$method == "None") NULL else input$method
      map <- as.logical(input$map)
      dominance <- as.logical(input$dominance)
      
      # Set working directory for H.Matrix_Optimized
      original_wd <- getwd()
      setwd(workdir())
      on.exit(setwd(original_wd), add = TRUE)
      
      tuned <- H.Matrix_Optimized(
        pheno.data = input$pheno_file$datapath,
        geno.data = input$geno_file$datapath,
        pedigree.data = ped_data,
        trait = input$traits,
        Fixed = input$fixed_effects,
        Random = input$random_effects,
        blend.lower = input$blend_lower,
        blend.upper = input$blend_upper,
        ploidy = input$ploidy,
        max.iter = input$max.iter,
        map = map,
        dominance = dominance,
        workspace = input$workspace,
        pworkspace = input$pworkspace,
        min.minor.allele = input$min.minor.allele,
        fix.eff.marker = fix.eff.marker,
        covariates = covariates,
        mask = mask_data(),
        method = method,
        non.add = input$non.add,
        tol = input$tol
      )
      analysis_results(tuned)
      
      output$save_ui <- renderUI({
        actionButton(ns("save_results"), "Save Results")
      })
      
      output$results_summary <- renderPrint({
        cat("Analysis completed. Save the results to proceed.")
      })
      
      showNotification("Two-stage analysis completed successfully.", type = "message")
    }, error = function(e) {
      showNotification(paste("Error running two-stage analysis:", e$message), type = "error")
      output$results_summary <- renderPrint({
        cat("Error: ", e$message)
      })
      analysis_results(NULL)
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
  
  observeEvent(input$save_results, {
    req(workdir(), analysis_results())
    if (!dir.exists(workdir())) {
      showNotification("Selected working directory does not exist.", type = "error")
      return()
    }
    output_path_val <- file.path(workdir(), "two_stage_results.rds")
    tryCatch({
      saveRDS(analysis_results(), output_path_val)
      output_path(output_path_val)
      showNotification(paste("File saved to:", output_path_val), type = "message")
      analysis_results(NULL)
    }, error = function(e) {
      showNotification(paste("Error saving file:", e$message), type = "error")
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

shinyApp(twoStageUI, twoStageServer)
