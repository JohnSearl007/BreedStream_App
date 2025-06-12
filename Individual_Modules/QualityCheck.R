gc()

library(shiny)
library(data.table)
library(dplyr)
library(BreedStream)
library(StageWise)
library(ggplot2)

qualityUI <- function(id) {
  ns <- NS(id)
  fluidPage(
    titlePanel("Quality Check"),
    sidebarLayout(
      sidebarPanel(
        h4("Data Inputs"),
        fileInput(ns("pheno_file"), "Upload Formatted Phenotypic Data (CSV)", accept = ".csv"),
        h4("Model"),
        selectInput(ns("traits"), "Select Traits:", multiple = TRUE, choices = NULL),
        selectInput(ns("fixed_effects"), "Fixed Effects:", multiple = TRUE, choices = NULL),
        selectInput(ns("random_effects"), "Random Effects:", multiple = TRUE, choices = NULL),
        h4("Run Analysis"),
        actionButton(ns("run_analysis"), "Run Analysis")
      ),
      mainPanel(
        h4("Run the first stage of a Two-Stage analysis using SPATS not ASReml in order to generate diagnostic plots for selected traits. This allows you to check the residuals, heritability, etc. before proceeding."),
        textOutput(ns("status_message")),
        uiOutput(ns("plots_ui"))
      )
    )
  )
}

qualityServer <- function(input, output, session, workdir) {
  ns <- session$ns
  # Reactive value to store the uploaded data
  data <- reactive({
    req(input$pheno_file)
    fread(input$pheno_file$datapath)
  })
  
  # Update select inputs when phenotypic data is uploaded
  observe({
    req(data())
    all_cols <- names(data())
    numeric_cols <- names(data())[sapply(data(), is.numeric)]
    updateSelectInput(session, "traits", choices = numeric_cols)
    updateSelectInput(session, "fixed_effects", choices = all_cols)
    updateSelectInput(session, "random_effects", choices = all_cols)
  })
  
  # Compute plots and save to files when run button is clicked
  plots_list <- eventReactive(input$run_analysis, {
    req(input$pheno_file, input$traits)
    withProgress(message = "Running analysis", value = 0, {
      n <- length(input$traits)
      plots_files <- list()
      for (i in 1:n) {
        trait <- input$traits[i]
        incProgress(1/n, detail = paste("Processing trait:", trait))
        plots <- Stage1_plots(
          filename = input$pheno_file$datapath,
          trait = trait,
          Fixed = input$fixed_effects,
          Random = input$random_effects
        )
        # Save each plot to a temporary PNG file
        boxplot_file <- tempfile(fileext = ".png")
        png(boxplot_file, width = 600, height = 400)
        print(plots$boxplot)
        dev.off()
        
        qqplot_file <- tempfile(fileext = ".png")
        png(qqplot_file, width = 600, height = 400)
        print(plots$qqplot)
        dev.off()
        
        heritability_file <- tempfile(fileext = ".png")
        png(heritability_file, width = 600, height = 400)
        print(plots$heritability)
        dev.off()
        
        # Handle multiple spatial plots
        spatial_files <- list()
        if (is.list(plots$spatial)) {
          for (env in names(plots$spatial)) {
            spatial_file <- tempfile(fileext = ".png")
            png(spatial_file, width = 600, height = 400)
            print(plots$spatial[[env]])
            dev.off()
            spatial_files[[env]] <- spatial_file
          }
        } else {
          spatial_file <- tempfile(fileext = ".png")
          png(spatial_file, width = 600, height = 400)
          print(plots$spatial)
          dev.off()
          spatial_files[["default"]] <- spatial_file
        }
        
        plots_files[[trait]] <- list(
          boxplot = boxplot_file,
          qqplot = qqplot_file,
          heritability = heritability_file,
          spatial = spatial_files
        )
      }
      plots_files
    })
  })
  
  # Render status message
  output$status_message <- renderText({
    if (is.null(plots_list())) {
      "Please select traits, fixed and random effects, then click 'Run Analysis' to generate plots."
    } else {
      "Analysis completed. Plots are displayed below."
    }
  })
  
  # Dynamically generate UI for image outputs
  output$plots_ui <- renderUI({
    req(plots_list())
    plot_output_list <- lapply(names(plots_list()), function(trait) {
      trait_id <- gsub("[^a-zA-Z0-9]", "_", trait)
      base_list <- list(
        h3(paste("Trait:", trait)),
        imageOutput(ns(paste0("boxplot_", trait_id))),
        imageOutput(ns(paste0("qqplot_", trait_id))),
        imageOutput(ns(paste0("heritability_", trait_id)))
      )
      # Add spatial plots for each environment
      spatial_outputs <- lapply(names(plots_list()[[trait]]$spatial), function(env) {
        imageOutput(ns(paste0("spatial_", trait_id, "_", gsub("[^a-zA-Z0-9]", "_", env))))
      })
      c(base_list, spatial_outputs)
    })
    do.call(tagList, unlist(plot_output_list, recursive = FALSE))
  })
  
  # Render images for each plot
  observe({
    req(plots_list())
    for (trait in names(plots_list())) {
      trait_id <- gsub("[^a-zA-Z0-9]", "_", trait)
      local({
        trait <- trait
        trait_id <- trait_id
        output[[paste0("boxplot_", trait_id)]] <- renderImage({
          list(src = plots_list()[[trait]]$boxplot, contentType = "image/png", width = 600, height = 400)
        }, deleteFile = FALSE)
        output[[paste0("qqplot_", trait_id)]] <- renderImage({
          list(src = plots_list()[[trait]]$qqplot, contentType = "image/png", width = 600, height = 400)
        }, deleteFile = FALSE)
        output[[paste0("heritability_", trait_id)]] <- renderImage({
          list(src = plots_list()[[trait]]$heritability, contentType = "image/png", width = 600, height = 400)
        }, deleteFile = FALSE)
        # Render each spatial plot with a local block for each environment
        for (env in names(plots_list()[[trait]]$spatial)) {
          env_id <- gsub("[^a-zA-Z0-9]", "_", env)
          local({
            current_env <- env
            output[[paste0("spatial_", trait_id, "_", env_id)]] <- renderImage({
              list(src = plots_list()[[trait]]$spatial[[current_env]], contentType = "image/png", width = 600, height = 400)
            }, deleteFile = FALSE)
          })
        }
      })
    }
  })
  
  # Ensure garbage collection when server function exits
  on.exit(gc())
}

shinyApp(qualityUI, qualityServer)
