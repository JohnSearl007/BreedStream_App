gc()

library(shiny)
library(shinyWidgets)
library(StageWise)

indexUI <- function(id) {
  ns <- NS(id)
  fluidPage(
    titlePanel("Selection Index"),
    sidebarLayout(
      sidebarPanel(
        fileInput(ns("rdsFile"), "Upload Two Stage Results .rds file", accept = ".rds"),
        radioButtons(ns("restricted"), "Restricted Index", choices = c("TRUE" = TRUE, "FALSE" = FALSE), selected = TRUE),
        uiOutput(ns("traitInputs")),
        actionButton(ns("processBtn"), "Calculate Index Coefficients"),
        uiOutput(ns("save_ui"))
      ),
      mainPanel(
        h4("Only for use with multi-trait selection indices. Skip for single trait selection."),
        tableOutput(ns("index_table")),
        uiOutput(ns("file_status_ui"))
      )
    )
  )
}

indexServer <- function(input, output, session, workdir) {
  ns <- session$ns
  tuned <- reactiveVal(NULL)
  index_coeff <- reactiveVal(NULL)
  output_path <- reactiveVal(NULL)
  
  observeEvent(input$rdsFile, {
    req(input$rdsFile)
    tryCatch({
      tuned_data <- readRDS(input$rdsFile$datapath)
      tuned(tuned_data)
      showNotification("Two Stage Results loaded successfully.", type = "message")
    }, error = function(e) {
      showNotification(paste("Error loading .rds file:", e$message), type = "error")
      tuned(NULL)
    })
  })
  
  traits <- reactive({
    req(tuned())
    names(tuned()@avg.env)
  })
  
  output$traitInputs <- renderUI({
    req(traits())
    trait_list <- lapply(traits(), function(trait) {
      fluidRow(
        column(12,
               h4(trait),
               numericInput(ns(paste0("coeff_", trait)), "Coefficient", value = 1, min = 0),
               if (input$restricted == "TRUE") {
                 textInput(ns(paste0("sign_", trait)), "Sign (>, =, <)", value = ">")
               }
        )
      )
    })
    do.call(tagList, trait_list)
  })
  
  observeEvent(input$processBtn, {
    req(traits())
    tryCatch({
      merit.coeff <- sapply(traits(), function(trait) input[[paste0("coeff_", trait)]])
      names(merit.coeff) <- traits()
      if (input$restricted == "TRUE") {
        restricted_df <- data.frame(
          trait = traits(),
          sign = sapply(traits(), function(trait) input[[paste0("sign_", trait)]])
        )
        gain_result <- gain(input = tuned(), merit = merit.coeff, restricted = restricted_df)
      } else {
        gain_result <- gain(input = tuned(), merit = merit.coeff)
      }
      result <- setNames(gain_result$table$index, gain_result$table$trait)
      index_coeff(result)
      output$save_ui <- renderUI({
        actionButton(ns("save_index"), "Save Index Coefficients")
      })
      output$index_table <- renderTable({
        as.data.frame(index_coeff())
      })
      showNotification("Index coefficients calculated successfully.", type = "message")
    }, error = function(e) {
      showNotification(paste("Error calculating index coefficients:", e$message), type = "error")
      index_coeff(NULL)
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
  
  observeEvent(input$save_index, {
    req(workdir(), index_coeff())
    if (!dir.exists(workdir())) {
      showNotification("Selected working directory does not exist.", type = "error")
      return()
    }
    output_path_val <- file.path(workdir(), "index_coeff.rds")
    tryCatch({
      saveRDS(index_coeff(), output_path_val)
      output_path(output_path_val)
      showNotification(paste("File saved to:", output_path_val), type = "message")
      index_coeff(NULL)
      tuned(NULL)
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

shinyApp(indexUI, indexServer)
