gc()

library(shiny)
library(dplyr)
library(BreedStream)

omaUI <- function(id) {
  ns <- NS(id)
  fluidPage(
    titlePanel("OMA"),
    sidebarLayout(
      sidebarPanel(
        fileInput(ns("uc_file"), "Usefulness Criterion .rds File", accept = ".rds"),
        textInput(ns("group_names"), "Group Names (comma-separated)", placeholder = "e.g., SS,NSS"),
        textInput(ns("dF"), "dF Values (comma-separated)", value = "0.005,0.01", placeholder = "e.g., 0.005,0.01"),
        numericInput(ns("ploidy"), "Ploidy Level:", value = 2, min = 2, step = 2),
        numericInput(ns("dF_adapt_step"), "dF.adapt Step:", value = 0.005, min = 0, step = 0.001),
        numericInput(ns("dF_adapt_max"), "dF.adapt Max:", value = 0.1, min = 0, step = 0.01),
        radioButtons(ns("solver"), "CVXR Solver:", choices = list("ECOS" = "ECOS", "SCS" = "SCS"), selected = "ECOS"),
        radioButtons(ns("base"), "Base Method:", choices = list("RM" = "RM", "Current" = "current"), selected = "RM"),
        actionButton(ns("processBtn"), "Process Data"),
        uiOutput(ns("save_ui"))
      ),
      mainPanel(
        h4("Using the Usefulness Criterion as merit, Optimal Mate Allocation is performed."),
        uiOutput(ns("group_tables")),
        uiOutput(ns("file_status_ui"))
      )
    )
  )
}

omaServer <- function(input, output, session, workdir) {
  ns <- session$ns
  oma_data <- reactiveVal(NULL)
  output_paths <- reactiveVal(NULL)
  
  uc_data <- reactive({
    req(input$uc_file)
    tryCatch({
      readRDS(input$uc_file$datapath)
    }, error = function(e) {
      showNotification(paste("Error loading UC file:", e$message), type = "error")
      stop(e$message)
    })
  })
  
  observeEvent(input$processBtn, {
    req(input$uc_file, input$group_names, input$dF)
    tryCatch({
      UC <- uc_data()
      group_names <- trimws(strsplit(input$group_names, ",")[[1]])
      dF_values <- as.numeric(unlist(strsplit(input$dF, ",")))
      if (any(is.na(dF_values))) stop("dF must contain valid comma-separated numbers")
      dF_adapt <- list(step = input$dF_adapt_step, max = input$dF_adapt_max)
      data_list <- list()
      for (group in group_names) {
        if (!group %in% names(UC)) stop(paste("Group", group, "not found in UC data"))
        matings <- UC[[group]][["usef_add"]][[2]][, -4]
        colnames(matings) <- c("parent1", "parent2", "merit")
        parents <- data.frame(id = unique(c(unique(matings$parent1), unique(matings$parent2))), min = 0, max = 1)
        matings <- data.frame(matings, min = 0, max = 1)
        ans <- oma_optimized(
          dF = dF_values,
          parents = parents,
          matings = matings,
          ploidy = input$ploidy,
          K = UC[[group]][["relMat"]],
          dF.adapt = dF_adapt,
          base = input$base,
          solver = input$solver
        )
        data_list[[group]] <- ans
      }
      oma_data(data_list)
      
      output$group_tables <- renderUI({
        lapply(names(data_list), function(group) {
          list(
            h4(paste("Group:", group)),
            tableOutput(ns(paste0("table_", group)))
          )
        })
      })
      
      output$save_ui <- renderUI({
        actionButton(ns("save_all"), "Save All Results")
      })
      
      lapply(names(data_list), function(group) {
        output[[paste0("table_", group)]] <- renderTable({
          # Filter rows where the 'value' column is not equal to 0
          as.data.frame(data_list[[group]]$om) %>% filter(value != 0)
        })
      })
      
      showNotification("OMA data processed successfully.", type = "message")
    }, error = function(e) {
      showNotification(paste("Error processing OMA data:", e$message), type = "error")
      oma_data(NULL)
    })
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
    req(workdir(), oma_data())
    if (!dir.exists(workdir())) {
      showNotification("Selected working directory does not exist.", type = "error")
      return()
    }
    paths <- list()
    tryCatch({
      # Save OMA.rds
      rds_path <- file.path(workdir(), "OMA.rds")
      saveRDS(oma_data(), rds_path)
      paths[["OMA.rds"]] <- rds_path
      
      # Save per-group om.csv files
      for (group in names(oma_data())) {
        group_path <- file.path(workdir(), paste0(group, "_om.csv"))
        data.table::fwrite(oma_data()[[group]]$om, group_path)
        paths[[paste0(group, "_om.csv")]] <- group_path
      }
      
      output_paths(paths)
      showNotification("All files saved successfully.", type = "message")
      oma_data(NULL)
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

shinyApp(omaUI, omaServer)
