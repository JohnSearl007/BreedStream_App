gc()

library(shiny)
library(dplyr)
library(BreedStream)

formatUI <- function(id) {
  ns <- NS(id)
  fluidPage(
    titlePanel("Cleaned Data Formatting"),
    sidebarLayout(
      sidebarPanel(
        h4("Environments"),
        textInput(ns("env_names"), "Enter environment names (comma-separated)", placeholder = "e.g., WM_24,ARL_24"),
        h4("Processing Parameters"),
        textInput(ns("fixed"), "Fixed terms (comma-separated)", value = "id_1"),
        textInput(ns("random"), "Random terms (comma-separated)", value = "rep,range,row"),
        numericInput(ns("moisture"), "Moisture (%)", value = 15.5),
        numericInput(ns("bushel"), "Bushel weight (lb)", value = 56),
        numericInput(ns("plot_length"), "Plot length (ft)", value = 22.5),
        numericInput(ns("row_spacing"), "Row spacing (in)", value = 30),
        numericInput(ns("plot_rows"), "Number of plot rows", value = 2),
        textInput(ns("nest_terms"), "Nest terms (e.g., term1,term2;term3,term4)", value = ""),
        textInput(ns("nest_name"), "Nest names (e.g., nested1,nested2)", value = ""),
        h4("Plot Filtering"),
        textInput(ns("na_plots"), "Plot designations to set to NA (comma-separated)", value = "B"),
        h4("File Uploads"),
        uiOutput(ns("file_inputs")),
        actionButton(ns("process"), "Format Data"),
        uiOutput(ns("save_ui"))
      ),
      mainPanel(
        h4("Final data formatting, calculation of yield from plot combine values and plot dimensions, and removal of any Border/Filler plots."),
        tableOutput(ns("combined_table")),
        uiOutput(ns("file_status_ui"))
      )
    )
  )
}

formatServer <- function(input, output, session, workdir) {
  ns <- session$ns
  combined_data <- reactiveVal(NULL)
  output_path <- reactiveVal(NULL)
  
  envs <- reactive({
    req(input$env_names)
    strsplit(input$env_names, ",")[[1]] %>% trimws()
  })
  
  output$file_inputs <- renderUI({
    env_list <- envs()
    do.call(tagList, lapply(env_list, function(env) {
      list(
        h3(env),
        fileInput(ns(paste0("design_", env)), "Design File"),
        fileInput(ns(paste0("combine_", env)), "Clean File")
      )
    }))
  })
  
  observeEvent(input$process, {
    req(input$env_names)
    env_list <- envs()
    
    tryCatch({
      # Validation: Ensure random terms are provided (required parameter)
      if (input$random == "") {
        stop("Please enter random terms.")
      }
      
      # Parse fixed terms
      fixed <- if (input$fixed == "") {
        c("id_1")
      } else {
        strsplit(input$fixed, ",")[[1]] %>% trimws()
      }
      
      # Parse random terms
      random <- strsplit(input$random, ",")[[1]] %>% trimws()
      
      # Parse nest.terms (defaults to NULL if empty)
      nest_terms <- if (input$nest_terms == "") {
        NULL
      } else {
        groups <- strsplit(input$nest_terms, ";")[[1]] %>% trimws()
        lapply(groups, function(g) strsplit(g, ",")[[1]] %>% trimws())
      }
      
      # Parse nest.name (defaults to NULL if empty)
      nest_name <- if (input$nest_name == "") {
        NULL
      } else {
        strsplit(input$nest_name, ",")[[1]] %>% trimws()
      }
      
      # Parse na_plots (designations to set to NA, default "B")
      na_plots <- strsplit(input$na_plots, ",")[[1]] %>% trimws()
      
      # Process data for each environment
      dfs <- lapply(env_list, function(env) {
        design_input <- input[[paste0("design_", env)]]
        combine_input <- input[[paste0("combine_", env)]]
        req(design_input, combine_input)
        harvest.master(
          Location = env,
          Design.File = design_input$datapath,
          Combine.File = combine_input$datapath,
          Fixed = fixed,
          Random = random,
          Moisture = input$moisture,
          Bushel = input$bushel,
          Plot.Length = input$plot_length,
          Row.Spacing = input$row_spacing,
          Plot.Rows = input$plot_rows,
          nest.terms = nest_terms,
          nest.name = nest_name
        )
      })
      
      # Combine, rename, filter plots, and order the data
      formatted <- do.call(rbind, dfs) %>%
        dplyr::rename(id = id_1) %>%
        dplyr::mutate(across(!dplyr::any_of(c("env", "range", "row")),
                             ~ ifelse(plot %in% na_plots, NA, .))) %>%
        dplyr::arrange(env, range, row)
      
      # Store the result
      combined_data(formatted)
      
      # Display preview
      output$combined_table <- renderTable({
        head(combined_data())
      })
      
      # Enable save button
      output$save_ui <- renderUI({
        actionButton(ns("save_formatted"), "Save Formatted Data")
      })
      
      showNotification("Data formatted successfully.", type = "message")
    }, error = function(e) {
      showNotification(paste("Error formatting data:", e$message), type = "error")
      combined_data(NULL)
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
  
  observeEvent(input$save_formatted, {
    req(workdir(), combined_data())
    if (!dir.exists(workdir())) {
      showNotification("Selected working directory does not exist.", type = "error")
      return()
    }
    output_path_val <- file.path(workdir(), "Pheno.csv")
    tryCatch({
      write.csv(combined_data(), output_path_val, row.names = FALSE)
      output_path(output_path_val)
      showNotification(paste("File saved to:", output_path_val), type = "message")
      combined_data(NULL) # Clear to free memory
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

shinyApp(formatUI, formatServer)