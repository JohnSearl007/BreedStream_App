gc()

library(shiny)
library(data.table)
library(dplyr)

cleanUI <- function(id) {
  ns <- NS(id)
  fluidPage(
    titlePanel("Phenotypic Data Cleaning"),
    sidebarLayout(
      sidebarPanel(
        fileInput(ns("pheno_file"), "Upload Phenotypic Data (CSV)", accept = ".csv"),
        textInput(ns("location"), "Enter Location Name:", value = "Unknown"),
        uiOutput(ns("trait_selector")),
        uiOutput(ns("cleaning_ui")),
        uiOutput(ns("post_cleaning_ui"))
      ),
      mainPanel(
        h4("Cleaning of data on a single environment basis using visual inspection for outliers."),
        imageOutput(ns("histogram"), height = "400px", width = "600px"),
        uiOutput(ns("file_status_ui"))
      )
    )
  )
}

cleanServer <- function(input, output, session, workdir) {
  ns <- session$ns
  data <- reactiveVal(NULL)
  selected_traits <- reactiveVal(NULL)
  current_trait_index <- reactiveVal(1)
  bounds_list <- reactiveValues()
  cleaned_data <- reactiveVal(NULL)
  final_data <- reactiveVal(NULL)
  output_path <- reactiveVal(NULL)
  
  observeEvent(input$pheno_file, {
    req(input$pheno_file)
    df <- fread(input$pheno_file$datapath)
    data(df)
    cleaned_data(NULL)
    final_data(NULL)
    numeric_cols <- names(df)[sapply(df, is.numeric)]
    output$trait_selector <- renderUI({
      checkboxGroupInput(ns("traits"), "Select traits to clean:", choices = numeric_cols)
    })
  })
  
  observeEvent(input$traits, {
    req(input$traits)
    selected_traits(input$traits)
    cleaned_data(NULL)
    final_data(NULL)
    for (trait in selected_traits()) bounds_list[[trait]] <- NULL
    current_trait_index(1)
    update_cleaning_ui()
  })
  
  update_cleaning_ui <- function() {
    output$cleaning_ui <- renderUI({
      req(selected_traits(), current_trait_index() <= length(selected_traits()))
      trait <- selected_traits()[current_trait_index()]
      min_val <- min(data()[[trait]], na.rm = TRUE)
      max_val <- max(data()[[trait]], na.rm = TRUE)
      list(
        h4(paste("Cleaning trait:", trait)),
        sliderInput(ns("lower_bound"), "Lower bound:", min = min_val, max = max_val, value = min_val),
        sliderInput(ns("upper_bound"), "Upper bound:", min = min_val, max = max_val, value = max_val),
        actionButton(ns("confirm"), "Confirm bounds")
      )
    })
  }
  
  output$histogram <- renderImage({
    req(data(), selected_traits(), current_trait_index() <= length(selected_traits()))
    trait <- selected_traits()[current_trait_index()]
    lower <- input$lower_bound
    upper <- input$upper_bound
    if (is.null(lower) || is.null(upper)) return(list(src = ""))
    temp <- data()[[trait]][data()[[trait]] >= lower & data()[[trait]] <= upper & !is.na(data()[[trait]])]
    outfile <- tempfile(fileext = ".png")
    png(outfile, width = 600, height = 400, res = 96)
    if (length(temp) == 0) {
      plot.new()
      text(0.5, 0.5, "No data within selected bounds", cex = 1.5)
    } else {
      hist(temp, main = paste("Histogram of", trait), xlab = trait, col = "lightblue", border = "black")
    }
    dev.off()
    list(src = outfile, contentType = "image/png", width = 600, height = 400, alt = "Histogram")
  }, deleteFile = TRUE)
  
  observeEvent(input$confirm, {
    trait <- selected_traits()[current_trait_index()]
    bounds_list[[trait]] <- c(input$lower_bound, input$upper_bound)
    if (current_trait_index() < length(selected_traits())) {
      current_trait_index(current_trait_index() + 1)
      update_cleaning_ui()
    } else {
      df <- data()
      for (trait in selected_traits()) {
        bounds <- bounds_list[[trait]]
        df[[trait]][df[[trait]] < bounds[1] | df[[trait]] > bounds[2]] <- NA
      }
      cleaned_data(df)
    }
  })
  
  output$post_cleaning_ui <- renderUI({
    if (is.null(cleaned_data())) return(NULL)
    list(
      h4("Manual NA Assignment"),
      fileInput(ns("na_file"), "Upload NA list (CSV: Id 1, trait)", accept = ".csv"),
      textAreaInput(ns("na_manual"), "Or enter manually (Id 1,trait per line):", rows = 3),
      actionButton(ns("apply_na"), "Apply NA"),
      h4("Select Columns to Retain"),
      checkboxGroupInput(ns("retain_columns"), "Select columns:",
                         choices = names(cleaned_data()),
                         selected = intersect(c("Date/Time", "Range", "Row", "Id 1", selected_traits()), names(cleaned_data()))),
      actionButton(ns("confirm_columns"), "Confirm Columns"),
      uiOutput(ns("save_ui"))
    )
  })
  
  output$file_status_ui <- renderUI({
    if (is.null(output_path())) return(NULL)
    tagList(
      h4("File Status"),
      p(paste("File saved to:", output_path())),
      actionButton(ns("open_dir"), "Open Working Directory")
    )
  })
  
  output$save_ui <- renderUI({
    if (!is.null(final_data())) {
      actionButton(ns("save_cleaned"), "Save Cleaned Data")
    }
  })
  
  observeEvent(input$apply_na, {
    req(cleaned_data())
    df <- cleaned_data()
    if (!is.null(input$na_file)) {
      na_list <- fread(input$na_file$datapath)
      req(all(c("Id 1", "trait") %in% names(na_list)))
      for (i in 1:nrow(na_list)) {
        id <- na_list$`Id 1`[i]
        trait <- na_list$trait[i]
        if (trait %in% names(df)) df[[trait]][df$`Id 1` == id] <- NA
      }
    }
    if (!is.null(input$na_manual) && nchar(input$na_manual) > 0) {
      na_lines <- strsplit(input$na_manual, "\n")[[1]]
      for (line in na_lines) {
        parts <- strsplit(trimws(line), ",")[[1]]
        if (length(parts) == 2) {
          id <- parts[1]
          trait <- parts[2]
          if (trait %in% names(df)) df[[trait]][df$`Id 1` == id] <- NA
        }
      }
    }
    cleaned_data(df)
    final_data(NULL)
  })
  
  observeEvent(input$confirm_columns, {
    req(cleaned_data(), input$retain_columns)
    if (length(input$retain_columns) == 0) {
      showNotification("Please select at least one column to retain.", type = "error")
      return()
    }
    final_data(cleaned_data() %>% dplyr::select(all_of(input$retain_columns)))
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
  
  observeEvent(input$save_cleaned, {
    req(workdir(), final_data())
    if (!dir.exists(workdir())) {
      showNotification("Selected working directory does not exist.", type = "error")
      return()
    }
    output_path_val <- file.path(workdir(), paste0(gsub("[^a-zA-Z0-9]", "_", input$location), "_cleaned_data.csv"))
    tryCatch({
      fwrite(final_data(), output_path_val, quote = FALSE, row.names = FALSE)
      output_path(output_path_val)
      showNotification(paste("File saved to:", output_path_val), type = "message")
      final_data(NULL)
      cleaned_data(NULL)
      data(NULL)
      selected_traits(NULL)
      current_trait_index(1)
      for (trait in names(bounds_list)) bounds_list[[trait]] <- NULL
    }, error = function(e) {
      showNotification(paste("Error saving file:", e$message), type = "error")
    })
  })
  
  # Ensure garbage collection when server function exits
  on.exit(gc())
}

shinyApp(cleanUI, cleanServer)
