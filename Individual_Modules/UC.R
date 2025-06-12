gc()

library(shiny)
library(dplyr)
library(StageWise)
library(BreedStream)
library(ibdsim2)

ucUI <- function(id) {
  ns <- NS(id)
  fluidPage(
    titlePanel("Usefulness Criterion"),
    sidebarLayout(
      sidebarPanel(
        numericInput(ns("num_chrom"), "Number of Chromosomes", value = 10, min = 1, step = 1),
        uiOutput(ns("chrom_map_inputs")),
        textInput(ns("group_names"), "Group Names (comma-separated)", placeholder = "e.g., SS,NSS"),
        uiOutput(ns("geno_file_inputs")),
        fileInput(ns("marker_effects_file"), "Marker Effects .rds File", accept = ".rds"),
        numericInput(ns("ploidy"), "Ploidy Level:", value = 2, min = 2, step = 2),
        radioButtons(ns("map"), "Map Included:", choices = list("Yes" = "TRUE", "No" = "FALSE"), selected = "TRUE"),
        numericInput(ns("min.minor.allele"), "Minor Allele Threshold:", value = 1, min = 1, step = 1),
        radioButtons(ns("dominance"), "Include Dominance:", choices = list("Yes" = "TRUE", "No" = "FALSE"), selected = "TRUE"),
        actionButton(ns("processBtn"), "Calculate Marker Effects"),
        uiOutput(ns("save_ui"))
      ),
      mainPanel(
        h4("Using an input for each chromosome that contains the physical position and the genetic position (in centimorgans), recombination is simulated between markers for the calculation of a usefullness criterion (the expected performance of the superior progeny). This usefullness criterion is used as the merit for Optimal Mate Allocation. For maize, physical and genetic position information can be navigated to on a chromosome basis from this webpage: https://www.maizegdb.org/data_center/map"),
        uiOutput(ns("group_tables")),
        uiOutput(ns("file_status_ui"))
      )
    )
  )
}

ucServer <- function(input, output, session, workdir) {
  ns <- session$ns
  group_data <- reactiveVal(NULL)
  output_paths <- reactiveVal(NULL)
  
  output$chrom_map_inputs <- renderUI({
    num_chrom <- input$num_chrom
    lapply(1:num_chrom, function(i) {
      fileInput(ns(paste0("chrom_map_", i)), paste("Chromosome", i, "Map File"), accept = ".txt")
    })
  })
  
  output$geno_file_inputs <- renderUI({
    req(input$group_names)
    group_names <- trimws(strsplit(input$group_names, ",")[[1]])
    lapply(group_names, function(group) {
      fileInput(ns(paste0("geno_file_", group)), paste(group, "Genotype File"), accept = ".csv")
    })
  })
  
  chrom_maps <- reactive({
    num_chrom <- input$num_chrom
    maps <- list()
    tryCatch({
      for (i in 1:num_chrom) {
        file <- input[[paste0("chrom_map_", i)]]
        if (!is.null(file)) {
          maps[[i]] <- ConvertToMap(read.delim(file$datapath, skip = 1))
        } else {
          maps[[i]] <- NULL
        }
      }
      maps
    }, error = function(e) {
      showNotification(paste("Error loading chromosome map files:", e$message), type = "error")
      stop(e$message)
    })
  })
  
  add_marker_effects <- reactive({
    req(input$marker_effects_file)
    tryCatch({
      marker_effects <- readRDS(input$marker_effects_file$datapath)
      marker_effects[["add_effects"]]
    }, error = function(e) {
      showNotification(paste("Error loading marker effects file:", e$message), type = "error")
      stop(e$message)
    })
  })
  
  observeEvent(input$processBtn, {
    req(input$group_names, input$marker_effects_file)
    tryCatch({
      group_names <- trimws(strsplit(input$group_names, ",")[[1]])
      chrom_maps_list <- chrom_maps()
      add_eff <- add_marker_effects()
      map_logical <- as.logical(input$map)
      dominance_logical <- as.logical(input$dominance)
      data_list <- list()
      
      for (group in group_names) {
        geno_file <- input[[paste0("geno_file_", group)]]
        req(geno_file)
        
        # Step 1: Read and filter genotype data
        genoI <- StageWise::read_geno(
          filename = geno_file$datapath,
          ploidy = input$ploidy,
          map = map_logical,
          min.minor.allele = input$min.minor.allele,
          dominance = dominance_logical
        )
        
        # Step 2: Extract map and convert positions to cM
        Map <- genoI@map
        Map$chrom <- gsub("chr", "", Map$chrom)
        Map$chrom <- as.integer(Map$chrom)
        
        for (row in 1:nrow(Map)) {
          chrom_num <- Map[row, "chrom"]
          if (!is.null(chrom_maps_list[[chrom_num]])) {
            Map[row, "cM"] <- convertPos(
              Mb = Map[row, "position"] / 1000000,
              map = chrom_maps_list[[chrom_num]]
            )
          } else {
            Map[row, "cM"] <- Map[row, "position"] / 1000000 # Fallback
          }
        }
        
        Map <- Map %>%
          dplyr::select(all_of(c("chrom", "cM", "marker"))) %>%
          dplyr::rename(chr = chrom, pos = cM, mkr = marker)
        if (any(is.na(Map$pos))) {
          warning("Markers with missing cM values for group ", group, ": ", paste(head(Map$mkr[is.na(Map$pos)], 5), collapse = ", "))
          Map <- Map[!is.na(Map$pos), ]
        }
        
        # Step 3: Read and format genotype matrix
        markers_full <- BreedStream::format_markers_usef_add(geno_file = geno_file$datapath)
        
        # Step 4: Find common markers
        common_markers <- intersect(colnames(markers_full), Map$mkr)
        common_markers <- intersect(common_markers, add_eff$marker)
        
        if (length(common_markers) == 0) {
          stop(paste("No common markers found for group", group))
        }
        
        # Step 5: Subset Map to common markers
        Map <- Map[Map$mkr %in% common_markers, ]
        
        # Step 6: Sort Map by chromosome and position
        Map <- Map[order(Map$chr, Map$pos), ]
        
        # Step 7: Subset and reorder markers to match the sorted Map$mkr
        markers <- markers_full[, Map$mkr, drop = FALSE]
        
        # Step 8: Subset and reorder add_eff to match the sorted Map$mkr
        add_eff_group <- add_eff[add_eff$marker %in% Map$mkr, ]
        add_eff_group <- add_eff_group[match(Map$mkr, add_eff_group$marker), ]
        if (any(is.na(add_eff_group$marker))) {
          stop("Some markers in Map are not in add_eff for group ", group)
        }
        
        # Name addEff to match Map$mkr
        addEff <- add_eff_group$effect
        names(addEff) <- add_eff_group$marker
        
        # Step 9: Validate marker alignment
        if (!all(colnames(markers) == Map$mkr) || !all(names(addEff) == Map$mkr)) {
          stop("Marker ordering mismatch for group ", group)
        }
        
        # Step 10: Generate relationship matrix and crossing plan
        relMat <- as.matrix(genoI@G)
        Cross_plan <- BreedStream::planCross(
          TargetPop = rownames(markers),
          MateDesign = "half"
        )
        Cross_plan <- Cross_plan[
          Cross_plan$Parent1 %in% rownames(markers) &
            Cross_plan$Parent2 %in% rownames(markers),
        ]
        
        # Step 11: Calculate usefulness criterion
        usef_add <- BreedStream::getUsefA(
          MatePlan = Cross_plan,
          Markers = markers,
          addEff = addEff,
          Map.In = Map,
          K = relMat,
          propSel = 0.05,
          Type = "DH",
          Generation = 1
        )
        
        data_list[[group]] <- list(relMat = relMat, usef_add = usef_add)
      }
      
      group_data(data_list)
      
      output$save_ui <- renderUI({
        actionButton(ns("save_all"), "Save All Results")
      })
      
      output$group_tables <- renderUI({
        lapply(names(data_list), function(group) {
          list(
            h4(paste("Group:", group)),
            tableOutput(ns(paste0("table_", group)))
          )
        })
      })
      
      for (group in names(data_list)) {
        output[[paste0("table_", group)]] <- renderTable({
          head(group_data()[[group]]$usef_add[[2]])
        })
      }
      
      showNotification("Usefulness Criterion calculated successfully.", type = "message")
    }, error = function(e) {
      showNotification(paste("Error calculating Usefulness Criterion:", e$message), type = "error")
      group_data(NULL)
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
    req(workdir(), group_data())
    if (!dir.exists(workdir())) {
      showNotification("Selected working directory does not exist.", type = "error")
      return()
    }
    paths <- list()
    tryCatch({
      # Save UC.rds
      rds_path <- file.path(workdir(), "UC.rds")
      saveRDS(group_data(), rds_path)
      paths[["UC.rds"]] <- rds_path
      
      # Save per-group usef_add.csv files
      for (group in names(group_data())) {
        group_path <- file.path(workdir(), paste0(group, "_usef_add.csv"))
        data.table::fwrite(group_data()[[group]]$usef_add[[2]], group_path)
        paths[[paste0(group, "_usef_add.csv")]] <- group_path
      }
      
      output_paths(paths)
      showNotification("All files saved successfully.", type = "message")
      group_data(NULL)
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

shinyApp(ucUI, ucServer)
