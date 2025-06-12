library(shiny)
library(data.table)
library(dplyr)
library(bigmemory)
library(BreedStream)
library(StageWise)
library(ibdsim2)
library(shinyWidgets)
library(shinyFiles)

# Source module files (adjust paths as needed)
source("~/Documents/BreedStream_App/Individual_Modules/Cleaning.R")
source("~/Documents/BreedStream_App/Individual_Modules/Format.R")
source("~/Documents/BreedStream_App/Individual_Modules/insilico.R")
source("~/Documents/BreedStream_App/Individual_Modules/TwoStage.R")
source("~/Documents/BreedStream_App/Individual_Modules/Index.R")
source("~/Documents/BreedStream_App/Individual_Modules/Effects.R")
source("~/Documents/BreedStream_App/Individual_Modules/UC.R")
source("~/Documents/BreedStream_App/Individual_Modules/OMA.R")
source("~/Documents/BreedStream_App/Individual_Modules/TransposedPredict.R")
source("~/Documents/BreedStream_App/Individual_Modules/Transposedinsilico.R")
source("~/Documents/BreedStream_App/Individual_Modules/QualityCheck.R")

# Set max upload size to 60GB
options(shiny.maxRequestSize = 60000 * 1024^2)

# Main UI
ui <- navbarPage(
  title = "BreedStream Shiny App",
  tabPanel(
    "Getting Started",
    sidebarLayout(
      sidebarPanel(
        shinyDirButton("workdir", "Select Working Directory", "Choose a directory for output files"),
        verbatimTextOutput("workdir_path")
      ),
      mainPanel(
        h4("Select a working directory to store output files, such as prediction results.")
      )
    )
  ),
  navbarMenu(
    "Clean & Format",
    tabPanel("Clean", cleanUI("clean")),
    tabPanel("Format", formatUI("format")),
    tabPanel("Quality Check", qualityUI("quality"))
  ),
  tabPanel("in silico Hybrids", insilicoUI("insilico")),
  navbarMenu(
    "Analysis",
    tabPanel("Model", twoStageUI("model")),
    tabPanel("Index", indexUI("index")),
    tabPanel("Marker Effects", effectsUI("effects"))
  ),
  navbarMenu(
    "Breeding Crosses",
    tabPanel("Usefulness Criterion", ucUI("uc")),
    tabPanel("Optimal Mate Allocation", omaUI("oma"))
  ),
  navbarMenu(
    "RR-BLUP Prediction",
    tabPanel("Transpose", insilico_transposedUI("insilico_transposed")),
    tabPanel("Prediction", predictUI("predict"))
  )
)

# Main Server
server <- function(input, output, session) {
  # Define directory chooser
  volumes <- c(Home = fs::path_home(), getVolumes()())
  shinyDirChoose(input, "workdir", roots = volumes, session = session)
  
  # Reactive value to store selected directory
  workdir <- reactive({
    if (is.integer(input$workdir)) {
      return(getwd())
    }
    parseDirPath(volumes, input$workdir)
  })
  
  # Display selected directory
  output$workdir_path <- renderText({
    workdir()
  })
  
  # Call modules, passing workdir to all modules
  callModule(cleanServer, "clean", workdir = workdir)
  callModule(formatServer, "format", workdir = workdir)
  callModule(insilicoServer, "insilico", workdir = workdir)
  callModule(twoStageServer, "model", workdir = workdir)
  callModule(indexServer, "index", workdir = workdir)
  callModule(effectsServer, "effects", workdir = workdir)
  callModule(ucServer, "uc", workdir = workdir)
  callModule(omaServer, "oma", workdir = workdir)
  callModule(predictServer, "predict", workdir = workdir)
  callModule(insilico_transposedServer, "insilico_transposed", workdir = workdir)
  callModule(qualityServer, "quality", workdir = workdir)
  
  # Clean up temporary files on session end
  session$onSessionEnded(function() {
    temp_files <- list.files(tempdir(), full.names = TRUE)
    if (length(temp_files) > 0) unlink(temp_files)
  })
}

# Run the app
shinyApp(ui, server)
