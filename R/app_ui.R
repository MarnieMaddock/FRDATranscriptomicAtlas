#' Application User Interface
#'
#' Defines the UI layout for the FRDA Transcriptomic Atlas app.
#' @keywords internal
#' @import shiny
#' @import bslib
#' @importFrom shinyjs useShinyjs
#' @importFrom fontawesome fa

get_css_path <- function() {
  local_path <- "inst/www/style.css"
  pkg_path   <- system.file("www", "style.css", package = "FRDATranscriptomicAtlas")

  if (file.exists(local_path)) local_path else pkg_path
}

get_logo_path <- function() {
  if (file.exists("inst/www/dottori_lab_pentagon.svg")) {
    return("inst/www/dottori_lab_pentagon.svg") #shinyapps.io
  } else {
    return(system.file("www", "dottori_lab_pentagon.svg", package = "ProntoPCR")) #github
  }
}

get_UOW_path <- function() {
  if (file.exists("inst/www/UOW.png")) {
    return("inst/www/UOW.png") #shinyapps.io
  } else {
    return(system.file("www", "UOW.png", package = "ProntoPCR")) #github
  }
}

app_ui <- function() {

  fluidPage(
    theme = bslib::bs_theme(version = 4, bootswatch = "pulse"),

    tags$head(includeCSS(get_css_path())),

    # --- Logos / header branding -----------------------------------------
    div(id = "logo", bslib::card_image(file = get_logo_path(), fill = FALSE, width = "70px")),
    div(id = "logo2", bslib::card_image(file = get_UOW_path(), fill = FALSE, width = "220px")),

    # --- Main Layout -----------------------------------------------------
    # Application title
    div(tags$h1("FRDA Transcriptomic Atlas", style = "margin-left: 65px;")),
    shinyjs::useShinyjs(), # Use shinyjs to hide and show elements

    #sidebar options
    sidebarLayout(
      sidebarPanel(
        style = "height: 85vh; overflow-y: auto;", # Set the sidebar height and add a scroll bar
        conditionalPanel(condition = "input.tabselected==2 && input.subInput == 2.1",
                         inputFileUI("file")
        ), #insert csv file and check that it meets the required formatting, enter housekeeper names and save them
        #Calculations tab: #Delta Cq tab
        conditionalPanel(condition = "input.tabselected == 3 && input.subPanel == 3.1",
                         wrangleDataSidebar("wrangleDataModule") #display delta Cq data, average housekeepers module
        ),
        conditionalPanel(condition = "input.tabselected == 3 && input.subPanel == 3.2 && input.subCalc2 == 3",
                         ddcqSidebar("ddcqModule") #display ddcq data module
        ),
        # Statistics tab
        # Suggested workflow tab
        conditionalPanel(condition = "input.tabselected == 4",
                         statsSidebar("statsModule") #display stats options module
        ),
        conditionalPanel(condition = "input.tabselected == 5",
                         graphsSidebar("graphsModule") #display graphs options
        ),
      ), #sidebarPanel closing bracket

      mainPanel(
        tabsetPanel(
          type = "tabs",
          id = "tabselected",
          selected = 1, # Default tab selected is 1
          tabPanel("About", icon = icon("home", lib = "font-awesome"), #display home icon in the tab
                   textOutput("about"), value = 1,
                   about_text
          ), #display about text from source("about.R")
          tabPanel("Input Data", textOutput("inputdata"), value = 2,
                   tabsetPanel(
                     id = "subInput",
                     selected = 2.1, #display inserted data by the user
                     tabPanel("Data", value = 2.1,
                              #Display uploaded data using DataTable (module_deltaCq.R)
                              inputDataUI("inputDataModule"),
                     ),
                     tabPanel("Example Data", value = 2.2,
                              exampleDataUI("exampleData")
                     ) #demonstrates an example file
                   )
          ),
          #calculations tab
          tabPanel("Calculations", value = 3,
                   tabsetPanel(
                     id = "subPanel",
                     selected = 3.1,
                     tabPanel(HTML("2<sup>-(∆Cq)</sup>"), value = 3.1, #dcq tab
                              tabsetPanel(
                                id = "subCalc",
                                selected = 1,
                                tabPanel("All Data", value = 1, #dcq tab
                                         wrangleDataUI("wrangleDataModule"), # Display delta Cq data here, and filtered data (module_deltaCq.R)
                                         tags$br(),
                                         tags$br()),
                                #display biological replicate averages
                                tabPanel("Biological Replicate", value = 2,
                                         repDataUI("rep_data"),
                                )
                              )
                     ),
                     #delta delta cq tab
                     tabPanel(HTML("2<sup>-(∆∆Cq)</sup>"), value = 3.2,
                              tabsetPanel(
                                id = "subCalc2",
                                selected = 3,
                                tabPanel("All Data", value = 3,
                                         ddcqMain("ddcqModule")
                                ), #display processed ddcq data
                                tabPanel("Biological Replicate", value = 4,
                                         DDCQrepMain("ddcqRep"),
                                )
                              ),
                     ),
                   ),
          ),
          #statistics tab
          tabPanel("Statistics", value = 4,
                   statsMain("statsModule")
          ),
          tabPanel("Graphs", value = 5,
                   graphsMain("graphsModule")
          )
        )
      ) #main panel close bracket
    ) #sidebarLayout close bracket
  ) #fluidPage close bracket
