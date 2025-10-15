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
        conditionalPanel(condition = "input.tabselected==2 && input.degs_tabs == 2.1",
                         degTablesSidebarUI("deg_tables")
        ),
        conditionalPanel(condition = "input.tabselected==2 && input.degs_tabs == 2.2",
                         degVennUI("deg_venn")
        ),
        conditionalPanel(
          condition = "input.tabselected == 6 && input.gene_plots == 6.1",
          genePlotsSidebarUI("gene_plots")
        ),
        conditionalPanel(
          condition = "input.tabselected == 6 && input.gene_plots == 6.2",
          forestPlotsUI("forest")
        )


      ), #sidebarPanel closing bracket

      mainPanel(
        tabsetPanel(
          type = "tabs",
          id = "tabselected",
          selected = 1, # Default tab selected is 1
          tabPanel("About", icon = icon("home", lib = "font-awesome"), #display home icon in the tab
                   #textOutput("about"),
                   value = 1,
                   #include daatsets info here
                   #about_text
          ), #display about text from source("about.R")
          tabPanel("DEGs",
                   value = 2,
                   tabsetPanel(
                     id = "degs_tabs",
                     type = "tabs",
                     tabPanel("Explore by Dataset",
                              value = 2.1,
                              degTablesMainUI("deg_tables")
                              ),
                     tabPanel("Compare Datasets",
                              value = 2.2,
                              degVennMainUI("deg_venn")
                              #upsetDegUI("deg_compare")
                              )
                   )
          ),
          #PCA graphs
          tabPanel("PCA", value = 3,
                   ),
          #statistics tab
          tabPanel("Functional Enrichment", value = 4,
          ),
          tabPanel("Volcano Plots", value = 5,
          ),
          #Boxplots or violin plots per gene using {ggplotly} after user selects a gene.
          tabPanel("Gene Plots", value = 6,
                   #sub tab with forrest plots to compare across studies
                   tabsetPanel(
                     id = "gene_plots",
                     type = "tabs",
                     tabPanel("Explore by Dataset",
                              value = 6.1,
                              genePlotsMainUI("gene_plots")
                     ),
                     tabPanel("Compare Datasets",
                              value = 6.2,
                              forestPlotMainUI("forest")
                     )
                   )

          ),
          tabPanel("DTU", value = 7,
                   tabsetPanel(
                     id = "DTU_tables",
                     type = "tabs",
                     tabPanel("Explore by Dataset",
                              value = 7.1,
                              #degTablesUI("deg_tables")
                     ),
                     tabPanel("Compare Datasets",
                              value = 7.2,
                              #degCompareUI("deg_compare")
                     )
                   )
          ),
          tabPanel("SwitchPlots", value = 8,
          ),
        )
      ) #main panel close bracket
    ) #sidebarLayout close bracket
  ) #fluidPage close bracket
}
