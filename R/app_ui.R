#' Application User Interface
#'
#' Defines the UI layout for the FRDA Transcriptomic Atlas app.
#' @keywords internal
#' @import shiny
#' @import bslib
#' @importFrom shinyjs useShinyjs
#' @importFrom fontawesome fa

# logos/css
# Define helper functions for resource paths ----------------------

get_logo_path <- function() {
  if (file.exists("inst/www/dottori_lab_pentagon.svg")) {
    return("inst/www/dottori_lab_pentagon.svg")  # shinyapps.io
  } else {
    return(system.file("www", "dottori_lab_pentagon.svg", package = "FRDATranscriptomicAtlas"))  # GitHub / package
  }
}

get_UOW_path <- function() {
  if (file.exists("inst/www/UOW.png")) {
    return("inst/www/UOW.png")  # shinyapps.io
  } else {
    return(system.file("www", "UOW.png", package = "FRDATranscriptomicAtlas"))
  }
}

get_css_path <- function() {
  if (file.exists("inst/www/style.css")) {
    return("inst/www/style.css")  # shinyapps.io
  } else {
    return(system.file("www", "style.css", package = "FRDATranscriptomicAtlas"))
  }
}

get_DTU_dir <- function() {
  if (dir.exists("inst/extdata/DTU")) {
    return("inst/extdata/DTU")  # shinyapps.io / local project
  } else {
    return(system.file("extdata", "DTU", package = "FRDATranscriptomicAtlas"))
  }
}

get_switch_src <- function() {
  # Prefer project file during development
  if (file.exists("inst/www/switchplot_example.svg")) {
    return("projwww/switchplot_example.svg")   # <-- URL, not a disk path
  }
  # Fallback to installed package asset (served at /pkgwww)
  if (nzchar(system.file("www", "switchplot_example.svg",
                         package = "FRDATranscriptomicAtlas"))) {
    return("pkgwww/switchplot_example.svg")    # <-- URL, not a disk path
  }
  NULL  # nothing found
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
        id = "sidebar",
        conditionalPanel(
          condition = "input.tabselected == 1",
          PCASidebarUI("pca")
        ),
        conditionalPanel(condition = "input.tabselected==2 && input.degs_tabs == 2.1",
                         degTablesSidebarUI("deg_tables")
        ),
        conditionalPanel(condition = "input.tabselected==2 && input.degs_tabs == 2.2",
                         degVennUI("deg_venn")
        ),
        conditionalPanel(
          condition = "input.tabselected == 3",
          volcanoSidebarUI("volc", pretty_map = pretty_map)
        ),
        conditionalPanel(
          condition = "input.tabselected == 4 && input.fe_tabs === 'explore'",
          GSEASidebarUI("gsea")
        ),
        conditionalPanel(
          condition = "input.tabselected == 4 && input.fe_tabs === 'compare'",
          gseaCompareUI("gsea_compare")
        ),
        conditionalPanel(
          condition = "input.tabselected == 5",
          tpmHeatmapSidebarUI("tpm_hm")
        ),
        conditionalPanel(
          condition = "input.tabselected == 6 && input.gene_plots == 6.1",
          genePlotsSidebarUI("gene_plots")
        ),
        conditionalPanel(
          condition = "input.tabselected == 6 && input.gene_plots == 6.2",
          forestPlotsUI("forest")
        ),
        conditionalPanel(
          condition = "input.tabselected == 7 && input.DTU_tables == 7.1",
          dtuResultsSidebarUI("dtu")
        ),
        conditionalPanel(
          condition = "input.tabselected == 7 && input.DTU_tables == 7.2",
          dtuVennUI("dtuVenn")
        ),
        conditionalPanel(
          condition = "input.tabselected == 7 && input.DTU_tables == 7.3",
          consequencesSidebarUI("dtu_func_cons")
        )

      ), #sidebarPanel closing bracket

      mainPanel(id = "main_wrap",
        tabsetPanel(
          type = "tabs",
          id = "tabselected",
          selected = 0, # Default tab selected is 1
          tabPanel("Home", icon = icon("home", lib = "font-awesome"), #display home icon in the tab
                   value = 0,
                   tabsetPanel(
                     id = "degs_tabs",
                     type = "tabs",
                     tabPanel("About",
                              value = 0.1,
                              aboutUI("about"),
                     ),
                     tabPanel("Datasets",
                              value = 0.2,
                              datasetsUI("datasets"),
                     ),
                     tabPanel("Contact",
                              value = 0.3,
                     ),
                     tabPanel("Help",
                              value = 0.4,

                     )
                   )
          ),
          tabPanel("PCA", value = 1,
                   PCAMainUI("pca")
          ),
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
                              )
                   )
          ),
          tabPanel("Volcano Plots", value = 3,
                   volcanoMainUI("volc")
          ),
          tabPanel("Heatmaps", value = 5,
                   tpmHeatmapMainUI("tpm_hm")
          ),
          tabPanel(
            "Functional Enrichment", value = 4,
            tabsetPanel(
              id = "fe_tabs",
              tabPanel(
                title = "Explore by Dataset",
                value = "explore",
                GSEAMainUI("gsea")
              ),
              tabPanel(
                title = "Compare Datasets",
                value = "compare",
                gseaCompareMainUI("gsea_compare")
              )
            )
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
                              dtuResultsMainUI("dtu")
                     ),
                     tabPanel("Compare Datasets",
                              value = 7.2,
                              dtuVennMainUI("dtuVenn")
                     ),
                     tabPanel("Functional Consequences",
                              value = 7.3,
                              consequencesMainUI("dtu_func_cons")
                              )
                   )
          ),
          tabPanel(
            "SwitchPlots", value = 8,
            fluidRow(
              column(
                width = 12,
                switchplotsHelpUI("switchplots_help")   # your full-width module,
              )
            )
          )

        )
      ) #main panel close bracket
    ) #sidebarLayout close bracket
  ) #fluidPage close bracket
}
