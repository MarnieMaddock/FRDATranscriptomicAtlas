#' Application User Interface
#'
#' Defines the UI layout for the FRDA Transcriptomic Atlas app.
#' @import shiny
#' @import bslib
#' @importFrom shinyjs useShinyjs
#' @importFrom fontawesome fa
#' @noRd

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


app_ui <- function(data_mode = "cloud") {
  data_mode <- match.arg(data_mode, c("cloud", "local"))

  fluidPage(
    theme = bslib::bs_theme(version = 4, bootswatch = "pulse"),

    tags$head(
      includeCSS(get_css_path()),

      tags$script(
        async = NA,
        src = "https://www.googletagmanager.com/gtag/js?id=G-G7LP80W47V"
      ),

      tags$script(HTML("
        window.dataLayer = window.dataLayer || [];
        function gtag(){dataLayer.push(arguments);}
        gtag('js', new Date());
        gtag('config', 'G-G7LP80W47V', {
          'page_title': 'FRDA Transcriptomic Atlas'
        });
      "))
    ),

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
        conditionalPanel(condition = "input.tabselected==2 && input.degs_tabs == 2.3",
                         queryGeneAcrossDatasetsSidebarUI("gene_query")
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
          condition = "input.tabselected == 9",
          biomarkerUI("biomarkers")
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
                     tabPanel("Sequencing Metrics",
                              value = 0.3,
                              isoformConfidenceUI("isoform_confidence")
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
                              ),
                     tabPanel("Query Gene Across Datasets",
                              value = 2.3,
                              queryGeneAcrossDatasetsMainUI("gene_query")
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
          tabPanel(
            "Biomarker Discovery", value = 9,
            biomarkerMainUI("biomarkers")
          )

        )
      ) #main panel close bracket
    ) #sidebarLayout close bracket
  ) #fluidPage close bracket
}
