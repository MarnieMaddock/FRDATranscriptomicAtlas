#' About module UI
#' @param id Shiny module id
#' @param title Title text
#' @param github_url Project URL
#' @param issues_url Issues URL
#' @param contact_email Contact email
#' @noRd
#' @importFrom shiny NS tagList div h2 p h3 uiOutput span moduleServer renderUI tags HTML
aboutUI <- function(id,
                    title = "About the FRDA Transcriptomic Atlas",
                    github_url = "https://github.com/MarnieMaddock/FRDATranscriptomicAtlas",
                    issues_url = "https://github.com/MarnieMaddock/FRDATranscriptomicAtlas/issues",
                    contact_email = "mlm715@uowmail.edu.au") {
  ns <- NS(id)
  tagList(
    div(
      class = "container-fluid",
      h2(title, class = "mb-3"),
      uiOutput(ns("pkg_version")),
      uiOutput(ns("data_updates")),
      actionButton(ns("check_updates"), "Check for data updates"),

      p(
        "The Friedreich's Ataxia (FRDA) Transcriptomic Atlas is an interactive Shiny application ",
        "that provides open access to harmonised transcriptomic analyses across publicly available FRDA RNA-seq datasets. ",
        "This resource was developed to enable researchers, clinicians, and students to explore and compare genes, differential expression, isoforms etc across a broad range of human FRDA models. The Atlas serves as a central reference for transcriptome-wide alterations in FRDA,
        supporting hypothesis generation, candidate gene prioritisation, and data reuse for downstream integrative analyses. ",
        "All analyses were performed using a consistent bioinformatic framework to ensure comparability across studies. Details of the datasets, quality control, and analytical pipelines are available in the associated publication <GIVE LINK HERE>.
        Users are encouraged to consult the paper for methodological details and to cite it in any work derived from or informed by this resource."
      ),

      h3("Features"),
      tags$ul(
        tags$li(
          strong("Principal Component Analysis (PCA)"), " - Inspect global sample structure and separation by condition, tissue/cell type, or other metadata."
        ),
        tags$li(
          strong("Differential Expression (DEG/DEI) Tables"), " - Gene- and isoform-level results ",
          "with log2 fold-change and adjusted p-values, filterable by direction and significance."
        ),
        tags$li(
          strong("Volcano Plots"), " - Fold-change vs. significance visualisation with options to highlight genes and adjust thresholds."
        ),
        tags$li(
          strong("Venn Diagrams"), " - Overlap of differentially expressed genes or isoforms across datasets."
        ),
        tags$li(
          strong("Expression Heatmaps"), " - Compare transcript abundance patterns across conditions and datasets."
        ),
        tags$li(
          strong("Gene TPM Plots (Gene Plots)"), " - Visualise expression of selected genes (TPM) by condition, replicate, and dataset for direct comparison."
        ),
        tags$li(
          strong("Cross-Dataset Forest Plots"), " - Compare per-dataset effect sizes for selected genes to assess concordance and heterogeneity across studies."
        ),
        tags$li(
          strong("Download and Export"), " - Tables and plots can be exported for reporting and reuse."
        )
      ),

      h3("How to Use"),
      p("All instructions on how to download the app and use its features can be found here: ",
        tags$a(
          href = "https://marniemaddock.github.io/FRDATranscriptomicAtlas/",
          "Instructions"
        )),

      h3("Authorship and Acknowledgements"),
      p(HTML(
        paste0(
          "<b>Developed by:</b> Marnie Maddock, Prof. Mirella Dottori, XXXX, XXXXX University of Wollongong, Australia<br/>",
          "<b>Funding:</b> Friedreich's Ataxia Research Alliance (FARA)<br/>",
          "<b>We thank all researchers who generated and made publicly available the original datasets used in this meta-analysis.</b> "
        )
      )),

      h3("Citation and Use"),
      p("If you use this Atlas or any underlying data/analyses in your work, ",
        "please cite the associated publication:"),
      tags$blockquote(
        HTML(
          "Maddock, M. <i>et&nbsp;al.</i> (2025). <i>XXXXXXXXXXXXXXXXXX</i>. [Journal details forthcoming]."
        )
      ),
      p("Users must cite the publication in any derivative analyses, figures, or reports generated from this resource. Users must also cite the original studies from which the data were derived. ",
        "Please refer to the 'Datasets' section within the app for details on each dataset and its original publication."
      ),

      h3("Future Updates and Contributions"),
      p(
        "This resource may be updated as new FRDA RNA-seq datasets become available. ",
        "Researchers who wish to contribute data or suggest additions are encouraged to contact the authors ",
        tags$a(
          href = "mailto:mlm715@uowmail.edu.au,mmaddock@uow.edu.au,mdottori@uow.edu.au?subject=FRDA%20RNA-seq%20Atlas%20Contribution",
          "via email"
        ),
        ". "
      ),

      h3("Contact and Bug Reports"),
      tags$ul(
        tags$li(HTML(paste0("<b>Contact:</b> ", contact_email))),
        tags$li(HTML(paste0("<b>Project GitHub:</b> <a href='", github_url, "' target='_blank'>", github_url, "</a>"))),
        tags$li(HTML(paste0("<b>Report an issue or request features:</b> <a href='", issues_url, "' target='_blank'>", issues_url, "</a>")))
      ),
      br(),
      #add footer.svg
      div(
        class = "text-center",
        bslib::card_image(
          file  = system.file("www", "footer.svg", package = "FRDATranscriptomicAtlas"),
          fill  = FALSE,
          width = "1200px"
        )
      ),
      br(),
      br()
    )
  )
}

aboutServer <- function(id, package_name = "FRDATranscriptomicAtlas") {
  moduleServer(id, function(input, output, session) {

    output$pkg_version <- renderUI({
      v <- tryCatch(utils::packageDescription(package_name)$Version, error = function(e) NULL)
      if (is.null(v)) return(NULL)
      shiny::tags$span(
        class = "text-muted",
        style = "display:block;margin-bottom:10px;",
        paste0("Version: ", v)
      )
    })

    # Cache results so we don't keep re-checking
    updates <- shiny::reactiveVal(NULL)
    check_err <- shiny::reactiveVal(NULL)
    checked_manually <- shiny::reactiveVal(FALSE)

    # Run ONCE when About tab module initializes
    shiny::observeEvent(TRUE, {
      out <- tryCatch(
        check_atlas_updates(package = package_name),
        error = function(e) { check_err(conditionMessage(e)); NULL }
      )
      updates(out)
    }, once = TRUE)

    # Optional: manual re-check button
    shiny::observeEvent(input$check_updates, {
      check_err(NULL)
      checked_manually(TRUE)
      out <- tryCatch(
        check_atlas_updates(package = package_name),
        error = function(e) {
          check_err(conditionMessage(e))
          NULL
          }
      )
      updates(out)
    })


    output$data_updates <- shiny::renderUI({

      # Error case
      if (!is.null(check_err())) {
        return(shiny::tags$div(
          class = "text-warning",
          paste0("Could not check for data updates: ", check_err())
        ))
      }

      up <- updates()

      # Manual check, no updates
      if (checked_manually() && (is.null(up) || !nrow(up))) {
        return(shiny::tags$div(
          class = "text-success",
          "Atlas data is up to date."
        ))
      }

      # Updates available
      if (!is.null(up) && nrow(up)) {
        return(shiny::tags$div(
          class = "text-warning",
          shiny::tags$strong("Atlas data updates available: "),
          paste(up$key, collapse = ", ")
        ))
      }

      # Default / initial state
      shiny::tags$div(class = "text-muted", "Atlas data: up to date.")
    })


  })
}

