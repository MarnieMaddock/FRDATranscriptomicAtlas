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

      p(
        "The Friedreich's Ataxia (FRDA) Transcriptomic Atlas is an interactive Shiny application ",
        "that provides open access to harmonised transcriptomic analyses across publicly available FRDA RNA-seq datasets. ",
        "This resource was developed to enable researchers, clinicians, and students to explore and compare genes, differential expression, isoforms etc across a broad range of human FRDA models. The Atlas serves as a central reference for transcriptome-wide alterations in FRDA,
        supporting hypothesis generation, candidate gene prioritisation, and data reuse for downstream integrative analyses. ",
        "All analyses were performed using a consistent bioinformatic framework to ensure comparability across studies.",
          HTML(
    "Details of the datasets, quality control, and analytical pipelines are available in the associated publication <a href='https://doi.org/10.64898/2026.03.18.712785' target='_blank'>here</a>. "
  ),

        "Users are encouraged to consult the paper for methodological details and to cite it in any work derived from or informed by this resource."
      ),

  div(
    style = "
      background-color: #e8f4f8;
      border-left: 4px solid #2a776d;
      padding: 10px;
      margin-bottom: 15px;
      font-size: 13px;
    ",
    strong("Data loading notice"),
    br(),
    "Datasets are retrieved from external servers and may take 20 seconds to 3 minutes to load, with longer delays possible during periods of high demand."
  ),

      h3("How to Use"),
      p("All instructions on how to use the app and its features can be found here: ",
        tags$a(
          href = "https://marniemaddock.github.io/FRDATranscriptomicAtlas/",
          "Instructions"
        )),

      h3("Citation and Use"),
      p("If you use this Atlas or any underlying data/analyses in your work, ",
        "please cite the associated publication:"),
      tags$blockquote(
        HTML(
          "Maddock, M. L. <i>et&nbsp;al.</i> (2026). <i>Friedreich ataxia transcriptomic dysregulation and identification of cell type-specific biomarkers: A systematic review and meta-analysis</i>. bioRxiv. https://doi.org/10.64898/2026.03.18.712785"
        )
      ),
      p("Users must cite the publication in any derivative analyses, figures, or reports generated from this resource. Users must also cite the original studies from which the data were derived. ",
        "Please refer to the 'Datasets' section within the app for details on each dataset and its original publication."
      ),

      h3("Future Updates and Contributions"),
      p(
        "If you are interested in contributing data to the FRDA Transcriptomic Atlas, please get in touch.",
        "For contribution guidelines and further information, see the ",
        tags$a(
          href = "https://marniemaddock.github.io/FRDATranscriptomicAtlas/contribute.html",
          "contribution page",
          target = "_blank"
        )
      ),

      h3("Contact and Bug Reports"),
      tags$ul(
        tags$li(HTML(paste0("<b>Contact:</b> ", contact_email))),
        tags$li(HTML(paste0("<b>Project GitHub:</b> <a href='", github_url, "' target='_blank'>", github_url, "</a>"))),
        tags$li(HTML(paste0("<b>Report an issue or request features:</b> <a href='", issues_url, "' target='_blank'>", issues_url, "</a>")))
      ),

      h3("Authorship and Acknowledgements"),
      p(HTML(
        paste0(
          "Developed by: Marnie Maddock, Prof. Mirella Dottori, University of Wollongong, Australia<br/>",
          "<b>We thank all researchers who generated and made publicly available the original datasets used in this atlas.</b> "
        )
      )),

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


