#' Application Server Logic
#'
#' Defines the server-side logic for the FRDA Transcriptomic Atlas app.
#' @keywords internal
#' @import shiny
#' @importFrom ggplot2 ggplot aes geom_point
#' @noRd
app_server <- function(input, output, session, data_mode = "local") {

  # --- make `pkg` robust for both project + installed package modes ----
  pkg <- tryCatch(utils::packageName(), error = function(e) "")
  if (!length(pkg) || !is.character(pkg) || !nzchar(pkg)) pkg <- "FRDATranscriptomicAtlas"
  pkg <- pkg[[1L]]  # ensure length 1

  pkg_www <- system.file("www", package = pkg, mustWork = FALSE)
  if (nzchar(pkg_www) && dir.exists(pkg_www)) {
    shiny::addResourcePath("pkgwww", pkg_www)  # /pkgwww -> <package>/inst/www
  }

  # Also serve the project copy of inst/www at /projwww (for dev / source tree runs)
  proj_www <- file.path("inst", "www")
  if (dir.exists(proj_www)) {
    shiny::addResourcePath("projwww", proj_www)  # /projwww -> <project>/inst/www
  }

  # --- ensure pretty_map exists (fallback is harmless) ---
  if (!exists("pretty_map", inherits = TRUE)) {
    pretty_map <- stats::setNames(character(0), character(0))
  }

  # --- Call modules (pass the same `pkg`) -------------------------------
  aboutServer("about", package_name = pkg)
  datasetsServer("datasets")
  #tpmHeatmapServer("tpm_hm", pkg = pkg)
  volcanoServer("volc", level = "genes", pkg = pkg, data_mode = data_mode)
  #GSEAServer("gsea",  pkg = pkg)
  #gseaCompareServer("gsea_compare",  pkg = pkg)
  genePlotsServer("gene_plots", pkg = pkg, data_mode = data_mode)
  biomarkerServer("biomarkers", data_mode = data_mode)
  degTablesServer("deg_tables", pkg = pkg, data_mode = data_mode)
  degVennServer("deg_venn", pkg = pkg, data_mode = data_mode)
  forestPlotsServer("forest", pkg = pkg, data_mode = data_mode)
  queryGeneAcrossDatasetsServer("gene_query", data_mode = data_mode)
  pcaServer("pca", data_mode = data_mode)
  isoformConfidenceServer("isoform_confidence")
  if (data_mode == "local") {
    GSEAServer("gsea", pkg = pkg)
    gseaCompareServer("gsea_compare", pkg = pkg)
  } else {
    message("[Atlas] Cloud mode: skipping GSEA modules until cloud backend is implemented.")
  }


  # --- sidebar behavior for SwitchPlots full-width tab ------------------
  observe({
    shinyjs::toggle(
      id = "sidebar",
      condition = !(input$tabselected %in% c(0, 8))
    )
  })

  observeEvent(input$tabselected, {
    if (input$tabselected %in% c(0, 8)) {
      shinyjs::addClass(id = "main_wrap", class = "fullwidth")
    } else {
      shinyjs::removeClass(id = "main_wrap", class = "fullwidth")
    }
  })

}
