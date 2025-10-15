#' Application Server Logic
#'
#' Defines the server-side logic for the FRDA Transcriptomic Atlas app.
#' @keywords internal
#' @import shiny
#' @importFrom ggplot2 ggplot aes geom_point
app_server <- function(input, output, session) {

  addResourcePath("pkgwww", system.file("www", package = "FRDATranscriptomicAtlas"))

  # Build long-format DEGs from your packaged results (padj ≤ 0.05 by default)
  deg_long_df <- build_deg_long(pkg = utils::packageName(), p_thr = 0.05, lfc_min = 0, level = "genes")

  # --- Call modules -----------------------------------------------------
  degTablesServer("deg_tables")
  genePlotsServer("gene_plots")
  degVennServer("deg_venn", pkg = utils::packageName())
  volcanoServer("volc", level = "genes", pkg = "FRDATranscriptomicAtlas", custom_loader = deseq_loader)
  forestPlotsServer("forest")

  # ---- In server.R or inside your app_server() ----
  observe({
    shinyjs::toggle(id = "sidebar", condition = !(input$tabselected == 8))
  })

  observeEvent(input$tabselected, {
    # tab values arrive as strings
    if (identical(input$tabselected, "8")) {
      shinyjs::addClass(id = "main_wrap", class = "fullwidth")
    } else {
      shinyjs::removeClass(id = "main_wrap", class = "fullwidth")
    }
  })


}
