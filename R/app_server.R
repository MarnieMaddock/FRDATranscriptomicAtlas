#' Application Server Logic
#'
#' Defines the server-side logic for the FRDA Transcriptomic Atlas app.
#' @keywords internal
#' @import shiny
#' @importFrom ggplot2 ggplot aes geom_point
app_server <- function(input, output, session) {

  # Build long-format DEGs from your packaged results (padj ≤ 0.05 by default)
  deg_long_df <- build_deg_long(pkg = utils::packageName(), p_thr = 0.05, lfc_min = 0, level = "genes")

  # --- Call modules -----------------------------------------------------
  #data_r <- module_data_server("data_input")
  #module_plot_server("plot_output", data = data_r)
  degTablesServer("deg_tables")
  genePlotsServer("gene_plots")
  upsetDegServer("deg_compare", deg_long = deg_long_df)
  forestPlotsServer("forest")
}
