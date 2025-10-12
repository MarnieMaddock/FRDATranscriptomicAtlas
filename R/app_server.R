#' Application Server Logic
#'
#' Defines the server-side logic for the FRDA Transcriptomic Atlas app.
#' @keywords internal
#' @import shiny
#' @importFrom ggplot2 ggplot aes geom_point
app_server <- function(input, output, session) {

  # --- Call modules -----------------------------------------------------
  #data_r <- module_data_server("data_input")
  #module_plot_server("plot_output", data = data_r)
  degTablesServer("deg_tables")
  #genePlotsServer("gene_plots")


}
