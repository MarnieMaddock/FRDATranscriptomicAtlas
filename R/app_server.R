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

  # --- Example reactive pipeline ---------------------------------------
  observeEvent(input$run_analysis, {
    shiny::showNotification("Running analysis...", type = "message", duration = 2)
    # add actual analysis steps here later
  })

  # --- Example outputs --------------------------------------------------
  output$summary_table <- DT::renderDataTable({
    req(data_r())
    head(data_r(), 10)
  })

  output$plot_display <- renderPlot({
    req(data_r())
    ggplot2::ggplot(data_r(), ggplot2::aes(x = 1, y = 1)) +
      ggplot2::geom_point()
  })

  output$results_text <- renderPrint({
    "Results will appear here once analysis is implemented."
  })
}
