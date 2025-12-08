#' Launch the FRDA Transcriptomic Atlas app
#'
#' @return A shiny.appobj that runs the app
#' @export
run_app <- function() {
  shiny::shinyApp(ui = app_ui(), server = app_server)
}
