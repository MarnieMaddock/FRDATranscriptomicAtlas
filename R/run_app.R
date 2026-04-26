#' Launch the FRDA Transcriptomic Atlas app
#'
#' @return A shiny.appobj that runs the app
#' @export
run_app <- function(data_mode = c("local", "cloud")) {
  data_mode <- match.arg(data_mode)

  options(FRDATranscriptomicAtlas.data_mode = data_mode)

  shiny::shinyApp(
    ui = app_ui(),
    server = function(input, output, session) {
      app_server(input, output, session, data_mode = data_mode)
    }
  )
}
