#' Launch the FRDA Transcriptomic Atlas app
#'
#' @param data_mode Data access mode. Defaults to "cloud". Use "local" for local access of data for Developers
#'
#' @return A shiny.appobj that runs the app
#' @export
run_app <- function(data_mode = "cloud") {
  data_mode <- match.arg(data_mode, c("cloud", "local"))

  options(FRDATranscriptomicAtlas.data_mode = data_mode)

  message("[Atlas] Running in ", data_mode, " mode")

  shiny::shinyApp(
    ui = app_ui(),
    server = function(input, output, session) {
      app_server(input, output, session, data_mode = data_mode)
    }
  )
}
