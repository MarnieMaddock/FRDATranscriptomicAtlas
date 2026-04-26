#' Launch the FRDA Transcriptomic Atlas app
#'
#' @return A shiny.appobj that runs the app
#' @export
run_app <- function(data_mode = c("local", "cloud")) {
data_mode <- match.arg(data_mode)

options(FRDATranscriptomicAtlas.data_mode = data_mode)

message("[startup] before app_ui: ", Sys.time())
ui <- app_ui(data_mode = data_mode)
message("[startup] after app_ui: ", Sys.time())

shiny::shinyApp(
  ui = ui,
  server = function(input, output, session) {
    message("[startup] before app_server: ", Sys.time())
    app_server(input, output, session, data_mode = data_mode)
    message("[startup] after app_server: ", Sys.time())
  }
)
}
# run_app <- function(data_mode = c("local", "cloud")) {
#   data_mode <- match.arg(data_mode)
#
#   options(FRDATranscriptomicAtlas.data_mode = data_mode)
#
#   shiny::shinyApp(
#     ui = app_ui(),
#     server = function(input, output, session) {
#       app_server(input, output, session, data_mode = data_mode)
#     }
#   )
# }
