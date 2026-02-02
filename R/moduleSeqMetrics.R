# R/modules/isoform_confidence_module.R
# -------------------------------------

isoformConfidenceUI <- function(id) {
  ns <- shiny::NS(id)

  tagList(
    h3("Dataset Sequencing Metrics and Isoform Confidence Scores"),

    div(
      class = "text-muted",
      style = "margin-bottom:10px; font-size:0.9em;",
      p(
        "Isoform inference confidence was classified as low, medium, or high using a composite score based on mean sequencing depth, mean read length, and library layout (paired-end vs single-end)."
      ),
      p(
        "Mean sequencing depth was scored as ",
        tags$strong("0"), " (<10 million reads), ",
        tags$strong("1"), " (10-30 million reads), or ",
        tags$strong("2"), " (>30 million reads); ",
        "mean read length was scored as ",
        tags$strong("0"), " (<75 bp), ",
        tags$strong("1"), " (75-99 bp), or ",
        tags$strong("2"), " (>=100 bp); ",
        "and paired-end libraries were assigned ",
        tags$strong("+1"), " point relative to single-end libraries. ",
        "Component scores were summed, with higher scores indicating greater expected isoform resolution. Composite scores of 0-2 were classified as low confidence, scores of 3-4 as medium confidence, and scores of 5 as high confidence."
      )
    ),

    shinycssloaders::withSpinner(
      DT::DTOutput(ns("isoform_conf_tbl")),
      type = 4, color = "#005249"
    )
  )
}


isoformConfidenceServer <- function(id) {
  shiny::moduleServer(id, function(input, output, session) {

    # ---- master table ----
    df <- tibble::tribble(
      ~Study,                   ~`Library Type`, ~`Mean Sequencing Depth (Million)`, ~`Mean Read Length (base pairs)`, ~`Isoform Confidence Score`, ~`Isoform Confidence`,
      "Chutake YK (2014)",      "Paired-End",     25.49,                              100.00,                          4,                          "Medium",
      "Erwin GS (2017)",        "Single-End",      4.72,                               64.00,                          0,                          "Low",
      "Indelicato E (2023)",    "Single-End",      5.24,                              166.31,                          2,                          "Low",
      "Lai J (2018)",           "Single-End",     34.70,                               83.75,                          3,                          "Medium",
      "Lees JG (2025*)",        "Paired-End",     30.25,                              110.00,                          5,                          "High",
      "Li J (2019)",            "Single-End",     32.19,                               50.00,                          2,                          "Low",
      "Maddock ML (*)",         "Paired-End",    124.20,                              100.00,                          5,                          "High",
      "Mishra P (2024)",        "Paired-End",     31.06,                              100.00,                          5,                          "High",
      "Napierala JS (2017)",    "Paired-End",     32.58,                               76.00,                          4,                          "Medium",
      "Vilema-Enriquez G (2020)","Paired-End",    30.25,                               74.00,                          3,                          "Medium",
      "Wang F (2022)",          "Single-End",     24.86,                               76.00,                          2,                          "Low"
    )

    output$isoform_conf_tbl <- DT::renderDT({
      DT::datatable(
        df,
        rownames = FALSE,
        escape = TRUE,
        options = list(
          dom = "Bfrtip",
          pageLength = 11,
          buttons = c("copy", "excel"),
          scrollX = TRUE
        ),
        extensions = "Buttons"
      )
    })
  })
}
