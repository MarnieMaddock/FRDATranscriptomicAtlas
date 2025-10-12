# R/moduleGenePlots.R

#' Gene Plots – sidebar
#' @noRd
genePlotsSidebarUI <- function(id) {
  ns <- NS(id)
  tagList(
    h4("Gene plots"),
    selectizeInput(ns("gp_dataset"), "Dataset", choices = NULL, options = list(placeholder = "Select a dataset…")),
    selectizeInput(ns("gp_gene"), "Gene", choices = NULL, selected = "FXN"),
    checkboxInput(ns("gp_logy"), "Log10 Y-axis", value = FALSE),
    tags$hr(),
    strong("Download"),
    br(),
    downloadButton(ns("dl_points"), "Replicates (CSV)"),
    downloadButton(ns("dl_summary"), "Summary (CSV)"),
    downloadButton(ns("dl_plot"), "Plot (PNG)")
  )
}

#' Gene Plots – main
#' @noRd
genePlotsMainUI <- function(id) {
  ns <- NS(id)
  tagList(
    shinycssloaders::withSpinner(plotOutput(ns("gp_plot"), height = 420), type = 4),
    br(),
    DT::dataTableOutput(ns("gp_table"))
  )
}

#' Gene Plots server
#' @noRd
genePlotsServer <- function(id, pkg = utils::packageName()) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # ---- load TPM wide CSV once ----
    tpm_path <- system.file("extdata/tpm/Combined_gene_tpm.csv",
                            package = pkg, mustWork = TRUE)
    tpms_wide <- readr::read_csv(tpm_path, show_col_types = FALSE)

    # ---- reshape to long + parse sample names ----
    # columns: gene_id, gene_name, then many sample columns "Dataset_COND_[extra]_REP#(.dup)?"
    tpms_long <- tpms_wide |>
      tidyr::pivot_longer(
        cols = -c(gene_id, gene_name),
        names_to = "sample",
        values_to = "TPM"
      ) |>
      # extract dataset, condition (first token after dataset), replicate number
      tidyr::extract(
        sample,
        into  = c("dataset", "condition", "rep"),
        regex = "^(.+?)_([A-Za-z]+)(?:_[A-Za-z0-9]+)*_REP(\\d+)(?:\\..+)?$",
        remove = FALSE
      ) |>
      dplyr::mutate(
        dataset  = as.character(dataset),
        condition = as.character(condition),
        rep = as.integer(rep)
      )

    # ---- populate dataset choices ----
    observe({
      ds <- tpms_long$dataset |> unique() |> sort()
      updateSelectizeInput(session, "gp_dataset", choices = ds, server = TRUE)
    })

    # ---- populate gene choices (default FXN if present) ----
    observe({
      genes <- tpms_wide$gene_name |> unique() |> sort()
      sel <- if ("FXN" %in% genes) "FXN" else genes[1]
      updateSelectizeInput(session, "gp_gene", choices = genes, selected = sel, server = TRUE)
    })

    # ---- filter to selected dataset + gene ----
    dat_points <- reactive({
      req(input$gp_dataset, input$gp_gene)
      tpms_long |>
        dplyr::filter(dataset == input$gp_dataset, gene_name == input$gp_gene)
    })

    # ---- summary (mean ± SD per condition) ----
    dat_summary <- reactive({
      dat_points() |>
        dplyr::group_by(condition) |>
        dplyr::summarise(
          n   = dplyr::n(),
          mean = mean(TPM, na.rm = TRUE),
          sd   = sd(TPM, na.rm = TRUE),
          .groups = "drop"
        )
    })

    # ---- plot ----
    output$gp_plot <- renderPlot({
      pts <- dat_points(); sms <- dat_summary()
      validate(need(nrow(pts) > 0, "No TPM values for this selection."))

      p <- ggplot2::ggplot() +
        # mean ± SD as crossbar + errorbar
        ggplot2::geom_crossbar(
          data = sms,
          ggplot2::aes(x = condition, y = mean, ymin = mean, ymax = mean),
          width = 0.6, alpha = 0.3
        ) +
        ggplot2::geom_errorbar(
          data = sms,
          ggplot2::aes(x = condition, ymin = mean - sd, ymax = mean + sd),
          width = 0.2
        ) +
        # replicate points (jittered)
        ggplot2::geom_point(
          data = pts,
          ggplot2::aes(x = condition, y = TPM),
          position = ggplot2::position_jitter(width = 0.15, height = 0, seed = 1),
          size = 2, alpha = 0.8
        ) +
        ggplot2::labs(
          title = paste(input$gp_gene, "|", input$gp_dataset),
          x = "Condition", y = "TPM"
        ) +
        ggplot2::theme_bw(base_size = 12)

      if (isTRUE(input$gp_logy)) {
        p <- p + ggplot2::scale_y_continuous(trans = "log10")
      }
      p
    })

    # ---- show replicates table ----
    output$gp_table <- DT::renderDataTable({
      DT::datatable(
        dat_points() |>
          dplyr::arrange(condition, rep) |>
          dplyr::select(dataset, condition, rep, gene_id, gene_name, TPM),
        rownames = FALSE,
        options = list(pageLength = 10, scrollX = TRUE)
      )
    }, server = TRUE)

    # ---- downloads ----
    output$dl_points <- downloadHandler(
      filename = function()
        sprintf("TPM_%s_%s_replicates.csv", input$gp_dataset, input$gp_gene),
      content = function(file) readr::write_csv(dat_points(), file)
    )

    output$dl_summary <- downloadHandler(
      filename = function()
        sprintf("TPM_%s_%s_summary.csv", input$gp_dataset, input$gp_gene),
      content = function(file) readr::write_csv(dat_summary(), file)
    )

    output$dl_plot <- downloadHandler(
      filename = function()
        sprintf("TPM_%s_%s.png", input$gp_dataset, input$gp_gene),
      content = function(file) {
        gr <- grDevices::png(file, width = 1800, height = 1200, res = 200)
        print(isolate(output$gp_plot()))  # uses current plot code
        grDevices::dev.off()
      }
    )
  })
}
