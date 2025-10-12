# R/moduleGenePlots.R

#' Gene Plots – sidebar
#' @noRd
genePlotsSidebarUI <- function(id) {
  ns <- NS(id)
  tagList(
    h4("Gene plots"),
    selectizeInput(ns("gp_dataset"), "Dataset", choices = c("Erwin (Lymphoblastoid Cells)", "Indelicato (Skeletal Muscle)", "Lai (iPSCs)", "Lai (CNS neurons)", "Lai (PNS neurons)", "Lees (Cardiomyocytes)", "Maddock (Lower Motor Neurons)", "Maddock (Sensory Neurons)", "Maddock (Neural Crest Cells)", "Mishra (Neurons)", "Napierala (Fibroblasts)", "Vilema-Enriquez (Fibroblasts)" ), options = list(placeholder = "Select a dataset…")),
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

#' Gene Plots – main UI
#' @noRd
genePlotsMainUI <- function(id) {
  ns <- NS(id)
  tagList(
    # if you don't have shinycssloaders installed yet, use plotOutput directly
    if (requireNamespace("shinycssloaders", quietly = TRUE)) {
      shinycssloaders::withSpinner(plotOutput(ns("gp_plot"), height = 420), type = 4)
    } else {
      plotOutput(ns("gp_plot"), height = 420)
    },
    br(),
    DT::dataTableOutput(ns("gp_table"))
  )
}

#' Gene Plots server (RDS first, CSV fallback; seeded dataset list)
#' @noRd
genePlotsServer <- function(id, pkg = utils::packageName()) {
  moduleServer(id, function(input, output, session) {

    # ---- load TPM (prefer RDS; fallback to CSV) ----
    tpm_rds <- system.file("extdata/tpm/Combined_gene_tpm.rds", package = pkg, mustWork = FALSE)
    tpms_wide <- readRDS(tpm_rds)
    tpms_wide <- as.data.frame(tpms_wide)

    # ---- seed dataset dropdown immediately (so UI isn't empty) ----
    datasets_expected <- c(
      "Erwin", "Indelicato",
      "Lai_CNS", "Lai_iPSC", "Lai_PNS",
      "Lees",
      "Mishra",
      "Maddock_NCC", "Maddock_SN", "Maddock_LMN",
      "Napierala", "Vilema"
    )
    updateSelectizeInput(session, "gp_dataset",
                         choices  = datasets_expected,
                         selected = "Lees",
                         server   = TRUE
    )

    # ---- long format ----
    tpms_long0 <- tpms_wide |>
      tidyr::pivot_longer(
        cols = -c(gene_id, gene_name),
        names_to = "sample",
        values_to = "TPM"
      )

    # ---- parse dataset + condition with your groupings ----
    parse_sample <- function(s) {
      if (grepl("^Erwin_", s))      return(c("Erwin",      sub("^Erwin_([A-Za-z]+).*", "\\1", s)))
      if (grepl("^Indelicato_", s)) return(c("Indelicato", sub("^Indelicato_([A-Za-z]+).*", "\\1", s)))

      if (grepl("^Lai_", s)) {
        tissue <- sub("^Lai_(?:FRDA2|FRDA|IC)_([A-Za-z]+).*", "\\1", s)
        status <- sub("^Lai_((?:FRDA2|FRDA|IC)).*", "\\1", s)
        status <- ifelse(status == "FRDA2", "FRDA", status)
        return(c(paste0("Lai_", tissue), status))
      }

      if (grepl("^Lees_", s))       return(c("Lees",       sub("^Lees_(FA[0-9]+).*", "\\1", s)))
      if (grepl("^Mishra_", s))     return(c("Mishra",     sub("^Mishra_(?:FRDA|IC)_([A-Za-z0-9]+).*", "\\1", s)))
      if (grepl("^Napierala_", s))  return(c("Napierala",  sub("^Napierala_([A-Za-z0-9]+).*", "\\1", s)))
      if (grepl("^Vilema_", s))     return(c("Vilema",     sub("^Vilema_([A-Za-z0-9]+).*", "\\1", s)))

      if (grepl("^Maddock_FA[12].*NCC", s)) {
        return(c("Maddock_NCC", ifelse(grepl("icNCC", s), "icNCC", "NCC")))
      }
      if (grepl("^Maddock_SN_FA[12]", s)) {
        return(c("Maddock_SN", sub("^Maddock_SN_(FA[0-9]+).*", "\\1", s)))
      }
      if (grepl("^Maddock_FA2", s) && !grepl("NCC|SN", s)) {
        cond <- dplyr::case_when(
          grepl("icNIL", s) ~ "icNIL",
          grepl("NIL", s)   ~ "NIL",
          grepl("icN", s)   ~ "icN",
          grepl("N", s)     ~ "N",
          TRUE              ~ "Other"
        )
        return(c("Maddock_LMN", cond))
      }

      parts <- strsplit(s, "_", fixed = TRUE)[[1]]
      c(parts[1], ifelse(length(parts) > 1, parts[2], NA_character_))
    }

    parsed <- t(vapply(tpms_long0$sample, parse_sample, FUN.VALUE = c("", "")))

    tpms_long <- tpms_long0 |>
      dplyr::mutate(
        dataset   = parsed[,1],
        condition = parsed[,2],
        rep       = as.integer(sub(".*_REP([0-9]+).*", "\\1", sample))
      )

    # ---- replace dataset choices with the ones actually found ----
    observe({
      ds_found <- tpms_long$dataset |> unique() |> sort()
      if (!length(ds_found)) {
        shiny::showNotification("No datasets parsed from TPM; showing default list.", type = "warning", duration = 6)
        ds_found <- datasets_expected
      }
      updateSelectizeInput(session, "gp_dataset",
                           choices  = ds_found,
                           selected = if ("Lees" %in% ds_found) "Lees" else ds_found[1],
                           server   = TRUE
      )
    })

    # ---- genes dropdown ----
    observe({
      genes <- tryCatch(unique(tpms_wide$gene_name), error = function(e) character(0))
      genes <- sort(genes)
      if (!length(genes)) {
        genes <- c("FXN")
        shiny::showNotification("No gene_name column found; defaulting to FXN.", type = "error", duration = 6)
      }
      sel <- if ("FXN" %in% genes) "FXN" else genes[1]
      updateSelectizeInput(session, "gp_gene", choices = genes, selected = sel, server = TRUE)
    })

    # ---- filtered data ----
    dat_points <- reactive({
      req(input$gp_dataset, input$gp_gene)
      tpms_long |>
        dplyr::filter(dataset == input$gp_dataset, gene_name == input$gp_gene)
    })

    dat_summary <- reactive({
      dat_points() |>
        dplyr::group_by(condition) |>
        dplyr::summarise(
          n    = dplyr::n(),
          mean = mean(TPM, na.rm = TRUE),
          sd   = sd(TPM,   na.rm = TRUE),
          .groups = "drop"
        )
    })

    # ---- plot ----
    build_plot <- function(points, summary, logy = FALSE, title_txt = "") {
      p <- ggplot2::ggplot() +
        ggplot2::geom_crossbar(
          data = summary,
          ggplot2::aes(x = condition, y = mean, ymin = mean, ymax = mean),
          width = 0.55, alpha = 0.25
        ) +
        ggplot2::geom_errorbar(
          data = summary,
          ggplot2::aes(x = condition, ymin = mean - sd, ymax = mean + sd),
          width = 0.18
        ) +
        ggplot2::geom_point(
          data = points,
          ggplot2::aes(x = condition, y = TPM),
          position = ggplot2::position_jitter(width = 0.15, height = 0, seed = 1),
          size = 2, alpha = 0.85
        ) +
        ggplot2::labs(title = title_txt, x = "Condition", y = "TPM") +
        ggplot2::theme_bw(base_size = 12)
      if (isTRUE(logy)) p <- p + ggplot2::scale_y_continuous(trans = "log10")
      p
    }

    output$gp_plot <- renderPlot({
      pts <- dat_points(); sms <- dat_summary()
      validate(need(nrow(pts) > 0, "No TPM values for this selection."))
      build_plot(pts, sms, input$gp_logy, paste(input$gp_gene, "—", input$gp_dataset))
    })

    # ---- table ----
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
      filename = function() sprintf("TPM_%s_%s_replicates.csv", input$gp_dataset, input$gp_gene),
      content  = function(file) readr::write_csv(dat_points(), file)
    )
    output$dl_summary <- downloadHandler(
      filename = function() sprintf("TPM_%s_%s_summary.csv", input$gp_dataset, input$gp_gene),
      content  = function(file) readr::write_csv(dat_summary(), file)
    )
    output$dl_plot <- downloadHandler(
      filename = function() sprintf("TPM_%s_%s.png", input$gp_dataset, input$gp_gene),
      content  = function(file) {
        pts <- isolate(dat_points()); sms <- isolate(dat_summary())
        grDevices::png(file, width = 1800, height = 1200, res = 200)
        print(build_plot(pts, sms, isolate(input$gp_logy),
                         paste(isolate(input$gp_gene), "—", isolate(input$gp_dataset))))
        grDevices::dev.off()
      }
    )
  })
}
