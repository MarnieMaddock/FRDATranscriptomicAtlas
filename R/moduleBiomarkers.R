#' @import ggplot2
NULL

biomarkerUI <- function(id) {
  ns <- NS(id)
  tagList(
    h4("Biomarker Discovery"),
    div(
      style = "display:flex; gap:8px; margin-bottom:8px;",
      actionButton(ns("datasets_all"),  "Select all", class = "btn btn-sm btn-default"),
      actionButton(ns("datasets_none"), "Clear",      class = "btn btn-sm btn-default")
    ),

    checkboxGroupInput(
      ns("datasets"),
      label = "Datasets",
      choices = character(0)
    ),
    br(),
    h4("Gene selection mode"),
    radioButtons(
      ns("gene_mode"),
      label = "",
      choices = c(
        "Discover top transcriptomic biomarkers" = "discover",
        "Plot specified genes" = "query"
      ),
      selected = "discover",
      inline = TRUE
    ),

    conditionalPanel(
      condition = sprintf("input['%s'] === 'query'", ns("gene_mode")),
      textAreaInput(
        ns("feature_query"),
        label = "Specify genes",
        placeholder = "FXN, PIP5K1B, RPS29 ...",
        rows = 4
      )
    ),

    hr(),
    sliderInput(
      ns("min_study_count"),
      "Dataset significance threshold:",
      min = 6, max = 12,
      value = 8, step = 1
    ),
    helpText("A gene is included only if it is significant in more than this number of datasets."),
    br(),
    radioButtons(
      ns("alpha"),
      "Adjusted p-value threshold:",
      choices = c("0.10", "0.05", "0.01", "0.001"),
      selected = "0.05"
    ),
    br(),
    helpText("Comma-separated gene symbols or Ensembl IDs. Leave blank to include all ranked genes."),
    br(),
    downloadButton(ns("download_data"), "Download Ranked Genes"),
    br(),
    h4("Export Figure"),
    fluidRow(
      column(
        width = 4,
        numericInput(
          ns("export_width"),
          "Width (cm)",
          value = 20,
          min = 5,
          max = 100,
          step = 1
        )
      ),
      column(
        width = 4,
        numericInput(
          ns("export_height"),
          "Height (cm)",
          value = 20,
          min = 5,
          max = 100,
          step = 1
        )
      ),
      column(
        width = 4,
        numericInput(
          ns("export_dpi"),
          "PNG DPI",
          value = 300,
          min = 72,
          max = 600,
          step = 50
        )
      )
    ),

    br(),

    fluidRow(
      column(
        width = 6,
        downloadButton(ns("download_svg"), "Download SVG")
      ),
      column(
        width = 6,
        downloadButton(ns("download_png"), "Download PNG")
      )
    )
  )
}




biomarkerMainUI <- function(id, pkg = utils::packageName()) {
  ns <- NS(id)
  tagList(
    h3("Common DEGs Across FRDA Studies"),
    shinycssloaders::withSpinner(
      plotOutput(ns("combined_plot")),
      type = 4, color = "#005249"
    ),
    br()
  )
}

biomarkerServer <- function(id,
                            pkg = utils::packageName(),
                            data_dir = NULL,
                            data_mode = "local") {

  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    `%||%` <- function(x, y) if (is.null(x)) y else x

    baseline_long <- if (identical(data_mode, "cloud")) {

      get_biomarker_baseline_cloud_cached()

    } else {

      ensure_atlas_data(
        keys = "biomarker",
        package = pkg,
        data_mode = data_mode
      )

      biomarker_dir <- function(package) {
        file.path(
          tools::R_user_dir(package, "cache"),
          "biomarker"
        )
      }

      baseline_path <- if (!is.null(data_dir)) {
        file.path(data_dir, "baseline_long.rds")
      } else {
        file.path(biomarker_dir(pkg), "baseline_long.rds")
      }

      validate(
        need(
          file.exists(baseline_path),
          "Biomarker baseline data not found. Please ensure data download completed."
        )
      )

      readRDS(baseline_path)
    }


    alpha <- reactive(as.numeric(input$alpha %||% 0.05))

    # ============================================================
    # pretty_map
    # ============================================================
    pretty_map <- tryCatch(
      get("pretty_map", envir = asNamespace(pkg)),
      error = function(e) character(0)
    )

    pretty_label <- function(ids) {
      out <- pretty_map[ids]
      out[is.na(out)] <- ids[is.na(out)]
      out
    }

    # ============================================================
    # Dataset selector
    # ============================================================
    observeEvent(input$datasets_all, {
      avail <- sort(unique(baseline_long$study))
      updateCheckboxGroupInput(session, "datasets", selected = avail)
    })

    observeEvent(input$datasets_none, {
      updateCheckboxGroupInput(session, "datasets", selected = character(0))
    })

    observe({
      avail <- sort(unique(baseline_long$study))
      updateCheckboxGroupInput(
        session, "datasets",
        choices  = stats::setNames(avail, pretty_label(avail)),
        selected = avail
      )
    })

    filtered_baseline <- reactive({
      req(input$datasets)
      baseline_long %>% dplyr::filter(study %in% input$datasets)
    })

    threshold_valid <- reactive({
      req(input$datasets)

      n_ds <- length(input$datasets)

      validate(
        need(
          input$min_study_count <= n_ds,
          paste0(
            "Dataset significance threshold (", input$min_study_count,
            ") exceeds the number of selected datasets (", n_ds, ")."
          )
        )
      )

      TRUE
    })

    # ============================================================
    # Parse comma-separated gene list (QUERY MODE)
    # ============================================================
    feature_query <- reactive({
      if (input$gene_mode != "query") return(NULL)

      q <- input$feature_query
      if (is.null(q) || q == "") return(character(0))

      feats <- unlist(strsplit(q, "[,\\n]+"))
      feats <- trimws(feats)
      feats[feats != ""]
    })

    # ============================================================
    # DISCOVERY MODE: ranked genes
    # ============================================================
    ranked_genes <- reactive({
      req(input$gene_mode == "discover")
      req(threshold_valid())
      df <- filtered_baseline()
      min_n <- input$min_study_count %||% 8


      top_up <- df %>%
        dplyr::filter(!is.na(padj), padj <= alpha(), log2FoldChange > 0) %>%
        dplyr::distinct(gene_id, gene_name, study) %>%
        dplyr::count(gene_id, gene_name, name = "n_studies") %>%
        dplyr::filter(n_studies > min_n)

      top_down <- df %>%
        dplyr::filter(!is.na(padj), padj <= alpha(), log2FoldChange < 0) %>%
        dplyr::distinct(gene_id, gene_name, study) %>%
        dplyr::count(gene_id, gene_name, name = "n_studies") %>%
        dplyr::filter(n_studies > min_n)

      dplyr::bind_rows(top_up, top_down) %>%
        dplyr::distinct(gene_id, gene_name)
    })

    # ============================================================
    # QUERY MODE: exact user genes
    # ============================================================
    queried_genes <- reactive({
      req(input$gene_mode == "query")

      feats <- feature_query()
      if (length(feats) == 0) return(NULL)

      filtered_baseline() %>%
        dplyr::filter(
          toupper(gene_name) %in% toupper(feats) |
            toupper(gene_id)   %in% toupper(feats)
        ) %>%
        dplyr::distinct(gene_id, gene_name)
    })

    # ============================================================
    # Unified gene set
    # ============================================================
    selected_genes <- reactive({
      if (input$gene_mode == "discover") ranked_genes()
      else queried_genes()
    })

    # ============================================================
    # Build heat_df
    # ============================================================
    heat_df <- reactive({
      req(threshold_valid())
      genes <- selected_genes()
      req(genes)

      df <- filtered_baseline() %>%
        dplyr::semi_join(genes, by = "gene_id") %>%
        dplyr::mutate(
          direction = dplyr::case_when(
            is.na(padj) ~ "Filtered (no adjusted p-value)",
            padj < alpha() & log2FoldChange > 0 ~ "Upregulated in FRDA (FDR < 0.05, log2FC > 0)",
            padj < alpha() & log2FoldChange < 0 ~ "Downregulated in FRDA (FDR < 0.05, log2FC < 0)",
            TRUE ~ "Not significant (FDR >= 0.05)"
          ),
          gene_lab = factor(gene_name),
          outline_class = dplyr::case_when(
            log2FoldChange > 0    ~ "Higher in FRDA",
            log2FoldChange < 0    ~ "Lower in FRDA",
            TRUE   ~ "No change"
          ),
          study        = factor(study, levels = input$datasets),
          study_pretty = factor(pretty_label(study),
                                levels = pretty_label(input$datasets)),
          abs_lfc = ifelse(is.finite(abs(log2FoldChange)),
                           abs(log2FoldChange), 0)
        )

      cap <- stats::quantile(df$abs_lfc, 0.99, na.rm = TRUE)
      df$size_scaled <- pmin(df$abs_lfc, cap)
      df
    })

    # ============================================================
    # Gene ordering
    # ============================================================
    heat_df_reordered <- reactive({
      req(threshold_valid())
      df <- heat_df()

      if (input$gene_mode == "query") {

        ord <- feature_query()

      } else {

        ord <- df %>%
          dplyr::summarise(
            up = sum(
              direction == "Upregulated in FRDA (FDR < 0.05, log2FC > 0)",
              na.rm = TRUE
            ),
            down = sum(
              direction == "Downregulated in FRDA (FDR < 0.05, log2FC < 0)",
              na.rm = TRUE
            ),
            .by = gene_lab
          ) %>%
          dplyr::mutate(score = down - up) %>%
          dplyr::arrange(dplyr::desc(score), gene_lab) %>%
          dplyr::pull(gene_lab)
      }

      df %>% dplyr::mutate(gene_lab = factor(gene_lab, levels = ord))
    })



    # ============================================================
    # Plot sizing
    # ============================================================
    plot_height <- reactive({
      req(threshold_valid())
      df <- heat_df_reordered()
      req(nrow(df) > 0)

      n_genes <- length(unique(df$gene_lab))

      per_gene <- if (n_genes <= 10) 26 else
          if (n_genes <= 80) 25 else 23

      top_padding    <- 20
      bottom_padding <- 80   # <- important for 45° labels

      h <- (n_genes * per_gene) + top_padding + bottom_padding

      max(700, min(h, 8000))
    })

    plot_width <- reactive({
      n_ds <- length(input$datasets)

      # base width for readability
      base_width <- 800

      # add extra space only when datasets exceed threshold
      extra_per_dataset <- 50

      extra <- if (n_ds > 8) {
        (n_ds - 8) * extra_per_dataset
      } else {
        0
      }

      base_width + extra
    })

    # ============================================================
    # Plot
    # ============================================================
    combined_plot_obj <- reactive({
      req(threshold_valid())
      df <- heat_df_reordered()
      req(nrow(df) > 0)

      # ------------------------------------------------------------
      # Colours
      # ------------------------------------------------------------
      dir_cols <- c(
        "Downregulated in FRDA (FDR < 0.05, log2FC < 0)"    = "#00B7C7",
        "Upregulated in FRDA (FDR < 0.05, log2FC > 0)"       = "#DC267F",
        "Not significant (FDR >= 0.05)"     = "gray60",
        "Filtered (no adjusted p-value)"         = "gray80"
      )
      outline_cols <- c(
        "Higher in FRDA" = "black",
        "Lower in FRDA" = "transparent",
        "No change"     = "grey50"
      )

      # ------------------------------------------------------------
      # Heatmap
      # ------------------------------------------------------------
      p_heat <- ggplot(df, aes(x = study, y = gene_lab)) +
        ggplot2::geom_point(
          ggplot2::aes(fill = direction, size = size_scaled, colour = outline_class),
          shape = 21,
          stroke = 0.7
        ) +
        ggplot2::scale_fill_manual(
          values = dir_cols,
          name   = "Differential expression",
          guide  = ggplot2::guide_legend(
            override.aes = list(shape = 21, size = 6, stroke = 0.3)
          )
        ) +
        ggplot2::scale_colour_manual(
          values = outline_cols,
          name   = expression("Direction of log"[2]*"FC change"),
          labels = c(
            "Up (FDR < 0.05)"              = expression("Up (FDR " < " 0.05)"),
            "Down (FDR < 0.05)"            = expression("Down (FDR " < " 0.05)"),
            "Not significant (FDR >= 0.05)" = expression("Not significant (FDR " >= " 0.05)")
          ),
          guide  = ggplot2::guide_legend(
            override.aes = list(fill = "white", shape = 21, size = 6, stroke = 1)
          )
        ) +
        ggplot2::scale_size_continuous(
          name = expression("Effect size ("* "|" *"log"[2]*"FC|)"),
          range = c(2, 10),
          guide = ggplot2::guide_legend(
            override.aes = list(shape = 21, fill = "grey60")
          )
        ) +
        ggplot2::labs(x = NULL, y = NULL) +
        ggplot2::theme_minimal(base_size = 18) +
        ggplot2::theme(
          panel.grid = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
          axis.text.y = element_text(face = "italic", color = "black"),
          legend.title = element_text(color = "black"),
          legend.text  = element_text(color = "black"),
          plot.margin = margin(t = 20, r = 30, b = 120, l = 120)
        ) +
        ggplot2:: scale_y_discrete(expand = expansion(add = c(0.7, 1))) +
        ggplot2::scale_x_discrete(
          labels = pretty_label
        )

      # ------------------------------------------------------------
      # Bar counts for each gene
      # ------------------------------------------------------------
      direction_levels <- c(
        "Upregulated in FRDA (FDR < 0.05, log2FC > 0)",
        "Not significant (FDR >= 0.05)",
        "Filtered (no adjusted p-value)",
        "Downregulated in FRDA (FDR < 0.05, log2FC < 0)"
      )

      counts_long <- df %>%
        count(gene_lab, direction, name = "n") %>%
        tidyr::complete(
          gene_lab,
          direction = direction_levels,
          fill = list(n = 0)
        ) %>%
        dplyr::mutate(
          direction = factor(direction, levels = direction_levels),
          gene_lab  = factor(gene_lab, levels = levels(df$gene_lab))
        ) %>%
        dplyr::arrange(gene_lab, direction)

      n_datasets <- length(input$datasets)

      p_bar <- ggplot(
        counts_long,
        ggplot2::aes(y = gene_lab, x = n, fill = direction)
      ) +
        ggplot2::geom_col(width = 0.85) +
        ggplot2::scale_fill_manual(values = dir_cols) +
        ggplot2::guides(fill = "none") +
        ggplot2::scale_x_continuous(
          limits = c(0, n_datasets),
          expand = expansion(mult = c(0, 0.02))
        ) +
        ggplot2::labs(x = "Number of datasets", y = NULL) +
        ggplot2::theme_minimal(base_size = 18) +
        ggplot2::theme(
          panel.grid.major.y = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_text(color = "black"),
          axis.title.x = element_text(color = "black"),
          plot.margin = margin(t = 20, r = 30, b = 120, l = 10)
        ) +
        ggplot2::scale_y_discrete(expand = expansion(add = c(0.7, 1))) +
        ggplot2::geom_text(
          data = dplyr::filter(counts_long, n > 0),
          aes(label = n),
          position = position_stack(vjust = 0.5),
          size = 4,
          color = "black"
        )

      # ------------------------------------------------------------
      # Combined figure
      # ------------------------------------------------------------
      patchwork::wrap_plots(
        p_heat,
        p_bar,
        widths = c(2, 0.5),
        guides = "collect"
      )

    })


    output$combined_plot <- renderPlot({

      req(input$datasets)

      n_ds <- length(input$datasets)

      validate(
        need(
          input$min_study_count <= n_ds,
          paste0(
            "Dataset significance threshold (", input$min_study_count,
            ") exceeds the number of selected datasets (", n_ds, ").\n\n",
            "Lower the slider or select more datasets."
          )
        )
      )

      combined_plot_obj()

    }, height = plot_height, width = plot_width)


    # ============================================================
    # Downloads
    # ============================================================
    output$download_data <- downloadHandler(
      filename = function() "biomarker_genes.csv",
      content = function(file) {
        utils::write.csv(selected_genes(), file, row.names = FALSE)
      }
    )

    output$download_svg <- downloadHandler(
      filename = function() {
        "biomarker_plot.svg"
      },
      content = function(file) {

        w_cm <- input$export_width
        h_cm <- input$export_height

        w_in <- w_cm / 2.54
        h_in <- h_cm / 2.54

        svglite::svglite(
          file = file,
          width = w_in,
          height = h_in
        )

        print(combined_plot_obj())
        grDevices::dev.off()
      }
    )

    output$download_png <- downloadHandler(
      filename = function() {
        "biomarker_plot.png"
      },
      content = function(file) {

        w_cm <- input$export_width
        h_cm <- input$export_height
        dpi  <- input$export_dpi

        w_in <- w_cm / 2.54
        h_in <- h_cm / 2.54

        ggplot2::ggsave(
          filename = file,
          plot     = combined_plot_obj(),
          width    = w_in,
          height   = h_in,
          dpi      = dpi,
          units    = "in"
        )
      }
    )

    px_to_cm <- function(px, dpi = 96) (px / dpi) * 2.54
    cm_to_px <- function(cm, dpi = 96) (cm / 2.54) * dpi

    observeEvent(
      list(plot_width(), plot_height(), input$export_dpi),
      {
        dpi <- input$export_dpi %||% 300

        # Convert the current on-screen plot px -> cm for the UI defaults
        updateNumericInput(session, "export_width",
                           value = round(plot_width() / 37, 1)
        )
        updateNumericInput(session, "export_height",
                           value = round(plot_height() / 37, 1)
        )
      },
      ignoreInit = FALSE
    )

    observeEvent(input$gene_mode == "discover", {
      updateTextAreaInput(session, "feature_query", value = "")
    })



  })
}
