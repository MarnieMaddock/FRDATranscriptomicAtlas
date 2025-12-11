biomarkerUI <- function(id) {
  ns <- NS(id)
  tagList(
    h4("Biomarker Discovery"),

    checkboxGroupInput(
      ns("datasets"),
      label = "Datasets (select any number)",
      choices = character(0)   # populated dynamically
    ),

    hr(),
    sliderInput(
      ns("min_study_count"),
      "Dataset significance threshold:",
      min = 1, max = 12,
      value = 8, step = 1
    ),
    helpText("A gene is included only if it is significant in more than this number of datasets."),
    selectInput(
      ns("alpha"),
      "Adjusted p-value threshold:",
      choices = c("0.10", "0.05", "0.01", "0.001"),
      selected = "0.05"
    ),
    br(),
    downloadButton(ns("download_data"), "Download Ranked Genes"),
    downloadButton(ns("download_plot"), "Download Figure")
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
                            baseline_long_path = system.file(
                              "extdata/biomarker/baseline_long.rds",
                              package = "FRDATranscriptomicAtlas"
                            ),
                            pkg = utils::packageName()) {

  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # ============================================================
    # Load data
    # ============================================================
    baseline_long <- readRDS(baseline_long_path)
    baseline_long <- baseline_long %>%
      mutate(study = gsub("_batchcorrection$", "", study))


    alpha <- reactive({
      as.numeric(input$alpha %||% 0.05)
    })


    # ============================================================
    # Load pretty_map safely
    # ============================================================
    `%||%` <- function(x, y) if (is.null(x)) y else x

    pretty_map <- tryCatch(
      get("pretty_map", envir = asNamespace(pkg)),
      error = function(e) {
        warning("pretty_map not found; using empty vector.")
        character(0)
      }
    )

    # vectorised label function
    pretty_label <- function(ids) {
      out <- pretty_map[ids]
      out[is.na(out)] <- ids[is.na(out)]
      return(out)
    }

    # ============================================================
    # Populate dataset choices (checkboxGroupInput)
    # ============================================================
    observe({
      avail <- sort(unique(baseline_long$study))
      labs  <- pretty_label(avail)

      updateCheckboxGroupInput(
        session, "datasets",
        choices  = stats::setNames(avail, labs),
        selected = avail   # default: all datasets included
      )
    })

    # ============================================================
    # Filter baseline_long by selected datasets
    # ============================================================
    filtered_baseline <- reactive({
      req(input$datasets)
      baseline_long %>% filter(study %in% input$datasets)
    })


    # ============================================================
    # 1. Gene ranking
    # ============================================================
    rank_df <- reactive({
      df <- filtered_baseline()
      min_n <- input$min_study_count %||% 8

      top_up <- df %>%
        filter(!is.na(padj), padj <= alpha(), log2FoldChange > 0) %>%
        distinct(gene_id, gene_name, study) %>%
        count(gene_id, gene_name, name = "n_studies") %>%
        filter(n_studies > min_n)

      top_down <- df %>%
        filter(!is.na(padj), padj <= alpha(), log2FoldChange < 0) %>%
        distinct(gene_id, gene_name, study) %>%
        count(gene_id, gene_name, name = "n_studies") %>%
        filter(n_studies > min_n)


      bind_rows(top_up, top_down) %>%
        distinct(gene_id, .keep_all = TRUE) %>%
        mutate(direction = if_else(gene_id %in% top_up$gene_id, "Up", "Down")) %>%
        relocate(gene_id, gene_name)
    })


    # ============================================================
    # 2. Build heat_df
    # ============================================================
    heat_df <- reactive({

      df  <- filtered_baseline()
      rnk <- rank_df()

      df <- df %>%
        semi_join(rnk, by = "gene_id") %>%
        mutate(
          direction = case_when(
            is.na(padj) ~ "No p-value (filtered)",
            padj < alpha() & log2FoldChange > 0 ~ "p < 0.05, log2FC > 0 (Up)",
            padj < alpha() & log2FoldChange < 0 ~ "p < 0.05, log2FC < 0 (Down)",
            TRUE ~ "p ≥ 0.05 (Not significant)"
          ),
          gene_lab = factor(gene_name, levels = rnk$gene_name),
          outline_class = case_when(
            is.na(log2FoldChange) ~ "NA",
            log2FoldChange > 0    ~ "Positive",
            log2FoldChange < 0    ~ "Negative",
            log2FoldChange == 0   ~ "Zero"
          ),
          study        = factor(study, levels = input$datasets),
          study_pretty = factor(pretty_label(study),
                                levels = pretty_label(input$datasets)),
          abs_lfc = ifelse(is.finite(abs(log2FoldChange)),
                           abs(log2FoldChange), 0)
        )

      lfc_cap <- quantile(df$abs_lfc, 0.99, na.rm = TRUE)
      df$size_scaled <- pmin(df$abs_lfc, lfc_cap)

      df
    })


    # ============================================================
    # 3. Reorder genes by (down - up) score
    # ============================================================
    heat_df_reordered <- reactive({

      df <- heat_df()

      gene_levels <- df %>%
        filter(direction %in% c(
          "p < 0.05, log2FC > 0 (Up)",
          "p < 0.05, log2FC < 0 (Down)"
        )) %>%
        count(gene_lab, direction, name = "n") %>%
        tidyr::complete(
          gene_lab,
          direction = c(
            "p < 0.05, log2FC < 0 (Down)",
            "p < 0.05, log2FC > 0 (Up)"
          ),
          fill = list(n = 0)
        ) %>%
        tidyr::pivot_wider(names_from = direction, values_from = n) %>%
        mutate(score = `p < 0.05, log2FC < 0 (Down)` -
                 `p < 0.05, log2FC > 0 (Up)`) %>%
        arrange(desc(score), gene_lab) %>%
        pull(gene_lab)

      df %>% mutate(gene_lab = factor(gene_lab, levels = gene_levels))
    })



    # ============================================================
    # 4. Combined plot (heatmap + barplot)
    # ============================================================
    plot_height <- reactive({
      df <- heat_df_reordered()
      n_genes <- length(unique(df$gene_lab))

      # pixels per gene row
      per_gene <- 20

      # compute base height
      h <- n_genes * per_gene

      # enforce min/max bounds
      h <- max(400, min(h, 3000))

      return(h)
    })

    plot_width <- reactive({
      n_ds <- length(input$datasets)
      w <- max(800, n_ds * 70)
      return(w)
    })

    output$combined_plot <- renderPlot({

      df <- heat_df_reordered()

      # ------------------------------------------------------------
      # Colours
      # ------------------------------------------------------------
      dir_cols <- c(
        "p < 0.05, log2FC < 0 (Down)"     = "#00B7C7",
        "p < 0.05, log2FC > 0 (Up)"       = "#DC267F",
        "p ≥ 0.05 (Not significant)"      = "gray60",
        "No p-value (filtered)"           = "gray80"
      )

      outline_cols <- c(
        "Positive" = "black",
        "Negative" = "transparent",
        "Zero"     = "grey50",
        "NA"       = "grey80"
      )

      # ------------------------------------------------------------
      # Heatmap
      # ------------------------------------------------------------
      p_heat <- ggplot(df, aes(x = study_pretty, y = gene_lab)) +
        geom_point(
          aes(fill = direction, size = size_scaled, colour = outline_class),
          shape = 21,
          stroke = 0.7
        ) +
        scale_fill_manual(
          values = dir_cols,
          name   = "DE Category",
          guide  = guide_legend(
            override.aes = list(shape = 21, size = 6, stroke = 0.3)
          )
        ) +
        scale_colour_manual(
          values = outline_cols,
          name   = "Direction (log2FC)",
          guide  = guide_legend(
            override.aes = list(fill = "white", shape = 21, size = 6, stroke = 1)
          )
        ) +
        scale_size_continuous(
          name = "abs(log2FC)",
          range = c(2, 10),
          guide = guide_legend(
            override.aes = list(shape = 21, fill = "grey60")
          )
        ) +
        labs(x = NULL, y = NULL) +
        theme_minimal(base_size = 18) +
        theme(
          panel.grid = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
          axis.text.y = element_text(face = "italic", color = "black"),
          legend.title = element_text(color = "black"),
          legend.text  = element_text(color = "black")
        )


      # ------------------------------------------------------------
      # Bar counts for each gene
      # ------------------------------------------------------------
      direction_levels <- c(
        "p < 0.05, log2FC > 0 (Up)",
        "p ≥ 0.05 (Not significant)",
        "No p-value (filtered)",
        "p < 0.05, log2FC < 0 (Down)"
      )

      counts_long <- df %>%
        count(gene_lab, direction, name = "n") %>%
        tidyr::complete(
          gene_lab,
          direction = direction_levels,
          fill = list(n = 0)
        ) %>%
        mutate(
          direction = factor(direction, levels = direction_levels),
          gene_lab  = factor(gene_lab, levels = levels(df$gene_lab))
        ) %>%
        arrange(gene_lab, direction)

      n_datasets <- length(input$datasets)

      p_bar <- ggplot(
        counts_long,
        aes(y = gene_lab, x = n, fill = direction)
      ) +
        geom_col(width = 0.85) +
        scale_fill_manual(values = dir_cols) +
        guides(fill = "none") +
        scale_x_continuous(
          limits = c(0, n_datasets),
          expand = expansion(mult = c(0, 0.02))
        ) +
        labs(x = "Number of datasets", y = NULL) +
        theme_minimal(base_size = 18) +
        theme(
          panel.grid.major.y = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_text(color = "black"),
          axis.title.x = element_text(color = "black")
        ) +
        geom_text(
          data = dplyr::filter(counts_long, n > 0),
          aes(label = n),
          position = position_stack(vjust = 0.5),
          size = 4,
          color = "black"
        )


      # ------------------------------------------------------------
      # Combined figure
      # ------------------------------------------------------------
      p_heat + p_bar + patchwork::plot_layout(
        widths = c(2, 0.5),
        guides = "collect"
      )
    }, height = plot_height, width = plot_width)




    # ============================================================
    # 6. Downloads
    # ============================================================
    output$download_data <- downloadHandler(
      filename = function() "biomarker_gene_rank.csv",
      content = function(file) {
        write.csv(rank_df(), file, row.names = FALSE)
      }
    )

    output$download_plot <- downloadHandler(
      filename = function() "biomarker_plot.svg",
      content = function(file) {
        svg(file, width = 14, height = 14)
        print(output$combined_plot())
        dev.off()
      }
    )

  })
}

