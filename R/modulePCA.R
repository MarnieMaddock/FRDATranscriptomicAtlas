# ---- PCA module ------------------------------------------------------------
# Expects files saved by your pipeline: <label>_pca_input.rds
# Each file is a list(vsd_mat, meta, percent_var)

pcaUI <- function(id, title = "PCA (VST)") {
  ns <- NS(id)
  tagList(
    h4(title),
    fluidRow(
      column(
        width = 3,
        helpText("Select a saved comparison and customise aesthetics."),
        uiOutput(ns("pick_file_ui")),
        hr(),
        uiOutput(ns("color_var_ui")),
        uiOutput(ns("shape_var_ui")),
        checkboxInput(ns("label_points"), "Show sample labels", FALSE),
        checkboxInput(ns("draw_ellipses"), "Group ellipses (95%)", FALSE),
        sliderInput(ns("pt_size"), "Point size", min = 1, max = 6, value = 3, step = 0.5),
        radioButtons(ns("engine"), "Plot engine", inline = TRUE,
                     choices = c("ggplot" = "ggplot", "plotly (lasso)" = "plotly"),
                     selected = "ggplot"),
        actionButton(ns("update"), "Update plot", class = "btn-primary"),
        hr(),
        downloadButton(ns("download_png"), "Download PNG"),
        downloadButton(ns("download_svg"), "Download SVG"),
        downloadButton(ns("download_scores_csv"), "Download PC scores CSV")
      ),
      column(
        width = 9,
        uiOutput(ns("var_explained_text")),
        conditionalPanel(
          sprintf("input['%s'] === 'plotly'", ns("engine")),
          plotly::plotlyOutput(ns("pca_plotly"), height = "560px")
        ),
        conditionalPanel(
          sprintf("input['%s'] === 'ggplot'", ns("engine")),
          plotOutput(ns("pca_plot"), height = "560px")
        ),
        br(),
        h5("Selected samples"),
        DT::DTOutput(ns("selected_table"))
      )
    )
  )
}

pcaServer <- function(id,
                      data_dir = system.file("extdata/deseq_objects", package = "FRDATranscriptomicAtlas"),
                      fallback_dir = "inst/extdata/deseq_objects") {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # ---- locate directory with *_pca_input.rds
    dir_use <- reactive({
      if (dir.exists(data_dir)) return(data_dir)
      if (dir.exists(fallback_dir)) return(fallback_dir)
      stop("No PCA directory found. Checked:\n  - ", data_dir, "\n  - ", fallback_dir)
    })

    # ---- list available files
    pca_files <- reactive({
      paths <- list.files(dir_use(), pattern = "_pca_input\\.rds$", full.names = TRUE)
      tibble::tibble(
        label = sub("_pca_input\\.rds$", "", basename(paths)),
        path  = paths
      ) %>% dplyr::arrange(label)
    })

    output$pick_file_ui <- renderUI({
      req(nrow(pca_files()) > 0)
      selectInput(ns("picked"), "Comparison",
                  choices = setNames(pca_files()$path, pca_files()$label))
    })

    # ---- load chosen pca_input
    pca_input <- reactive({
      req(input$picked)
      readRDS(input$picked)
    })

    # ---- meta columns for aesthetics
    meta_cols <- reactive({
      req(pca_input()$meta)
      # offer factors/characters with <= 30 unique values (avoid IDs)
      m <- pca_input()$meta
      keep <- vapply(m, function(x) {
        is.factor(x) || is.character(x) || is.logical(x)
      }, logical(1))
      m <- m[, keep, drop = FALSE]
      small <- vapply(m, function(x) length(unique(x)) <= 30, logical(1))
      names(m[, small, drop = FALSE])
    })

    output$color_var_ui <- renderUI({
      req(meta_cols())
      # prefer these if present
      preferred <- c("group", "FRDA_CTRL", "cell_type", "Study_Alias")
      default <- intersect(preferred, meta_cols())
      selectInput(ns("color_var"), "Colour by",
                  choices = meta_cols(),
                  selected = if (length(default)) default[1] else meta_cols()[1])
    })

    output$shape_var_ui <- renderUI({
      req(meta_cols())
      choices <- c("None", meta_cols())
      selectInput(ns("shape_var"), "Shape by", choices = choices, selected = "None")
    })

    # ---- compute PCA
    pr_obj <- reactive({
      req(pca_input()$vsd_mat)
      X <- pca_input()$vsd_mat
      # prcomp expects samples as rows
      pr <- prcomp(t(X), center = TRUE, scale. = FALSE)
      pr
    })

    percent_var <- reactive({
      pr <- pr_obj()
      round(100 * (pr$sdev^2) / sum(pr$sdev^2), 2)
    })

    scores_df <- reactive({
      req(pr_obj(), pca_input()$meta)
      df <- as.data.frame(pr_obj()$x)
      df$sample_id <- rownames(df)
      meta <- pca_input()$meta
      meta$sample_id <- rownames(meta)
      dplyr::left_join(df, meta, by = "sample_id")
    })

    output$var_explained_text <- renderUI({
      pv <- percent_var()
      HTML(sprintf("<p><em>PC1: %0.2f%% variance • PC2: %0.2f%%</em></p>", pv[1], pv[2]))
    })

    # ---- ggplot renderer
    plot_obj <- reactive({
      req(scores_df(), input$color_var, input$pt_size)
      df <- scores_df()
      aes_color <- rlang::sym(input$color_var)
      p <- ggplot2::ggplot(df, ggplot2::aes(PC1, PC2, colour = !!aes_color)) +
        ggplot2::geom_point(size = input$pt_size, alpha = 0.9)

      # optional shape
      if (!is.null(input$shape_var) && input$shape_var != "None") {
        aes_shape <- rlang::sym(input$shape_var)
        p <- ggplot2::ggplot(df, ggplot2::aes(PC1, PC2, colour = !!aes_color, shape = !!aes_shape)) +
          ggplot2::geom_point(size = input$pt_size, alpha = 0.9)
      }

      # labels
      if (isTRUE(input$label_points)) {
        p <- p + ggrepel::geom_text_repel(ggplot2::aes(label = sample_id), size = 3, max.overlaps = 30)
      }

      # ellipses
      if (isTRUE(input$draw_ellipses)) {
        p <- p + ggplot2::stat_ellipse(ggplot2::aes(group = !!aes_color), level = 0.95, linetype = 2)
      }

      # axis titles with % variance
      pv <- percent_var()
      p <- p +
        ggplot2::labs(x = sprintf("PC1 (%.2f%%)", pv[1]),
                      y = sprintf("PC2 (%.2f%%)", pv[2]),
                      color = input$color_var,
                      shape = if (!is.null(input$shape_var) && input$shape_var != "None") input$shape_var else NULL) +
        ggplot2::theme_classic(base_size = 14)

      # Apply user's theme if available
      if (exists("theme_Marnie", mode = "function") || exists("theme_Marnie")) {
        p <- p + theme_Marnie
      }
      p
    })

    # re-render on Update button
    plot_trig <- reactive(input$update)

    output$pca_plot <- renderPlot({
      plot_trig()  # depend
      plot_obj()
    })

    # ---- plotly with lasso selection
    selected_ids <- reactiveVal(character(0))

    output$pca_plotly <- plotly::renderPlotly({
      req(plot_obj())
      plot_trig()
      gg <- plot_obj()
      plt <- plotly::ggplotly(gg, tooltip = c("sample_id", input$color_var))
      plt <- plotly::layout(plt, dragmode = "lasso")
      plt <- plotly::event_register(plt, "plotly_selected")
      plt
    })

    # capture lasso selections
    observeEvent(plotly::event_data("plotly_selected", source = NULL), {
      ed <- plotly::event_data("plotly_selected", source = NULL)
      if (is.null(ed) || !nrow(ed)) {
        selected_ids(character(0))
      } else {
        # ggplotly maps key if present; fallback to near points using PC1/PC2
        df <- scores_df()
        # approximate match
        hits <- dplyr::semi_join(df,
                                 dplyr::distinct(ed, x = round(x, 6), y = round(y, 6)),
                                 by = c("PC1" = "x", "PC2" = "y"))
        selected_ids(hits$sample_id)
      }
    }, ignoreInit = TRUE)

    # show selected samples (works for either engine; empty for ggplot)
    output$selected_table <- DT::renderDT({
      df <- scores_df()
      sel <- selected_ids()
      out <- if (length(sel)) dplyr::filter(df, sample_id %in% sel) else df[0, ]
      DT::datatable(out[, c("sample_id", "PC1", "PC2", input$color_var), drop = FALSE],
                    rownames = FALSE, options = list(pageLength = 8))
    })

    # ---- downloads
    output$download_png <- downloadHandler(
      filename = function() {
        paste0(basename(sub("_pca_input\\.rds$", "", input$picked)), "_PCA.png")
      },
      content = function(file) {
        ggplot2::ggsave(file, plot = plot_obj(), width = 7, height = 5.5, dpi = 300)
      }
    )

    output$download_svg <- downloadHandler(
      filename = function() {
        paste0(basename(sub("_pca_input\\.rds$", "", input$picked)), "_PCA.svg")
      },
      content = function(file) {
        ggplot2::ggsave(file, plot = plot_obj(), width = 7, height = 5.5, device = "svg")
      }
    )

    output$download_scores_csv <- downloadHandler(
      filename = function() {
        paste0(basename(sub("_pca_input\\.rds$", "", input$picked)), "_PC_scores.csv")
      },
      content = function(file) {
        readr::write_csv(scores_df(), file)
      }
    )

    # return selected IDs if you want to use downstream
    return(list(selected_samples = selected_ids))
  })
}
