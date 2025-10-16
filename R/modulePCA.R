# ---- PCA module (multi-select) ---------------------------------------------
# Files expected: <label>_pca_input.rds with list(vsd_mat, meta, percent_var)

PCASidebarUI <- function(id, title = "PCA (VST)") {
  ns <- NS(id)
  tagList(
    h4(title),
    helpText("Select one or more saved comparisons. If multiple are selected, a joint PCA is computed on the intersected gene set."),
    uiOutput(ns("pick_files_ui")),
    tags$hr(),
    uiOutput(ns("color_var_ui")),
    uiOutput(ns("shape_var_ui")),
    checkboxInput(ns("label_points"), "Show sample labels", FALSE),
    checkboxInput(ns("draw_ellipses"), "Group ellipses (95%)", FALSE),
    sliderInput(ns("pt_size"), "Point size", min = 1, max = 6, value = 3, step = 0.5),
    radioButtons(ns("engine"), "Plot engine", inline = TRUE,
                 choices = c("Static Plot" = "ggplot", "Interactive Plot" = "plotly"),
                 selected = "ggplot")
  )
}

PCAMainUI <- function(id) {
  ns <- NS(id)
  tagList(
    conditionalPanel(
      sprintf("input['%s'] === 'plotly'", ns("engine")),
      plotly::plotlyOutput(ns("pca_plotly"), height = "600px")
    ),
    conditionalPanel(
      sprintf("input['%s'] === 'ggplot'", ns("engine")),
      plotOutput(ns("pca_plot"), height = "600px")
    ),
    tags$br(),
    fluidRow(
      column(4, downloadButton(ns("download_png"), "Download PNG")),
      column(4, downloadButton(ns("download_svg"), "Download SVG")),
      column(4, downloadButton(ns("download_scores_csv"), "Download PC scores CSV"))
    )
  )
}

pcaServer <- function(id,
                      data_dir = system.file("extdata/deseq_objects", package = "FRDATranscriptomicAtlas"),
                      fallback_dir = "inst/extdata/deseq_objects") {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Locate directory with *_pca_input.rds
    dir_use <- reactive({
      if (dir.exists(data_dir)) return(data_dir)
      if (dir.exists(fallback_dir)) return(fallback_dir)
      stop("No PCA directory found. Checked:\n  - ", data_dir, "\n  - ", fallback_dir)
    })

    # List available PCA inputs
    pca_files <- reactive({
      paths <- list.files(dir_use(), pattern = "_pca_input\\.rds$", full.names = TRUE)
      tibble::tibble(
        label = sub("_pca_input\\.rds$", "", basename(paths)),
        path  = paths
      ) |>
        dplyr::arrange(label)
    })

    output$pick_files_ui <- renderUI({
      req(nrow(pca_files()) > 0)
      pf <- pca_files()   # tibble(label, path)

      wanted <- c("Lai_iPSC_FRDA_vs_IC", "Lai_CNS_FRDA_vs_IC", "Lai_PNS_FRDA_vs_IC")
      sel <- pf$path[match(wanted, pf$label, nomatch = 0)]   # keep those that exist

      if (length(sel) == 0) sel <- pf$path[1]  # fallback

      selectizeInput(
        ns("picked"), "Comparison(s)",
        choices  = setNames(pf$path, pf$label),
        selected = sel,
        multiple = TRUE,
        options  = list(plugins = list("remove_button"))
      )
    })


    # Load and merge (intersect genes across selections)
    merged_input <- reactive({
      req(input$picked)
      paths  <- as.character(input$picked)
      objs   <- lapply(paths, readRDS)
      labels <- sub("_pca_input\\.rds$", "", basename(paths))

      if (length(objs) == 1L) {
        X <- objs[[1]]$vsd_mat
        M <- as.data.frame(objs[[1]]$meta)
        M$.dataset <- labels[1]
        if (anyDuplicated(rownames(M))) {
          rn <- make.unique(rownames(M))
          rownames(M) <- rn
          colnames(X) <- rn
        }
        return(list(vsd_mat = X, meta = M, files = labels))
      }

      genes_common <- Reduce(intersect, lapply(objs, \(o) rownames(o$vsd_mat)))
      if (length(genes_common) < 200) {
        warning("Very small gene intersection (", length(genes_common), "). PCA may be unstable.")
      }

      mats  <- list()
      metas <- list()
      for (i in seq_along(objs)) {
        Xi <- objs[[i]]$vsd_mat[genes_common, , drop = FALSE]
        Mi <- as.data.frame(objs[[i]]$meta)
        new_ids <- paste0(labels[i], "_", colnames(Xi))  # keep IDs unique
        colnames(Xi) <- new_ids
        rownames(Mi) <- new_ids
        Mi$.dataset  <- labels[i]
        mats[[i]]  <- Xi
        metas[[i]] <- Mi
      }
      Xall <- do.call(cbind, mats)
      Mall <- dplyr::bind_rows(metas)
      Mall <- Mall[colnames(Xall), , drop = FALSE]

      list(vsd_mat = Xall, meta = Mall, files = labels)
    })

    # Aesthetic fields (categorical-ish)
    meta_cols <- reactive({
      M <- merged_input()$meta
      keep  <- vapply(M, \(x) is.factor(x) || is.character(x) || is.logical(x), logical(1))
      M2    <- M[, keep, drop = FALSE]
      small <- vapply(M2, \(x) length(unique(x)) <= 50, logical(1))
      nm    <- names(M2[, small, drop = FALSE])
      unique(c(".dataset", nm))
    })

    output$color_var_ui <- renderUI({
      req(meta_cols())
      preferred <- c(".dataset", "group", "FRDA_CTRL", "cell_type", "Study_Alias")
      default   <- intersect(preferred, meta_cols())
      selectInput(ns("color_var"), "Colour by",
                  choices = meta_cols(),
                  selected = if (length(default)) default[1] else meta_cols()[1])
    })

    output$shape_var_ui <- renderUI({
      req(meta_cols())
      choices <- c("None", meta_cols())
      selectInput(ns("shape_var"), "Shape by", choices = choices, selected = "None")
    })

    # PCA
    pr_obj <- reactive({
      X <- merged_input()$vsd_mat
      validate(need(ncol(X) >= 2, "Need at least two samples to compute PCA."))
      set.seed(1234)
      prcomp(t(X), center = TRUE, scale. = FALSE)
    })

    percent_var <- reactive({
      pr <- pr_obj()
      round(100 * (pr$sdev^2) / sum(pr$sdev^2), 2)
    })

    scores_df <- reactive({
      df <- as.data.frame(pr_obj()$x)
      df$sample_id <- rownames(df)
      M  <- merged_input()$meta
      M$sample_id <- rownames(M)
      dplyr::left_join(df, M, by = "sample_id")
    })

    # Build ggplot
    plot_obj <- reactive({
      req(scores_df(), input$color_var, input$pt_size)
      df <- scores_df()
      aes_color <- rlang::sym(input$color_var)

      if (!is.null(input$shape_var) && input$shape_var != "None") {
        aes_shape <- rlang::sym(input$shape_var)
        p <- ggplot2::ggplot(df, ggplot2::aes(PC1, PC2, colour = !!aes_color, shape = !!aes_shape)) +
          ggplot2::geom_point(size = input$pt_size, alpha = 0.9)
      } else {
        p <- ggplot2::ggplot(df, ggplot2::aes(PC1, PC2, colour = !!aes_color)) +
          ggplot2::geom_point(size = input$pt_size, alpha = 0.9)
      }

      if (isTRUE(input$label_points)) {
        p <- p + ggrepel::geom_text_repel(ggplot2::aes(label = sample_id), size = 3, max.overlaps = 50)
      }
      if (isTRUE(input$draw_ellipses)) {
        p <- p + ggplot2::stat_ellipse(ggplot2::aes(group = !!aes_color), level = 0.95, linetype = 2)
      }

      pv <- percent_var()
      p <- p +
        ggplot2::labs(
          x = sprintf("PC1 (%.2f%%)", pv[1]),
          y = sprintf("PC2 (%.2f%%)", pv[2]),
          color = input$color_var,
          shape = if (!is.null(input$shape_var) && input$shape_var != "None") input$shape_var else NULL
        ) +
        ggplot2::theme_classic(base_size = 14)

      if (exists("theme_Marnie", mode = "function") || exists("theme_Marnie")) p <- p + theme_Marnie
      p
    })


    output$pca_plot <- renderPlot({
      plot_obj()
    })

    # Plotly (lasso)
    output$pca_plotly <- plotly::renderPlotly({
      gg  <- plot_obj()
      plt <- plotly::ggplotly(gg, tooltip = c("sample_id", input$color_var))
      plotly::layout(plt, dragmode = "lasso")
    })

    # Downloads
    filename_stub <- reactive({
      labs <- merged_input()$files
      if (length(labs) == 1L) labs else paste0("MULTI_", paste(labs, collapse = "_"))
    })

    output$download_png <- downloadHandler(
      filename = function() sprintf("%s_PCA.png", filename_stub()),
      content  = function(file) ggplot2::ggsave(file, plot = plot_obj(), width = 7, height = 5.5, dpi = 300)
    )

    output$download_svg <- downloadHandler(
      filename = function() sprintf("%s_PCA.svg", filename_stub()),
      content  = function(file) ggplot2::ggsave(file, plot = plot_obj(), width = 7, height = 5.5, device = "svg")
    )

    output$download_scores_csv <- downloadHandler(
      filename = function() sprintf("%s_PC_scores.csv", filename_stub()),
      content  = function(file) readr::write_csv(scores_df(), file)
    )
  })
}
