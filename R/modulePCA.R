#' @importFrom stats prcomp
PCASidebarUI <- function(id, title = "Principal Component Analysis") {
  ns <- NS(id)
  tagList(
    h4(title),
    helpText("Select one or more saved comparisons. If multiple are selected, a joint PCA is computed on the intersected gene set."),
    uiOutput(ns("pick_files_ui")),
    tags$hr(),
    uiOutput(ns("color_var_ui")),
    uiOutput(ns("shape_var_ui")),
    checkboxInput(ns("draw_ellipses"), "Group ellipses (95%)", FALSE),
    sliderInput(ns("pt_size"), "Point size", min = 1, max = 6, value = 3, step = 0.5),
    radioButtons(ns("engine"), "Plot engine", inline = TRUE,
                 choices = c("Static Plot" = "ggplot", "Interactive Plot" = "plotly"),
                 selected = "ggplot")
  )
}

PCAMainUI <- function(id) {
  ns <- NS(id)
  spin <- function(x) shinycssloaders::withSpinner(x, type = 4, color = "#005249")
  tagList(
    conditionalPanel(
      sprintf("input['%s'] === 'plotly'", ns("engine")),
      spin(plotly::plotlyOutput(ns("pca_plotly"), height = "600px"))
    ),
    conditionalPanel(
      sprintf("input['%s'] === 'ggplot'", ns("engine")),
      spin(plotOutput(ns("pca_plot"), height = "600px"))
    ),
    tags$br(),
    fluidRow(
      column(4, downloadButton(ns("download_png"), "Download PNG")),
      column(4, downloadButton(ns("download_svg"), "Download SVG")),
      column(4, downloadButton(ns("download_scores_csv"), "Download PC scores CSV"))
    )
  )
}

pcaServer <- function(id, pkg = utils::packageName(), data_mode = c("cloud", "local")) {
  data_mode <- match.arg(data_mode)

  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    `%||%` <- function(a, b) if (is.null(a)) b else a

    # Set your package name here if you want it explicit
    pkg <- "FRDATranscriptomicAtlas"

    pretty_map <- tryCatch(
      get("pretty_map", envir = asNamespace(pkg)),
      error = function(e) {
        warning("pretty_map could not be found; using empty vector.")
        character(0)
      }
    )
    aesthetic_map <- c(
      "Study"                  = "Study",
      ".dataset_pretty"        = "Dataset",
      "case_diff_controls"     = "Disease Status",
      "cell_type"              = "Cell Type"
    )
    # tolerant pretty-name resolver
    pretty_lookup <- function(study_key, pretty_map) {
      # 1) exact key
      hit <- unname(pretty_map[study_key])

      # 2) try batchcorrection alias
      if (is.na(hit) || !nzchar(hit)) {
        hit <- unname(pretty_map[paste0(study_key, "_batchcorrection")])
      }

      # 3) fallback
      if (is.na(hit) || !nzchar(hit)) hit <- study_key

      hit
    }


    # ---- ensure PCA (DESeq2 objects) are available ----
    manifest <- reactiveVal(NULL)

    observeEvent(TRUE, {
      m <- if (identical(data_mode, "cloud")) {
        get_pca_manifest_cloud()
      } else {
        # local fallback
        cache_root <- tools::R_user_dir(pkg, "cache")
        paths <- list.files(
          file.path(cache_root, "pca_input"),
          pattern = "_pca_input\\.rds$",
          full.names = TRUE
        )

        tibble::tibble(
          label = sub("_pca_input\\.rds$", "", basename(paths)),
          path = paths,
          filename = basename(paths)
        )
      }

      manifest(m)
    }, once = TRUE)

    pca_files <- reactive({
      req(manifest())

      m <- manifest()

      if (identical(data_mode, "cloud")) {
        tibble::tibble(
          label = m$dataset,
          value = m$filename
        )
      } else {
        tibble::tibble(
          label = m$label,
          value = m$path
        )
      }
    })

    output$pick_files_ui <- renderUI({
      req(nrow(pca_files()) > 0)
      pf <- pca_files()

      # Extract study key (e.g. Lai_PNS, Maddock_SN_FA1, Li)
      study_key <- sub("(_FRDA.*$)", "", pf$label)

      # Default mapping via pretty_map
      pretty_label <- vapply(
        study_key,
        pretty_lookup,
        character(1),
        pretty_map = pretty_map
      )

      # ---- SPECIAL CASE: Li PCA ----
      is_li <- study_key == "Li"
      pretty_label[is_li] <- "Li (Cardiomyocytes)"



      # Named vector: values = file paths, names = pretty labels
      choices_named <- stats::setNames(pf$value, pretty_label)

      # wanted <- c("Indelicato_FRDA_vs_CTRL")
      # sel <- pf$value[match(wanted, pf$label, nomatch = 0)]
      # if (length(sel) == 0) sel <- pf$value[1]

      selectizeInput(
        ns("picked"), "Comparison(s)",
        choices  = choices_named,
        selected = character(0),
        multiple = TRUE,
        options  = list(
          plugins = list("remove_button"),
          placeholder = "Choose one or more PCA datasets..."
        )
      )


    })


    # Load and merge (intersect genes across selections)
    merged_input <- reactive({
      validate(
        need(length(input$picked) >= 1, "Select at least one PCA dataset.")
      )

      paths  <- as.character(input$picked)
      objs <- lapply(paths, function(p) {

        if (identical(data_mode, "cloud")) {
          get_pca_data_cloud_cached(p)
        } else {
          readRDS(p)
        }

      })
      labels <- basename(paths)
      labels <- sub("_pca_input\\.rds$", "", labels)
      labels <- sub("\\.rds$", "", labels)

      # parent study inferred from dataset label
      parent_studies <- sub("_.*$", "", labels)


      # Z-score only if more than one study is present
      zscore_required <- length(unique(parent_studies)) > 1
      ## --------------------------------------------------
      ## SINGLE DATASET -> standard VST PCA
      ## --------------------------------------------------
      if (length(objs) == 1L) {

        X <- objs[[1]]$vsd_mat
        M <- as.data.frame(objs[[1]]$meta)

        M$.dataset <- labels[1]

        study_key <- sub("(_FRDA.*$)", "", labels[1])

        pretty_label <- pretty_lookup(study_key, pretty_map)


        # SPECIAL CASE: Li PCA
        if (identical(study_key, "Li")) {
          pretty_label <- "Li (Cardiomyocytes)"
        }

        M$.dataset_pretty <- pretty_label

        # ---- CRITICAL FIX: align samples ----
        common_ids <- intersect(colnames(X), rownames(M))
        X <- X[, common_ids, drop = FALSE]
        M <- M[common_ids, , drop = FALSE]

        return(list(
          vsd_mat = X,
          meta    = M,
          files   = labels,
          zscore  = FALSE
        ))

      }


      ## --------------------------------------------------
      ## MULTIPLE DATASETS -> intersect genes + Z-score
      ## --------------------------------------------------
      genes_common <- Reduce(
        intersect,
        lapply(objs, function(o) rownames(o$vsd_mat))
      )

      if (length(genes_common) < 200) {
        warning(
          "Very small gene intersection (", length(genes_common),
          "). PCA may be unstable."
        )
      }

      mats  <- list()
      metas <- list()

      for (i in seq_along(objs)) {

        ## ---- subset to shared genes ----
        Xi_raw <- objs[[i]]$vsd_mat[genes_common, , drop = FALSE]

        if (zscore_required) {
          Xi <- t(scale(t(Xi_raw)))
          Xi[is.na(Xi)] <- 0
        } else {
          Xi <- Xi_raw
        }


        ## ---- metadata ----
        Mi <- as.data.frame(objs[[i]]$meta)

        # enforce sample_id from expression matrix ONLY
        Mi$sample_id <- colnames(Xi)

        # make sample IDs globally unique
        new_ids <- paste0(labels[i], "_", Mi$sample_id)

        # apply consistently
        Mi$sample_id <- new_ids
        colnames(Xi) <- new_ids


        Mi$.dataset <- labels[i]

        # ---- dataset pretty name (PCA-specific) ----
        study_key <- sub("(_FRDA.*$)", "", labels[i])

        pretty_label <- pretty_lookup(study_key, pretty_map)


        # SPECIAL CASE: Li PCA
        if (identical(study_key, "Li")) {
          pretty_label <- "Li (Cardiomyocytes)"
        }

        Mi$.dataset_pretty <- pretty_label


        mats[[i]]  <- Xi
        metas[[i]] <- Mi
      }

      ## ---- combine matrices + metadata ----
      Xall <- do.call(cbind, mats)
      Mall <- dplyr::bind_rows(metas)

      # final hard alignment (authoritative)
      Mall <- Mall[match(colnames(Xall), Mall$sample_id), , drop = FALSE]

      stopifnot(identical(colnames(Xall), Mall$sample_id))

      list(
        vsd_mat = Xall,
        meta    = Mall,
        files   = labels,
        zscore  = zscore_required
      )
    })

    ## Decide PCA input matrix (VST vs Z-scored)
    pca_matrix <- reactive({
      X <- merged_input()$vsd_mat
      req(ncol(X) >= 2)
      X
    })


    # Aesthetic fields (categorical-ish)
    meta_cols <- reactive({
      M <- merged_input()$meta

      # only categorical-ish fields
      keep <- vapply(
        M,
        function(x) is.factor(x) || is.character(x) || is.logical(x),
        logical(1)
      )

      M2 <- M[, keep, drop = FALSE]

      # drop high-cardinality junk
      small <- vapply(
        M2,
        function(x) length(unique(x)) <= 50,
        logical(1)
      )

      available <- names(M2[, small, drop = FALSE])

      # intersect with allowed aesthetic variables
      allowed <- intersect(names(aesthetic_map), available)

      # named vector: values = column names, names = pretty labels
      stats::setNames(allowed, aesthetic_map[allowed])
    })

    output$color_var_ui <- renderUI({
      choices <- meta_cols()
      req(length(choices) > 0)

      selectInput(
        ns("color_var"),
        "Colour by",
        choices  = choices,
        selected = ".dataset_pretty"
      )
    })

    output$shape_var_ui <- renderUI({
      choices <- meta_cols()
      req(length(choices) > 0)

      selectInput(
        ns("shape_var"),
        "Shape by",
        choices  = c("None" = "None", choices),
        selected = "case_diff_controls"
      )
    })


    # PCA
    pr_obj <- reactive({
      set.seed(1234)
      prcomp(t(pca_matrix()), center = TRUE, scale. = FALSE)
    })

    percent_var <- reactive({
      pr <- pr_obj()
      round(100 * (pr$sdev^2) / sum(pr$sdev^2), 2)
    })

    scores_df <- reactive({
      df <- as.data.frame(pr_obj()$x)
      df$sample_id <- rownames(df)
      M  <- merged_input()$meta

      # Ensure M has sample_id; in single-dataset branch it doesn't yet, so add it safely
      if (!"sample_id" %in% names(M)) {
        M$sample_id <- rownames(M)
      }

      dplyr::left_join(df, M, by = "sample_id")
    })

    # Build ggplot
    plot_obj <- reactive({
      req(scores_df(), input$color_var, input$pt_size)
      df <- scores_df()
      aes_color <- rlang::sym(input$color_var)

      p <- ggplot2::ggplot(df, ggplot2::aes(PC1, PC2, colour = !!aes_color))

      if (!is.null(input$shape_var) && input$shape_var != "None") {
        aes_shape <- rlang::sym(input$shape_var)
        p <- p + ggplot2::aes(shape = !!aes_shape)
      }

      p <- p + ggplot2::geom_point(size = input$pt_size, alpha = 0.9)


      if (isTRUE(input$draw_ellipses)) {
        p <- p + ggplot2::stat_ellipse(
          ggplot2::aes(group = !!aes_color),
          level = 0.95,
          linetype = 2
        )
      }

      pv <- percent_var()
      title_txt <- if (!merged_input()$zscore) {
        "PCA (VST)"
      } else {
        "PCA (Z-scored across studies)"
      }


      color_label <- aesthetic_map[[input$color_var]] %||% input$color_var
      shape_label <- if (!is.null(input$shape_var) && input$shape_var != "None") {
        aesthetic_map[[input$shape_var]] %||% input$shape_var
      } else {
        NULL
      }

      p <- p +
        ggplot2::labs(
          title = title_txt,
          x = sprintf("PC1 (%.2f%%)", pv[1]),
          y = sprintf("PC2 (%.2f%%)", pv[2]),
          color = color_label,
          shape = shape_label
        )

      if (exists("theme_Marnie", mode = "function", inherits = TRUE)) {
        p <- p + theme_Marnie()
      }

      p
    })

    output$pca_plot <- renderPlot({
      validate(
        need(nrow(scores_df()) > 0, "Not enough samples to compute PCA.")
      )
      plot_obj()
    })

    output$pca_plotly <- plotly::renderPlotly({
      validate(
        need(nrow(scores_df()) > 0, "Not enough samples to compute PCA.")
      )
      plt <- plotly::ggplotly(
        plot_obj(),
        tooltip = c("sample_id", input$color_var)
      )
      plotly::layout(plt, dragmode = "lasso")
    })

    # Downloads
    filename_stub <- reactive({
      paste0(
        "PCA_",
        format(Sys.time(), "%Y%m%d_%H%M%S")
      )
    })

    output$download_png <- downloadHandler(
      filename = function() sprintf("%s_PCA.png", filename_stub()),
      content  = function(file) {
        ggplot2::ggsave(file, plot = plot_obj(), width = 7, height = 5.5, dpi = 300)
      }
    )

    output$download_svg <- downloadHandler(
      filename = function() sprintf("%s_PCA.svg", filename_stub()),
      content  = function(file) {
        ggplot2::ggsave(file, plot = plot_obj(), width = 7, height = 5.5, device = "svg")
      }
    )

    output$download_scores_csv <- downloadHandler(
      filename = function() sprintf("%s_PC_scores.csv", filename_stub()),
      content  = function(file) {
        readr::write_csv(scores_df(), file)
      }
    )
  })
}
