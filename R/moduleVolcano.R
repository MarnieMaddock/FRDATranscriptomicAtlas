#moduleVolcano.R
#' Volcano Plot — Sidebar UI
#' @importFrom rlang %||%
#' @param id module id
#' @param pretty_map named character vector dataset_key = "Pretty Label"
#' @noRd
volcanoSidebarUI <- function(id, pretty_map) {
  ns <- shiny::NS(id)

  choices_keyed <- stats::setNames(names(pretty_map), pretty_map)


  shiny::tagList(
    shiny::h4("Volcano plot"),
    shiny::selectizeInput(
      ns("dataset"),
      label = "Dataset",
      choices = choices_keyed,
      selected = names(pretty_map)[1],
      multiple = FALSE,
      options = list(placeholder = "Choose a dataset…")
    ),
    shiny::hr(),
    shiny::radioButtons(
      ns("deg_level"),
      label = "Differential level",
      choices = c(
        "Genes"      = "genes",
        "Isoforms"   = "transcripts"
      ),
      selected = "genes",
      inline = TRUE
    ),
    shiny::hr(),
    shiny::sliderInput(
      ns("lfc_thresh"),
      label = "Absolute log2FC threshold",
      min = 0, max = 4, value = 1, step = 0.1
    ),
    shiny::sliderInput(
      ns("padj_thresh"),
      label = "Adjusted p-value (FDR) threshold",
      min = 1e-6, max = 0.25, value = 0.05, step = 0.005
    ),
    shiny::numericInput(
      ns("label_topn"),
      label = "Max labels for tops (by |log2FC|)",
      value = 20, min = 0, max = 200, step = 1
    ),
    shiny::textInput(
      ns("highlight_genes"),
      label = "Highlight genes (comma-separated symbols)",
      placeholder = "e.g. FXN, TP53, NFE2L2"
    ),
    shiny::checkboxInput(
      ns("show_ns"),
      label = "Show non-significant points",
      value = TRUE
    ),
    tags$br(),
    shiny::helpText(
      shiny::HTML(
        "<b>Notes</b><ul style='margin-top:4px'>
          <li>Y-axis uses −log10(FDR) and is <b>capped at 50</b> to avoid extreme values blowing out the scale.</li>
          <li>Use the Plotly toolbar to <b>Box Select</b> or <b>Lasso Select</b> points. This is in the top right corner of the plot. Selected genes appear in the table.</li>
        </ul>"
      )
    ),
  )
}

#' Volcano Plot — Main UI
#' @param id module id
#' @noRd
volcanoMainUI <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shinycssloaders::withSpinner(
      plotly::plotlyOutput(ns("volcano_plot"), height = "650px"),
      type = 4, color = "#005249"
    ),
    shiny::br(),
    shiny::fluidRow(
      shiny::column(
        width = 6,
        shiny::h5("Points you lasso/box-select"),
        DT::DTOutput(ns("selected_table"))
      ),
      shiny::column(
        width = 6,
        shiny::h5("Downloads"),
        shiny::div(
          style = "display:flex; gap:8px; flex-wrap:wrap;",
          shiny::downloadButton(ns("dl_plot_svg"), "Download SVG"),
          shiny::downloadButton(ns("dl_plot_png"), "Download PNG")
        ),
        shiny::br(),
        shiny::verbatimTextOutput(ns("summary_text"))
      )
    )
  )
}

#' Volcano Plot — Server
#' @param id module id
#' @param level "genes" or "transcripts" (default "genes")
#' @param pkg package name for system.file lookup (defaults to calling package)
#' @param custom_loader optional function(dataset, level, pkg) -> data.frame
#'        Return columns: gene, log2FC, padj (or pvalue). You may ignore padj if not used.
#' @noRd
# Volcano module (drop-in)
volcanoServer <- function(
    id,
    level = "genes",                    # "genes" or "transcripts"
    pkg   = utils::packageName()
) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns
    `%||%` <- function(a, b) if (is.null(a)) b else a
    cur_level <- function() {
      if (!is.null(input$deg_level) && nzchar(input$deg_level)) {
        input$deg_level
      } else if (is.function(level)) {
        level()
      } else {
        level
      }
    }


    .safe_num <- function(x) suppressWarnings(as.numeric(as.character(x)))
    .norm_filename <- function(x) gsub("[^A-Za-z0-9_\\-]+", "_", x)

    # -------- tx2gene cache (optional) --------
    .tx2_cache <- NULL
    .get_tx2 <- function(pkg) {
      if (!is.null(.tx2_cache)) return(.tx2_cache)
      p <- system.file("extdata", "maps", "tx2gene.tsv", package = pkg, mustWork = FALSE)
      if (nzchar(p) && file.exists(p)) {
        .tx2_cache <<- readr::read_tsv(p, col_types = "ccc") |>
          dplyr::mutate(
            transcript_id = sub("\\.\\d+$", "", transcript_id),
            gene_id       = sub("\\.\\d+$", "", gene_id)
          )
      } else {
        .tx2_cache <<- NULL
      }
      .tx2_cache
    }

    .add_display_label <- function(df, lvl, pkg) {

      # Ensure gene column exists
      if (!"gene" %in% names(df)) {
        if ("gene_id" %in% names(df)) df$gene <- df$gene_id
        else if (!is.null(rownames(df))) df$gene <- rownames(df)
        else df$gene <- NA_character_
      }

      # Always strip version suffixes for display
      df$gene <- sub("\\.\\d+$", "", df$gene)

      # ---- transcript mode: show transcript ID ----
      if (identical(lvl, "transcripts")) {
        df$display_label <- df$gene
        return(df)
      }

      # ---- gene mode: map to symbol if possible ----
      m <- .get_tx2(pkg)

      if (is.null(m) || !"gene_name" %in% names(m)) {
        df$display_label <- df$gene
        return(df)
      }

      sym <- m$gene_name[match(df$gene, m$gene_id)]
      df$display_label <- dplyr::coalesce(sym, df$gene)

      df
    }

    .deg_data_ready <- reactiveValues(
      genes = FALSE,
      transcripts = FALSE
    )

    # -------- robust loader (data frame files) --------
    # -------- robust loader (Zenodo-backed DEG files) --------
    load_deg <- function(dataset, level, pkg) {

      level <- match.arg(level, c("genes", "transcripts"))

      # ---- map level -> manifest key + cache subdir ----
      manifest_key <- if (level == "genes") "deg_genes" else "deg_transcripts"
      level_dir    <- if (level == "genes") "genes"     else "txs"

      # ---- ensure data are present locally (download if needed) ----
      if (!.deg_data_ready[[level]]) {
        ensure_atlas_data(
          keys    = manifest_key,
          package = pkg
        )
        .deg_data_ready[[level]] <- TRUE
      }


      # ---- cached DEG directory ----
      cache_root <- tools::R_user_dir(pkg, which = "cache")
      deg_dir    <- file.path(cache_root, "deg", level_dir)

      if (!dir.exists(deg_dir)) {
        stop(
          "Cached DEG directory not found after download:\n  ",
          deg_dir,
          call. = FALSE
        )
      }

      # ---- filename patterns ----
      pat_exact <- paste0("^DESEQ2_res_", dataset, "_0\\.05_all_", level, "\\.rds$")
      pat_loose <- paste0("^DESEQ2_res_", dataset, "_0\\.[0-9]+_all_", level, "\\.rds$")
      simple    <- paste0("^", dataset, "\\.rds$")

      files <- list.files(deg_dir, full.names = TRUE)

      f <- files[grepl(pat_exact, basename(files))]
      if (!length(f)) f <- files[grepl(pat_loose, basename(files))]
      if (!length(f)) f <- files[grepl(simple, basename(files))]

      if (!length(f)) {
        stop(
          "DEG file not found for dataset '", dataset, "'.\n",
          "Looked in:\n  ", deg_dir, "\n",
          "Expected one of:\n",
          "  • DESEQ2_res_", dataset, "_0.05_all_", level, ".rds\n",
          "  • DESEQ2_res_", dataset, "_0.xx_all_", level, ".rds\n",
          "  • ", dataset, ".rds",
          call. = FALSE
        )
      }

      # ---- load file ----
      x <- readRDS(f[1])
      if (!is.data.frame(x)) x <- as.data.frame(x)
      nm <- names(x)

      # ---- standardise ID column ----
      if (!"gene" %in% nm) {
        if ("symbol" %in% nm)                  x <- dplyr::rename(x, gene = symbol)
        else if ("Gene" %in% nm)               x <- dplyr::rename(x, gene = Gene)
        else if ("external_gene_name" %in% nm) x <- dplyr::rename(x, gene = external_gene_name)
        else if ("ensembl_gene_id" %in% nm)    x <- dplyr::rename(x, gene = ensembl_gene_id)
        else if ("gene_id" %in% nm)            x <- dplyr::rename(x, gene = gene_id)
        else if (!is.null(rownames(x)))        x$gene <- rownames(x)
        else                                   x$gene <- NA_character_
      }

      # ---- standardise log2FC ----
      if (!"log2FC" %in% nm) {
        if ("log2FoldChange" %in% nm) x <- dplyr::rename(x, log2FC = log2FoldChange)
        else if ("beta" %in% nm)       x <- dplyr::rename(x, log2FC = beta)
      }

      # ---- sanity check for p-values ----
      if (!("padj" %in% names(x) || "pvalue" %in% names(x))) {
        stop(
          "Need 'padj' or 'pvalue' column in DEG file:\n  ",
          basename(f[1]),
          call. = FALSE
        )
      }

      x
    }



    # -------- main reactive table --------
    deg_tbl <- shiny::reactive({
      req(input$dataset)
      df <- load_deg(input$dataset, cur_level(), pkg)

      # numeric & choose p column
      if ("padj" %in% names(df)) {
        df$P <- .safe_num(df$padj)
      } else {
        df$P <- .safe_num(df$pvalue)
      }
      df$log2FC <- .safe_num(df$log2FC)

      # derived values
      df$negLog10P <- -log10(pmax(df$P, .Machine$double.eps))
      df$negLog10P <- pmin(df$negLog10P, 50)  # cap at 50 (UI note)

      df$row_id <- seq_len(nrow(df))

      # ensure display gene (symbol if available)
      df <- .add_display_label(df, cur_level(), pkg)

      # status
      lfc_th  <- input$lfc_thresh %||% 1
      padj_th <- input$padj_thresh %||% 0.05
      df$status <- dplyr::case_when(
        df$P <= padj_th & df$log2FC >=  lfc_th ~ "Up",
        df$P <= padj_th & df$log2FC <= -lfc_th ~ "Down",
        TRUE ~ "NS"
      )
      df$status <- factor(df$status, levels = c("Down","NS","Up"))

      # highlights (by symbol/display)
      hi <- trimws(unlist(strsplit(input$highlight_genes %||% "", ",")))
      hi <- hi[nzchar(hi)]
      df$highlight <- if (length(hi)) tolower(df$display_label) %in% tolower(hi) else FALSE

      df$.pcol <- if ("padj" %in% names(df)) "padj" else "pvalue"
      df
    })

    # -------- summary text --------
    output$summary_text <- shiny::renderText({
      df <- deg_tbl(); req(nrow(df))
      padj_th <- input$padj_thresh; lfc_th <- input$lfc_thresh
      ns_n <- sum(df$status == "NS", na.rm = TRUE)
      up_n <- sum(df$status == "Up", na.rm = TRUE)
      dn_n <- sum(df$status == "Down", na.rm = TRUE)
      paste0(
        "n = ", nrow(df), " total; Up = ", up_n, ", Down = ", dn_n,
        ", NS = ", ns_n, "  |  thresholds: |log2FC| ≥ ", lfc_th,
        ", ", df$.pcol[1], " ≤ ", signif(padj_th, 3)
      )
    })

    # -------- gg builder (used for both plotly and downloads) --------
    make_gg <- function(df, show_ns = TRUE) {
      dplot <- if (isTRUE(show_ns)) df else dplyr::filter(df, status != "NS")
      dplot <- dplyr::filter(dplot, is.finite(log2FC), is.finite(negLog10P))

      col_map <- c("Down" = "#1f77b4", "NS" = "grey80", "Up" = "#d62728")

      p <- ggplot2::ggplot(
        dplot,
        ggplot2::aes(x = log2FC, y = negLog10P, color = status, key = row_id)
      ) +
        ggplot2::geom_point(alpha = 0.8, size = 1.6, stroke = 0) +
        ggplot2::scale_color_manual(values = col_map, drop = FALSE) +
        ggplot2::geom_vline(xintercept = c(-input$lfc_thresh, input$lfc_thresh),
                            linetype = "dashed", linewidth = 0.4) +
        ggplot2::geom_hline(yintercept = -log10(input$padj_thresh),
                            linetype = "dashed", linewidth = 0.4) +
        ggplot2::labs(x = "log2 Fold Change", y = "-log10(p.adj)",
                      subtitle = if (cur_level() == "transcripts") "Isoform-level differential expression"
                      else "Gene-level differential expression",
                      color = NULL)

      if (exists("theme_Marnie", inherits = TRUE)) p <- p + get("theme_Marnie", inherits = TRUE) else p <- p + ggplot2::theme_bw()
      p + ggplot2::theme(legend.position = "top", panel.grid.minor = ggplot2::element_blank())

      # ring-highlight requested genes
      if (any(dplot$highlight, na.rm = TRUE)) {
        p <- p + ggplot2::geom_point(
          data = dplot[dplot$highlight %in% TRUE, ],
          ggplot2::aes(x = log2FC, y = negLog10P),
          inherit.aes = FALSE, size = 3.2, shape = 21, fill = NA, color = "black", stroke = 0.6
        )
      }
      p
    }


    # -------- interactive plot --------
    output$volcano_plot <- plotly::renderPlotly({
      df <- deg_tbl(); req(nrow(df))
      df$hover_txt <- paste0(
        "<b>", df$display_label, "</b>",
        "<br>log2FC: ", sprintf("%.3f", df$log2FC),
        "<br>-log10(p): ", sprintf("%.3f", df$negLog10P),
        "<br>", toupper(df$.pcol[1]), ": ", ifelse(is.na(df$P), "NA", signif(df$P, 3)),
        "<br>Status: ", df$status
      )
      gp <- make_gg(df, show_ns = isTRUE(input$show_ns)) + ggplot2::aes(text = hover_txt)
      plt <- plotly::ggplotly(gp, tooltip = "text", dynamicTicks = TRUE)
      plt$x$source <- "volc"
      plt <- plotly::config(plt, modeBarButtonsToAdd = c("select2d", "lasso2d"))
      plotly::event_register(plt, "plotly_selected")
    })

    # -------- selection table --------
    output$selected_table <- DT::renderDT({
      df <- deg_tbl(); req(nrow(df))
      ev <- tryCatch(plotly::event_data("plotly_selected", source = "volc"), error = function(e) NULL)
      if (is.null(ev) || !nrow(ev) || is.null(ev$key)) {
        sel <- df[0, ]
      } else {
        sel <- df[df$row_id %in% as.integer(unique(ev$key)), , drop = FALSE]
      }

      id_col <- if (cur_level() == "transcripts") "Transcript ID" else "Gene"

      DT::datatable(
        sel |>
          dplyr::transmute(
            !!id_col := display_label,
            log2FC,
            padj = P,
            status
          ) |>
          dplyr::arrange(padj, dplyr::desc(abs(log2FC))),
        rownames = FALSE,
        options = list(pageLength = 10, scrollX = TRUE, dom = "tip"),
        filter = "top"
      )
    })

    # -------- downloads --------
    make_base_name <- function(suffix) {
      ds <- input$dataset %||% "dataset"
      paste0(.norm_filename(ds), "_volcano_F", input$lfc_thresh,
             "_", (deg_tbl()$.pcol[1]), input$padj_thresh, "_", suffix)
    }

    output$dl_plot_svg <- shiny::downloadHandler(
      filename = function() paste0(make_base_name("plot"), ".svg"),
      content = function(file) {
        df <- deg_tbl()
        gp <- make_gg(df, show_ns = isTRUE(input$show_ns))
        svglite::svglite(file, width = 8, height = 6)
        on.exit(grDevices::dev.off(), add = TRUE)
        print(gp)
      }
    )

    output$dl_plot_png <- shiny::downloadHandler(
      filename = function() paste0(make_base_name("plot"), ".png"),
      content = function(file) {
        df <- deg_tbl()
        gp <- make_gg(df, show_ns = isTRUE(input$show_ns))
        ggplot2::ggsave(filename = file, plot = gp, width = 8, height = 6, dpi = 300)
      }
    )
  })
}
