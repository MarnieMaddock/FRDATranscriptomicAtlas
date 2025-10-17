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
    shiny::actionButton(ns("update_plot"), "Update plot"),
    tags$br(),
    shiny::helpText(
      shiny::HTML(
        "<b>Notes</b><ul style='margin-top:4px'>
          <li>Y-axis uses −log10(FDR) and is <b>capped at 50</b> to avoid extreme values blowing out the scale.</li>
          <li>Use the Plotly toolbar to <b>Box Select</b> or <b>Lasso Select</b> points. This is in the top right corner of the plot. Selected genes appear in the table.</li>
          <li>After changing thresholds, click <b>Update plot</b> to refresh point colours.</li>
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
volcanoServer <- function(
    id,
    level = "genes",
    pkg = utils::packageName(),
    custom_loader = deseq_loader
) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns
    # helper: allow either string or reactive
    cur_level <- function() if (is.function(level)) level() else level

    # ---- pretty helpers ----
    .safe_num <- function(x) as.numeric(as.character(x))
    .norm_filename <- function(x) gsub("[^A-Za-z0-9_\\-]+", "_", x)

    # ---- tx2gene mapping helpers ----
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

    .add_display_gene <- function(df, lvl, pkg) {
      m <- .get_tx2(pkg)
      if (is.null(m) || !all(c("transcript_id","gene_id","gene_name") %in% names(m))) {
        df$display_gene <- df$gene
        return(df)
      }
      ids <- sub("\\.\\d+$", "", df$gene)
      if (identical(lvl, "transcripts")) {
        sym <- m$gene_name[ match(ids, m$transcript_id) ]
      } else {
        sym <- m$gene_name[ match(ids, m$gene_id) ]
      }
      df$display_gene <- dplyr::coalesce(sym, df$gene)
      df
    }


    # Prefer a custom loader if you pass one; else read from inst/extdata/deg/<level>/
    load_deg <- function(dataset, level, pkg) {
      if (is.function(custom_loader)) {
        out <- custom_loader(dataset, level, pkg)
      } else {
        f <- system.file("extdata", "deg", level, paste0(dataset, ".rds"), package = pkg)
        if (!nzchar(f) || !file.exists(f)) {
          stop("DEG file not found for dataset '", dataset, "' at: ", f)
        }
        out <- readRDS(f)
      }

      # allow level as string or reactive
      .cur_level <- function() if (is.function(level)) level() else level


      # ---- Standardise column names here if needed ----
      # EXPECTED: gene (symbol), log2FC, padj or pvalue
      # Try to coerce if common variants exist:
      nm <- names(out)
      if (!"gene" %in% nm) {
        if ("symbol" %in% nm) out <- dplyr::rename(out, gene = symbol)
        else if ("Gene" %in% nm) out <- dplyr::rename(out, gene = Gene)
      }
      if (!"log2FC" %in% nm) {
        if ("log2FoldChange" %in% nm) out <- dplyr::rename(out, log2FC = log2FoldChange)
      }
      if (!"padj" %in% names(out) && !"pvalue" %in% names(out)) {
        stop("Need either 'padj' or 'pvalue' column in DEG table for dataset ", dataset)
      }

      out
    }

    # Cache/load per dataset
    deg_tbl <- shiny::eventReactive(
      input$update_plot,
      {
        req(input$dataset)
        df <- load_deg(input$dataset, cur_level(), pkg)

        # Ensure numeric
        if ("padj" %in% names(df)) df$padj <- .safe_num(df$padj)
        if ("pvalue" %in% names(df)) df$pvalue <- .safe_num(df$pvalue)
        df$log2FC <- .safe_num(df$log2FC)

        # Pick p column
        validate(need("padj" %in% names(df), "padj column is required for this plot."))
        df$P <- df$padj
        pcol <- "padj"

        # always use padj (required)
        df$P <- df$padj
        df$negLog10P <- -log10(pmax(df$P, .Machine$double.eps))
        if (isTRUE(input$clip_extremes)) df$negLog10P <- pmin(df$negLog10P, 50)

        df$row_id <- seq_len(nrow(df))

        # map Ensembl → symbol for hover/labels  👈 ADD THIS
        df <- .add_display_gene(df, cur_level(), pkg)


        # status factor
        lfc_th <- input$lfc_thresh %||% 1
        padj_th <- input$padj_thresh %||% 0.05
        df$status <- dplyr::case_when(
          df$P <= padj_th & df$log2FC >=  lfc_th ~ "Up",
          df$P <= padj_th & df$log2FC <= -lfc_th ~ "Down",
          TRUE ~ "NS"
        )
        df$status <- factor(df$status, levels = c("Down","NS","Up"))

        # highlight vector (match against symbols)
        hi <- trimws(unlist(strsplit(input$highlight_genes %||% "", ",")))
        hi <- hi[nzchar(hi)]

        # add display symbols first
        df <- .add_display_gene(df, cur_level(), pkg)

        if (length(hi)) {
          df$highlight <- tolower(df$display_gene) %in% tolower(hi)
        } else {
          df$highlight <- FALSE
        }



        df$.pcol <- pcol
        df
      },
      ignoreInit = FALSE
    )

    # Text summary
    output$summary_text <- shiny::renderText({
      df <- deg_tbl()
      req(nrow(df))
      padj_th <- input$padj_thresh
      lfc_th  <- input$lfc_thresh
      ns_n  <- sum(df$status == "NS", na.rm = TRUE)
      up_n  <- sum(df$status == "Up", na.rm = TRUE)
      dn_n  <- sum(df$status == "Down", na.rm = TRUE)
      paste0(
        "n = ", nrow(df), " total; Up = ", up_n, ", Down = ", dn_n,
        ", NS = ", ns_n, "  |  thresholds: |log2FC| ≥ ", lfc_th,
        ", ", df$.pcol[1], " ≤ ", signif(padj_th, 3)
      )
    })

    # Build ggplot volcano
    make_gg <- function(df, show_ns = TRUE, label_sig = TRUE, label_topn = 20) {
      dplot <- if (isTRUE(show_ns)) df else dplyr::filter(df, status != "NS")
      dplot <- dplyr::filter(dplot, !is.na(log2FC), !is.na(negLog10P))  # drop NA points

      col_map <- c("Down" = "#1f77b4", "NS" = "grey80", "Up" = "#d62728")

      p <- ggplot2::ggplot(
        dplot,
        ggplot2::aes(x = log2FC, y = negLog10P, color = status, key = row_id)  # << key!
      ) +
        ggplot2::geom_point(alpha = 0.8, size = 1.6, stroke = 0) +
        ggplot2::scale_color_manual(values = col_map, drop = FALSE) +
        ggplot2::geom_vline(xintercept = c(-input$lfc_thresh, input$lfc_thresh),
                            linetype = "dashed", linewidth = 0.4) +
        ggplot2::geom_hline(yintercept = -log10(input$padj_thresh),
                            linetype = "dashed", linewidth = 0.4) +
        ggplot2::labs(x = "log2 Fold Change", y = "-log10(p.adj)", color = NULL)

      if (exists("theme_Marnie", inherits = TRUE)) p <- p + get("theme_Marnie", inherits = TRUE) else p <- p + ggplot2::theme_bw()
      p + ggplot2::theme(legend.position = "top", panel.grid.minor = ggplot2::element_blank())

      # Highlight ring for user-specified genes
      if (any(dplot$highlight, na.rm = TRUE)) {
        p <- p + ggplot2::geom_point(
          data = dplot[dplot$highlight %in% TRUE, ],
          ggplot2::aes(x = log2FC, y = negLog10P),
          inherit.aes = FALSE, size = 3.2, shape = 21, fill = NA, color = "black", stroke = 0.6
        )
      }

      p
    }


    # --- PLOT ---
    output$volcano_plot <- plotly::renderPlotly({
      df <- deg_tbl(); req(nrow(df))
      # force re-evaluation of these so the lines and colours refresh
      lfc_th <- input$lfc_thresh
      padj_th <- input$padj_thresh

      # one clean hover field (no duplicates)
      df$hover_txt <- paste0(
        "<b>", df$display_gene, "</b>",
        "<br>log2FC: ", sprintf("%.3f", df$log2FC),
        "<br>-log10(p): ", sprintf("%.3f", df$negLog10P),
        "<br>FDR (padj): ", ifelse(is.na(df$P), "NA", signif(df$P, 3)),
        "<br>Status: ", df$status
      )

      gp <- make_gg(
        df,
        show_ns   = isTRUE(input$show_ns),
        label_sig = FALSE,                        # no ggrepel in interactive plotly
        label_topn = 0
      ) + ggplot2::aes(text = hover_txt)     # attach custom hover

      plt <- plotly::ggplotly(gp, tooltip = "text", dynamicTicks = TRUE)
      plt$x$source <- "volc"                      # plain, non-namespaced id
      plt <- plotly::config(plt, modeBarButtonsToAdd = c("select2d", "lasso2d"))
      plt <- plotly::event_register(plt, "plotly_selected")
      plt
    })


    # --- SELECTED POINTS TABLE ---
    output$selected_table <- DT::renderDT({
      df <- deg_tbl(); req(nrow(df))

      # Be tolerant while the plot initializes
      ev <- tryCatch(plotly::event_data("plotly_selected", source = "volc"),
                     error = function(e) NULL)

      if (is.null(ev) || !nrow(ev) || is.null(ev$key)) {
        sel <- df[0, ]
      } else {
        keys <- unique(ev$key)        # character or numeric
        # keys come from plotted data; row_id was assigned in df
        sel <- df[df$row_id %in% as.integer(keys), , drop = FALSE]
      }

      DT::datatable(
        sel |>
          dplyr::transmute(
            gene = display_gene,
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





    # Downloads ---------------------------------------------------------------

    # Filenames
    make_base_name <- function(suffix) {
      ds <- input$dataset %||% "dataset"
      paste0(
        .norm_filename(ds), "_volcano_F", input$lfc_thresh,
        "_", (deg_tbl()$.pcol[1]), input$padj_thresh, "_", suffix
      )
    }

    # SVG
    output$dl_plot_svg <- shiny::downloadHandler(
      filename = function() paste0(make_base_name("plot"), ".svg"),
      content = function(file) {
        df <- deg_tbl()
        gp <- make_gg(
          df,
          show_ns   = isTRUE(input$show_ns),
          label_sig = isTRUE(input$label_sig),
          label_topn = input$label_topn %||% 20
        )
        svglite::svglite(file, width = 8, height = 6)
        on.exit(grDevices::dev.off(), add = TRUE)
        print(gp)
      }
    )

    # PNG
    output$dl_plot_png <- shiny::downloadHandler(
      filename = function() paste0(make_base_name("plot"), ".png"),
      content = function(file) {
        df <- deg_tbl()
        gp <- make_gg(
          df,
          show_ns   = isTRUE(input$show_ns),
          label_sig = isTRUE(input$label_sig),
          label_topn = input$label_topn %||% 20
        )
        ggplot2::ggsave(
          filename = file, plot = gp, width = 8, height = 6, dpi = 300
        )
      }
    )

  })
}
