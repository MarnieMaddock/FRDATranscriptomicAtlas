#moduleGSEA.R
# R/modules/goGSEA_module.R
# -------------------------

# Sidebar controls
GSEASidebarUI <- function(id) {
  ns <- shiny::NS(id)
  tagList(
    selectInput(ns("dataset"), "Dataset", choices = NULL),
    radioButtons(ns("ont"), "Ontology",
                 choices = c("BP", "CC", "MF"), inline = TRUE),
    sliderInput(ns("ncat"), "Show top categories",
                min = 2, max = 500, value = 10, step = 5),
    tags$hr(),
    br(),
  )
}
GSEAMainUI <- function(id) {
  ns <- shiny::NS(id)
  tagList(
    h4("Dotplot (GeneRatio)"),
    uiOutput(ns("plot_gr_ui")),
    div(class = "text-center mt-2",
        downloadButton(ns("dl_plot_gr_png"), "Download PNG", class = "btn-sm"),
        downloadButton(ns("dl_plot_gr_svg"), "Download SVG", class = "btn-sm")
    ),
    hr(),
    h4("Dotplot (NES)"),
    uiOutput(ns("plot_nes_ui")),
    div(class = "text-center mt-2",
        downloadButton(ns("dl_plot_nes_png"), "Download PNG", class = "btn-sm"),
        downloadButton(ns("dl_plot_nes_svg"), "Download SVG", class = "btn-sm")
    ),
    hr(),
    h4("All Results (FDR < 0.05)"),
    shinycssloaders::withSpinner(
      DT::DTOutput(ns("tbl")),
      type = 4, color = "#005249"
    ),
    div(class = "text-center mt-2",
        downloadButton(ns("dl_table_csv"), "Download table (CSV)", class = "btn-sm"),
        div(class = "text-muted small",
            "Tip: click a row in the table, then use “Show genes” or “Download genes (CSV)”."),

        actionButton(ns("show_genes"), "Show genes for selected term", class = "btn-sm"),
        downloadButton(ns("dl_genes_csv"), "Download genes (CSV)", class = "btn-sm")
    ),
    br()
  )
}


# R/modules/goGSEA_module.R  (replace just the server)
# R/modules/goGSEA_module.R
# -------------------------
# Server for GSEA (CSV/RDS-on-disk, on-demand)
GSEAServer <- function(id, base_dir = NULL, pkg = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # --- maps (tx2gene) ---
    tx2_path <- system.file("extdata/maps/tx2gene.tsv", package = pkg, mustWork = FALSE)
    if (!nzchar(tx2_path)) tx2_path <- file.path("inst", "extdata", "maps", "tx2gene.tsv")

    tx2gene <- if (nzchar(tx2_path) && file.exists(tx2_path)) {
      readr::read_tsv(tx2_path, col_types = "ccc") |>
        dplyr::mutate(
          transcript_id = sub("\\.\\d+$","", transcript_id),
          gene_id       = sub("\\.\\d+$","", gene_id)
        )
    } else NULL

    norm_id <- function(x) sub("\\.\\d+$","", as.character(x))

    # ENSEMBL -> SYMBOL lookup (vector)
    ensg_to_symbol <- if (!is.null(tx2gene) && all(c("gene_id","gene_name") %in% names(tx2gene))) {
      ids  <- norm_id(tx2gene$gene_id)
      syms <- as.character(tx2gene$gene_name)
      keep <- !duplicated(ids)
      stats::setNames(syms[keep], ids[keep])
    } else stats::setNames(character(0), character(0))

    # Vectorized mapper (returns symbol; if missing, returns cleaned id)
    map_ids_to_symbol <- function(ids) {
      if (!length(ids)) return(character())
      ids2 <- norm_id(ids)
      hit  <- ensg_to_symbol[ids2]
      out  <- ifelse(is.na(hit) | !nzchar(hit), ids2, unname(hit))
      out
    }


    # ---- make `pkg` safe here too ----
    if (is.null(pkg) || !is.character(pkg) || !nzchar(pkg)) {
      pkg <- tryCatch(utils::packageName(), error = function(e) "")
      if (!length(pkg) || !nzchar(pkg)) pkg <- "FRDATranscriptomicAtlas"
      pkg <- pkg[[1L]]  # ensure length 1
    }

    # optional: guard for pretty_map if it isn't globally defined
    if (!exists("pretty_map", inherits = TRUE)) {
      pretty_map <- setNames(character(0), character(0))
    }

    # ---- deps ----
    requireNamespace("DT", quietly = TRUE)
    requireNamespace("ggplot2", quietly = TRUE)
    requireNamespace("svglite", quietly = TRUE)
    requireNamespace("enrichplot", quietly = TRUE)

    `%||%` <- function(x,y) if (is.null(x) || (is.character(x) && !nzchar(x))) y else x

    # --- helper: extract core genes vector from a one-row data.frame
    .extract_genes_raw <- function(rowdf) {
      if (nrow(rowdf) < 1) return(character())
      col <- if ("core_enrichment" %in% names(rowdf)) "core_enrichment"
      else if ("geneID" %in% names(rowdf)) "geneID" else return(character())
      unique(norm_id(strsplit(as.character(rowdf[[col]][[1]]), "/", fixed = TRUE)[[1]]))
    }

    .add_core_cols <- function(df) {
      if (!nrow(df)) return(df)
      col <- if ("core_enrichment" %in% names(df)) "core_enrichment"
      else if ("geneID" %in% names(df)) "geneID" else return(df)

      sp <- strsplit(as.character(df[[col]]), "/", fixed = TRUE)
      genes_list   <- lapply(sp, function(v) unique(norm_id(v[nzchar(v)])))
      symbols_list <- lapply(genes_list, map_ids_to_symbol)

      df$core_count   <- vapply(genes_list, length, integer(1))
      df$core_preview <- vapply(symbols_list, function(v) {
        if (!length(v)) return("")
        paste0(paste(utils::head(v, 8), collapse = ", "),
               if (length(v) > 8) sprintf(" … (+%d)", length(v) - 8) else "")
      }, character(1))
      df
    }

    # ---- locate extdata/GSEA_results no matter how the app is run ----
    if (is.null(base_dir) || !nzchar(base_dir)) {
      # installed package path
      base_dir <- system.file("extdata", "GSEA_results", package = pkg, mustWork = FALSE)
      # project (dev) fallback
      if (!nzchar(base_dir) || !dir.exists(base_dir)) {
        base_dir <- file.path("inst", "extdata", "GSEA_results")
      }
    }

    # ---- index: dataset folder -> {BP,CC,MF} -> file path ----
    build_index <- function(root) {
      if (!dir.exists(root)) return(list())
      ds_dirs <- list.dirs(root, recursive = FALSE, full.names = TRUE)
      idx <- lapply(ds_dirs, function(d) {
        # prefer .rds; fallback to .csv (if you kept CSVs)
        map <- setNames(vector("list", 3), c("BP","CC","MF"))
        for (ont in c("BP","CC","MF")) {
          p_rds <- file.path(d, sprintf("gsea_GO_%s.rds", ont))
          p_csv <- file.path(d, sprintf("gsea_GO_%s.csv", ont))
          if (file.exists(p_rds)) {
            map[[ont]] <- p_rds
          } else if (file.exists(p_csv)) {
            map[[ont]] <- p_csv
          } else {
            map[[ont]] <- ""
          }
        }
        map
      })
      names(idx) <- basename(ds_dirs)
      # keep datasets that have at least one ontology file
      idx[vapply(idx, function(x) any(nzchar(unlist(x))), logical(1))]
    }

    idx <- build_index(base_dir)
    datasets_vec <- names(idx)

    # ---- pretty labels for the dataset select ----
    # expects a global/internal object `pretty_map` (named character vector)
    # ---- pretty labels for the dataset select (dataset_key3) ----
    # keep up to first 3 tokens (e.g., "Lees_FA1", "Maddock_SN_FA2")
    dataset_key3 <- function(x) {
      b <- basename(x)
      b <- sub("\\.rds$", "", b)
      b <- sub("_batchcorrection", "", b, ignore.case = TRUE)
      b <- sub("(_0\\.[0-9].*)$", "", b)  # strip trailing "_0.05_all_genes..." if present
      parts <- strsplit(b, "_", fixed = TRUE)[[1]]
      if (length(parts) >= 3) paste(parts[1:3], collapse = "_")
      else if (length(parts) >= 2) paste(parts[1:2], collapse = "_")
      else parts[1]
    }

    labels <- vapply(datasets_vec, function(d) {
      key <- dataset_key3(d)
      pretty_map[[key]] %||% gsub("_", " ", key)   # fallback: "Lees FA1"
    }, character(1))

    choices_named <- stats::setNames(datasets_vec, labels)

    label_from_value <- stats::setNames(labels, datasets_vec)

    # ---- reactive: are we safe to render? ----
    ready <- reactive({
      ds_ok  <- !is.null(input$dataset) && input$dataset %in% datasets_vec
      ont_ok <- !is.null(input$ont) && input$ont %in% c("BP","CC","MF")
      file_ok <- if (ds_ok && ont_ok) nzchar(idx[[input$dataset]][[input$ont]]) else FALSE
      ds_ok && ont_ok && file_ok
    })

    # ---- populate controls ----
    observe({
      validate(need(length(datasets_vec) > 0, "No GSEA files found under extdata/GSEA_results."))
      updateSelectInput(session, "dataset", choices = choices_named, selected = datasets_vec[[1]])
    })
    observeEvent(input$dataset, {
      updateRadioButtons(session, "ont", choices = c("BP","CC","MF"), selected = "BP")
    }, ignoreInit = TRUE)

    # ---- current file path & object ----
    r_path <- reactive({
      req(ready())
      idx[[input$dataset]][[input$ont]]
    })

    read_any <- function(p) {
      if (grepl("\\.rds$", p, ignore.case = TRUE)) {
        readRDS(p)
      } else {
        # CSV fallback support (returns a df, not a gseaResult). We convert.
        df <- utils::read.csv(p, check.names = FALSE)
        df
      }
    }

    r_obj <- reactive({
      req(ready())
      read_any(r_path())
    })

    # ---- title ----
    title_txt <- reactive({
      ds_label <- label_from_value[[input$dataset %||% ""]] %||% (input$dataset %||% "")
      paste0(ds_label, " - GO ", input$ont %||% "", " GSEA")
    })
    # ---- df accessor (works for gseaResult or data.frame fallback) ----
    r_df <- reactive({
      x <- r_obj()
      df <- if (inherits(x, "gseaResult")) as.data.frame(x) else x
      if ("p.adjust" %in% names(df)) {
        df <- df[order(df$p.adjust, -abs(df$NES %||% 0)), , drop = FALSE]
      }
      .add_core_cols(df)   # uses symbols now
    })



    # ---- plot helpers ----
    # ---- plot helpers ----
    make_dotplot <- function(metric = c("GeneRatio","NES")) {
      metric <- match.arg(metric)
      x <- r_obj()

      # ----------------------------
      # Case 1: object is a gseaResult from clusterProfiler
      # ----------------------------
      if (inherits(x, "gseaResult")) {

        # enrichplot NOT installed → fallback message plot
        if (!requireNamespace("enrichplot", quietly = TRUE)) {
          return(
            ggplot2::ggplot() +
              ggplot2::annotate(
                "text",
                x = 0.5, y = 0.5,
                hjust = 0.5,
                label = "Dotplot requires the Bioconductor package 'enrichplot',\nwhich is not available for R ≥ 4.5.\nShowing message instead.",
                size = 5
              ) +
              ggplot2::theme_void()
          )
        }

        # enrichplot available → use standard dotplot
        return(
          enrichplot::dotplot(
            x,
            x = metric,
            showCategory = input$ncat %||% 10
          ) +
            ggplot2::ggtitle(title_txt())
        )
      }

      # ----------------------------
      # Case 2: fallback for data.frame CSV results
      # ----------------------------
      df <- r_df()
      n  <- max(1, as.integer(input$ncat %||% 10))
      d  <- head(df, n)

      # Compute GeneRatio if missing
      if (!"GeneRatio" %in% names(d)) {
        if ("core_enrichment" %in% names(d) && "setSize" %in% names(d)) {
          core_n <- vapply(
            strsplit(as.character(d$core_enrichment), "/", fixed = TRUE),
            function(v) sum(nzchar(v)),
            integer(1)
          )
          d$GeneRatio <- core_n / as.numeric(d$setSize)
        }
      }

      metric_col <- switch(metric,
                           "GeneRatio" = "GeneRatio",
                           "NES"       = "NES")

      # When metric column is missing → safe fallback
      if (!metric_col %in% names(d)) {
        return(
          ggplot2::ggplot() +
            ggplot2::annotate(
              "text",
              x = 0.5, y = 0.5,
              hjust = 0.5,
              label = sprintf("Metric '%s' not available in results.", metric),
              size = 5
            ) +
            ggplot2::theme_void()
        )
      }

      # ----------------------------
      # Clean fallback dotplot (ggplot2 only)
      # ----------------------------
      ggplot2::ggplot(
        d,
        ggplot2::aes(
          x = reorder(Description, !!rlang::sym(metric_col)),
          y = !!rlang::sym(metric_col)
        )
      ) +
        ggplot2::geom_point(size = 3, color = "#3366AA") +
        ggplot2::coord_flip() +
        ggplot2::labs(
          title = title_txt(),
          x     = "Term",
          y     = metric
        ) +
        ggplot2::theme_bw(base_size = 12)
    }

    # make_dotplot <- function(metric = c("GeneRatio","NES")) {
    #   metric <- match.arg(metric)
    #   x <- r_obj()
    #   if (inherits(x, "gseaResult")) {
    #     enrichplot::dotplot(x, x = metric, showCategory = input$ncat %||% 10) +
    #       ggplot2::ggtitle(title_txt())
    #   } else {
    #     # data.frame fallback (if you kept only CSVs) — manual dotplot
    #     df <- r_df()
    #     n  <- max(1, as.integer(input$ncat %||% 10))
    #     d  <- head(df, n)
    #     if (!"GeneRatio" %in% names(d)) {
    #       # compute GeneRatio from core_enrichment/setSize if present
    #       if ("core_enrichment" %in% names(d) && "setSize" %in% names(d)) {
    #         core_n <- vapply(strsplit(as.character(d$core_enrichment), "/", fixed = TRUE),
    #                          function(v) sum(nzchar(v)), integer(1))
    #         d$GeneRatio <- core_n / as.numeric(d$setSize)
    #       }
    #     }
    #     aes_x <- if (metric == "NES") ggplot2::aes(x = NES) else ggplot2::aes(x = GeneRatio)
    #     ggplot2::ggplot(d, aes_x + ggplot2::aes(y = stats::reorder(Description, !!rlang::sym(names(d)[match(metric, c("GeneRatio","NES"))])))) +
    #       ggplot2::geom_point(ggplot2::aes(size = setSize, alpha = -log10(p.adjust))) +
    #       { if (metric == "NES") ggplot2::geom_vline(xintercept = 0, linetype = 2) else NULL } +
    #       ggplot2::labs(x = if (metric == "NES") "Normalized Enrichment Score (NES)" else "GeneRatio (core / set size)",
    #                     y = NULL, title = title_txt(), alpha = "-log10(adj.P)") +
    #       ggplot2::theme_minimal(base_size = 12)
    #   }
    # }

    # helper to compute how many rows will be drawn (guard against short tables)
    n_rows_to_plot <- reactive({
      req(ready())
      min(as.integer(input$ncat), nrow(r_df()))
    })

    row_px   <- 50   # pixels per category row (tweak 28–34 as you like)
    marginpx <- 130  # fixed padding for title/axes/legend
    minpx    <- 380  # minimum so tiny plots still look decent

    dyn_height <- reactive({
      h <- row_px * n_rows_to_plot() + marginpx
      max(h, minpx)
    })


    # dynamic left margin based on longest label (keeps full, unwrapped text)
    left_margin_px <- reactive({
      req(ready())
      d <- head(r_df()[, "Description", drop = TRUE], n_rows_to_plot())
      # ~6 px per character is a reasonable default for 12–13 pt fonts
      px <- 6 * max(nchar(d), na.rm = TRUE)
      px <- max(140, min(px, 320))  # clamp to sensible bounds
      px
    })

    # dyn_height() already defined in your server
    output$plot_gr_ui <- renderUI({
      req(ready())
      shinycssloaders::withSpinner(
        plotOutput(ns("plot_gr"), height = dyn_height()),   # height set in UI
        type = 4, color = "#005249"
      )
    })

    output$plot_nes_ui <- renderUI({
      req(ready())
      shinycssloaders::withSpinner(
        plotOutput(ns("plot_nes"), height = dyn_height()),   # height set in UI
        type = 4, color = "#005249"
      )
    })

    # keep your existing renderPlot() calls:
    # output$plot_gr  <- renderPlot({ ... }, res = 96)
    # output$plot_nes <- renderPlot({ ... }, res = 96)

    output$plot_gr <- renderPlot(
      { req(ready()); make_dotplot("GeneRatio") },
      res    = 96,
      height = function() dyn_height()
    )


    output$plot_nes <- renderPlot(
      { req(ready()); make_dotplot("NES") },
      res    = 96,
      height = function() dyn_height()
    )


    output$tbl <- DT::renderDT({
      req(ready())
      df <- r_df()
      # choose a sensible subset/order if you like:
      keep <- intersect(c("ID","Description","NES","p.adjust","setSize","core_count","core_preview"),
                        names(df))
      DT::datatable(
        df[, keep, drop = FALSE],
        rownames = FALSE,
        filter   = "top",
        selection = "single",
        options = list(pageLength = 10, scrollX = TRUE,
                       columnDefs = list(
                         list(targets = which(colnames(df[, keep]) == "core_preview") - 1L,
                              render = DT::JS(
                                "function(data,type,row,meta){",
                                " if(type==='display' && data && data.length>120){",
                                "   return '<span title=\"'+data+'\">'+data.slice(0,120)+'…</span>';",
                                " } return data; }"
                              ))
                       ))
      )
    })
    # ---- show genes modal ----
    observeEvent(input$show_genes, {
      sel <- input$tbl_rows_selected
      validate(need(length(sel) == 1, "Select one term in the table first."))
      df  <- r_df()
      row <- df[sel, , drop = FALSE]

      # get IDs, then map to symbols for display
      ids   <- .extract_genes_raw(row)
      genes <- map_ids_to_symbol(ids)   # pretty names in the modal

      term_title <- if ("Description" %in% names(row)) row$Description[[1]] else "Selected term"

      shiny::showModal(
        modalDialog(
          title = term_title,
          size = "l",
          easyClose = TRUE,
          footer = modalButton("Close"),
          tagList(
            p(sprintf("Genes in leading edge/core set: %d", length(genes))),
            tags$div(
              style = "max-height: 50vh; overflow:auto; font-family: monospace; white-space: pre-wrap;",
              paste(genes, collapse = ", ")
            )
          )
        )
      )
    })
    # --- download genes CSV ---
    output$dl_genes_csv <- downloadHandler(
      filename = function() {
        sel  <- input$tbl_rows_selected
        base <- if (length(sel) == 1 && "ID" %in% names(r_df()))
          paste0("genes_", r_df()$ID[[sel]]) else "genes_selected_term"
        paste0(gsub("[^A-Za-z0-9_]+", "_", base), ".csv")
      },
      contentType = "text/csv; charset=UTF-8",
      content = function(file) {
        sel <- input$tbl_rows_selected
        validate(need(length(sel) == 1, "Select one term in the table first."))

        df  <- r_df()
        row <- df[sel, , drop = FALSE]

        ids  <- .extract_genes_raw(row)      # cleaned ENSG IDs
        syms <- map_ids_to_symbol(ids)       # HGNC symbols (or ID fallback)

        out <- data.frame(
          term_id   = row$ID          %||% NA_character_,
          term_name = row$Description %||% NA_character_,
          NES       = row$NES         %||% NA_real_,
          padj      = row$p.adjust    %||% NA_real_,
          gene_id   = ids,
          symbol    = syms,
          stringsAsFactors = FALSE
        )

        if (requireNamespace("readr", quietly = TRUE)) {
          readr::write_excel_csv(out, file)
        } else {
          utils::write.csv(out, file, row.names = FALSE)
        }
      }
    )


    # ---- downloads ----
    output$dl_plot_gr_png <- downloadHandler(
      filename = function() sprintf("dotplot_GeneRatio_%s_%s.png", input$dataset, input$ont),
      content  = function(file) ggplot2::ggsave(file, make_dotplot("GeneRatio"), width = 7, height = 5, dpi = 300)
    )
    output$dl_plot_gr_svg <- downloadHandler(
      filename = function() sprintf("dotplot_GeneRatio_%s_%s.svg", input$dataset, input$ont),
      content  = function(file) { svglite::svglite(file, 7, 5); print(make_dotplot("GeneRatio")); dev.off() }
    )
    output$dl_plot_nes_png <- downloadHandler(
      filename = function() sprintf("dotplot_NES_%s_%s.png", input$dataset, input$ont),
      content  = function(file) ggplot2::ggsave(file, make_dotplot("NES"), width = 7, height = 5, dpi = 300)
    )
    output$dl_plot_nes_svg <- downloadHandler(
      filename = function() sprintf("dotplot_NES_%s_%s.svg", input$dataset, input$ont),
      content  = function(file) { svglite::svglite(file, 7, 5); print(make_dotplot("NES")); dev.off() }
    )
    output$dl_table_csv <- downloadHandler(
      filename = function() sprintf("gsea_table_%s_%s.csv", input$dataset, input$ont),
      content  = function(file) utils::write.csv(r_df(), file, row.names = FALSE)
    )
  })
}
