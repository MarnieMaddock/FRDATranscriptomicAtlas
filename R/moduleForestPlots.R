#' @importFrom graphics par
# UI ----
forestPlotsUI <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(
        width = 5,
        # in forestPlotsUI()
        textInput(ns("gene_text"),
                  "Gene (Symbol):",
                  value = "FXN",
                  placeholder = "e.g., FXN or ENSG00000165060"),
        actionButton(ns("plot_gene"), "Plot", class = "btn-primary")
      ),
      column(
        width = 3,
        br(),
        downloadButton(ns("download_svg"), "Download SVG"),
        downloadButton(ns("download_tiff"), "Download TIFF")
      ),
      tags$br(),
      tags$div(
        class = "alert alert-info",
        style = "margin-top:.5rem; margin-bottom:.5rem;",
        HTML("
              <strong>Interpretation</strong><br/>
              Squares are study estimates, sized by weight; horizontal bars are 95% CIs. The diamond is the random-effects summary (95% CI).<br/><br/>
              <strong>Direction</strong><br/>
              • <em>log2FC &gt; 0</em> → higher expression in FRDA.<br/>
              • <em>log2FC &lt; 0</em> → lower expression in FRDA.")
      ),
      tags$br(),
      shinycssloaders::withSpinner(
        verbatimTextOutput(ns("model_stats")),
        type = 4, color = "#005249", size = 0.6
      )
    )
  )
}


forestPlotMainUI <- function(id) {
  ns <- NS(id)
  # Spinner shows while renderPlot is busy
  shinycssloaders::withSpinner(
    plotOutput(ns("forest"), width = "1500px", height = "680px"),
    type = 4,
    color = "#005249",  # match theme
    size = 0.9,         # relative size
    proxy.height = "680px"  # reserve space to prevent layout jump
  )
}


# SERVER ----
forestPlotsServer <- function(id, pkg = utils::packageName(), data_dir = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # ---------- Helpers ----------
    # --- symbol map (tx2gene) ---
    tx2_path <- system.file("extdata/maps/tx2gene.tsv", package = pkg, mustWork = FALSE)
    gene_map <- if (nzchar(tx2_path) && file.exists(tx2_path)) {
      readr::read_tsv(tx2_path, col_types = "ccc") |>
        dplyr::distinct(gene_id, gene_name) |>
        dplyr::rename(external_gene_name = gene_name) |>
        dplyr::mutate(gene_id = sub("\\.\\d+$","", gene_id))
    } else NULL


    # Convert filename to a compact study key, e.g. "DESEQ2_res_Lees_FA2_0.05_all_genes.rds" -> "Lees_FA2"
    .study_key_from_path <- function(p) {
      nm <- basename(p)
      nm <- sub("\\.rds$", "", nm)
      nm <- sub("^DESEQ2_res_", "", nm)
      nm <- sub("_0[.]?0*5?_all_genes$", "", nm)  # tolerant of "0.05_all_genes"
      nm
    }

    # Choose load path
    data_path <- if (is.null(data_dir)) {
      system.file("extdata", "deg", "genes", package = pkg)
    } else {
      data_dir
    }
    validate(need(dir.exists(data_path), paste0("Data folder not found: ", data_path)))

    # List RDS files
    rds_files <- list.files(data_path, pattern = "\\.rds$", full.names = TRUE)
    validate(need(length(rds_files) > 0, paste0("No .rds files found in ", data_path)))

    # Read and harmonize each dataset into a standard shape
    read_one <- function(f) {
      df <- readRDS(f)

      study <- .study_key_from_path(f)
      # rename effect/se columns with study suffix
      ren <- setNames(
        paste0(c("log2FoldChange","lfcSE","padj")[c("log2FoldChange","lfcSE","padj") %in% names(df)],
               "_", study),
        c("log2FoldChange","lfcSE","padj")[c("log2FoldChange","lfcSE","padj") %in% names(df)]
      )
      dplyr::rename(df, !!!ren)
    }

    # helpers
    resolve_gene <- function(q, df) {
      q <- trimws(q)
      if (!nzchar(q)) return(integer(0))
      # Ensembl?
      if (grepl("^ENSG", q)) {
        ensg <- sub("\\.\\d+$","", q)
        return(which(df$gene_id == ensg))
      }
      # Otherwise treat as symbol
      if ("external_gene_name" %in% names(df)) {
        return(which(df$external_gene_name == q))
      }
      integer(0)
    }

    gene_query <- eventReactive(input$plot_gene, {
      input$gene_text
    }, ignoreInit = TRUE)

    # meta_df() uses the resolver
    meta_df <- reactive({
      df <- merged_all()
      q  <- gene_query()
      idx <- resolve_gene(q, df)
      validate(need(length(idx) >= 1, paste0("Gene not found: '", q, "'. ",
                                             if (!"external_gene_name" %in% names(df))
                                               "No symbol map loaded; try an Ensembl ID."
                                             else "")))
      row_idx <- idx[1]

      lfc_cols <- grep("^log2FoldChange_", names(df), value = TRUE)
      se_cols  <- grep("^lfcSE_",          names(df), value = TRUE)
      studies  <- intersect(sub("^log2FoldChange_", "", lfc_cols),
                            sub("^lfcSE_",          "", se_cols))
      validate(need(length(studies) > 1, "Need at least two studies with LFC and SE."))

      tibble::tibble(
        study = studies,
        yi    = as.numeric(df[row_idx, paste0("log2FoldChange_", studies), drop = TRUE]),
        sei   = as.numeric(df[row_idx, paste0("lfcSE_",          studies), drop = TRUE]),
        padj  = {
          pc <- paste0("padj_", studies)
          have <- pc %in% names(df)
          v <- rep(NA_real_, length(studies))
          if (any(have)) v[have] <- as.numeric(df[row_idx, pc[have], drop = TRUE])
          v
        }
      ) |>
        dplyr::filter(is.finite(yi), is.finite(sei)) |>
        dplyr::arrange(tolower(study))
    })


    # Merge all datasets by both keys when available
    merged_all <- reactive({

      rds_paths <- list.files(data_path, pattern = "all_genes\\.rds$", full.names = TRUE)
      validate(need(length(rds_paths) > 0, paste0("No *all_genes.rds files in ", data_path)))

      # 1) read all (log class) – only drop NULLs
      read_ok <- list()
      for (f in rds_paths) {
        nm <- sub("\\.rds$", "", basename(f))
        obj <- tryCatch(readRDS(f), error = function(e) {
          message("[read] ERROR: ", nm, " -> ", e$message); NULL
        })
        if (!is.null(obj)) {
          read_ok[[nm]] <- obj
        }
      }
      validate(need(length(read_ok) > 0, "No readable DEG RDS files."))

      # 2) map source object name -> study key
      study_key <- function(objname) {
        x <- objname
        x <- sub("^DESEQ2_res_", "", x)
        x <- sub("^DESEQ2_results_", "", x)   # Erwin naming
        x <- sub("_0[.]?0*5?_all_genes$", "", x)
        x
      }

      # 3) prep per-study df (coerce to data.frame here)
      prep_df <- function(df, study_name) {
        df <- as.data.frame(df)  # <— critical for DESeqResults

        # ensure gene_id
        if (!"gene_id" %in% names(df)) {
          if ("gene" %in% names(df)) {
            df <- dplyr::rename(df, gene_id = gene)
          } else if (!is.null(rownames(df))) {
            df <- tibble::rownames_to_column(df, var = "gene_id")
          } else {
            stop("'", study_name, "' has no gene_id/gene column or rownames.")
          }
        }
        df$gene_id <- sub("\\.\\d+$", "", df$gene_id)

        # keep typical DESeq2 columns (if present)
        keep <- intersect(c("log2FoldChange","lfcSE","padj","external_gene_name"), names(df))
        if (!all(c("log2FoldChange","lfcSE") %in% keep)) {
          stop("'", study_name, "' missing log2FoldChange or lfcSE.")
        }

        out <- dplyr::select(df, gene_id, dplyr::all_of(keep))

        # suffix non-id, non-symbol with study
        non_id <- setdiff(names(out), c("gene_id","external_gene_name"))
        names(out)[match(non_id, names(out))] <- paste0(non_id, "_", study_name)

        out
      }

      # 4) build named list of prepared dfs
      all_dfs <- list()
      for (nm in names(read_ok)) {
        st <- study_key(nm)
        all_dfs[[st]] <- tryCatch(
          prep_df(read_ok[[nm]], st),
          error = function(e) { message("[prep] SKIP ", nm, " (", st, "): ", e$message); NULL }
        )
      }
      all_dfs <- Filter(Negate(is.null), all_dfs)
      validate(need(length(all_dfs) > 0, "No usable DEG tables after prep."))

      # 5) reduce join
      merged <- purrr::reduce(all_dfs, dplyr::full_join, by = "gene_id")


      # log which studies have BOTH stats
      nms <- names(merged)
      studies <- intersect(sub("^log2FoldChange_", "", grep("^log2FoldChange_", nms, value=TRUE)),
                           sub("^lfcSE_",          "", grep("^lfcSE_",          nms, value=TRUE)))

      # attach symbols if available
      if (!is.null(gene_map)) {
        merged <- merged |>
          dplyr::left_join(gene_map, by = "gene_id")
      }


      # coalesce any stray symbol cols from individual tables
      sym_cols <- grep("^external_gene_name($|[_\\.])", names(merged), value = TRUE)
      if (length(sym_cols) > 1) {
        merged$external_gene_name <- dplyr::coalesce(!!!rlang::syms(sym_cols))
        merged <- dplyr::select(merged, -dplyr::all_of(setdiff(sym_cols, "external_gene_name")))
      }

      merged
    })




    # Map pretty labels for studies (fallback to the raw key if not listed)
    study_labels <- c(
      Erwin           = "Erwin GS (2017) - LCL",
      Indelicato      = "Indelicato E (2023) - SM",
      Lai_CNS         = "Lai J (2018) - CNS",
      Lai_PNS         = "Lai J (2018) - PNS",
      Lai_iPSC        = "Lai J (2018) - iPSC",
      Lees_FA1        = "Lees JG (2025*) - CM FA1",
      Lees_FA2        = "Lees JG (2025*) - CM FA2",
      Lees_FA3        = "Lees JG (2025*) - CM FA3",
      Maddock_LMN_FA2 = "Maddock ML (*) - LMN FA2",
      Maddock_NCC_FA1 = "Maddock ML (*) - NCC FA1",
      Maddock_NCC_FA2 = "Maddock ML (*) - NCC FA2",
      Maddock_SN_FA1  = "Maddock ML (*) - SN FA1",
      Maddock_SN_FA2  = "Maddock ML (*) - SN FA2",
      Mishra_223      = "Mishra P (2024) - N 223",
      Mishra_850      = "Mishra P (2024) - N 850",
      Mishra_FF1      = "Mishra P (2024) - N FF1",
      Mishra_FF2      = "Mishra P (2024) - N FF2",
      Napierala       = "Napierala JS (2017) - FB",
      Vilema          = "Vilema-Enríquez G (2020) - FB"
    )

    # Parse user's selection into either symbol or Ensembl ID
    parse_query <- function(q) {
      if (is.null(q) || !nzchar(q)) return(list(symbol = NULL, ensg = NULL))
      # If user picked "SYMBOL (ENSG…)" split it
      m <- regmatches(q, regexec("^(.*?)\\s*\\((ENSG\\w+)\\)$", q))[[1]]
      if (length(m) == 3) {
        list(symbol = if (m[2] == m[3]) NULL else m[2], ensg = m[3])
      } else if (grepl("^ENSG", q)) {
        list(symbol = NULL, ensg = sub("\\.\\d+$","", q))
      } else {
        list(symbol = q, ensg = NULL)
      }
    }


    # Fit models + stats
    models <- reactive({
      md <- meta_df()
      validate(need(nrow(md) >= 2, "Need at least two studies with LFC and SE to run meta-analysis."))
      res_re <- metafor::rma(yi = md$yi, sei = md$sei, data = md, method = "REML")
      res_fe <- metafor::rma(yi = md$yi, sei = md$sei, data = md, method = "FE")
      list(md = md, re = res_re, fe = res_fe)
    })

    # Text stats
    output$model_stats <- renderText({
      m <- models()
      re <- m$re
      df <- re$k - re$p
      I2 <- max(0, (re$QE - df) / re$QE) * 100
      paste0(
        "Random-effects (REML)\n",
        sprintf("Summary LFC = %.3f  [%.3f, %.3f]\n", re$b[1], re$ci.lb, re$ci.ub),
        sprintf("Tau^2 = %.4f; Q(df=%d) = %.2f, p = %.4g; I^2 = %.1f%%",
                re$tau2, df, re$QE, re$QEp, I2)
      )
    })

    # Forest plot
    draw_forest <- function(md, re, svg_path = NULL) {
      # compute wide-enough x limits from data + model (no clipping)
      # use ±2.5*SE for study CIs to be safe
      ci_lo <- md$yi - 2.5 * md$sei
      ci_hi <- md$yi + 2.5 * md$sei
      xmin  <- min(ci_lo, re$ci.lb, na.rm = TRUE)
      xmax  <- max(ci_hi, re$ci.ub, na.rm = TRUE)
      pad   <- max(0.2, 0.4 * (xmax - xmin))  # small padding
      alim  <- c(xmin - pad, xmax + pad)

      # a bit of extra room for labels on the right
      xlim  <- c(alim[1] - 0.5*(alim[2]-alim[1])*0.00,  # adjust if you add ilab later
                 alim[2] + 0.10*(alim[2]-alim[1]))

      # ticks
      at_ticks <- pretty(alim, n = 6)

      # I^2 label
      df  <- re$k - re$p
      I2  <- max(0, (re$QE - df) / re$QE) * 100
      mlab_txt <- sprintf("Summary (I^2 = %.1f%%)", I2)

      # open device if saving
      if (!is.null(svg_path)) {
        svglite::svglite(svg_path, width = 7.5, height = 6)  # taller, nicer aspect
        on.exit(grDevices::dev.off(), add = TRUE)
      }

      oldpar <- par(no.readonly = TRUE); on.exit(par(oldpar), add = TRUE)

      # margins: bottom, left (labels), top, right (CI text)
      par(mar = c(4.2, 8.5, 2, 6))

      # pretty study labels if you have your map
      slabs <- if (exists("study_labels", inherits = TRUE)) {
        unname(ifelse(md$study %in% names(study_labels), study_labels[md$study], md$study))
      } else md$study

      metafor::forest(
        re,
        slab    = slabs,
        xlab    = "log2 fold change",
        alim    = alim,              # <- dynamic, no clipping, so no arrows
        xlim    = xlim,              # <- room for label column on right
        at      = at_ticks,
        refline = 0,
        cex     = 1.5,              # slightly larger text
        # psize = 1,                  # uncomment for uniform box sizes (otherwise weight-scaled)
        mlab    = mlab_txt
      )
    }


    output$forest <- renderPlot({
      m <- models()
      draw_forest(m$md, m$re, svg_path = NULL)
    },
    width  = 1500,  # pixels
    height = 680   # pixels
    )

    px_to_in <- function(px, dpi = 96) px / dpi

    # Download handler (SVG)
    output$download_svg <- downloadHandler(
      filename = function() {
        g <- gene_query()
        if (is.null(g) || !nzchar(g)) g <- "gene"
        paste0(gsub("[^A-Za-z0-9_\\-]+", "_", g), "_forest.svg")
      },
      content = function(file) {
        m <- models()
        draw_forest(m$md, m$re, svg_path = file)
      }
    )

    safe_name <- function(x) gsub("[^A-Za-z0-9_\\-]+", "_", x)

    output$download_tiff <- downloadHandler(
      filename = function() paste0(safe_name(gene_query()), "_forest.tiff"),
      contentType = "image/tiff",
      content = function(file) {
        m <- models()
        # open a TIFF device (match UI size or pick print size)
        grDevices::tiff(
          filename = file,
          width = 1500, height = 680, units = "px",  # or: width=8,height=6, units="in", res=300
          res = 96,
          compression = "lzw",
          pointsize = 11
        )
        on.exit(grDevices::dev.off(), add = TRUE)
        draw_forest(m$md, m$re)
      }
    )


  })
}
