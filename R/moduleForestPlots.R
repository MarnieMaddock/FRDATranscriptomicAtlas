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


## SERVER ----
forestPlotsServer <- function(id, pkg = utils::packageName(), data_dir = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # --- make pkg safe (length-1 string) ---------------------------------
    pkg <- tryCatch(pkg, error = function(e) "")
    if (!length(pkg) || !is.character(pkg) || !nzchar(pkg)) pkg <- "FRDATranscriptomicAtlas"
    pkg <- pkg[[1L]]

    # --- helpers to resolve package-or-project paths ----------------------
    resolve_dir <- function(subpath) {
      d <- system.file(subpath, package = pkg, mustWork = FALSE)
      if (!nzchar(d)) d <- file.path("inst", subpath)
      d
    }
    resolve_file <- function(subpath, fname) file.path(resolve_dir(subpath), fname)

    # ---------- symbol map (tx2gene) -------------------------------------
    tx2_path <- resolve_file(file.path("extdata", "maps"), "tx2gene.tsv")
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

    # ---------- choose DEG folder (genes) --------------------------------
    data_path <- if (!is.null(data_dir)) {
      data_dir
    } else {
      resolve_dir(file.path("extdata", "deg", "genes"))
    }
    validate(need(nzchar(data_path) && dir.exists(data_path),
                  paste0("Data folder not found: ", data_path)))

    # ---------- list RDS files -------------------------------------------
    rds_files <- list.files(data_path, pattern = "\\.rds$", full.names = TRUE)
    validate(need(length(rds_files) > 0,
                  paste0("No .rds files found in ", data_path)))

    # ---------- read + harmonize per study --------------------------------
    read_one <- function(f) {
      df <- readRDS(f)
      study <- .study_key_from_path(f)
      # rename effect/se columns with study suffix
      have <- c("log2FoldChange","lfcSE","padj")
      have <- have[have %in% names(df)]
      ren  <- setNames(paste0(have, "_", study), have)
      dplyr::rename(df, !!!ren)
    }

    # ---------- gene resolver ---------------------------------------------
    resolve_gene <- function(q, df) {
      q <- trimws(q)
      if (!nzchar(q)) return(integer(0))
      if (grepl("^ENSG", q)) {
        ensg <- sub("\\.\\d+$","", q)
        return(which(df$gene_id == ensg))
      }
      if ("external_gene_name" %in% names(df)) {
        return(which(df$external_gene_name == q))
      }
      integer(0)
    }

    gene_query <- eventReactive(input$plot_gene, input$gene_text, ignoreInit = TRUE)

    # ---------- merge all studies ----------------------------------------
    merged_all <- reactive({
      rds_paths <- list.files(data_path, pattern = "all_genes\\.rds$", full.names = TRUE)
      validate(need(length(rds_paths) > 0,
                    paste0("No *all_genes.rds files in ", data_path)))

      read_ok <- list()
      for (f in rds_paths) {
        nm <- sub("\\.rds$", "", basename(f))
        obj <- tryCatch(readRDS(f), error = function(e) NULL)
        if (!is.null(obj)) read_ok[[nm]] <- obj
      }
      validate(need(length(read_ok) > 0, "No readable DEG RDS files."))

      study_key <- function(objname) {
        x <- objname
        x <- sub("^DESEQ2_res_", "", x)
        x <- sub("^DESEQ2_results_", "", x)   # Erwin naming
        x <- sub("_0[.]?0*5?_all_genes$", "", x)
        x
      }

      prep_df <- function(df, study_name) {
        df <- as.data.frame(df)

        if (!"gene_id" %in% names(df)) {
          if ("gene" %in% names(df)) df <- dplyr::rename(df, gene_id = gene)
          else if (!is.null(rownames(df))) df <- tibble::rownames_to_column(df, var = "gene_id")
          else stop("'", study_name, "' has no gene_id/gene column or rownames.")
        }
        df$gene_id <- sub("\\.\\d+$", "", df$gene_id)

        keep <- intersect(c("log2FoldChange","lfcSE","padj","external_gene_name"), names(df))
        if (!all(c("log2FoldChange","lfcSE") %in% keep))
          stop("'", study_name, "' missing log2FoldChange or lfcSE.")

        out <- dplyr::select(df, gene_id, dplyr::all_of(keep))
        non_id <- setdiff(names(out), c("gene_id","external_gene_name"))
        names(out)[match(non_id, names(out))] <- paste0(non_id, "_", study_name)
        out
      }

      all_dfs <- list()
      for (nm in names(read_ok)) {
        st <- study_key(nm)
        all_dfs[[st]] <- tryCatch(prep_df(read_ok[[nm]], st),
                                  error = function(e) NULL)
      }
      all_dfs <- Filter(Negate(is.null), all_dfs)
      validate(need(length(all_dfs) > 0, "No usable DEG tables after prep."))

      merged <- purrr::reduce(all_dfs, dplyr::full_join, by = "gene_id")

      if (!is.null(gene_map)) {
        merged <- dplyr::left_join(merged, gene_map, by = "gene_id")
      }

      # coalesce stray symbol columns if present
      sym_cols <- grep("^external_gene_name($|[_\\.])", names(merged), value = TRUE)
      if (length(sym_cols) > 1) {
        merged$external_gene_name <- dplyr::coalesce(!!!rlang::syms(sym_cols))
        merged <- dplyr::select(merged, -dplyr::all_of(setdiff(sym_cols, "external_gene_name")))
      }

      merged
    })

    # ---------- build meta-data frame for metafor -------------------------
    meta_df <- reactive({
      df <- merged_all()
      q  <- gene_query()
      idx <- resolve_gene(q, df)
      validate(need(length(idx) >= 1, paste0(
        "Gene not found: '", q, "'. ",
        if (!"external_gene_name" %in% names(df)) "No symbol map loaded; try an Ensembl ID." else ""
      )))
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

    # ---------- pretty study labels (unchanged) ---------------------------
    study_labels <- c(
      Chutake="Chutake YK (2014) - LCL",
      Erwin="Erwin GS (2017) - LCL", Indelicato="Indelicato E (2023) - SM",
      Lai_CNS="Lai J (2018) - CNS", Lai_PNS="Lai J (2018) - PNS", Lai_iPSC="Lai J (2018) - iPSC",

      Lees_FA1="Lees JG (2025*) - CM FA1", Lees_FA2="Lees JG (2025*) - CM FA2", Lees_FA3="Lees JG (2025*) - CM FA3",
      Li_FRDA_CTRL="Li J (2019) - CM^", Li_FRDA_IC="Li J (2019) - CM#",
      Maddock_LMN_FA2="Maddock ML (*) - LMN FA2", Maddock_NCC_FA1="Maddock ML (*) - NCC FA1",
      Maddock_NCC_FA2="Maddock ML (*) - NCC FA2", Maddock_SN_FA1="Maddock ML (*) - SN FA1",
      Maddock_SN_FA2="Maddock ML (*) - SN FA2", Mishra_223="Mishra P (2024) - N 223",
      Mishra_850="Mishra P (2024) - N 850", Mishra_FF1="Mishra P (2024) - N FF1", Mishra_FF2="Mishra P (2024) - N FF2",
      Napierala="Napierala JS (2017) - FB", Vilema="Vilema-Enríquez G (2020) - FB", Wang="Wang F (2022) - FB"
    )

    # ---------- models + outputs (unchanged) ------------------------------
    models <- reactive({
      md <- meta_df()
      validate(need(nrow(md) >= 2, "Need at least two studies with LFC and SE to run meta-analysis."))
      res_re <- metafor::rma(yi = md$yi, sei = md$sei, data = md, method = "REML")
      res_fe <- metafor::rma(yi = md$yi, sei = md$sei, data = md, method = "FE")
      list(md = md, re = res_re, fe = res_fe)
    })

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

    draw_forest <- function(md, re, svg_path = NULL) {
      ci_lo <- md$yi - 2.5 * md$sei
      ci_hi <- md$yi + 2.5 * md$sei
      xmin  <- min(ci_lo, re$ci.lb, na.rm = TRUE)
      xmax  <- max(ci_hi, re$ci.ub, na.rm = TRUE)
      pad   <- max(0.2, 0.4 * (xmax - xmin))
      alim  <- c(xmin - pad, xmax + pad)
      xlim  <- c(alim[1], alim[2] + 0.10*(alim[2]-alim[1]))
      at_ticks <- pretty(alim, n = 6)

      if (!is.null(svg_path)) {
        svglite::svglite(svg_path, width = 7.5, height = 6)
        on.exit(grDevices::dev.off(), add = TRUE)
      }
      oldpar <- par(no.readonly = TRUE); on.exit(par(oldpar), add = TRUE)
      par(mar = c(4.2, 8.5, 2, 6))

      slabs <- if (exists("study_labels", inherits = TRUE)) {
        unname(ifelse(md$study %in% names(study_labels), study_labels[md$study], md$study))
      } else md$study

      metafor::forest(
        re, slab = slabs, xlab = "log2 fold change",
        alim = alim, xlim = xlim, at = at_ticks, refline = 0, cex = 1.5,
        mlab = sprintf("Summary (I^2 = %.1f%%)", {
          df <- re$k - re$p; max(0, (re$QE - df) / re$QE) * 100
        })
      )
    }

    output$forest <- renderPlot({
      m <- models(); draw_forest(m$md, m$re, svg_path = NULL)
    }, width = 1500, height = 680)

    safe_name <- function(x) gsub("[^A-Za-z0-9_\\-]+", "_", x)

    output$download_svg <- downloadHandler(
      filename = function() paste0(safe_name(gene_query() %||% "gene"), "_forest.svg"),
      content  = function(file) { m <- models(); draw_forest(m$md, m$re, svg_path = file) }
    )

    output$download_tiff <- downloadHandler(
      filename = function() paste0(safe_name(gene_query() %||% "gene"), "_forest.tiff"),
      contentType = "image/tiff",
      content  = function(file) {
        m <- models()
        grDevices::tiff(file, width = 1500, height = 680, units = "px", res = 96,
                        compression = "lzw", pointsize = 11)
        on.exit(grDevices::dev.off(), add = TRUE)
        draw_forest(m$md, m$re)
      }
    )
  })
}
