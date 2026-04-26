#' @importFrom graphics plot.new text
degVennUI <- function(id) {
  ns <- NS(id)
  tagList(
    h4("Venn of shared DEGs / DETs"),
    fluidRow(
      column(
        width = 4,
        radioButtons(
          ns("feature_level"), "Level",
          choices = c("Genes" = "genes", "Isoforms (transcripts)" = "transcripts"),
          selected = "genes"
        ),
        radioButtons(
          ns("direction"), "Direction", inline = TRUE,
          choices = c("Both" = "both", "Up" = "up", "Down" = "down"),
          selected = "both"
        ),
        radioButtons(
          ns("p_filter_mode"),
          "Adjusted P-value threshold",
          inline = FALSE,
          choiceNames = list(
            "None",
            HTML("&le; 0.10"),
            HTML("&le; 0.05"),
            HTML("&le; 0.01"),
            HTML("&le; 0.001")
          ),
          choiceValues = list(
            NA,
            0.10,
            0.05,
            0.01,
            0.001
          ),
          selected = 0.05
        ),
        numericInput(
          ns("lfc_min"),
          HTML("Minimum |log<sub>2</sub>FC|"),
          value = 0, min = 0, max = 10, step = 0.1
        ),
        tags$small(
          em("Venn diagrams are accurate up to 6 datasets. When > 6 are selected, shared-gene tables are shown instead.")
        )
      ),
      column(
        width = 8,
        div(
          style = "display:flex; gap:8px; margin-bottom:8px;",
          actionButton(ns("datasets_all"),  "Select all", class = "btn btn-sm btn-default"),
          actionButton(ns("datasets_none"), "Clear",      class = "btn btn-sm btn-default")
        ),

        checkboxGroupInput(
          ns("datasets"),
          label = "Datasets (select any number)",
          choices = character(0)
        ),

        div(class = "alert alert-info", uiOutput(ns("selection_info")))
      )
    )
  )
}

degVennMainUI <- function(id) {
  ns <- NS(id)
  tagList(
    uiOutput(ns("mode_msg")),
    shinycssloaders::withSpinner(
      plotOutput(ns("venn_plot"), height = 720, width = 1100),
      type = 4, color = "#005249"
    ),
    fluidRow(
      column(
        width = 12,
        div(class = "mb-3",
            downloadButton(ns("dl_venn_svg"), "Download Venn (SVG)"),
            tags$span(" "),
            downloadButton(ns("dl_venn_png"), "Download Venn (PNG)")
        )
      )
    ),
    br(),
    h5("Total filtered genes/isoforms per dataset"),
    DT::dataTableOutput(ns("venn_totals")),
    div(class = "mb-2", downloadButton(ns("dl_totals_csv"), "Download totals (CSV)")),
    br(),
    h5("Shared genes/isoforms across datasets"),
    DT::dataTableOutput(ns("venn_overlaps")),
    div(class = "mb-2", downloadButton(ns("dl_overlaps_csv"), "Download overlaps (CSV)")),
    br(),
    DT::dataTableOutput(ns("overlap_items")),
    div(class = "mb-2", downloadButton(ns("dl_overlap_items_csv"), "Download item list (CSV)"))
  )
}

# --- SERVER -------------------------------------------------------------------

degVennServer <- function(
    id,
    pkg = utils::packageName(),
    data_mode = "local"
) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # --- packages quietly ---
    requireNamespace("ggVennDiagram", quietly = TRUE)
    requireNamespace("ggplot2", quietly = TRUE)
    requireNamespace("stringr", quietly = TRUE)
    requireNamespace("svglite", quietly = TRUE)

    ensure_atlas_data(
      keys = c("deg_genes", "deg_transcripts"),
      package = pkg,
      data_mode = data_mode
    )

    # --- paths + cached reader ---
    # cache_root <- tools::R_user_dir(pkg, which = "cache")
    #
    # deg_dir_genes <- file.path(cache_root, "genes")
    # deg_dir_transcripts <- file.path(cache_root, "txs")
    #
    # validate(
    #   need(dir.exists(deg_dir_genes), "Gene-level DEG data not available."),
    #   need(dir.exists(deg_dir_transcripts), "Transcript-level DEG data not available.")
    # )
    # # --- cached RDS reader ---
    # read_cached <- memoise::memoise(readRDS)


    # --- maps (tx2gene) ---
    tx2_path <- system.file("extdata/maps/tx2gene.tsv", package = pkg, mustWork = FALSE)
    if (!nzchar(tx2_path)) tx2_path <- file.path("inst", "extdata", "maps", "tx2gene.tsv")
    tx2gene <- if (nzchar(tx2_path) && file.exists(tx2_path)) {
      readr::read_tsv(tx2_path, col_types = "ccc") |>
        # expected cols: transcript_id, gene_id, gene_name
        dplyr::mutate(
          transcript_id = sub("\\.\\d+$","", transcript_id),
          gene_id       = sub("\\.\\d+$","", gene_id)
        )
    } else NULL

    norm_id <- function(x) sub("\\.\\d+$","", as.character(x))

    # --- manifest ---
    # manifest <- reactive({
    #   files <- c(
    #     if (nzchar(deg_dir_genes))       list.files(deg_dir_genes,       full.names = TRUE) else character(0),
    #     if (nzchar(deg_dir_transcripts)) list.files(deg_dir_transcripts, full.names = TRUE) else character(0)
    #   )
    #   if (!length(files)) return(tibble::tibble())
    #   rx <- "^.*/DESEQ2_res_(.+)_(0\\.[0-9]+)_all_(genes|transcripts)\\.rds$"
    #   tibble::tibble(path = files) |>
    #     tidyr::extract(path, into = c("dataset","p_str","level"), regex = rx, remove = FALSE) |>
    #     dplyr::mutate(p = suppressWarnings(as.numeric(p_str)))
    # })
    manifest <- reactiveVal(NULL)

    observeEvent(TRUE, {
      m <- get_deg_manifest(pkg = pkg, data_mode = data_mode)
      manifest(m)
    }, once = TRUE)

    # --- pretty names ---
    `%||%` <- function(x, y) if (is.null(x)) y else x

    pretty_map <- tryCatch(
      get("pretty_map", envir = asNamespace(pkg)),
      error = function(e) {
        warning("pretty_map not found in namespace; using empty vector.")
        character(0)
      }
    )

    pretty_label <- function(id) pretty_map[[id]] %||% id

    # --- populate dataset list per level ---
    observeEvent(input$feature_level, {
      # m <- manifest()
      m <- req(manifest())
      lvl <- input$feature_level %||% "genes"
      avail <- sort(unique(m$dataset[m$level == lvl]))
      labs <- pretty_map[avail]
      labs[is.na(labs)] <- avail[is.na(labs)]
      updateCheckboxGroupInput(session, "datasets",
                               choices = stats::setNames(avail, labs),
                               selected = intersect(input$datasets %||% character(0), avail))
    }, ignoreInit = FALSE)

    # Select all (for current feature level)
    observeEvent(input$datasets_all, {
      #m   <- manifest()
      m <- req(manifest())
      lvl <- input$feature_level %||% "genes"

      avail <- sort(unique(m$dataset[m$level == lvl]))
      updateCheckboxGroupInput(session, "datasets", selected = avail)
    })

    # Clear
    observeEvent(input$datasets_none, {
      updateCheckboxGroupInput(session, "datasets", selected = character(0))
    })

    # --- helper: read + filter IDs ---
    # Return the filtered ID vector for one dataset (genes or transcripts)
    # one_set_ids <- function(dataset_id, lvl = c("genes","transcripts"),
    #                         thr, lfc_min = 0, direction = c("both","up","down")) {
    #   lvl       <- match.arg(lvl)
    #   direction <- match.arg(direction)
    #
    #   # normalise Ensembl-style IDs: ENSGxxxx.xx -> ENSGxxxx
    #   norm_id <- function(x) sub("\\.\\d+$", "", as.character(x))
    #
    #   #m <- manifest()
    #   m <- req(manifest())
    #   f <- dplyr::filter(m, dataset == dataset_id, level == lvl)
    #   if (!nrow(f)) return(character(0))
    #
    #   # if multiple p-threshold files exist, pick the most permissive (min p)
    #   f <- f[which.min(f$p), , drop = FALSE]
    #
    #   x <- read_cached(f$path[1]) |> as.data.frame()
    #
    #   # --- harmonise common column names ---
    #   # log2FC
    #   if (!"log2FoldChange" %in% names(x)) {
    #     if ("log2FC" %in% names(x))        x$log2FoldChange <- x$log2FC
    #     else if ("beta" %in% names(x))     x$log2FoldChange <- x$beta
    #   }
    #   # padj
    #   if (!"padj" %in% names(x)) {
    #     if ("qvalue" %in% names(x))        x$padj <- x$qvalue
    #     else if ("adj.P.Val" %in% names(x)) x$padj <- x$adj.P.Val
    #   }
    #
    #   # ID column per level (accept several variants; fallback to rownames)
    #   id_candidates <- if (lvl == "genes") {
    #     c("ensembl_gene_id", "gene_id", "EnsemblGeneID")
    #   } else {
    #     c("transcript_id", "ensembl_transcript_id", "tx_id")
    #   }
    #   id_col <- id_candidates[id_candidates %in% names(x)][1]
    #   if (is.na(id_col) || is.null(id_col)) {
    #     if (!is.null(rownames(x))) {
    #       x[["__tmp_id"]] <- rownames(x)
    #       id_col <- "__tmp_id"
    #     } else {
    #       return(character(0))
    #     }
    #   }
    #
    #   # --- apply filters ---
    #   # p-value / FDR
    #   thr_num <- suppressWarnings(as.numeric(thr))
    #   if (!is.na(thr_num) && "padj" %in% names(x)) {
    #     x <- subset(x, !is.na(padj) & padj <= thr_num)
    #   }
    #
    #   # LFC / direction
    #   if ("log2FoldChange" %in% names(x) && is.finite(lfc_min) && lfc_min >= 0) {
    #     if (direction == "up") {
    #       x <- subset(x, !is.na(log2FoldChange) & log2FoldChange >=  lfc_min)
    #     } else if (direction == "down") {
    #       x <- subset(x, !is.na(log2FoldChange) & log2FoldChange <= -lfc_min)
    #     } else { # both
    #       x <- subset(x, !is.na(log2FoldChange) & abs(log2FoldChange) >= lfc_min)
    #     }
    #   }
    #
    #   ids_raw <- x[[id_col]]
    #   ids <- ids_raw[!is.na(ids_raw) & nzchar(ids_raw)]
    #   unique(norm_id(ids))
    # }
    one_set_ids <- function(dataset_id,
                            lvl = c("genes", "transcripts"),
                            thr,
                            lfc_min = 0,
                            direction = c("both", "up", "down")) {
      lvl <- match.arg(lvl)
      direction <- match.arg(direction)

      m <- req(manifest())

      f <- dplyr::filter(m, dataset == dataset_id, level == lvl)
      if (!nrow(f)) return(character(0))

      thr_num <- suppressWarnings(as.numeric(thr))

      # if (identical(data_mode, "local")) {
      #   if (is.na(thr_num)) {
      #     f <- f[which.max(f$p), , drop = FALSE]
      #   } else {
      #     f <- dplyr::filter(f, p == thr_num)
      #     if (!nrow(f)) return(character(0))
      #     f <- f[1, , drop = FALSE]
      #   }
      #
      #   file_path <- f$path[1]
      # } else {
      #   file_path <- NULL
      # }
      if (identical(data_mode, "local")) {

        # Local files may only exist for one threshold, often p = 0.05.
        # Use the most permissive available local file, then filter inside get_deg_data().
        f2 <- f[!is.na(f$p), , drop = FALSE]

        if (nrow(f2)) {
          f <- f2[which.max(f2$p), , drop = FALSE]
        } else {
          f <- f[1, , drop = FALSE]
        }

        file_path <- f$path[1]

      } else {
        file_path <- NULL
      }

      x <- get_deg_data(
        dataset = dataset_id,
        level = lvl,
        file_path = file_path,
        data_mode = data_mode,
        padj_max = thr_num,
        lfc_min = lfc_min,
        direction = direction
      )

      if (!is.data.frame(x)) x <- as.data.frame(x)

      id_candidates <- if (lvl == "genes") {
        c("ensembl_gene_id", "gene_id", "EnsemblGeneID")
      } else {
        c("transcript_id", "ensembl_transcript_id", "tx_id")
      }

      id_col <- id_candidates[id_candidates %in% names(x)][1]

      if (is.na(id_col) || is.null(id_col)) {
        if (!is.null(rownames(x))) {
          x[["__tmp_id"]] <- rownames(x)
          id_col <- "__tmp_id"
        } else {
          return(character(0))
        }
      }

      ids <- x[[id_col]]
      ids <- ids[!is.na(ids) & nzchar(ids)]

      unique(sub("\\.\\d+$", "", as.character(ids)))
    }

    # # --- NEW: build ID -> symbol map for current level ---
    id_symbol_map <- reactive({
      lvl <- input$feature_level %||% "genes"

      if (!is.null(tx2gene)) {
        if (lvl == "genes") {
          return(
            tx2gene |>
              dplyr::distinct(gene_id, gene_name) |>
              dplyr::rename(id = gene_id, symbol = gene_name)
          )
        } else {
          return(
            tx2gene |>
              dplyr::distinct(transcript_id, gene_name) |>
              dplyr::rename(id = transcript_id, symbol = gene_name)
          )
        }
      }

      tibble::tibble(id = character(0), symbol = character(0))
    })
    # id_symbol_map <- reactive({
    #   lvl <- input$feature_level %||% "genes"
    #
    #   # If we have tx2gene, build clean maps up front
    #   if (!is.null(tx2gene)) {
    #     if (lvl == "genes") {
    #       # gene_id -> gene_name (HGNC symbol)
    #       gm <- tx2gene |>
    #         dplyr::distinct(gene_id, gene_name) |>
    #         dplyr::rename(id = gene_id, symbol = gene_name)
    #       return(gm)
    #     } else {
    #       # transcripts: map transcript_id -> gene_name (show gene symbol for each transcript)
    #       tm <- tx2gene |>
    #         dplyr::distinct(transcript_id, gene_name) |>
    #         dplyr::rename(id = transcript_id, symbol = gene_name)
    #       return(tm)
    #     }
    #   }
    #
    #   # ---- Fallback (if no tx2gene file available) ----
    #   #m <- manifest() |> dplyr::filter(level == lvl)
    #   m <- req(manifest()) |> dplyr::filter(level == lvl)
    #   if (!nrow(m)) return(tibble::tibble(id = character(0), symbol = character(0)))
    #
    #   prefer_gene_syms <- c("gene_name","symbol","external_gene_name","hgnc_symbol")
    #   prefer_tx_syms   <- c("transcript_name","tx_name")
    #
    #   maps <- lapply(m$path, function(p) {
    #     df <- read_cached(p) |> as.data.frame()
    #     if (lvl == "genes") {
    #       id_col <- if ("ensembl_gene_id" %in% names(df)) "ensembl_gene_id"
    #       else if (!is.null(rownames(df))) { df[["__tmp_id"]] <- rownames(df); "__tmp_id" }
    #       else NA_character_
    #       sym_col <- intersect(prefer_gene_syms, names(df))
    #     } else {
    #       id_col <- if ("transcript_id" %in% names(df)) "transcript_id"
    #       else if (!is.null(rownames(df))) { df[["__tmp_id"]] <- rownames(df); "__tmp_id" }
    #       else NA_character_
    #       sym_col <- intersect(prefer_tx_syms, names(df))
    #     }
    #     if (is.na(id_col)) return(NULL)
    #     ids <- norm_id(df[[id_col]])
    #     sym <- if (length(sym_col)) df[[sym_col[1]]] else ids
    #     tibble::tibble(id = ids, symbol = as.character(sym))
    #   })
    #
    #   d <- dplyr::bind_rows(Filter(Negate(is.null), maps)) |>
    #     dplyr::distinct(id, .keep_all = TRUE)
    #   d[!is.na(d$id) & nzchar(d$id), , drop = FALSE]
    # })



    # --- reactive list of sets ---
    sets_list <- reactive({
      req(length(input$datasets) >= 2)
      lvl <- input$feature_level %||% "genes"
      thr <- suppressWarnings(as.numeric(input$p_filter_mode))
      lfc <- input$lfc_min %||% 0
      dir <- input$direction %||% "both"
      ids <- lapply(input$datasets, one_set_ids,
                    lvl = lvl, thr = thr, lfc_min = lfc, direction = dir)
      names(ids) <- vapply(input$datasets, pretty_label, "", USE.NAMES = FALSE)
      ids <- ids[vapply(ids, length, 1L) > 0L]
      validate(need(length(ids) >= 2, "Need at least two non-empty sets."))
      ids
    })

    # --- overlap table helper (non-exclusive counts) ---
    venn_overlap_tbl <- function(s) {
      Universe <- unique(unlist(s, use.names = FALSE))
      M <- vapply(s, function(v) Universe %in% v, logical(length(Universe)))
      colnames(M) <- names(s)
      which_sets <- apply(M, 1, function(r) names(s)[r])
      which_sets <- which_sets[lengths(which_sets) >= 2]
      if (!length(which_sets)) {
        return(data.frame(
          `Dataset Combination`              = character(0),
          `Number of Datasets`               = integer(0),
          `Number of Shared Genes/Isoforms`  = integer(0),
          check.names = FALSE
        ))
      }
      keys <- vapply(which_sets, function(v) paste(sort(v), collapse = " & "), character(1))
      tt <- sort(table(keys), decreasing = TRUE)
      out <- data.frame(
        `Dataset Combination`             = names(tt),
        `Number of Shared Genes/Isoforms` = as.integer(tt),
        check.names = FALSE, row.names = NULL
      )
      out$`Number of Datasets` <- 1 + stringr::str_count(out$`Dataset Combination`, " & ")
      out[order(out$`Number of Datasets`, -out$`Number of Shared Genes/Isoforms`), ]
    }

    # --- mode message ---
    output$mode_msg <- renderUI({
      req(input$datasets)
      n <- length(input$datasets)
      if (n <= 6) {
        div(class = "alert alert-success",
            sprintf("Showing Venn diagram for %d datasets. Counts = overlapping genes/isoforms.", n))
      } else {
        div(class = "alert alert-warning",
            sprintf("You selected %d datasets (> 6). Displaying shared-gene tables instead of a Venn diagram.", n))
      }
    })

    # --- filenames ---
    fname_prefix <- reactive({
      thr <- suppressWarnings(as.numeric(input$p_filter_mode))
      thr_txt <- if (is.na(thr)) "pnone" else sprintf("p%g", thr)
      paste0(
        "venn_",
        (input$feature_level %||% "genes"), "_",
        (input$direction %||% "both"), "_",
        thr_txt, "_lfc", (input$lfc_min %||% 0), "_",
        length(input$datasets), "sets"
      )
    })

    # --- tables reused by DT + download ---
    totals_tbl <- reactive({
      s <- sets_list()
      data.frame(
        Dataset = names(s),
        `Total Filtered Genes/Isoforms` = as.integer(vapply(s, length, 1L)),
        check.names = FALSE
      )
    })

    overlaps_tbl <- reactive({
      s <- sets_list()
      tbl <- venn_overlap_tbl(s)
      num_cols <- c("Number of Datasets", "Number of Shared Genes/Isoforms")
      if (nrow(tbl)) tbl[num_cols] <- lapply(tbl[num_cols], as.integer)
      tbl
    })

    # --- Venn ggplot object (NULL if > 6) ---
    venn_plot_obj <- reactive({
      s <- sets_list()
      if (length(s) > 6) return(NULL)
      labs <- names(s)
      labs <- stringr::str_wrap(labs, width = 24)
      names(s) <- labs
      suppressWarnings(suppressMessages(
      ggVennDiagram::ggVennDiagram(s, label = "count", label_size = 8, set_size = 8) +
        ggplot2::scale_fill_gradient(low = "#ccdcda", high = "#005249") +
        ggplot2::theme_void(base_size = 30) +
        ggplot2::theme(legend.position = "right",
                       plot.margin = ggplot2::margin(60, 120, 60, 120)) +
        ggplot2::coord_cartesian(clip = "off")
      ))
    })

    # --- plot render ---
    output$venn_plot <- renderPlot({
      p <- venn_plot_obj()
      if (is.null(p)) {
        plot.new(); text(0.5, 0.5, "Overlap tables shown below", cex = 1.6)
      } else {
        suppressWarnings(print(p))
      }
    })

    # --- enable/disable plot download buttons when >6 datasets ---
    observe({
      have_plot <- !is.null(venn_plot_obj())
      if (requireNamespace("shinyjs", quietly = TRUE)) {
        shinyjs::toggleState(ns("dl_venn_svg"), condition = have_plot)
        shinyjs::toggleState(ns("dl_venn_png"), condition = have_plot)
      }
    })

    # --- selection info ---
    output$selection_info <- renderUI({
      req(input$datasets)
      tags$span(sprintf("Selected: %d dataset(s).", length(input$datasets)))
    })

    # --- Totals table ---
    output$venn_totals <- DT::renderDataTable({
      DT::datatable(
        totals_tbl(),
        rownames = FALSE,
        options = list(
          dom = "tip",
          pageLength = 10,
          lengthMenu = list(c(10, 25, 50, 100), c("10", "25", "50", "100")),
          deferRender = TRUE,
          scrollX = TRUE
        )
      )
    }, server = TRUE)

    # --- Overlaps table ---
    output$venn_overlaps <- DT::renderDataTable({
      DT::datatable(
        overlaps_tbl(),
        rownames = FALSE,
        selection = "single",
        options = list(
          dom = "tip",
          pageLength = 10,
          lengthMenu = list(c(10, 25, 50, 100), c("10", "25", "50", "100")),
          deferRender = TRUE,
          scrollX = TRUE
        )
      )
    }, server = TRUE)

    # keep combo list in sync with overlaps table
    # observe({
    #   tbl <- overlaps_tbl()
    #   combos <- tbl$`Dataset Combination`
    #
    #   # preserve current selection if still valid; else fall back to first
    #   sel <- isolate(input$combo_pick)
    #   if (is.null(sel) || !nzchar(sel) || !(sel %in% combos)) {
    #     sel <- combos[[1]] %||% ""
    #   }
    #
    #   updateSelectizeInput(
    #     session,
    #     "combo_pick",
    #     choices  = combos,
    #     selected = sel,
    #     server   = TRUE
    #   )
    # })
    observe({
      tbl <- overlaps_tbl()
      combos <- tbl$`Dataset Combination`

      if (!length(combos)) {
        updateSelectizeInput(
          session,
          "combo_pick",
          choices = character(0),
          selected = character(0),
          server = TRUE
        )
        return()
      }

      sel <- isolate(input$combo_pick)

      if (is.null(sel) || !nzchar(sel) || !(sel %in% combos)) {
        sel <- combos[1]
      }

      updateSelectizeInput(
        session,
        "combo_pick",
        choices = combos,
        selected = sel,
        server = TRUE
      )
    })

    # react to row click in overlaps DT
    observeEvent(input$venn_overlaps_rows_selected, {
      idx <- input$venn_overlaps_rows_selected
      tbl <- overlaps_tbl()
      if (length(idx) && nrow(tbl) >= idx) {
        updateSelectizeInput(
          session,
          "combo_pick",
          selected = tbl$`Dataset Combination`[idx],
          server   = TRUE
        )
      }
    })


    combo_items <- reactive({
      s <- sets_list()
      req(length(s) >= 2)

      # Universe + membership matrix
      Universe <- unique(unlist(s, use.names = FALSE))
      if (!length(Universe)) {
        return(tibble::tibble(
          ID = character(0),
          Symbol = character(0),
          Sum = integer(0)
        ))
      }

      M <- vapply(s, function(v) Universe %in% v, logical(length(Universe)))
      colnames(M) <- names(s)

      # ID -> symbol mapping
      map <- id_symbol_map()
      sym <- map$symbol[match(Universe, map$id)]
      sym[is.na(sym) | !nzchar(sym)] <- Universe

      presence <- as.data.frame(M, stringsAsFactors = FALSE)
      presence[] <- lapply(presence, as.integer)
      presence$Sum <- rowSums(presence)

      out <- data.frame(
        ID = Universe,
        Symbol = sym,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
      out <- dplyr::bind_cols(out, presence)
      out <- out[, c("ID", "Symbol", "Sum", setdiff(names(out), c("ID","Symbol","Sum")))]
      out <- dplyr::arrange(out, dplyr::desc(Sum))


      out
    })


    # --- render items table ---
    output$overlap_items <- DT::renderDataTable({
      dat <- combo_items()
      DT::datatable(
        dat,
        rownames = FALSE,
        options = list(
          dom = "tip",
          pageLength = 25,
          lengthMenu = list(c(25, 50, 100, 200), c("25", "50", "100", "200")),
          deferRender = TRUE,
          scrollX = TRUE
        )
      )
    }, server = TRUE)

    # --- downloads: plot + tables ---
    output$dl_venn_svg <- downloadHandler(
      filename = function() paste0(fname_prefix(), ".svg"),
      content = function(file) {
        p <- venn_plot_obj(); req(p)
        svglite::svglite(file, width = 14, height = 9)
        on.exit(grDevices::dev.off(), add = TRUE)
        print(p)
      }
    )

    output$dl_venn_png <- downloadHandler(
      filename = function() paste0(fname_prefix(), ".png"),
      content = function(file) {
        p <- venn_plot_obj(); req(p)
        ggplot2::ggsave(filename = file, plot = p, width = 14, height = 9, dpi = 300)
      }
    )

    output$dl_totals_csv <- downloadHandler(
      filename = function() paste0(fname_prefix(), "_totals.csv"),
      content = function(file) utils::write.csv(totals_tbl(), file, row.names = FALSE)
    )

    output$dl_overlaps_csv <- downloadHandler(
      filename = function() paste0(fname_prefix(), "_overlaps.csv"),
      content = function(file) utils::write.csv(overlaps_tbl(), file, row.names = FALSE)
    )

    output$dl_overlap_items_csv <- downloadHandler(
      filename = function() {
        mode <- input$overlap_mode %||% "exact"
        paste0(fname_prefix(), "_items_", mode, ".csv")
      },
      content = function(file) utils::write.csv(combo_items(), file, row.names = FALSE)
    )
  })
}
