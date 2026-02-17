# mod_query_gene_across_datasets.R
#' Query Gene Across Datasets (DEG/DET lookup)
#' @noRd

queryGeneAcrossDatasetsSidebarUI <- function(id) {
  ns <- shiny::NS(id)

  shiny::tagList(
    shiny::h4("Query Gene Across Datasets"),

    shiny::radioButtons(
      ns("feature_level"),
      label = "Level",
      choices = c("Genes" = "genes",
                  "Isoforms (transcripts)" = "transcripts"),
      selected = "genes"
    ),

    shiny::conditionalPanel(
      condition = sprintf("input['%s'] === 'transcripts'", ns("feature_level")),
      shiny::div(
        class = "alert alert-warning",
        style = "margin-bottom:10px; font-size:0.9em;",
        shiny::strong("Important: "),
        "Isoform-level results are sensitive to sequencing depth, read length, and library layout. ",
        "Users are strongly encouraged to consult the dataset-level ",
        shiny::strong("isoform confidence scores in the About --> Sequencing Metrics tab"),
        " before interpreting transcript-level differential expression."
      )
    ),

    shiny::uiOutput(ns("datasets_ui")),


    shiny::textInput(
      ns("gene_query"),
      label = "Gene / ID (symbol or Ensembl ID)",
      value = "",
      placeholder = "e.g., FXN or ENSG00000165060 (or ENST... if transcripts)"
    ),


    shiny::hr(),

    shiny::radioButtons(
      ns("p_filter_mode"),
      "Adjusted P-value Threshold (for highlighting)",
      inline = FALSE,
      choiceNames = list(
        "None",
        shiny::HTML("&le; 0.10"),
        shiny::HTML("&le; 0.05"),
        shiny::HTML("&le; 0.01"),
        shiny::HTML("&le; 0.001")
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

    shiny::numericInput(
      ns("lfc_min"),
      label = "Minimum |log2FC| (for highlighting)",
      value = 0, min = 0, max = 10, step = 0.1
    ),

    shiny::radioButtons(
      ns("direction"),
      label = "Direction (for highlighting)",
      inline = TRUE,
      choices = c("Both" = "both", "Up" = "up", "Down" = "down"),
      selected = "both"
    ),

    shiny::strong("Download"),
    shiny::br(),
    shiny::downloadButton(ns("download_gene_table"), "Download table as CSV")
  )
}


queryGeneAcrossDatasetsMainUI <- function(id) {
  ns <- shiny::NS(id)

  shiny::tagList(
    shiny::br(),
    shiny::uiOutput(ns("summary_bar")),
    shinycssloaders::withSpinner(
      DT::dataTableOutput(ns("gene_table"), width = "100%"),
      type = 4, color = "#005249"
    )
  )
}


queryGeneAcrossDatasetsServer <- function(id, pkg = utils::packageName()) {
  shiny::moduleServer(id, function(input, output, session) {

    # ---- safe pkg string ----
    pkg <- tryCatch(pkg, error = function(e) "")
    if (!length(pkg) || !is.character(pkg) || !nzchar(pkg)) pkg <- "FRDATranscriptomicAtlas"
    pkg <- pkg[[1L]]

    # ---- helpers ----
    `%||%` <- function(a, b) if (is.null(a)) b else a
    norm_id <- function(x) sub("\\.\\d+$", "", as.character(x))

    # ---- ensure cached data ----
    ensure_atlas_data(
      keys    = c("deg_genes", "deg_transcripts"),
      package = pkg
    )

    cache_root <- tools::R_user_dir(pkg, which = "cache")
    deg_dir_genes <- file.path(cache_root, "genes")
    deg_dir_transcripts <- file.path(cache_root, "txs")

    if (!dir.exists(deg_dir_genes)) {
      stop("Gene-level DEG data not found: ", deg_dir_genes, call. = FALSE)
    }
    if (!dir.exists(deg_dir_transcripts)) {
      stop("Transcript-level DEG data not found: ", deg_dir_transcripts, call. = FALSE)
    }

    # ---- pretty names ----
    pretty_map <- tryCatch(
      get("pretty_map", envir = asNamespace(pkg)),
      error = function(e) character(0)
    )
    pretty_label <- function(id) pretty_map[[id]] %||% id


    # ---- maps (tx2gene) ----
    tx2_path <- system.file("extdata", "maps", "tx2gene.tsv", package = pkg, mustWork = FALSE)
    tx2 <- if (nzchar(tx2_path) && file.exists(tx2_path)) {
      readr::read_tsv(tx2_path, col_types = "ccc") |>
        dplyr::mutate(
          transcript_id = norm_id(transcript_id),
          gene_id       = norm_id(gene_id)
        ) |>
        dplyr::distinct()
    } else NULL

    gene_map <- if (!is.null(tx2)) {
      tx2 |>
        dplyr::select(gene_id, gene_name) |>
        dplyr::distinct() |>
        dplyr::rename(symbol = gene_name)
    } else NULL

    tx_map <- if (!is.null(tx2)) {
      tx2 |>
        dplyr::select(transcript_id, gene_id, gene_name) |>
        dplyr::distinct() |>
        dplyr::rename(symbol = gene_name)
    } else NULL

    read_cached <- memoise::memoise(readRDS)

    # ---- manifest of DEG files ----
    manifest <- shiny::reactive({
      files <- c(
        if (nzchar(deg_dir_genes)       && dir.exists(deg_dir_genes))
          list.files(deg_dir_genes,       full.names = TRUE) else character(0),
        if (nzchar(deg_dir_transcripts) && dir.exists(deg_dir_transcripts))
          list.files(deg_dir_transcripts, full.names = TRUE) else character(0)
      )
      if (!length(files)) return(tibble::tibble())

      rx <- "^.*[\\\\/]{1}DESEQ2_res_(.+)_(0\\.[0-9]+)_all_(genes|transcripts)\\.rds$"
      tibble::tibble(path = files) |>
        tidyr::extract(path, into = c("dataset", "p_str", "level"), regex = rx, remove = FALSE) |>
        dplyr::mutate(p = suppressWarnings(as.numeric(p_str)))
    })

    output$datasets_ui <- shiny::renderUI({
      m   <- manifest()
      lvl <- input$feature_level %||% "genes"

      avail <- sort(unique(m$dataset[m$level == lvl]))

      labs <- pretty_map[avail]
      labs[is.na(labs)] <- avail[is.na(labs)]

      shiny::checkboxGroupInput(
        session$ns("datasets"),
        label   = "Datasets (select multiple)",
        choices = stats::setNames(avail, labs),
        selected = intersect(input$datasets %||% character(0), avail)
      )
    })




    # ---- pick the “most permissive” file for each dataset/level ----
    # If p in filename corresponds to a pre-filter, we want the largest p (least strict).
    file_for_dataset <- function(dataset_id, lvl) {
      m <- manifest()
      cand <- dplyr::filter(m, dataset == dataset_id, level == lvl)
      if (!nrow(cand)) return(NA_character_)

      # choose max p; if p missing, just take first
      if (all(is.na(cand$p))) return(cand$path[1])
      cand <- cand[order(cand$p, decreasing = TRUE), , drop = FALSE]
      cand$path[1]
    }

    # ---- harmonise DESeq2 result columns ----
    standardise_res <- function(df, lvl) {
      if (!is.data.frame(df)) df <- as.data.frame(df)

      # ID column by level
      id_col <- if (lvl == "genes") "ensembl_gene_id" else "transcript_id"
      if (!(id_col %in% names(df)) && !is.null(rownames(df))) {
        df <- tibble::rownames_to_column(df, var = id_col)
      }

      # common alias fixes
      if (!"log2FoldChange" %in% names(df)) {
        if ("log2FC" %in% names(df)) df <- dplyr::rename(df, log2FoldChange = log2FC)
        else if ("beta" %in% names(df)) df <- dplyr::rename(df, log2FoldChange = beta)
      }
      if (!"padj" %in% names(df) && "qvalue" %in% names(df)) df <- dplyr::rename(df, padj = qvalue)

      # normalise ID versions
      if (id_col %in% names(df)) df[[id_col]] <- norm_id(df[[id_col]])

      df
    }

    # ---- main: build per-dataset row table for the selected gene ----
    gene_table_dat <- shiny::reactive({
      req(input$feature_level, input$gene_query)

      q_raw <- trimws(input$gene_query)
      shiny::validate(shiny::need(nzchar(q_raw), "Enter a gene symbol or Ensembl ID."))
      q <- norm_id(q_raw)

      req(input$datasets)
      shiny::validate(shiny::need(length(input$datasets) >= 1, "Select at least one dataset."))

      query_id <- NULL
      query_symbol <- NULL

      lvl <- input$feature_level %||% "genes"

      if (lvl == "genes") {
        if (grepl("^ENSG", q, ignore.case = FALSE)) {
          query_id <- q
          if (!is.null(gene_map)) {
            query_symbol <- gene_map$symbol[match(query_id, gene_map$gene_id)]
          }
        } else {
          # treat as symbol
          query_symbol <- q_raw
          if (!is.null(gene_map)) {
            hit <- gene_map$gene_id[match(toupper(query_symbol), toupper(gene_map$symbol))]
            if (!is.na(hit)) query_id <- hit
          }
        }
      } else {
        if (grepl("^ENST", q, ignore.case = FALSE)) {
          query_id <- q
          if (!is.null(tx_map)) {
            query_symbol <- tx_map$symbol[match(query_id, tx_map$transcript_id)]
          }
        } else {
          # treat as symbol (returns many transcripts; we can’t guess one)
          query_symbol <- q_raw
        }
      }

      # pull each dataset result and extract row(s)
      rows <- lapply(input$datasets, function(ds) {
        fp <- file_for_dataset(ds, lvl)
        if (!nzchar(fp) || is.na(fp) || !file.exists(fp)) {
          return(tibble::tibble(
            dataset_id = ds,
            dataset = pretty_label(ds),
            query = q_raw,
            found = FALSE
          ))
        }

        res <- standardise_res(read_cached(fp), lvl)

        id_col <- if (lvl == "genes") "ensembl_gene_id" else "transcript_id"

        # if user provided symbol at transcript level, we can only match via tx_map
        if (lvl == "transcripts" && (is.null(query_id) || !nzchar(query_id))) {
          if (is.null(tx_map)) {
            return(tibble::tibble(
              dataset_id = ds,
              dataset = pretty_label(ds),
              query = q_raw,
              found = FALSE
            ))
          }
          # select all transcripts whose symbol matches
          tx_hits <- tx_map$transcript_id[toupper(tx_map$symbol) == toupper(query_symbol)]
          tx_hits <- unique(tx_hits)
          if (!length(tx_hits)) {
            return(tibble::tibble(
              dataset_id = ds,
              dataset = pretty_label(ds),
              query = q_raw,
              found = FALSE
            ))
          }
          sub <- res[res[[id_col]] %in% tx_hits, , drop = FALSE]
          if (!nrow(sub)) {
            return(tibble::tibble(
              dataset_id = ds,
              dataset = pretty_label(ds),
              query = q_raw,
              found = FALSE
            ))
          }
          # keep all transcript hits for that symbol in that dataset
          sub <- sub |>
            dplyr::mutate(dataset_id = ds, dataset = pretty_label(ds), query = q_raw, found = TRUE)

          # add symbol/gene_id if missing
          if (!("symbol" %in% names(sub)) && !is.null(tx_map) && "transcript_id" %in% names(sub)) {
            sub <- dplyr::left_join(sub,
                                    tx_map |> dplyr::select(transcript_id, gene_id, symbol),
                                    by = "transcript_id")
          }
          return(sub)
        }

        # match by ID (genes or transcripts)
        if (is.null(query_id) || !nzchar(query_id)) {
          # typed symbol but no mapping available
          return(tibble::tibble(
            dataset_id = ds,
            dataset = pretty_label(ds),
            query = q_raw,
            found = FALSE
          ))
        }

        hit <- res[res[[id_col]] == query_id, , drop = FALSE]
        if (!nrow(hit)) {
          return(tibble::tibble(
            dataset_id = ds,
            dataset = pretty_label(ds),
            query = q_raw,
            found = FALSE
          ))
        }

        hit <- hit |>
          dplyr::mutate(dataset_id = ds, dataset = pretty_label(ds), query = q_raw, found = TRUE)

        # add symbol for genes if missing
        if (lvl == "genes" && !("symbol" %in% names(hit)) && !is.null(gene_map) && "ensembl_gene_id" %in% names(hit)) {
          hit <- dplyr::left_join(hit,
                                  gene_map |> dplyr::select(gene_id, symbol),
                                  by = c("ensembl_gene_id" = "gene_id"))
        }

        # add gene_id/symbol for transcripts if missing
        if (lvl == "transcripts" && !("symbol" %in% names(hit)) && !is.null(tx_map) && "transcript_id" %in% names(hit)) {
          hit <- dplyr::left_join(hit,
                                  tx_map |> dplyr::select(transcript_id, gene_id, symbol),
                                  by = "transcript_id")
        }

        hit
      })

      out <- dplyr::bind_rows(rows)

      # ---- highlight rules (optional) ----
      thr     <- suppressWarnings(as.numeric(input$p_filter_mode))
      lfc_min <- input$lfc_min %||% 0
      dir     <- input$direction %||% "both"

      n <- nrow(out)

      # pass_p: only meaningful if thr is set AND padj exists
      if (!is.na(thr) && "padj" %in% names(out)) {
        pass_p <- !is.na(out$padj) & (out$padj <= thr)
      } else {
        pass_p <- rep(NA, n)
      }

      # pass_dir: only meaningful if log2FoldChange exists
      if ("log2FoldChange" %in% names(out)) {
        lfc <- out$log2FoldChange
        pass_dir <- dplyr::case_when(
          is.na(lfc) ~ FALSE,
          dir == "up"   ~ lfc >=  lfc_min,
          dir == "down" ~ lfc <= -lfc_min,
          TRUE          ~ abs(lfc) >= lfc_min
        )
      } else {
        pass_dir <- rep(NA, n)
      }

      # meets_threshold logic
      if (is.na(thr)) {
        meets_threshold <- pass_dir
      } else {
        meets_threshold <- pass_p & pass_dir
      }

      # attach as columns
      out <- out |>
        dplyr::mutate(
          pass_p         = pass_p,
          pass_dir       = pass_dir,
          meets_threshold = meets_threshold
        ) |>
        dplyr::mutate(
          `Present in dataset`    = found,
          `FDR ≤ threshold`       = pass_p,
          `Meets log2FC filter`   = pass_dir,
          `Meets all filters`     = meets_threshold
        ) |>
        dplyr::select(-found, -pass_p, -pass_dir, -meets_threshold, -dataset_id)


      # Prefer a tidy, stable column order
      keep_first <- c("dataset", "query", "found",
                      "ensembl_gene_id", "transcript_id", "gene_id", "symbol",
                      "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj",
                      "meets_threshold")
      keep_first <- keep_first[keep_first %in% names(out)]
      out <- dplyr::relocate(out, dplyr::all_of(keep_first), .before = dplyr::everything())

      # round numeric columns (except pvalue/padj)
      out <- out |>
        dplyr::mutate(
          dplyr::across(
            .cols = dplyr::where(is.numeric) & !c(pvalue, padj),
            ~ round(.x, 4)
          )
        )

      out
    })

    output$summary_bar <- shiny::renderUI({
      req(input$feature_level, input$gene_query, input$datasets)
      x <- gene_table_dat()

      tags <- shiny::tags
      tags$div(
        class = "alert alert-info",
        tags$span(
          shiny::HTML(sprintf(
            "Level: %s | Gene/ID query: <b>%s</b> | Datasets selected: %d | Rows returned: %s",
            input$feature_level,
            htmltools::htmlEscape(input$gene_query),
            length(input$datasets),
            format(nrow(x), big.mark = ",")
          ))
        )
      )
    })

    output$gene_table <- DT::renderDataTable({
      req(gene_table_dat())
      DT::datatable(
        gene_table_dat(),
        filter = "top",
        rownames = FALSE,
        extensions = "Buttons",
        options = list(
          dom = "Bfrtip",
          buttons = c("copy"),
          pageLength = 25,
          deferRender = TRUE,
          scrollX = TRUE
        )
      )
    }, server = TRUE)

    output$download_gene_table <- shiny::downloadHandler(
      filename = function() {
        lvl <- input$feature_level %||% "genes"
        q <- gsub("[^A-Za-z0-9_\\-\\.]+", "_", input$gene_query %||% "query")
        sprintf("DEG_query_%s_%s_%ddatasets.csv", lvl, q, length(input$datasets %||% character(0)))
      },
      content = function(file) {
        readr::write_csv(gene_table_dat(), file)
      }
    )
  })
}
