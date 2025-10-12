#' DEGs-by-dataset UI module
#'
#' @param id Module id
#' @return A Shiny UI for exploring DEGs by dataset
#' @noRd
degTablesSidebarUI <- function(id) {
  ns <- NS(id)
  tagList(
    h4("DEG filters"),

    radioButtons(
      ns("feature_level"),
      label = "Level",
      choices = c("Genes" = "genes",
                  "Isoforms (transcripts)" = "transcripts"),
      selected = "genes"
    ),

    selectizeInput(
      ns("dataset"),
      "Dataset",
      choices = NULL, multiple = FALSE,
      options = list(placeholder = "Select a dataset…")
    ),
    radioButtons(
      ns("p_filter_mode"),
      "P-value threshold",
      inline = FALSE,
      choices = c("None" = NA,
                  "≤ 0.10" = 0.10,
                  "≤ 0.05"  = 0.05,
                  "≤ 0.01"  = 0.01,
                  "≤ 0.001" = 0.001),
      selected = "0.05"
    ),
    numericInput(
      ns("lfc_min"),
      label = "Minimum |log2FC|",
      value = 0, min = 0, max = 10, step = 0.1
    ),
    radioButtons(
      ns("direction"),
      label = "Direction",
      inline = TRUE,
      choices = c("Both" = "both", "Up" = "up", "Down" = "down"),
      selected = "both"
    ),
    strong("Download"),
    br(),
    downloadButton(ns("download_filtered"), "Download as CSV")
  )
}

#' Main area for DEG explorer (summary + table)
#' @noRd
degTablesMainUI <- function(id) {
  ns <- NS(id)
  tagList(
    br(),
    uiOutput(ns("summary_bar")),
    shinycssloaders::withSpinner(
      DT::dataTableOutput(ns("deg_table"), width = "100%"),
      type = 4
    )
  )
}

# R/mod_deg_tables_server.R

#' Server logic for DEG-by-dataset
#' @noRd
degTablesServer <- function(id, pkg = utils::packageName()) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # --- cache reads to speed repeat views
    read_cached <- memoise::memoise(readRDS)

    # --- Build a manifest of available files (once per session) ---
    manifest <- reactiveVal({
      roots <- list(
        genes       = system.file("extdata/deg/genes",       package = pkg, mustWork = FALSE),
        transcripts = system.file("extdata/deg/transcripts", package = pkg, mustWork = FALSE)
      )
      files <- unlist(lapply(roots, function(p) if (nzchar(p)) list.files(p, full.names = TRUE) else character()),
                      use.names = FALSE)

      if (!length(files)) return(dplyr::tibble())

      rx <- "^.*/DESEQ2_res_(.+)_(0\\.[0-9]+)_all_(genes|transcripts)\\.rds$"
      tibble::tibble(path = files) |>
        tidyr::extract(path, into = c("dataset","p_str","level"), regex = rx, remove = FALSE) |>
        dplyr::mutate(p = as.numeric(p_str))
    })

    # ---- Update dataset list when level changes ----
    observe({
      m <- manifest()
      lvl <- input$feature_level %||% "genes"
      ds <- sort(unique(m$dataset[m$level == lvl]))
      updateSelectizeInput(session, "dataset", choices = ds, selected = if (length(ds)) ds[1] else NULL, server = TRUE)
    })


    # Locate the correct file for the current selections
      # file selection: ignore p threshold; just take any for dataset+level
      file_sel <- reactive({
        req(input$dataset, input$feature_level)
        m <- manifest()
        cand <- dplyr::filter(m, dataset == input$dataset, level == input$feature_level)
        if (nrow(cand)) cand$path[1] else NULL   # or pick by max/min p, doesn’t matter
      })

    # --- load optional gene symbol map once ---
      # Load tx2gene once (ENST -> ENSG + symbol)
      tx2_path <- system.file("extdata/maps/tx2gene.tsv", package = pkg, mustWork = FALSE)
      tx2 <- if (nzchar(tx2_path) && file.exists(tx2_path)) {
        readr::read_tsv(tx2_path, col_types = "ccc") |>  # transcript_id, gene_id, gene_name
          dplyr::distinct()
      } else NULL

      # Build gene-level map from tx2gene (ENSG -> symbol)
      gene_map <- if (!is.null(tx2)) {
        dplyr::select(tx2, gene_id, gene_name) |>
          dplyr::distinct() |>
          dplyr::rename(symbol = gene_name)
      } else NULL

    # Load table
      dat <- reactive({
        fp <- file_sel()
        validate(need(!is.null(fp) && file.exists(fp), "No results file found."))

        x <- readRDS(fp)
        x <- as.data.frame(x)

        # --- decide which ID we need (genes vs transcripts) ---
        lvl <- input$feature_level %||% "genes"
        id_col <- if (identical(lvl, "genes")) "ensembl_gene_id" else "transcript_id"

        # --- promote rownames to the correct ID column (BEFORE any other renames/joins) ---
        if (!(id_col %in% names(x)) && !is.null(rownames(x))) {
          x <- tibble::rownames_to_column(x, var = id_col)
        }

        # If the wrong ID column contains the right IDs (e.g., ENST in ensembl_gene_id), fix it:
        if (identical(lvl, "transcripts") &&
            ("ensembl_gene_id" %in% names(x)) &&
            !("transcript_id" %in% names(x)) &&
            any(grepl("^ENST", head(x$ensembl_gene_id, 20)))) {
          x <- dplyr::rename(x, transcript_id = ensembl_gene_id)
        }
        if (identical(lvl, "genes") &&
            ("transcript_id" %in% names(x)) &&
            !("ensembl_gene_id" %in% names(x)) &&
            any(grepl("^ENSG", head(x$transcript_id, 20)))) {
          x <- dplyr::rename(x, ensembl_gene_id = transcript_id)
        }

        # --- standardize column names we filter on ---
        if (!"log2FoldChange" %in% names(x)) {
          if ("log2FC" %in% names(x))    x <- dplyr::rename(x, log2FoldChange = log2FC)
          else if ("beta" %in% names(x)) x <- dplyr::rename(x, log2FoldChange = beta)
        }
        if (!"padj" %in% names(x) && "qvalue" %in% names(x)) x <- dplyr::rename(x, padj = qvalue)

        # --- apply p-value + LFC filters ---
        thr <- suppressWarnings(as.numeric(input$p_filter_mode))
        if (!is.na(thr) && "padj" %in% names(x)) x <- x[!is.na(x$padj) & x$padj <= thr, , drop = FALSE]

        lfc_min <- input$lfc_min %||% 0
        dir     <- input$direction %||% "both"
        if ("log2FoldChange" %in% names(x)) {
          if (dir == "up")   x <- x[x$log2FoldChange >=  lfc_min, , drop = FALSE]
          if (dir == "down") x <- x[x$log2FoldChange <= -lfc_min, , drop = FALSE]
          if (dir == "both") x <- x[abs(x$log2FoldChange) >= lfc_min, , drop = FALSE]
        }

        # --- join symbols ---
        # If you have tx2gene.tsv (transcripts → gene_id + gene_name) as you showed:
        tx2_path <- system.file("extdata/maps/tx2gene.tsv", package = pkg, mustWork = FALSE)
        tx2 <- if (nzchar(tx2_path) && file.exists(tx2_path)) readr::read_tsv(tx2_path, col_types = "ccc") else NULL

        if (identical(lvl, "genes")) {
          # Build a gene map from tx2 if you don't also ship a separate gene map RDS
          gene_map <- if (!is.null(tx2)) {
            dplyr::select(tx2, gene_id, gene_name) |>
              dplyr::distinct() |>
              dplyr::rename(symbol = gene_name)
          } else NULL

          if (!is.null(gene_map) && "ensembl_gene_id" %in% names(x)) {
            x <- dplyr::left_join(x, gene_map, by = c("ensembl_gene_id" = "gene_id")) |>
              dplyr::relocate(ensembl_gene_id, symbol, .before = dplyr::everything())
          }

        } else { # transcripts
          if (!is.null(tx2) && "transcript_id" %in% names(x)) {
            x <- dplyr::left_join(x, tx2, by = "transcript_id") |>
              dplyr::rename(symbol = gene_name) |>
              dplyr::relocate(transcript_id, gene_id, symbol, .before = dplyr::everything())
          }
        }

        x
      })


    # Summary bar
    output$summary_bar <- renderUI({
      req(dat())
      x <- dat()
      tags$div(
        class = "alert alert-info",
        sprintf("Level: %s | Dataset: %s | p ≤ %s | |log2FC| ≥ %s | Rows: %s",
                input$feature_level, input$dataset,
                input$p_filter_mode, input$lfc_min, format(nrow(x), big.mark=",")))
    })

    # Render big table (DT with server-side processing)
    output$deg_table <- DT::renderDataTable({
      req(dat())
      DT::datatable(
        dat(),
        filter = "top",
        rownames = FALSE,
        extensions = "Buttons",
        options = list(
          dom = "Bfrtip",
          buttons = "copy",
          pageLength = 25,
          deferRender = TRUE,
          scrollX = TRUE
        )
      )
    }, server = TRUE)

    # Downloads
    output$download_filtered <- downloadHandler(
      filename = function() sprintf("DEG_%s_%s_p%s_filtered.csv",
                                    input$feature_level, input$dataset, input$p_filter_mode),
      content  = function(file) {
        x <- dat()
        readr::write_csv(x, file)
      }
    )

  })
}

