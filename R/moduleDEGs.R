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
      type = 4,  color = "#005249"
    )
  )
}

#' Server logic for DEG-by-dataset (robust for package + project)
#' @noRd
degTablesServer <- function(id, pkg = utils::packageName()) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # --- make pkg safe (length-1 string) ---------------------------------
    pkg <- tryCatch(pkg, error = function(e) "")
    if (!length(pkg) || !is.character(pkg) || !nzchar(pkg)) pkg <- "FRDATranscriptomicAtlas"
    pkg <- pkg[[1L]]

    # --- resolve paths with package-or-project fallback -------------------
    resolve_dir <- function(subpath) {
      d <- system.file(subpath, package = pkg, mustWork = FALSE)
      if (!nzchar(d)) d <- file.path("inst", subpath)
      d
    }
    resolve_file <- function(subpath, fname) {
      d <- resolve_dir(subpath)
      file.path(d, fname)
    }

    # one-time loads
    deg_dir_genes       <- resolve_dir(file.path("extdata", "deg", "genes"))
    deg_dir_transcripts <- resolve_dir(file.path("extdata", "deg", "transcripts"))

    tx2_path <- resolve_file(file.path("extdata", "maps"), "tx2gene.tsv")
    tx2 <- if (nzchar(tx2_path) && file.exists(tx2_path)) {
      readr::read_tsv(tx2_path, col_types = "ccc") |> dplyr::distinct()
    } else NULL

    gene_map <- if (!is.null(tx2)) {
      dplyr::select(tx2, gene_id, gene_name) |>
        dplyr::distinct() |>
        dplyr::rename(symbol = gene_name)
    } else NULL

    read_cached <- memoise::memoise(readRDS)

    # ---------- manifest of DEG files ----------
    manifest <- reactive({
      files <- c(
        if (nzchar(deg_dir_genes)       && dir.exists(deg_dir_genes))
          list.files(deg_dir_genes,       full.names = TRUE) else character(0),
        if (nzchar(deg_dir_transcripts) && dir.exists(deg_dir_transcripts))
          list.files(deg_dir_transcripts, full.names = TRUE) else character(0)
      )
      if (!length(files)) return(tibble::tibble())

      rx <- "^.*/DESEQ2_res_(.+)_(0\\.[0-9]+)_all_(genes|transcripts)\\.rds$"
      tibble::tibble(path = files) |>
        tidyr::extract(path, into = c("dataset","p_str","level"), regex = rx, remove = FALSE) |>
        dplyr::mutate(p = suppressWarnings(as.numeric(p_str)))
    })

    # ---------- dataset dropdown ----------
    `%||%` <- function(a,b) if (is.null(a)) b else a
    # --- load pretty_map from package namespace (internal object) ----------
    pretty_map <- tryCatch(
      get("pretty_map", envir = asNamespace(pkg)),
      error = function(e) {
        warning("pretty_map could not be found; using empty vector.")
        character(0)
      }
    )


    observe({
      m <- manifest()
      lvl <- input$feature_level %||% "genes"
      avail_ids <- sort(unique(m$dataset[m$level == lvl]))

      pm_sub <- pretty_map[avail_ids]
      pm_sub[is.na(pm_sub)] <- avail_ids[is.na(pm_sub)]
      labelled_choices <- stats::setNames(avail_ids, pm_sub)

      updateSelectizeInput(session, "dataset",
                           choices  = labelled_choices,
                           selected = if (length(avail_ids)) avail_ids[1] else NULL,
                           server   = TRUE
      )
    })

    # ---------- file selection ----------
    file_sel <- reactive({
      req(input$dataset, input$feature_level)
      m <- manifest()
      cand <- dplyr::filter(m, dataset == input$dataset, level == input$feature_level)
      if (nrow(cand)) cand$path[1] else NULL
    })

    # ---------- main data reactive ----------
    dat <- reactive({
      fp <- file_sel()
      validate(need(!is.null(fp) && file.exists(fp), "No results file found."))

      x <- read_cached(fp)
      if (!is.data.frame(x)) x <- as.data.frame(x)

      # ID column by level
      lvl <- input$feature_level %||% "genes"
      id_col <- if (identical(lvl, "genes")) "ensembl_gene_id" else "transcript_id"

      if (!(id_col %in% names(x)) && !is.null(rownames(x))) {
        x <- tibble::rownames_to_column(x, var = id_col)
      }

      # Fix occasional mislabeled IDs
      if (identical(lvl, "transcripts") &&
          ("ensembl_gene_id" %in% names(x)) && !("transcript_id" %in% names(x)) &&
          any(grepl("^ENST", head(x$ensembl_gene_id, 20)))) {
        x <- dplyr::rename(x, transcript_id = ensembl_gene_id)
      }
      if (identical(lvl, "genes") &&
          ("transcript_id" %in% names(x)) && !("ensembl_gene_id" %in% names(x)) &&
          any(grepl("^ENSG", head(x$transcript_id, 20)))) {
        x <- dplyr::rename(x, ensembl_gene_id = transcript_id)
      }

      # standardize columns
      if (!"log2FoldChange" %in% names(x)) {
        if ("log2FC" %in% names(x))    x <- dplyr::rename(x, log2FoldChange = log2FC)
        else if ("beta" %in% names(x)) x <- dplyr::rename(x, log2FoldChange = beta)
      }
      if (!"padj" %in% names(x) && "qvalue" %in% names(x)) x <- dplyr::rename(x, padj = qvalue)

      # p-value filter
      thr <- suppressWarnings(as.numeric(input$p_filter_mode))
      if (!is.na(thr) && "padj" %in% names(x)) {
        x <- x[!is.na(x$padj) & x$padj <= thr, , drop = FALSE]
      }

      # |log2FC| + direction
      lfc_min <- input$lfc_min %||% 0
      dir     <- input$direction %||% "both"
      if ("log2FoldChange" %in% names(x)) {
        if (identical(dir, "up"))   x <- x[x$log2FoldChange >=  lfc_min, , drop = FALSE]
        if (identical(dir, "down")) x <- x[x$log2FoldChange <= -lfc_min, , drop = FALSE]
        if (identical(dir, "both")) x <- x[abs(x$log2FoldChange) >= lfc_min, , drop = FALSE]
      }

      # symbols mapping
      if (identical(lvl, "genes")) {
        if (!is.null(gene_map) && "ensembl_gene_id" %in% names(x)) {
          x <- dplyr::left_join(x, gene_map, by = c("ensembl_gene_id" = "gene_id")) |>
            dplyr::relocate(ensembl_gene_id, symbol, .before = dplyr::everything())
        }
      } else {
        if (!is.null(tx2) && "transcript_id" %in% names(x)) {
          x <- dplyr::left_join(x, tx2, by = "transcript_id") |>
            dplyr::rename(symbol = gene_name) |>
            dplyr::relocate(transcript_id, gene_id, symbol, .before = dplyr::everything())
        }
      }

      x
    })

    # ---------- UI bits ----------
    pretty_label <- function(id) pretty_map[[id]] %||% id

    output$summary_bar <- renderUI({
      req(dat())
      x <- dat()
      tags$div(
        class = "alert alert-info",
        sprintf("Level: %s | Dataset: %s | p ≤ %s | |log2FC| ≥ %s | Number of Results: %s",
                input$feature_level, pretty_label(input$dataset),
                input$p_filter_mode, input$lfc_min, format(nrow(x), big.mark=",")))
    })

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

    output$download_filtered <- downloadHandler(
      filename = function() sprintf("DEG_%s_%s_p%s_filtered.csv",
                                    input$feature_level, input$dataset, input$p_filter_mode),
      content  = function(file) readr::write_csv(dat(), file)
    )
  })
}
