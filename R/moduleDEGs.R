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

#' Server logic for DEG-by-dataset (simplified & fast)
#' @noRd
degTablesServer <- function(id, pkg = utils::packageName()) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # ---------- one-time loads ----------
    deg_dir_genes       <- system.file("extdata/deg/genes",       package = pkg, mustWork = FALSE)
    deg_dir_transcripts <- system.file("extdata/deg/transcripts", package = pkg, mustWork = FALSE)

    tx2_path <- system.file("extdata/maps/tx2gene.tsv", package = pkg, mustWork = FALSE)
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
        if (nzchar(deg_dir_genes))       list.files(deg_dir_genes,       full.names = TRUE) else character(0),
        if (nzchar(deg_dir_transcripts)) list.files(deg_dir_transcripts, full.names = TRUE) else character(0)
      )
      if (!length(files)) return(dplyr::tibble())

      rx <- "^.*/DESEQ2_res_(.+)_(0\\.[0-9]+)_all_(genes|transcripts)\\.rds$"
      tibble::tibble(path = files) |>
        tidyr::extract(path, into = c("dataset","p_str","level"), regex = rx, remove = FALSE) |>
        dplyr::mutate(p = as.numeric(p_str))
    })

    # ---------- dataset dropdown ----------

    pretty_map <- c(
      "Erwin"             = "Erwin (Lymphoblastoid Cells)",
      "Indelicato"        = "Indelicato (Skeletal Muscle)",
      "Lai_iPSC"          = "Lai (iPSCs)",
      "Lai_CNS"           = "Lai (CNS neurons)",
      "Lai_PNS"           = "Lai (PNS neurons)",
      "Lees_FA1"          = "Lees (Cardiomyocytes) - FA1",
      "Lees_FA2"          = "Lees (Cardiomyocytes) - FA2",
      "Lees_FA3"          = "Lees (Cardiomyocytes) - FA3",
      "Maddock_LMN_FA2"   = "Maddock (Lower Motor Neurons) - FA2",
      "Maddock_SN_FA1"    = "Maddock (Sensory Neurons) - FA1",
      "Maddock_SN_FA2"    = "Maddock (Sensory Neurons) - FA2",
      "Maddock_NCC_FA1"   = "Maddock (Neural Crest Cells) - FA1",
      "Maddock_NCC_FA2"   = "Maddock (Neural Crest Cells) - FA2",
      "Mishra_223"        = "Mishra (Neurons) - 223",
      "Mishra_850"        = "Mishra (Neurons) - 850",
      "Mishra_FF1"        = "Mishra (Neurons) - FF1",
      "Mishra_FF2"        = "Mishra (Neurons) - FF2",
      "Napierala"         = "Napierala (Fibroblasts)",
      "Vilema"            = "Vilema-Enriquez (Fibroblasts)"
    )


    observe({
      m <- manifest()
      lvl <- input$feature_level %||% "genes"
      avail_ids <- sort(unique(m$dataset[m$level == lvl]))
      # keep only ids we actually have files for
      pm_sub <- pretty_map[avail_ids]                         # names = ids, values = labels
      # Build a named vector: names = labels (shown), values = ids (returned)
      labelled_choices <- stats::setNames(avail_ids, pm_sub)

      # Fallback for any ids that have no pretty name
      missing_ids <- setdiff(avail_ids, names(pretty_map))
      if (length(missing_ids)) {
        add <- stats::setNames(missing_ids, missing_ids)
        labelled_choices <- c(labelled_choices, add)
      }
      updateSelectizeInput(session, "dataset",
                           choices = labelled_choices,
                           selected = if (length(avail_ids)) avail_ids[1] else NULL,
                           server = TRUE
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

      # symbols
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

    # Summary bar
    pretty_label <- function(id) pretty_map[[id]] %||% id

    output$summary_bar <- renderUI({
      req(dat())
      x <- dat()
      tags$div(
        class = "alert alert-info",
        sprintf("Level: %s | Dataset: %s | p ≤ %s | |log2FC| ≥ %s | Numer of Results: %s",
                input$feature_level, pretty_label(input$dataset),
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

