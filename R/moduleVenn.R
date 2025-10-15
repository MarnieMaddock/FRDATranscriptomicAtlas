#' Venn of shared DEGs / DEIs across datasets — UI
#' @param id module id
#' @noRd
# Sidebar
degVennUI <- function(id) {
  ns <- NS(id)
  tagList(
    h4("Venn of shared DEGs / DEIs"),
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
          ns("p_filter_mode"), "P-value threshold",
          choices = c("None" = NA, "≤ 0.10" = 0.10, "≤ 0.05" = 0.05,
                      "≤ 0.01" = 0.01, "≤ 0.001" = 0.001),
          selected = "0.05"
        ),
        numericInput(
          ns("lfc_min"), "Minimum |log₂FC|",
          value = 0, min = 0, max = 10, step = 0.1
        ),
        tags$small(
          em("Venn diagrams are accurate up to 6 datasets. When > 6 are selected, shared-gene tables are shown instead.")
        )
      ),
      column(
        width = 8,
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
    # ---- new: plot download row ----
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
    div(
      class = "mb-2",
      # download button for totals
      downloadButton(ns("dl_totals_csv"), "Download totals (CSV)")
    ),
    br(),
    h5("Shared genes/isoforms across datasets"),
    DT::dataTableOutput(ns("venn_overlaps")),
    div(
      class = "mb-2",
      # download button for overlaps
      downloadButton(ns("dl_overlaps_csv"), "Download overlaps (CSV)")
    ),
  )
}



#' Venn of shared DEGs / DEIs across datasets — SERVER
#' @noRd
degVennServer <- function(id, pkg = utils::packageName()) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # --- packages quietly ---
    requireNamespace("ggVennDiagram", quietly = TRUE)
    requireNamespace("ggplot2", quietly = TRUE)
    requireNamespace("stringr", quietly = TRUE)
    requireNamespace("svglite", quietly = TRUE)

    # --- paths + cached reader ---
    deg_dir_genes       <- system.file("extdata/deg/genes",       package = pkg, mustWork = FALSE)
    deg_dir_transcripts <- system.file("extdata/deg/transcripts", package = pkg, mustWork = FALSE)
    read_cached <- memoise::memoise(readRDS)

    # --- manifest ---
    manifest <- reactive({
      files <- c(
        if (nzchar(deg_dir_genes))       list.files(deg_dir_genes,       full.names = TRUE) else character(0),
        if (nzchar(deg_dir_transcripts)) list.files(deg_dir_transcripts, full.names = TRUE) else character(0)
      )
      if (!length(files)) return(tibble::tibble())
      rx <- "^.*/DESEQ2_res_(.+)_(0\\.[0-9]+)_all_(genes|transcripts)\\.rds$"
      tibble::tibble(path = files) |>
        tidyr::extract(path, into = c("dataset","p_str","level"), regex = rx, remove = FALSE) |>
        dplyr::mutate(p = suppressWarnings(as.numeric(p_str)))
    })

    # --- pretty names ---
    pretty_map <- c(
      "Erwin"             = "Erwin (Lymphoblastoid Cells)",
      "Indelicato"        = "Indelicato (Skeletal Muscle)",
      "Lai_iPSC"          = "Lai (iPSCs)",
      "Lai_CNS"           = "Lai (CNS neurons)",
      "Lai_PNS"           = "Lai (PNS neurons)",
      "Lees_FA1"          = "Lees (Cardiomyocytes) – FA1",
      "Lees_FA2"          = "Lees (Cardiomyocytes) – FA2",
      "Lees_FA3"          = "Lees (Cardiomyocytes) – FA3",
      "Maddock_LMN_FA2"   = "Maddock (Lower Motor Neurons) – FA2",
      "Maddock_SN_FA1"    = "Maddock (Sensory Neurons) – FA1",
      "Maddock_SN_FA2"    = "Maddock (Sensory Neurons) – FA2",
      "Maddock_NCC_FA1"   = "Maddock (Neural Crest Cells) – FA1",
      "Maddock_NCC_FA2"   = "Maddock (Neural Crest Cells) – FA2",
      "Mishra_223"        = "Mishra (Neurons) – 223",
      "Mishra_850"        = "Mishra (Neurons) – 850",
      "Mishra_FF1"        = "Mishra (Neurons) – FF1",
      "Mishra_FF2"        = "Mishra (Neurons) – FF2",
      "Napierala"         = "Napierala (Fibroblasts)",
      "Vilema"            = "Vilema-Enriquez (Fibroblasts)"
    )
    `%||%` <- function(x, y) if (is.null(x)) y else x
    pretty_label <- function(id) pretty_map[[id]] %||% id

    # --- populate dataset list per level ---
    observeEvent(input$feature_level, {
      m <- manifest()
      lvl <- input$feature_level %||% "genes"
      avail <- sort(unique(m$dataset[m$level == lvl]))
      labs <- unname(pretty_map[avail]); labs[is.na(labs)] <- avail[is.na(labs)]
      updateCheckboxGroupInput(session, "datasets",
                               choices = stats::setNames(avail, labs))
    }, ignoreInit = FALSE)

    # --- helper: read + filter IDs ---
    one_set_ids <- function(dataset_id, lvl, thr, lfc_min, direction) {
      m <- manifest()
      f <- dplyr::filter(m, dataset == dataset_id, level == lvl)
      if (!nrow(f)) return(character(0))
      f <- f[which.min(f$p), , drop = FALSE]
      x <- read_cached(f$path[1]) |> as.data.frame()

      if (!"log2FoldChange" %in% names(x)) {
        if ("log2FC" %in% names(x)) x$log2FoldChange <- x$log2FC
        else if ("beta" %in% names(x)) x$log2FoldChange <- x$beta
      }
      if (!"padj" %in% names(x) && "qvalue" %in% names(x)) x$padj <- x$qvalue

      id_col <- if (lvl == "genes") "ensembl_gene_id" else "transcript_id"
      if (!(id_col %in% names(x)) && !is.null(rownames(x))) x[[id_col]] <- rownames(x)

      if (!is.na(thr) && "padj" %in% names(x))
        x <- subset(x, !is.na(padj) & padj <= thr)
      if ("log2FoldChange" %in% names(x)) {
        if (direction == "up")   x <- subset(x, log2FoldChange >=  lfc_min)
        if (direction == "down") x <- subset(x, log2FoldChange <= -lfc_min)
        if (direction == "both") x <- subset(x, abs(log2FoldChange) >= lfc_min)
      }
      unique(x[[id_col]][!is.na(x[[id_col]])])
    }

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

    # --- overlap table helper ---
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

    # ======================
    # New helpers + downloads
    # ======================

    # filenames
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

    # tables as reactives (reused by DT + download)
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

    # Venn ggplot object (NULL if > 6)
    venn_plot_obj <- reactive({
      s <- sets_list()
      if (length(s) > 6) return(NULL)

      labs <- names(s)
      labs <- stringr::str_wrap(labs, width = 24)
      names(s) <- labs
        ggVennDiagram::ggVennDiagram(s, label = "count", label_size = 8, set_size = 10) +
          ggplot2::scale_fill_gradient(low = "#ccdcda", high = "#005249") +
          ggplot2::theme_void(base_size = 30) +
          ggplot2::theme(
            legend.position = "right",
            plot.margin = ggplot2::margin(60, 120, 60, 120)
          ) +
          ggplot2::coord_cartesian(clip = "off")
    })

    # --- plot render ---
    output$venn_plot <- renderPlot({
      p <- venn_plot_obj()
      if (is.null(p)) {
        plot.new(); text(0.5, 0.5, "Overlap tables shown below", cex = 1.6)
      } else {
        print(p)
      }
    })

    # --- optionally disable plot download buttons when >6 datasets
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

    # --- Totals table (paginated; no DT copy/csv buttons) ---
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

    # --- Overlaps table (paginated; no DT copy/csv buttons) ---
    output$venn_overlaps <- DT::renderDataTable({
      DT::datatable(
        overlaps_tbl(),
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

    # --- downloads: Venn (SVG/PNG) ---
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

    # --- downloads: tables (CSV) ---
    output$dl_totals_csv <- downloadHandler(
      filename = function() paste0(fname_prefix(), "_totals.csv"),
      content = function(file) utils::write.csv(totals_tbl(), file, row.names = FALSE)
    )

    output$dl_overlaps_csv <- downloadHandler(
      filename = function() paste0(fname_prefix(), "_overlaps.csv"),
      content = function(file) utils::write.csv(overlaps_tbl(), file, row.names = FALSE)
    )
  })
}


