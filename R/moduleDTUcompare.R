# ---- DTU Venn module ----------------------------------------------------
# expects rds files named like "<DATASET>_significant_isoforms.rds"
# with columns at least:
#   gene_id, gene_name, isoform_id, dIF, isoform_switch_q_value, direction
# direction strings should be like "Higher in FRDA" / "Higher in Control"

# helper info
dtuOverlapSidebarHelp <- function(ns) {
  tags$div(
    class = "alert alert-secondary", style = "margin-top: .5rem;",
    tags$h5("About this analysis"),
    tags$p(HTML(
      "This panel compares <strong>differential transcript usage (DTU)</strong> across datasets.<br/>
       You can view overlaps at two levels:"
    )),
    tags$ul(
      tags$li(HTML("<strong>Switching genes</strong>: a gene is counted if <em>any</em> isoform switches in that dataset (higher recall, tolerant to which isoform switched). i.e. Do the same genes show isoform switching across datasets? Note: It can hide heterogeneity. Two datasets might overlap at gene FXN, but one switched ENST00000A and the other ENST00000B.")),
      tags$li(HTML("<strong>Switching isoforms</strong>: matches by exact transcript ID (ENST…); reveals shared mechanisms (higher specificity, lower recall). i.e. Is the same transcript switching across datasets?"))
    ),
    tags$h6("Filters"),
    tags$ul(
      tags$li(HTML("<code>q-value</code> (isoform_switch_q_value): FDR threshold for switching calls.")),
      tags$li(HTML("<code>|dIF|</code>: minimum absolute change in isoform fraction.")),
      tags$li(HTML("<strong>Direction</strong>: “Higher in FRDA”, “Higher in Control”, or Both."))
    ),
    tags$h6("Outputs"),
    tags$ul(
      tags$li(HTML("<strong>Venn diagram</strong> (≤ 6 datasets): visual overlap of the selected sets.")),
      tags$li(HTML("<strong>Totals</strong>: number of filtered switching items per dataset.")),
      tags$li(HTML("<strong>Shared combinations</strong>: counts for each dataset combination.")),
      tags$li(HTML("<strong>Presence matrix</strong>: list of IDs (+ symbols) with a 0/1 membership per dataset."))
    )
  )
}



# UI: sidebar
dtuVennUI <- function(id) {
  ns <- NS(id)
  tagList(
    h4("Venn of shared DTUs (isoform switching)"),
    fluidRow(
      column(
        width = 4,
        radioButtons(
          ns("dtu_level"), "Compare",
          choices = c("Switching genes" = "genes", "Switching isoforms" = "isoforms"),
          selected = "genes"
        ),
        radioButtons(
          ns("dtu_direction"), "Direction", inline = TRUE,
          choices = c("Both" = "both", "Higher in FRDA" = "frda", "Higher in Control" = "control"),
          selected = "both"
        ),
        radioButtons(
          ns("dtu_q_thr"), "q-value threshold (isoform_switch_q_value)",
          choices = c("≤ 0.10" = "0.10", "≤ 0.05" = "0.05",
                      "≤ 0.01" = "0.01", "≤ 0.001" = "0.001", "None" = "NA"),
          selected = "0.05"
        ),
        numericInput(
          ns("dtu_dif_min"), "Minimum |dIF|",
          value = 0.1, min = 0, max = 1, step = 0.01
        ),
        tags$small(
          em("Venn diagrams are accurate up to 6 datasets. When > 6 are selected, shared-item tables are shown instead.")
        )
      ),
      column(
        width = 8,
        checkboxGroupInput(
          ns("dtu_datasets"),
          label = "Datasets (select any number)",
          choices = character(0)
        ),
        div(class = "alert alert-info", uiOutput(ns("dtu_selection_info")))
      ),
      dtuOverlapSidebarHelp(ns)
    ),
    column(
      width = 8,
      checkboxGroupInput(ns("dtu_datasets"), "Datasets (select any number)", choices = character(0)),
      div(class = "alert alert-info", uiOutput(ns("dtu_selection_info")))
    )
  )
}


# UI: main panel
dtuVennMainUI <- function(id) {
  ns <- NS(id)
  tagList(
    uiOutput(ns("dtu_mode_msg")),
    shinycssloaders::withSpinner(
      plotOutput(ns("dtu_venn_plot"), height = 720, width = 1100),
      type = 4, color = "#005249"
    ),
    fluidRow(
      column(
        width = 12,
        div(class = "mb-3",
            downloadButton(ns("dl_dtu_venn_svg"), "Download Venn (SVG)"),
            tags$span(" "),
            downloadButton(ns("dl_dtu_venn_png"), "Download Venn (PNG)")
        )
      )
    ),
    br(),
    h5("Total filtered switching genes/isoforms per dataset"),
    DT::dataTableOutput(ns("dtu_totals")),
    div(class = "mb-2", downloadButton(ns("dl_dtu_totals_csv"), "Download totals (CSV)")),
    br(),
    h5("Shared switching genes/isoforms across datasets"),
    DT::dataTableOutput(ns("dtu_overlaps")),
    div(class = "mb-2", downloadButton(ns("dl_dtu_overlaps_csv"), "Download overlaps (CSV)")),
    br(),
    h5("Presence matrix (IDs + Symbols where available)"),
    DT::dataTableOutput(ns("dtu_overlap_items")),
    div(class = "mb-2", downloadButton(ns("dl_dtu_overlap_items_csv"), "Download item list (CSV)"))
  )
}

# SERVER
dtuVennServer <- function(id, pkg = utils::packageName()){
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    requireNamespace("ggVennDiagram", quietly = TRUE)
    requireNamespace("ggplot2", quietly = TRUE)
    requireNamespace("readr", quietly = TRUE)
    requireNamespace("dplyr", quietly = TRUE)
    requireNamespace("tidyr", quietly = TRUE)
    requireNamespace("DT", quietly = TRUE)
    requireNamespace("svglite", quietly = TRUE)
    requireNamespace("memoise", quietly = TRUE)

    `%||%` <- function(x, y) if (is.null(x)) y else x
    norm_id <- function(x) sub("\\.\\d+$","", as.character(x))

    # ---- ensure DTU data are available (Zenodo cache) ----
    ensure_atlas_data(
      keys    = "dtu_results",
      package = pkg
    )

    cache_root <- tools::R_user_dir(pkg, which = "cache")
    dtu_dir    <- file.path(cache_root, "dtu")

    if (!dir.exists(dtu_dir)) {
      stop(
        "DTU cache directory not found after download:\n  ",
        dtu_dir,
        call. = FALSE
      )
    }


    # ---- pretty dataset labels (optional – reuse/extend yours) ----
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
    pretty_label <- function(id) pretty_map[[id]] %||% id

    # ---- manifest of available DTU CSVs ----
    dtu_manifest <- reactive({
      if (!dir.exists(dtu_dir)) return(tibble::tibble())
      files <- list.files(dtu_dir, pattern = "_significant_isoforms\\.rds$", full.names = TRUE)
      if (!length(files)) return(tibble::tibble())
      tibble::tibble(
        path = files,
        dataset = sub("_significant_isoforms\\.rds$", "", basename(files))
      )
    })

    # populate dataset selector
    observe({
      m <- dtu_manifest()
      avail <- sort(unique(m$dataset))
      labs <- unname(pretty_map[avail]); labs[is.na(labs)] <- avail[is.na(labs)]
      updateCheckboxGroupInput(session, "dtu_datasets",
                               choices = stats::setNames(avail, labs))
    })

    # ---- read + filter one dataset to an ID vector ----
    # lvl: "genes" or "isoforms"
    # q_thr: numeric or NA
    # dIF_min: numeric (0..1)
    # dir_mode: "both" | "frda" | "control"
    dtu_one_set_ids <- function(dataset_id, lvl = c("genes","isoforms"),
                                q_thr, dIF_min, dir_mode) {
      lvl <- match.arg(lvl)
      m <- dtu_manifest()
      f <- dplyr::filter(m, dataset == dataset_id)
      if (!nrow(f)) return(character(0))

      df <- read_cached_rds(f$path[1]) |> as.data.frame()

      # Harmonise columns
      if (!("isoform_switch_q_value" %in% names(df))) return(character(0))
      if (!("dIF" %in% names(df))) return(character(0))

      # q-value filter
      thr_num <- suppressWarnings(as.numeric(q_thr))
      if (!is.na(thr_num)) {
        df <- subset(df, !is.na(isoform_switch_q_value) & isoform_switch_q_value <= thr_num)
      }

      # |dIF| filter
      if (is.finite(dIF_min) && dIF_min > 0) {
        df <- subset(df, !is.na(dIF) & abs(dIF) >= dIF_min)
      }

      # direction filter (expects strings like "Higher in FRDA" or "Higher in Control")
      if (dir_mode == "frda" && "direction" %in% names(df)) {
        df <- subset(df, tolower(direction) %in% "higher in frda")
      } else if (dir_mode == "control" && "direction" %in% names(df)) {
        df <- subset(df, tolower(direction) %in% "higher in control")
      }

      if (!nrow(df)) return(character(0))

      if (lvl == "genes") {
        # prefer gene_id if you want strict Ensembl IDs; gene_name is also available
        ids <- df$gene_id %||% df$gene_name
        ids <- ids[!is.na(ids) & nzchar(ids)]
        unique(norm_id(ids))
      } else {
        ids <- df$isoform_id
        ids <- ids[!is.na(ids) & nzchar(ids)]
        unique(norm_id(ids))
      }
    }

    # ---- ID -> display symbol (for presence matrix) ----
    # For genes: symbol = gene_name (fallback gene_id)
    # For isoforms: symbol = "<gene_name> (<isoform_id>)"
    dtu_id_symbol_map <- reactive({
      lvl <- input$dtu_level %||% "genes"
      m <- dtu_manifest()
      if (!nrow(m)) return(tibble::tibble(id = character(0), symbol = character(0)))

      # combine from all datasets (dedupe)
      maps <- lapply(m$path, function(p) {
        df <- read_cached_rds(p) |> as.data.frame()
        if (lvl == "genes") {
          gid <- norm_id(df$gene_id %||% df$gene_name)
          sym <- if ("gene_name" %in% names(df)) as.character(df$gene_name) else gid
          tibble::tibble(id = gid, symbol = sym)
        } else {
          tid <- norm_id(df$isoform_id)
          gname <- if ("gene_name" %in% names(df)) as.character(df$gene_name) else NA_character_
          sym <- ifelse(!is.na(gname) & nzchar(gname),
                        paste0(gname, " (", tid, ")"), tid)
          tibble::tibble(id = tid, symbol = sym)
        }
      })
      dplyr::bind_rows(maps) |>
        dplyr::filter(!is.na(id) & nzchar(id)) |>
        dplyr::distinct(id, .keep_all = TRUE)
    })

    # ---- build list of sets for the selected datasets ----
    dtu_sets_list <- reactive({
      req(length(input$dtu_datasets) >= 2)
      lvl <- input$dtu_level %||% "genes"
      q_thr <- suppressWarnings(as.numeric(input$dtu_q_thr))
      if (is.na(q_thr) && (input$dtu_q_thr %||% "NA") == "NA") q_thr <- NA_real_
      dif_min <- input$dtu_dif_min %||% 0
      dir_mode <- input$dtu_direction %||% "both"

      ids <- lapply(input$dtu_datasets, dtu_one_set_ids,
                    lvl = lvl, q_thr = q_thr, dIF_min = dif_min, dir_mode = dir_mode)
      names(ids) <- vapply(input$dtu_datasets, pretty_label, "", USE.NAMES = FALSE)
      ids <- ids[vapply(ids, length, 1L) > 0L]
      validate(need(length(ids) >= 2, "Need at least two non-empty sets after filtering."))
      ids
    })

    # ---- overlap summary table (non-exclusive) ----
    dtu_overlap_tbl <- function(s) {
      Universe <- unique(unlist(s, use.names = FALSE))
      M <- vapply(s, function(v) Universe %in% v, logical(length(Universe)))
      colnames(M) <- names(s)
      which_sets <- apply(M, 1, function(r) names(s)[r])
      which_sets <- which_sets[lengths(which_sets) >= 2]
      if (!length(which_sets)) {
        return(data.frame(
          `Dataset Combination`               = character(0),
          `Number of Datasets`                = integer(0),
          `Number of Shared Genes/Isoforms`   = integer(0),
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

    # ---- mode message ----
    output$dtu_mode_msg <- renderUI({
      req(input$dtu_datasets)
      n <- length(input$dtu_datasets)
      lbl <- if ((input$dtu_level %||% "genes") == "genes") "genes" else "isoforms"
      if (n <= 6) {
        div(class = "alert alert-success",
            sprintf("Showing Venn diagram for %d datasets. Counts = overlapping %s.", n, lbl))
      } else {
        div(class = "alert alert-warning",
            sprintf("You selected %d datasets (> 6). Displaying shared-item tables instead of a Venn diagram.", n))
      }
    })

    # ---- filename stem ----
    dtu_fname_prefix <- reactive({
      thr <- suppressWarnings(as.numeric(input$dtu_q_thr))
      thr_txt <- if (is.na(thr)) "qnone" else sprintf("q%g", thr)
      dir_txt <- input$dtu_direction %||% "both"
      lvl_txt <- input$dtu_level %||% "genes"
      paste0("dtu_venn_", lvl_txt, "_", dir_txt, "_", thr_txt,
             "_dif", (input$dtu_dif_min %||% 0), "_",
             length(input$dtu_datasets), "sets")
    })

    # ---- tables for DT + download ----
    dtu_totals_tbl <- reactive({
      s <- dtu_sets_list()
      data.frame(
        Dataset = names(s),
        `Total Filtered Switching Items` = as.integer(vapply(s, length, 1L)),
        check.names = FALSE
      )
    })

    dtu_overlaps_tbl <- reactive({
      s <- dtu_sets_list()
      tbl <- dtu_overlap_tbl(s)
      nc <- c("Number of Datasets", "Number of Shared Genes/Isoforms")
      if (nrow(tbl)) tbl[nc] <- lapply(tbl[nc], as.integer)
      tbl
    })

    # ---- Venn plot object (NULL if > 6 datasets) ----
    dtu_venn_plot_obj <- reactive({
      s <- dtu_sets_list()
      if (length(s) > 6) return(NULL)
      labs <- names(s)
      labs <- stringr::str_wrap(labs, width = 24)
      names(s) <- labs
      ggVennDiagram::ggVennDiagram(s, label = "count", label_size = 8, set_size = 10) +
        ggplot2::scale_fill_gradient(low = "#ccdcda", high = "#005249") +
        ggplot2::theme_void(base_size = 30) +
        ggplot2::theme(legend.position = "right",
                       plot.margin = ggplot2::margin(60, 120, 60, 120)) +
        ggplot2::coord_cartesian(clip = "off")
    })

    # ---- plot render ----
    output$dtu_venn_plot <- renderPlot({
      p <- dtu_venn_plot_obj()
      if (is.null(p)) {
        plot.new(); text(0.5, 0.5, "Overlap tables shown below", cex = 1.6)
      } else {
        suppressWarnings(print(p))
      }
    })

    observe({
      have_plot <- !is.null(dtu_venn_plot_obj())
      if (requireNamespace("shinyjs", quietly = TRUE)) {
        shinyjs::toggleState(ns("dl_dtu_venn_svg"), condition = have_plot)
        shinyjs::toggleState(ns("dl_dtu_venn_png"), condition = have_plot)
      }
    })

    # ---- selection info ----
    output$dtu_selection_info <- renderUI({
      req(input$dtu_datasets)
      tags$span(sprintf("Selected: %d dataset(s).", length(input$dtu_datasets)))
    })

    # ---- presence matrix (IDs + symbols) ----
    dtu_combo_items <- reactive({
      s <- dtu_sets_list(); req(length(s) >= 2)
      Universe <- unique(unlist(s, use.names = FALSE))
      M <- vapply(s, function(v) Universe %in% v, logical(length(Universe)))
      colnames(M) <- names(s)

      map <- dtu_id_symbol_map()
      sym <- map$symbol[match(Universe, map$id)]
      bad <- is.na(sym) | !nzchar(sym)
      sym[bad] <- Universe[bad]   # <- match lengths (no warning)


      presence <- as.data.frame(M, stringsAsFactors = FALSE)
      presence[] <- lapply(presence, as.integer)
      presence$Sum <- rowSums(presence)

      out <- data.frame(
        ID = Universe,
        Symbol = sym,
        Sum = presence$Sum,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
      out <- dplyr::bind_cols(out, presence[, setdiff(names(presence), "Sum"), drop = FALSE])
      dplyr::arrange(out, dplyr::desc(Sum))
    })

    # ---- tables render ----
    output$dtu_totals <- DT::renderDataTable({
      DT::datatable(
        dtu_totals_tbl(),
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

    output$dtu_overlaps <- DT::renderDataTable({
      DT::datatable(
        dtu_overlaps_tbl(),
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

    output$dtu_overlap_items <- DT::renderDataTable({
      DT::datatable(
        dtu_combo_items(),
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

    # ---- downloads ----
    output$dl_dtu_venn_svg <- downloadHandler(
      filename = function() paste0(dtu_fname_prefix(), ".svg"),
      content = function(file) {
        p <- dtu_venn_plot_obj(); req(p)
        svglite::svglite(file, width = 14, height = 9)
        on.exit(grDevices::dev.off(), add = TRUE)
        print(p)
      }
    )
    output$dl_dtu_venn_png <- downloadHandler(
      filename = function() paste0(dtu_fname_prefix(), ".png"),
      content = function(file) {
        p <- dtu_venn_plot_obj(); req(p)
        ggplot2::ggsave(filename = file, plot = p, width = 14, height = 9, dpi = 300)
      }
    )
    output$dl_dtu_totals_csv <- downloadHandler(
      filename = function() paste0(dtu_fname_prefix(), "_totals.csv"),
      content = function(file) utils::write.csv(dtu_totals_tbl(), file, row.names = FALSE)
    )
    output$dl_dtu_overlaps_csv <- downloadHandler(
      filename = function() paste0(dtu_fname_prefix(), "_overlaps.csv"),
      content = function(file) utils::write.csv(dtu_overlaps_tbl(), file, row.names = FALSE)
    )
    output$dl_dtu_overlap_items_csv <- downloadHandler(
      filename = function() paste0(dtu_fname_prefix(), "_items.csv"),
      content = function(file) utils::write.csv(dtu_combo_items(), file, row.names = FALSE)
    )
  })
}
