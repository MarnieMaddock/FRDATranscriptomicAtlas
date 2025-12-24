# =========================
# GSEA Compare (Venn) Module
# =========================

# ---- UI (sidebar) ----
gseaCompareUI <- function(id) {
  ns <- shiny::NS(id)
  tagList(
    h4("Venn of shared GSEA terms"),
    fluidRow(
      column(
        width = 4,
        radioButtons(
          ns("gsea_direction"), "Direction", inline = TRUE,
          choices = c("Both" = "both", "Up" = "up", "Down" = "down"),
          selected = "both"
        ),
        # --- explanatory info text ---
        div(
          class = "text-muted small mb-2",
          tags$p(
            HTML(
              "<strong>Matrix legend:</strong><br>
       <span style='color:#005249;'>+1</span> = pathway enriched (positive NES)<br>
       <span style='color:#a83232;'>−1</span> = pathway depleted (negative NES)<br>
       <span style='color:#666;'>0</span> = not significant or absent"
            )
          )
        ),

        tags$small(
          em("Venn diagrams are accurate up to 6 datasets. When > 6 are selected, shared-term tables are shown instead.")
        )
      ),
      column(
        width = 8,
        checkboxGroupInput(
          ns("gsea_datasets"),
          label = "Datasets (select any number)",
          choices = character(0)
        ),
        div(class = "alert alert-info", uiOutput(ns("gsea_selection_info")))
      )
    )
  )
}

# ---- UI (main panel) ----
gseaCompareMainUI <- function(id) {
  ns <- shiny::NS(id)
  tagList(
    uiOutput(ns("gsea_mode_msg")),
    shinycssloaders::withSpinner(
      plotOutput(ns("gsea_venn_plot"), height = 720, width = 1100),
      type = 4, color = "#005249"
    ),
    fluidRow(
      column(
        width = 12,
        div(class = "mb-3",
            downloadButton(ns("dl_gsea_venn_svg"), "Download Venn (SVG)"),
            tags$span(" "),
            downloadButton(ns("dl_gsea_venn_png"), "Download Venn (PNG)")
        )
      )
    ),
    br(),
    h5("Total filtered terms per dataset"),
    DT::dataTableOutput(ns("gsea_totals")),
    div(class = "mb-2", downloadButton(ns("dl_gsea_totals_csv"), "Download totals (CSV)")),
    br(),
    h5("Shared terms across datasets"),
    DT::dataTableOutput(ns("gsea_overlaps")),
    div(class = "mb-2", downloadButton(ns("dl_gsea_overlaps_csv"), "Download overlaps (CSV)")),
    br(),
    DT::dataTableOutput(ns("gsea_overlap_items")),
    div(class = "mb-2", downloadButton(ns("dl_gsea_overlap_items_csv"), "Download term list (CSV)"))
  )
}

# ---- SERVER ----
gseaCompareServer <- function(id, pkg = utils::packageName()) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # ---- deps ----
    requireNamespace("readr", quietly = TRUE)
    requireNamespace("dplyr", quietly = TRUE)
    requireNamespace("tidyr", quietly = TRUE)
    requireNamespace("stringr", quietly = TRUE)
    requireNamespace("ggplot2", quietly = TRUE)
    requireNamespace("ggVennDiagram", quietly = TRUE)
    requireNamespace("DT", quietly = TRUE)
    requireNamespace("svglite", quietly = TRUE)

    `%||%` <- function(x, y) if (is.null(x)) y else x

    # ---- ensure GSEA summary data are available (Zenodo-backed) ----
    ensure_atlas_data(
      keys    = "gsea_summary",
      package = pkg
    )

    # ---- resolve cached summary CSV ----
    cache_root <- tools::R_user_dir(pkg, which = "cache")

    csv_path <- file.path(
      cache_root,
      "GSEA_results", "summary",
      "GO_GSEA_direction_summary.csv"
    )

    validate(
      need(
        file.exists(csv_path),
        paste0(
          "GSEA summary CSV not found after download.\nExpected at:\n",
          csv_path
        )
      )
    )

    # ---- read summary table once ----
    gsea_tbl <- reactive({
      validate(need(file.exists(csv_path), paste0("GSEA summary CSV not found at: ", csv_path)))
      # Expect columns:
      # term, <dataset columns...> with values in {-1, 0, 1}, and optionally up_count, down_count, sig_count, consensus
      readr::read_csv(csv_path, show_col_types = FALSE)
    })

    # ---- identify dataset columns & dataset IDs ----
    dataset_cols <- reactive({
      df <- gsea_tbl()
      # known non-dataset cols
      drop_cols <- c("term", "up_count", "down_count", "sig_count", "consensus")
      setdiff(names(df), intersect(names(df), drop_cols))
    })

    # map column names like "Lees_FA1_0.05_all_genes.rds" -> dataset id "Lees_FA1"
    col_to_dataset_id <- function(x) {

      # Remove trailing pattern "_0.05_all_genes.rds"
      id <- sub("_0\\.[0-9]+_all_genes\\.rds$", "", x)

      # Remove unwanted batchcorrection suffix if present
      id <- sub("_batchcorrection$", "", id)

      # In case someone used plural or uppercase variants (optional but robust)
      id <- sub("_batchCorrection$", "", id, ignore.case = TRUE)

      # Ensure valid ID
      if (!nzchar(id)) return(NA_character_)

      id
    }


    dataset_ids <- reactive({
      ids <- vapply(dataset_cols(), col_to_dataset_id, "")
      unname(ids)
    })

    # ---- pretty labels (reuse your mapping from DEG module) ----
    pretty_map <- tryCatch(
      get("pretty_map", envir = asNamespace(pkg)),
      error = function(e) {
        warning("pretty_map not found in package namespace; using empty mapping.")
        character(0)
      }
    )

    pretty_label <- function(id) pretty_map[[id]] %||% id

    # ---- populate dataset choices ----
    observe({
      ids <- dataset_ids()
      labs <- vapply(ids, pretty_label, character(1))
      updateCheckboxGroupInput(
        session, "gsea_datasets",
        choices = stats::setNames(ids, labs)
      )
    })

    # ---- build term sets per dataset given direction ----
    # value semantics in CSV:  1 = NES up (sig), -1 = NES down (sig), 0 = not sig / not present
    get_terms_for <- function(df, col_name, direction = c("both","up","down")) {
      direction <- match.arg(direction)
      v <- df[[col_name]]
      if (is.null(v)) return(character(0))
      keep <- switch(
        direction,
        "up"   = v ==  1,
        "down" = v == -1,
        "both" = v !=  0
      )
      terms <- df$term[keep]
      terms <- terms[!is.na(terms) & nzchar(terms)]
      unique(terms)
    }

    # ---- reactive sets for selected datasets ----
    gsea_sets <- reactive({
      req(input$gsea_datasets); validate(need(length(input$gsea_datasets) >= 2, "Select at least two datasets."))
      df <- gsea_tbl()
      cols <- dataset_cols()

      # map dataset id -> corresponding CSV column(s)
      id2col <- stats::setNames(cols, vapply(cols, col_to_dataset_id, ""))
      # keep only requested dataset IDs that exist in the table
      want_ids <- intersect(input$gsea_datasets, names(id2col))
      validate(need(length(want_ids) >= 2, "Need at least two datasets with data."))

      dir <- input$gsea_direction %||% "both"
      sets <- lapply(want_ids, function(id) get_terms_for(df, id2col[[id]], dir))
      names(sets) <- vapply(want_ids, pretty_label, "", USE.NAMES = FALSE)
      # drop empties
      sets <- sets[vapply(sets, length, 1L) > 0L]
      validate(need(length(sets) >= 2, "Need at least two non-empty sets."))
      sets
    })

    # ---- helper: overlap table (non-exclusive) ----
    venn_overlap_tbl <- function(s) {
      Universe <- unique(unlist(s, use.names = FALSE))
      M <- vapply(s, function(v) Universe %in% v, logical(length(Universe)))
      colnames(M) <- names(s)
      which_sets <- apply(M, 1, function(r) names(s)[r])
      which_sets <- which_sets[lengths(which_sets) >= 2]
      if (!length(which_sets)) {
        return(data.frame(
          `Dataset Combination`        = character(0),
          `Number of Terms`            = integer(0),
          `Number of Datasets`         = integer(0),
          check.names = FALSE
        ))
      }
      keys <- vapply(which_sets, function(v) paste(sort(v), collapse = " & "), character(1))
      tt <- sort(table(keys), decreasing = TRUE)
      out <- data.frame(
        `Dataset Combination` = names(tt),
        `Number of Terms`     = as.integer(tt),
        check.names = FALSE, row.names = NULL
      )
      out$`Number of Datasets` <- 1 + stringr::str_count(out$`Dataset Combination`, " & ")
      out[order(out$`Number of Datasets`, -out$`Number of Terms`), ]
    }

    # ---- mode message ----
    output$gsea_mode_msg <- shiny::renderUI({
      req(input$gsea_datasets)
      n <- length(input$gsea_datasets)
      if (n <= 6) {
        div(class = "alert alert-success",
            sprintf("Showing Venn diagram for %d datasets. Counts = overlapping terms.", n))
      } else {
        div(class = "alert alert-warning",
            sprintf("You selected %d datasets (> 6). Displaying shared-term tables instead of a Venn diagram.", n))
      }
    })

    # ---- selection info ----
    output$gsea_selection_info <- shiny::renderUI({
      req(input$gsea_datasets)
      tags$span(sprintf("Selected: %d dataset(s). Direction: %s.", length(input$gsea_datasets),
                        input$gsea_direction %||% "both"))
    })

    # ---- filename prefix helper ----
    gsea_fname_prefix <- reactive({
      paste0(
        "gsea_venn_",
        (input$gsea_direction %||% "both"), "_",
        length(input$gsea_datasets), "sets"
      )
    })

    # ---- totals table ----
    gsea_totals_tbl <- reactive({
      s <- gsea_sets()
      data.frame(
        Dataset = names(s),
        `Total Terms` = as.integer(vapply(s, length, 1L)),
        check.names = FALSE
      )
    })

    # ---- overlaps (combination summary) ----
    gsea_overlaps_tbl <- reactive({
      s <- gsea_sets()
      venn_overlap_tbl(s)
    })

    # ---- presence matrix (items) with original signs ----
    # ---- presence matrix (items) with original signs + summary counts ----
    gsea_items_tbl <- reactive({
      s <- gsea_sets();                       # sets after applying current direction filter
      req(length(s) >= 2)

      # universe of displayed terms (respecting direction choice)
      Universe <- unique(unlist(s, use.names = FALSE))

      # original CSV for signs
      df <- gsea_tbl()

      # map selected dataset IDs -> CSV column names
      cols <- dataset_cols()
      id2col <- stats::setNames(cols, vapply(cols, col_to_dataset_id, ""))
      raw_ids   <- input$gsea_datasets
      want_raw  <- intersect(raw_ids, names(id2col))
      want_cols <- unname(id2col[want_raw])

      # subset df to selected terms & datasets (keep sign values 1 / -1 / 0)
      sub <- df[df$term %in% Universe, c("term", want_cols), drop = FALSE]
      for (cn in want_cols) if (!cn %in% names(sub)) sub[[cn]] <- 0L

      # coerce to integer sign matrix
      sign_mat <- as.data.frame(sub[want_cols], stringsAsFactors = FALSE)
      suppressWarnings({
        sign_mat[] <- lapply(sign_mat, function(x) as.integer(as.character(x)))
      })
      sign_mat[is.na(sign_mat)] <- 0L

      # counts across selected datasets
      up_count   <- rowSums(sign_mat ==  1L)
      down_count <- rowSums(sign_mat == -1L)
      sig_count  <- up_count + down_count
      consensus  <- ifelse(sig_count == 0L, "none",
                           ifelse(up_count > 0L & down_count == 0L, "up",
                                  ifelse(down_count > 0L & up_count == 0L, "down", "mixed")))

      # Sum that respects current direction
      dir <- input$gsea_direction %||% "both"
      Sum <- switch(dir,
                    "up"   = up_count,
                    "down" = down_count,
                    "both" = sig_count)

      # rename dataset columns to pretty labels for display
      pretty_names <- stats::setNames(vapply(want_raw, pretty_label, character(1)), want_cols)
      names(sub)[match(want_cols, names(sub))] <- unname(pretty_names)

      # assemble output
      ds_labels <- sort(unname(pretty_names))  # stable order for dataset columns
      out <- sub[, c("term", ds_labels), drop = FALSE]
      names(out)[1] <- "Term"

      out$up_count   <- up_count
      out$down_count <- down_count
      out$sig_count  <- sig_count
      out$consensus  <- consensus
      out$Sum        <- Sum

      # column order: Term, Sum, up/down/sig/consensus, then per-dataset signs
      out <- out[, c("Term", "Sum", "up_count", "down_count", "sig_count", "consensus", ds_labels), drop = FALSE]

      # order rows: Sum desc, then sig_count desc, then Term
      out[order(-out$Sum, -out$sig_count, out$Term), , drop = FALSE]
    })



    # ---- Venn plot (≤6 sets) ----
    gsea_venn_plot_obj <- reactive({
      s <- gsea_sets()
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

    output$gsea_venn_plot <- shiny::renderPlot({
      p <- gsea_venn_plot_obj()
      if (is.null(p)) {
        plot.new(); text(0.5, 0.5, "Overlap tables shown below", cex = 1.6)
      } else {
        suppressWarnings(print(p))
      }
    })

    # ---- enable/disable download buttons when >6 datasets ----
    observe({
      have_plot <- !is.null(gsea_venn_plot_obj())
      if (requireNamespace("shinyjs", quietly = TRUE)) {
        shinyjs::toggleState(ns("dl_gsea_venn_svg"), condition = have_plot)
        shinyjs::toggleState(ns("dl_gsea_venn_png"), condition = have_plot)
      }
    })

    # ---- DT renders ----
    output$gsea_totals <- DT::renderDataTable({
      DT::datatable(
        gsea_totals_tbl(),
        rownames = FALSE,
        options = list(dom = "tip", pageLength = 10, deferRender = TRUE, scrollX = TRUE)
      )
    }, server = TRUE)

    output$gsea_overlaps <- DT::renderDataTable({
      DT::datatable(
        gsea_overlaps_tbl(),
        rownames = FALSE,
        selection = "single",
        options = list(dom = "tip", pageLength = 10, deferRender = TRUE, scrollX = TRUE)
      )
    }, server = TRUE)

    output$gsea_overlap_items <- DT::renderDataTable({
      DT::datatable(
        gsea_items_tbl(),
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
    output$dl_gsea_venn_svg <- downloadHandler(
      filename = function() paste0(gsea_fname_prefix(), ".svg"),
      content = function(file) {
        p <- gsea_venn_plot_obj(); req(p)
        svglite::svglite(file, width = 14, height = 9)
        on.exit(grDevices::dev.off(), add = TRUE)
        print(p)
      }
    )
    output$dl_gsea_venn_png <- downloadHandler(
      filename = function() paste0(gsea_fname_prefix(), ".png"),
      content = function(file) {
        p <- gsea_venn_plot_obj(); req(p)
        ggplot2::ggsave(filename = file, plot = p, width = 14, height = 9, dpi = 300)
      }
    )
    output$dl_gsea_totals_csv <- downloadHandler(
      filename = function() paste0(gsea_fname_prefix(), "_totals.csv"),
      content = function(file) utils::write.csv(gsea_totals_tbl(), file, row.names = FALSE)
    )
    output$dl_gsea_overlaps_csv <- downloadHandler(
      filename = function() paste0(gsea_fname_prefix(), "_overlaps.csv"),
      content = function(file) utils::write.csv(gsea_overlaps_tbl(), file, row.names = FALSE)
    )
    output$dl_gsea_overlap_items_csv <- downloadHandler(
      filename = function() paste0(gsea_fname_prefix(), "_items.csv"),
      content = function(file) utils::write.csv(gsea_items_tbl(), file, row.names = FALSE)
    )
  })
}
