#' @importFrom grid gpar unit
tpmHeatmapSidebarUI <- function(id) {
  ns <- NS(id)
  tagList(
    h4("Heatmap options"),
    radioButtons(
      ns("feature_level"), "Level",
      choices = c("Genes" = "genes", "Isoforms (transcripts)" = "transcripts"),
      selected = "genes"
    ),
    checkboxGroupInput(
      ns("datasets"),
      label = "Datasets (select any number)",
      choices = character(0)
    ),
    uiOutput(ns("datasets_note")),
    textAreaInput(
      ns("feature_query"),
      label = "Genes / transcripts",
      placeholder = "FXN, PIP5K1B, RPS29 …",
      value = "FXN, PIP5K1B, RPS29",
      rows = 4
    ),
    uiOutput(ns("transform_ui")),   # dynamic transform panel
    radioButtons(
      ns("group_filter"), "Groups",
      choices = c("Both" = "both", "CTRL only" = "ctrl", "FRDA only" = "frda"),
      selected = "both"
    ),
    checkboxInput(ns("cluster_rows"), "Cluster rows", value = TRUE),
    checkboxInput(ns("cluster_cols"), "Cluster columns", value = TRUE),
    conditionalPanel(
      condition = sprintf("!input['%s'] && input['%s'].length > 1",
                          ns("cluster_cols"), ns("datasets")),
      selectInput(
        ns("column_order"),
        label = "Column ordering",
        choices = c(
          "Group by FRDA vs CTRL" = "group",
          "Group by dataset"       = "dataset"
        ),
        selected = "dataset"
      )
    ),
    br(),
    downloadButton(ns("dl_svg"), "Download SVG"),
    downloadButton(ns("dl_png"), "Download PNG"),
    uiOutput(ns("plot_notes"))
  )
}

# --- TPM Heatmap: Main UI (plot) ---
tpmHeatmapMainUI <- function(id) {
  ns <- NS(id)
  tagList(
    shinycssloaders::withSpinner(
      uiOutput(ns("heatmap_ui")),
      type = 4, color = "#005249"
    ),
    #add whitespace below the plot 100 px
    br(), br(), br()
  )
}


# --- TPM/VST Heatmap: Server (metadata-aware, single vs multi dataset) ---
tpmHeatmapServer <- function(
    id,
    pkg = utils::packageName(),
    sample_meta = NULL   # pass a data.frame or leave NULL to auto-load from extdata
) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # ---- make pkg robust ----
    pkg <- tryCatch(pkg, error = function(e) "")
    if (!length(pkg) || !is.character(pkg) || !nzchar(pkg)) pkg <- "FRDATranscriptomicAtlas"
    pkg <- pkg[[1L]]

    `%||%` <- function(a, b) if (is.null(a)) b else a

    # ---------------- helpers ----------------
    # -------- Pretty map (internal default) --------
    pretty_map <- tryCatch(
      get("pretty_map", envir = asNamespace(pkg)),
      error = function(e) {
        warning("pretty_map not found; using empty vector.")
        character(0)
      }
    )

    pretty_label <- function(id) pretty_map[[id]] %||% id

    sample_prefix_map <- list(
      Erwin        = "^Erwin_",
      Indelicato   = "^Indelicato_",
      Lai_CNS      = "^Lai_.*_CNS_",
      Lai_PNS      = "^Lai_.*_PNS_",
      Lai_iPSC     = "^Lai_.*_iPSC_",
      Lees         = "^Lees_",
      Maddock_LMN  = "^Maddock_.*LMN_",
      Maddock_SN   = "^Maddock_.*N_",
      Maddock_NCC  = "^Maddock_.*NCC_",
      Mishra       = "^Mishra_",
      Napierala    = "^Napierala_",
      Vilema       = "^Vilema_",
      Wang         = "^Wang_"
    )

    # ---- Define your dataset colour palette ----
    dataset_colors <- c(
      "#00B7C7", "#DC267F", "#FFB000", "#FE6100", "#785EF0",
      "#648FFF", "#00359C", "#009E73", "#E377C2", "#1F77B4",
      "#D62728", "#2CA02C", "#9467BD", "#F564E3", "#A6761D",
      "#BCBD22", "#8C564B"
    )

    harmonize_maddock_names <- function(vsd_names, tpm_names) {

      fixed <- vsd_names

      # --- Fix SN ---
      # N → SN, but only when surrounded by FA…_ … _REP pattern
      sn_pattern1 <- "^Maddock_FA([0-9]+)icN(_REP[0-9]+)$"
      sn_pattern2 <- "^Maddock_FA([0-9]+)N(_REP[0-9]+)$"

      fixed <- sub(sn_pattern1, "Maddock_FA\\1icSN\\2", fixed)
      fixed <- sub(sn_pattern2, "Maddock_FA\\1SN\\2", fixed)

      # --- Fix LMN (NIL → LMN) ---
      lmn_pattern1 <- "^Maddock_FA([0-9]+)icNIL(_REP[0-9]+)$"
      lmn_pattern2 <- "^Maddock_FA([0-9]+)NIL(_REP[0-9]+)$"

      fixed <- sub(lmn_pattern1, "Maddock_FA\\1icLMN\\2", fixed)
      fixed <- sub(lmn_pattern2, "Maddock_FA\\1LMN\\2", fixed)

      # --- ensure no change in vector length & no invalid names ---
      # Replace only when the corrected name actually exists in TPM
      fixed <- ifelse(fixed %in% tpm_names, fixed, vsd_names)

      return(fixed)
    }


    base_of <- function(dataset_id) sub("(_.*)$", "", dataset_id)

    .auto_downshift <- function(mat) {
      nr <- nrow(mat); nc <- ncol(mat)
      too_many_cells <- (nr * nc) > 1e6
      too_many_rows  <- nr > 1500
      too_many_cols  <- nc > 800

      if (too_many_cells || too_many_rows) {
        if (isTRUE(input$cluster_rows))
          updateCheckboxInput(session, "cluster_rows", value = FALSE)
      }
      if (too_many_cells || too_many_cols) {
        if (isTRUE(input$cluster_cols))
          updateCheckboxInput(session, "cluster_cols", value = FALSE)
      }
    }

    # extdata path resolver (works for installed pkg or dev tree)
    extdata_path <- function(..., package = utils::packageName()) {
      p <- system.file("extdata", package = package)
      if (nzchar(p)) return(file.path(p, ...))
      file.path(getwd(), "inst", "extdata", ...)
    }

    # ---- TPM loader (your original files) ----
    load_one_tpm <- function(dataset_id, level = c("genes","transcripts")) {
      level  <- match.arg(level)
      subdir <- if (level == "genes") "tpm" else "transcript_tpm"
      fname  <- if (level == "genes")
        paste0(dataset_id, "_gene_tpm.rds") else paste0(dataset_id, "_transcript_tpm.rds")
      path <- extdata_path(subdir, fname, package = pkg)
      if (!file.exists(path)) stop("Missing TPM file: ", path)
      readRDS(path)
    }

    # ---- VST loader from existing *_vsd.rds ----
    load_vsd <- function(dataset_id) {

      ddir <- extdata_path("deseq_objects", package = pkg)

      # Load ALL matching VSDs, not just one
      pat <- paste0("^", dataset_id, ".*_vsd\\.rds$")
      files <- list.files(ddir, pattern = pat, full.names = TRUE)

      if (!length(files))
        stop("No VST object (*.vsd.rds) found for dataset group: ", dataset_id)

      # Load each VST, return a list
      lapply(files, readRDS)
    }


    # feature query parser
    parse_feature_query <- function(txt) {
      if (is.null(txt) || !nzchar(txt)) return(character(0))
      toks <- unique(unlist(strsplit(txt, "[,\\s]+", perl = TRUE)))
      toks[nzchar(toks)]
    }

    # metadata loading / normalization
    load_sample_meta <- function(path = extdata_path("metadata", "samples.csv")) {
      stopifnot(file.exists(path))
      readr::read_csv(path, show_col_types = FALSE)
    }

    # ==== Group token mappers ====
    map_to_group_meta <- function(x) {
      x <- toupper(trimws(x))
      is_ctrl <- grepl("(?:^|[_-])(CTRL|IC)(?:[_-]|$)", x, perl = TRUE) ||
        grepl("FA\\d*IC", x, perl = TRUE) ||
        grepl("FAIC", x, perl = TRUE)
      is_frda <- grepl("FRDA\\d*", x, perl = TRUE) ||
        grepl("(?:(?:^|[_-])FA\\d*)(?!IC)", x, perl = TRUE) ||
        grepl("(?:(?:^|[_-])FA)(?!IC)", x, perl = TRUE)
      if (is_ctrl) "CTRL" else if (is_frda) "FRDA" else NA_character_
    }

    map_to_group_name <- function(x) {
      x <- toupper(x)
      is_ctrl <- grepl("(?:^|[_-])(CTRL|IC)(?:[_-]|$)", x, perl = TRUE) ||
        grepl("FA\\d*IC", x, perl = TRUE) ||
        grepl("FAIC", x, perl = TRUE)
      is_frda <- grepl("FRDA\\d*", x, perl = TRUE) ||
        grepl("(?:(?:^|[_-])FA\\d*)(?!IC)", x, perl = TRUE) ||
        grepl("(?:(?:^|[_-])FA)(?!IC)", x, perl = TRUE)
      if (is_ctrl) "CTRL" else if (is_frda) "FRDA" else NA_character_
    }

    normalise_meta <- function(df) {
      stopifnot(all(c("sample_id","Study") %in% names(df)))
      df$sample_id <- gsub("\\s+", "", trimws(df$sample_id))
      df$Study     <- trimws(as.character(df$Study))

      if ("FRDA_CTRL" %in% names(df)) {
        grp <- vapply(df$FRDA_CTRL, map_to_group_meta, character(1))
      } else if ("case_diff_controls" %in% names(df)) {
        grp <- vapply(df$case_diff_controls, map_to_group_meta, character(1))
      } else {
        grp <- NA_character_
      }

      df$Group <- factor(grp, levels = c("CTRL","FRDA"))
      df
    }

    group_from_name <- function(x) map_to_group_name(x)

    annot_from_meta <- function(cols, meta_norm) {
      out <- data.frame(sample = cols, stringsAsFactors = FALSE)

      if (!is.null(meta_norm)) {
        out <- dplyr::left_join(
          out,
          dplyr::select(meta_norm, sample_id, Study, Group),
          by = c("sample" = "sample_id")
        )
      }

      if (!"Study" %in% names(out) || anyNA(out$Study)) {
        out$Study <- if (!"Study" %in% names(out)) sub("_.*$", "", out$sample) else
          ifelse(is.na(out$Study), sub("_.*$", "", out$sample), out$Study)
      }

      if (!"Group" %in% names(out) || anyNA(out$Group)) {
        fallback <- vapply(out$sample, group_from_name, character(1))
        out$Group <- if (!"Group" %in% names(out)) fallback else
          ifelse(is.na(out$Group), fallback, as.character(out$Group))
      }

      out$Group <- factor(out$Group, levels = c("CTRL","FRDA"))
      out
    }

    # choose columns for a dataset with group filter (uses metadata if present)
    columns_for_dataset <- function(df, dataset_id, group_mode, meta_norm) {

      # 1. Detect Maddock subgroup datasets
      if (dataset_id == "Maddock_SN") {
        sample_cols <- grep("SN", names(df), value = TRUE)
      } else if (dataset_id == "Maddock_LMN") {
        sample_cols <- grep("LMN", names(df), value = TRUE)
      } else if (dataset_id == "Maddock_NCC") {
        sample_cols <- grep("NCC", names(df), value = TRUE)
      } else {
        # 2. Default behaviour: prefix match
        prefix <- sample_prefix_map[[dataset_id]]

        if (is.null(prefix)) {
          # fallback to old behaviour if needed
          sample_cols <- grep(paste0("^", dataset_id, "_"), names(df), value = TRUE)
        } else {
          sample_cols <- grep(prefix, names(df), value = TRUE)
        }

      }

      # remove non-sample columns
      id_cols <- intersect(c("tx","gene_id","gene_name"), names(df))
      sample_cols <- setdiff(sample_cols, id_cols)

      if (!length(sample_cols)) return(character(0))

      # metadata filtering stays the same
      ann <- annot_from_meta(sample_cols, meta_norm)
      keep <- switch(group_mode,
                     "ctrl" = ann$Group == "CTRL",
                     "frda" = ann$Group == "FRDA",
                     "both" = ann$Group %in% c("CTRL","FRDA"),
                     rep(FALSE, nrow(ann)))
      ann$sample[keep]
    }


    row_zscore <- function(mat) {
      m <- t(scale(t(mat)))
      m[!is.finite(m)] <- 0
      m
    }

    # Per-dataset row Z-score on VST matrix
    zscore_by_dataset <- function(mat, dataset_ids) {
      out <- mat
      for (ds in unique(dataset_ids)) {
        cols <- which(dataset_ids == ds)
        out[, cols] <- t(scale(t(out[, cols, drop = FALSE])))
      }
      out[!is.finite(out)] <- 0
      out
    }

    # dataset choices: show pretty labels but return ids
    observe({
      updateCheckboxGroupInput(
        session, "datasets",
        choices  = setNames(names(pretty_map), pretty_map),
        selected = character(0)
      )
    })

    output$datasets_note <- renderUI({
      req(input$datasets)
      labels <- unname(pretty_map[input$datasets])
      htmltools::HTML(
        sprintf("<small>Selected: <b>%s</b></small>", paste(labels, collapse = ", "))
      )
    })

    # normalize metadata (if provided or auto-loaded)
    if (is.null(sample_meta)) {
      sample_meta <- tryCatch(load_sample_meta(), error = function(e) NULL)
    }
    meta_norm <- if (!is.null(sample_meta)) normalise_meta(sample_meta) else NULL

    # ---- Transform UI: TPM options for 1 dataset, VST for >1 ----
    output$transform_ui <- renderUI({
      if (length(input$datasets) <= 1) {
        radioButtons(
          ns("transform_mode"), "Transform",
          choices  = c("log2(TPM+1)" = "log2p1",
                       "Row Z-score of log2(TPM+1)" = "zscore"),
          selected = "log2p1"
        )
      } else {
        radioButtons(
          ns("transform_mode"), "Transform",
          choices  = c("Z-score VST" = "zvst"),
          selected = "zvst"
        )
      }
    })

    # --- utility: infer dataset from sample name ---
    dataset_from_sample <- function(samples, selected_ids) {
      vapply(samples, function(s) {
        hit <- selected_ids[startsWith(s, paste0(selected_ids, "_"))]
        if (length(hit)) hit[1] else sub("_.*$", "", s)
      }, character(1))
    }

    # ---------------- build matrix ----------------
    unified_matrix <- reactive({
      req(input$datasets)
      order_mode <- input$column_order

      q <- parse_feature_query(input$feature_query)
      validate(need(length(q) > 0, "Enter at least one gene symbol or gene ID."))

      level <- input$feature_level %||% "genes"

      # ---- SINGLE DATASET: TPM mode (original behaviour) ----
      if (length(input$datasets) == 1) {

        ds <- input$datasets
        df <- load_one_tpm(ds, level = level)

        if (level == "genes") {
          id_col <- "gene_id"; name_col <- "gene_name"
        } else {
          id_col <- "tx";      name_col <- "gene_id"
        }

        sc <- columns_for_dataset(df, ds, input$group_filter %||% "both", meta_norm)
        if (!length(sc))
          validate("No samples available for the current group filter.")

        has_name <- name_col %in% names(df)
        by_id    <- df[[id_col]] %in% q
        by_sym   <- if (has_name) toupper(df[[name_col]]) %in% toupper(q) else FALSE
        keep_rows <- by_id | by_sym

        sub <- df[keep_rows, c(id_col, if (has_name) name_col, sc), drop = FALSE]
        validate(need(nrow(sub) > 0, "No matching genes in this dataset."))

        key <- if (has_name) {
          ifelse(nzchar(sub[[name_col]]),
                 sub[[name_col]],
                 sub[[id_col]])
        } else sub[[id_col]]

        mat <- as.matrix(sub[, sc, drop = FALSE])
        storage.mode(mat) <- "double"
        rownames(mat) <- key

        # transform
        if ((input$transform_mode %||% "log2p1") == "log2p1") {
          mat <- log2(mat + 1)
        } else {
          mat <- row_zscore(log2(mat + 1))
        }

        # when not clustering columns: order Study then CTRL→FRDA→UNKNOWN
        # optional ordered columns if no clustering
        if (!isTRUE(input$cluster_cols)) {

          ann <- annot_from_meta(colnames(mat), meta_norm)

          # Dataset rank: preserve user-selected order
          ds_order <- setNames(seq_along(input$datasets), input$datasets)
          ann$ds_rank <- ds_order[ann$Study] %||% (max(ds_order, na.rm = TRUE) + 1)

          # Group rank
          ann$grp_rank <- dplyr::recode(as.character(ann$Group),
                                        "CTRL" = 0L,
                                        "FRDA" = 1L,
                                        .default = 2L)

          # Apply selected ordering
          if (order_mode == "group") {
            # All CTRL across all datasets, then FRDA
            ann <- dplyr::arrange(ann, grp_rank, ds_rank, sample)
          } else {
            # Group by dataset block
            ann <- dplyr::arrange(ann, ds_rank, grp_rank, sample)
          }

          mat <- mat[, ann$sample, drop = FALSE]
        }


        return(mat)
      }

      # ---- MULTI-DATASET: VST + per-dataset Z-score ----
      validate(need(level == "genes",
                    "Multi-dataset mode currently supports genes only."))

      long_list <- list()
      all_keys  <- NULL

      for (ds in input$datasets) {

        # TPM for annotation + sample selection
        tpm_df <- load_one_tpm(ds, level = "genes")
        sc <- columns_for_dataset(tpm_df, ds, input$group_filter %||% "both", meta_norm)
        if (!length(sc)) next

        # VST values
        # VST mode: load ALL VSD objects for this umbrella dataset
        vsd_list <- load_vsd(ds)

        # Loop through each sub-dataset VST file (e.g., Lees_FA1, Lees_FA2, Lees_FA3)
        for (vsd_obj in vsd_list) {

          vsd_mat <- SummarizedExperiment::assay(vsd_obj)

          # Determine samples (columns) to keep
          # Attempt direct intersection first
          sc_use <- intersect(sc, colnames(vsd_mat))

          # If nothing matches, attempt safe harmonization
          if (!length(sc_use) && grepl("^Maddock", ds)) {

            new_names <- harmonize_maddock_names(colnames(vsd_mat), sc)

            # Reassign ONLY if the length matches and no duplicates
            if (length(new_names) == length(colnames(vsd_mat)) &&
                length(unique(new_names)) == length(new_names)) {

              colnames(vsd_mat) <- new_names
              sc_use <- intersect(sc, new_names)
            }
          }

          if (!length(sc_use)) next

          # Identify gene IDs of interest
          annot <- tpm_df[, c("gene_id", "gene_name")]
          wanted_ids <- unique(c(
            annot$gene_id[annot$gene_id %in% q],
            annot$gene_id[toupper(annot$gene_name) %in% toupper(q)]
          ))

          # Match rows
          keep_ids <- intersect(rownames(vsd_mat), wanted_ids)
          if (!length(keep_ids)) next

          sub_mat <- vsd_mat[keep_ids, sc_use, drop = FALSE]

          # Row names → pretty key (gene_name (gene_id))
          map_sub <- annot[match(keep_ids, annot$gene_id), ]
          key <- ifelse(
            !is.na(map_sub$gene_name) & nzchar(map_sub$gene_name),
            paste0(map_sub$gene_name, " (", map_sub$gene_id, ")"),
            map_sub$gene_id
          )

          # Long-format entry for this ONE vsd file
          lng <- as.data.frame(sub_mat)
          colnames(lng) <- sc_use
          lng$feature <- key

          lng <- tidyr::pivot_longer(
            lng,
            cols      = dplyr::all_of(sc_use),
            names_to  = "sample",
            values_to = "expr"
          )

          long_list[[length(long_list) + 1]] <- lng
          all_keys <- union(all_keys, unique(key))
        }

}
      validate(need(length(long_list) > 0,
                    "No matching genes/samples for the current filters."))

      all_long <- dplyr::bind_rows(long_list)
      all_long$feature <- factor(all_long$feature, levels = all_keys)

      wide <- tidyr::pivot_wider(
        all_long,
        id_cols    = "feature",
        names_from = "sample",
        values_from = "expr"
      )
      wide <- dplyr::arrange(wide, feature)

      mat <- as.matrix(wide[, -1, drop = FALSE])
      storage.mode(mat) <- "double"
      rownames(mat) <- wide$feature

      # per-dataset row Z-score
      ds_ids <- dataset_from_sample(colnames(mat), input$datasets)
      mat <- zscore_by_dataset(mat, ds_ids)

      # optional ordered columns if no clustering
      # optional ordered columns if no clustering
      if (!isTRUE(input$cluster_cols)) {
        if (!is.null(meta_norm)) {

          ann <- annot_from_meta(colnames(mat), meta_norm)

          # Determine dataset for each sample (Study is often correct,
          # but fallback to prefix-based detection if mismatched)
          ds_order <- setNames(seq_along(input$datasets), input$datasets)

          ann$Dataset <- ann$Study
          ann$Dataset[is.na(ann$Dataset)] <- sub("_.*$", "", ann$sample)

          ann$ds_rank <- ds_order[ann$Dataset]
          ann$grp_rank <- ifelse(ann$Group == "CTRL", 0L, 1L)

          # ---- ORDERING LOGIC ----
          if (input$column_order == "group") {
            # Option A: All CTRL first (across datasets), then FRDA
            ann <- ann[order(ann$grp_rank, ann$ds_rank, ann$sample), ]
          } else {
            # Option B: Dataset blocks first, then CTRL→FRDA inside each dataset
            ann <- ann[order(ann$ds_rank, ann$grp_rank, ann$sample), ]
          }

          mat <- mat[, ann$sample, drop = FALSE]
        }
      }


      mat
    })

    # Dynamically size the graphics device by matrix dimensions
    plot_dims <- reactive({
      mat <- unified_matrix()
      nr <- nrow(mat); nc <- ncol(mat)

      cell_h_pt <- 16
      extra_pad <- if (nr < 30) 700 else if (nr < 60) 800 else if (nr > 100) 1200 else 1000

      list(
        dev_h_px = max(500, round(nr * (cell_h_pt * 96/72)) + extra_pad),
        dev_w_px = NULL,
        cell_h_pt = cell_h_pt,
        show_row_names = nr <= 120,
        show_col_names = nc <= 120,
        row_cex = if (nr <= 80) 10 else if (nr <= 150) 8 else 6,
        col_cex = if (nc <= 60) 10 else if (nc <= 120) 8 else 6
      )
    })

    # keep last heatmap for downloads
    .last_ht <- reactiveVal(NULL)

    # ---------------- render plot ----------------
    output$heatmap_plot <- renderPlot({
      mat  <- unified_matrix()
      .auto_downshift(mat)
      dims <- plot_dims()

      ann_df <- {
        ds_id <- dataset_from_sample(colnames(mat), input$datasets)
        grp   <- vapply(colnames(mat), group_from_name, character(1))
        data.frame(Dataset = ds_id,
                   Group = factor(grp, levels = c("CTRL","FRDA")),
                   row.names = colnames(mat), check.names = FALSE)
      }

      ds_levels <- unique(ann_df$Dataset)
      ds_pal <- setNames(dataset_colors[seq_along(ds_levels) %% length(dataset_colors)],
                         ds_levels)

      present <- unique(as.character(ann_df$Group))
      lvls <- intersect(c("CTRL", "FRDA"), present)
      group_cols <- setNames(c("#a9a9a9ff", "#333333ff")[match(lvls, c("CTRL","FRDA"))], lvls)

      ha_top <- ComplexHeatmap::HeatmapAnnotation(
        Dataset = ann_df$Dataset,
        Group   = ann_df$Group,
        col = list(Dataset = ds_pal, Group = group_cols),
        annotation_legend_param = list(
          Dataset = list(title = "Dataset"),
          Group   = list(title = "Group")
        )
      )

      rng <- range(mat, na.rm = TRUE)
      if (length(input$datasets) == 1) {
        col_fun <- circlize::colorRamp2(c(rng[1], rng[2]), c("white", "#030058"))
        legend_title <- if ((input$transform_mode %||% "log2p1") == "log2p1")
          "log2(TPM+1)" else "Row Z-score"
      } else {
        col_fun <- circlize::colorRamp2(c(-3, 0, 3), c("#2166AC", "white", "#B2182B"))
        legend_title <- "Z-score (VST)"
      }

      cl_rows <- isTRUE(input$cluster_rows)
      cl_cols <- isTRUE(input$cluster_cols)

      ht <- ComplexHeatmap::Heatmap(
        mat,
        name = legend_title,
        col  = col_fun,
        top_annotation = ha_top,
        show_row_dend     = cl_rows,
        show_column_dend  = cl_cols,
        cluster_rows      = cl_rows,
        cluster_columns   = cl_cols,
        row_names_side    = "left",
        column_names_side = "top",
        show_row_names    = dims$show_row_names,
        show_column_names = dims$show_col_names,
        row_names_gp      = grid::gpar(fontsize = dims$row_cex),
        column_names_gp   = grid::gpar(fontsize = dims$col_cex),
        na_col            = "grey80",
        border            = TRUE,
        use_raster        = TRUE
      )

      ComplexHeatmap::draw(
        ht,
        heatmap_legend_side    = "right",
        annotation_legend_side = "right",
        padding = grid::unit(c(6, 20, 20, 6), "mm")
      )

      .last_ht(ht)
      gc()
    }, res = 120)

    output$heatmap_ui <- renderUI({
      dims <- plot_dims()
      h_px <- paste0(round(dims$dev_h_px), "px")
      div(
        style = "overflow: visible;",
        plotOutput(ns("heatmap_plot"), height = h_px, width = "100%")
      )
    })

    output$plot_notes <- renderUI({
      if (length(input$datasets) <= 1) {
        htmltools::HTML("<small>Using TPM (log2(TPM+1) or row Z-score). Interpretation is within-dataset only.</small>")
      } else {
        htmltools::HTML("<small>Using DESeq2 VST values with per-dataset row Z-score (appropriate for cross-dataset comparison).</small>")
      }
    })

    # ---------------- downloads ----------------
    output$dl_svg <- downloadHandler(
      filename = function() paste0("heatmap_", Sys.Date(), ".svg"),
      content = function(file) {
        dims <- plot_dims()
        w_in <- (session$clientData[[paste0("output_", ns("heatmap_plot"), "_width")]] %||% 1000) / 96
        h_in <- (dims$dev_h_px %||% 800) / 96
        svglite::svglite(file, width = w_in, height = h_in)
        on.exit(grDevices::dev.off(), add = TRUE)
        ComplexHeatmap::draw(
          .last_ht(),
          heatmap_legend_side    = "right",
          annotation_legend_side = "right",
          padding = grid::unit(c(6,10,16,6), "mm")
        )
      }
    )

    output$dl_png <- downloadHandler(
      filename = function() paste0("heatmap_", Sys.Date(), ".png"),
      content = function(file) {
        dims <- plot_dims()
        w_px <- session$clientData[[paste0("output_", ns("heatmap_plot"), "_width")]] %||% 1000
        h_px <- dims$dev_h_px %||% 800
        ragg::agg_png(file, width = w_px, height = h_px, units = "px", res = 150)
        on.exit(grDevices::dev.off(), add = TRUE)
        ComplexHeatmap::draw(
          .last_ht(),
          heatmap_legend_side    = "right",
          annotation_legend_side = "right",
          padding = grid::unit(c(6,10,16,6), "mm")
        )
      }
    )

  })
}
