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
      label = "Genes / transcripts (symbols or IDs; comma/space/newline separated)",
      placeholder = "e.g. FXN, PIP5K1B, RPS29 or ENST00000…",
      value = "FXN, PIP5K1B, RPS29",
      rows = 4
    ),
    radioButtons(
      ns("transform_mode"), "Transform",
      choices = c("log2(TPM + 1)" = "log2p1", "Z-score" = "zscore"),
      selected = "log2p1", inline = TRUE
    ),
    radioButtons(
      ns("group_filter"), "Groups to include",
      choices = c("Both" = "both", "CTRL only" = "ctrl", "FRDA only" = "frda"),
      selected = "both", inline = TRUE
    ),
    checkboxInput(ns("cluster_rows"), "Cluster rows (genes)", value = TRUE),
    checkboxInput(ns("cluster_cols"), "Cluster columns (samples)", value = TRUE),
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

# --- TPM Heatmap: Server module ---
# --- TPM Heatmap: Server (static, metadata-aware) ---
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

    # ---------------- helpers ----------------
    # -------- Pretty map (internal default) --------
    pretty_map  <- c(
      "Erwin"        = "Erwin (Lymphoblastoid Cells)",
      "Indelicato"   = "Indelicato (Skeletal Muscle)",
      "Lai_iPSC"     = "Lai (iPSCs)",
      "Lai_CNS"      = "Lai (CNS neurons)",
      "Lai_PNS"      = "Lai (PNS neurons)",
      "Lees"         = "Lees (Cardiomyocytes)",
      "Maddock_LMN"  = "Maddock (Lower Motor Neurons)",
      "Maddock_SN"   = "Maddock (Sensory Neurons)",
      "Maddock_NCC"  = "Maddock (Neural Crest Cells)",
      "Mishra"       = "Mishra (Neurons)",
      "Napierala"    = "Napierala (Fibroblasts)",
      "Vilema"       = "Vilema-Enriquez (Fibroblasts)"
    )
    pretty_map <- pretty_map %||% PRETTY_MAP_LOCAL

    # ---- Define your dataset colour palette ----
    dataset_colors <- c(
      # your original order
      "#00B7C7",  # teal–cyan
      "#DC267F",  # magenta–pink
      "#FFB000",  # golden yellow
      "#FE6100",  # orange
      "#785EF0",  # purple
      "#648FFF",  # blue
      "#00359C",  # navy blue
      "#009E73",  # emerald green
      "#E377C2",  # rose pink
      "#1F77B4",  # sky blue
      "#D62728",  # red
      "#2CA02C",  # forest green
      "#9467BD",  # muted violet
      "#F564E3",  # fuchsia
      "#A6761D",  # ochre–brown
      "#BCBD22",  # olive
      "#8C564B"   # taupe brown
    )


    `%||%` <- function(a, b) if (is.null(a)) b else a
    base_of <- function(dataset_id) sub("(_.*)$", "", dataset_id)

    .auto_downshift <- function(mat) {
      nr <- nrow(mat); nc <- ncol(mat)
      # thresholds you can tune
      too_many_cells <- (nr * nc) > 1e6
      too_many_rows  <- nr > 1500
      too_many_cols  <- nc > 800

      if (too_many_cells || too_many_rows) {
        # turn off row clustering
        if (isTRUE(input$cluster_rows)) updateCheckboxInput(session, "cluster_rows", value = FALSE)
      }
      if (too_many_cells || too_many_cols) {
        # turn off column clustering
        if (isTRUE(input$cluster_cols)) updateCheckboxInput(session, "cluster_cols", value = FALSE)
      }
    }


    # extdata path resolver (works for installed pkg or dev tree)
    extdata_path <- function(..., package = utils::packageName()) {
      p <- system.file("extdata", package = package)
      if (nzchar(p)) return(file.path(p, ...))
      file.path(getwd(), "inst", "extdata", ...)
    }

    load_one <- function(dataset_id, level = c("genes","transcripts")) {
      level  <- match.arg(level)
      subdir <- if (level == "genes") "tpm" else "transcript_tpm"
      fname  <- if (level == "genes")
        paste0(dataset_id, "_gene_tpm.rds") else paste0(dataset_id, "_transcript_tpm.rds")
      path <- extdata_path(subdir, fname, package = pkg)
      if (!file.exists(path)) stop("Missing file: ", path)
      readRDS(path)
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

    # 1) Metadata mapper (accepts messy tokens like "FA1CM", "FA2IC", "FRDA2", etc.)
    # Precedence: CTRL first so FA2IC/FAIC are CTRL, not FRDA.
    map_to_group_meta <- function(x) {
      x <- toupper(trimws(x))

      is_ctrl <- grepl("(?:^|[_-])(CTRL|IC)(?:[_-]|$)", x, perl = TRUE) ||
        grepl("FA\\d*IC", x, perl = TRUE) ||     # FA1IC, FA2IC, FAIC…
        grepl("FAIC", x, perl = TRUE)

      # FRDA if FRDA[digits] …or FA[digits] NOT followed by IC (supports suffixes like CM/LMN)
      is_frda <- grepl("FRDA\\d*", x, perl = TRUE) ||
        grepl("(?:(?:^|[_-])FA\\d*)(?!IC)", x, perl = TRUE) ||
        grepl("(?:(?:^|[_-])FA)(?!IC)", x, perl = TRUE)

      if (is_ctrl) "CTRL" else if (is_frda) "FRDA" else NA_character_
    }

    # 2) Sample-name fallback (lenient; finds tokens anywhere in the string)
    map_to_group_name <- function(x) {
      x <- toupper(x)

      is_ctrl <- grepl("(?:^|[_-])(CTRL|IC)(?:[_-]|$)", x, perl = TRUE) ||
        grepl("FA\\d*IC", x, perl = TRUE) ||     # catches FA1ICCM, FA2IC_LMN, etc.
        grepl("FAIC", x, perl = TRUE)

      is_frda <- grepl("FRDA\\d*", x, perl = TRUE) ||
        grepl("(?:(?:^|[_-])FA\\d*)(?!IC)", x, perl = TRUE) ||
        grepl("(?:(?:^|[_-])FA)(?!IC)", x, perl = TRUE)

      if (is_ctrl) "CTRL" else if (is_frda) "FRDA" else NA_character_
    }

    # ==== Normalise metadata ====

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

    # ==== Fallback from sample names ====

    group_from_name <- function(x) {
      map_to_group_name(x)
    }




    # --- ANNOTATE A VECTOR OF SAMPLE COLS USING META (WITH FALLBACKS) ---
    annot_from_meta <- function(cols, meta_norm) {
      out <- data.frame(sample = cols, stringsAsFactors = FALSE)

      if (!is.null(meta_norm)) {
        out <- dplyr::left_join(
          out,
          dplyr::select(meta_norm, sample_id, Study, Group),
          by = c("sample" = "sample_id")
        )
      }

      # Fill missing Study and Group
      if (!"Study" %in% names(out) || anyNA(out$Study))
        out$Study <- if (!"Study" %in% names(out)) sub("_.*$", "", out$sample) else
          ifelse(is.na(out$Study), sub("_.*$", "", out$sample), out$Study)

      if (!"Group" %in% names(out) || anyNA(out$Group)) {
        fallback <- vapply(out$sample, group_from_name, character(1))
        out$Group <- if (!"Group" %in% names(out)) fallback else ifelse(is.na(out$Group), fallback, as.character(out$Group))
      }

      # Keep only CTRL/FRDA; coerce to factor with those two levels
      out$Group <- factor(out$Group, levels = c("CTRL","FRDA"))
      out
    }


    # choose columns for a dataset with group filter (uses metadata if present)
    columns_for_dataset <- function(df, dataset_id, group_mode, meta_norm) {
      sample_cols <- grep(paste0("^", dataset_id, "_"), names(df), value = TRUE)
      if (!length(sample_cols)) {
        id_cols <- intersect(c("tx","gene_id","gene_name"), names(df))
        sample_cols <- setdiff(names(df), id_cols)
        sample_cols <- sample_cols[grepl("_", sample_cols)]
      }
      if (!length(sample_cols)) return(character(0))

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

    # dataset choices: show pretty labels but return ids
    dataset_choices <- setNames(names(pretty_map), pretty_map)
    observe({
      updateCheckboxGroupInput(session, "datasets",
                               choices = setNames(names(pretty_map), pretty_map),
                               selected = character(0)
      )
    })
    output$datasets_note <- renderUI({
      req(input$datasets)
      labels <- unname(pretty_map[input$datasets])
      htmltools::HTML(sprintf("<small>Selected: <b>%s</b></small>", paste(labels, collapse = ", ")))
    })

    # normalize metadata (if provided or auto-loaded)
    if (is.null(sample_meta)) {
      sample_meta <- tryCatch(load_sample_meta(), error = function(e) NULL)
    }
    meta_norm <- if (!is.null(sample_meta)) normalise_meta(sample_meta) else NULL

    # ---------------- build matrix ----------------
    unified_matrix <- reactive({
      req(input$datasets)
      q <- parse_feature_query(input$feature_query)
      validate(need(length(q) > 0, "Enter at least one gene symbol, gene ID, or transcript ID."))

      level <- input$feature_level %||% "genes"
      if (level == "genes") {
        id_col <- "gene_id"; name_col <- "gene_name"
      } else {
        id_col <- "tx";      name_col <- "gene_id"  # allow matching by gene_id too
      }

      long_list <- list()
      all_keys  <- NULL

      for (ds in input$datasets) {
        df <- load_one(ds, level = level)
        sc <- columns_for_dataset(df, ds, input$group_filter %||% "both", meta_norm)
        if (!length(sc)) next

        has_name <- name_col %in% names(df)
        by_id    <- df[[id_col]] %in% q
        by_sym   <- if (has_name) toupper(df[[name_col]]) %in% toupper(q) else rep(FALSE, nrow(df))
        keep_rows <- by_id | by_sym

        sub <- df[keep_rows, c(id_col, if (has_name) name_col, sc), drop = FALSE]
        if (!nrow(sub)) next

        key <- if (has_name) {
          ifelse(nzchar(sub[[name_col]]),
                 paste0(sub[[name_col]], " (", sub[[id_col]], ")"),
                 sub[[id_col]])
        } else sub[[id_col]]

        lng <- tidyr::pivot_longer(
          sub, cols = dplyr::all_of(sc),
          names_to = "sample", values_to = "TPM"
        )
        lng$feature <- rep(key, each = length(sc))
        long_list[[ds]] <- lng
        all_keys <- union(all_keys, unique(key))
      }

      validate(need(length(long_list) > 0, "No matching rows/samples for the current filters."))

      all_long <- dplyr::bind_rows(long_list)
      all_long$feature <- factor(all_long$feature, levels = all_keys)

      wide <- tidyr::pivot_wider(
        all_long, id_cols = "feature",
        names_from = "sample", values_from = "TPM"
      )
      wide <- dplyr::arrange(wide, feature)

      mat <- as.matrix(wide[, -1, drop = FALSE])
      storage.mode(mat) <- "double"
      rownames(mat) <- wide$feature

      if ((input$transform_mode %||% "log2p1") == "log2p1") {
        mat <- log2(mat + 1)
      } else {
        mat <- row_zscore(mat)
      }

      # when not clustering columns: order Study (per selected datasets) then CTRL→FRDA→UNKNOWN
      if (!isTRUE(input$cluster_cols)) {
        if (!is.null(meta_norm)) {
          ann <- annot_from_meta(colnames(mat), meta_norm)
          ds_order <- setNames(seq_along(input$datasets), input$datasets)
          ann$ds_rank  <- ds_order[ann$Study] %||% (max(ds_order, na.rm = TRUE) + 1)
          ann$grp_rank <- ifelse(ann$Group == "CTRL", 0L,
                                 ifelse(ann$Group == "FRDA", 1L, 2L))
          ann <- dplyr::arrange(ann, ds_rank, grp_rank, sample)
          mat <- mat[, ann$sample, drop = FALSE]
        } else {
          # regex fallback
          tmp <- data.frame(
            sample = colnames(mat),
            Study  = sub("_.*$", "", colnames(mat)),
            Group  = factor(vapply(colnames(mat), group_from_name, character(1)),
                            levels = c("CTRL","FRDA"))
          )
          ds_order <- setNames(seq_along(input$datasets), input$datasets)
          tmp$ds_rank  <- ds_order[tmp$Study] %||% (max(ds_order, na.rm = TRUE) + 1)
          tmp$grp_rank <- ifelse(tmp$Group == "CTRL", 0L,
                                 ifelse(tmp$Group == "FRDA", 1L, 2L))
          tmp <- dplyr::arrange(tmp, ds_rank, grp_rank, sample)
          mat <- mat[, tmp$sample, drop = FALSE]
        }
      }

      mat
    })

    # Dynamically size the graphics device by matrix dimensions
    plot_dims <- reactive({
      mat <- unified_matrix()
      nr <- nrow(mat); nc <- ncol(mat)

      # per-cell sizing (tweak to taste)
      cell_h_pt <- 16  # points per gene row
      cell_w_px <- 20  # pixels per sample col

      extra_pad <- if (nr < 30) {
        700
      } else if (nr < 60) {
        800
      } else if (nr > 100) {
        1200
      } else {
        1000
      }


      list(
        dev_h_px = max(500, round(nr * (cell_h_pt * 96/72)) + extra_pad), # device height (px); 72pt=96px
        dev_w_px = NULL,                                # let it fill width
        cell_h_pt = cell_h_pt,                          # pass to Heatmap 'height'
        show_row_names = nr <= 120,
        show_col_names = nc <= 120,
        row_cex = if (nr <= 80) 10 else if (nr <= 150) 8 else 6,
        col_cex = if (nc <= 60) 10 else if (nc <= 120) 8 else 6
      )
    })




    # keep last heatmap for downloads
    .last_ht <- reactiveVal(NULL)

    dataset_from_sample <- function(samples, selected_ids) {
      vapply(samples, function(s) {
        hit <- selected_ids[startsWith(s, paste0(selected_ids, "_"))]
        if (length(hit)) hit[1] else sub("_.*$", "", s)
      }, character(1))
    }

    # ---------------- render plot ----------------
    output$heatmap_plot <- renderPlot({
      mat  <- unified_matrix()
      .auto_downshift(mat)  # may toggle the cluster checkboxes
      dims <- plot_dims()

      ann_df <- {
        ds_id <- dataset_from_sample(colnames(mat), input$datasets)
        grp   <- vapply(colnames(mat), group_from_name, character(1))
        data.frame(Dataset = ds_id,
                   Group = factor(grp, levels = c("CTRL","FRDA")),
                   row.names = colnames(mat), check.names = FALSE)
      }

      # keep palettes as-is; they'll now key off proper character labels
      ds_levels <- unique(ann_df$Dataset)
      ds_pal <- setNames(dataset_colors[seq_along(ds_levels) %% length(dataset_colors)],
                         ds_levels)


      # Only include colors for groups that exist in the current plot
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
      if ((input$transform_mode %||% "log2p1") == "log2p1") {
        col_fun <- circlize::colorRamp2(c(rng[1], rng[2]), c("white", "#030058"))
        legend_title <- "log2(TPM+1)"
      } else {
        lim <- max(3, ceiling(max(abs(rng), na.rm = TRUE)))
        col_fun <- circlize::colorRamp2(c(-lim, 0, lim), c("#2166AC", "white", "#B2182B"))
        legend_title <- "Z-score"
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
        # Height matches your per-row intention (points -> ComplexHeatmap expects grid units)
        #height = grid::unit(max(1, nrow(mat)) * dims$cell_h_pt, "pt"),
        use_raster = TRUE
      )

      # Explicitly move legends to the right and add padding at the bottom to avoid clipping
      ComplexHeatmap::draw(
        ht,
        heatmap_legend_side      = "right",
        annotation_legend_side   = "right",
        padding = grid::unit(c(6, 20, 20, 6), "mm")  # top, right, bottom, left
      )

      .last_ht(ht)
      gc()
    }, res = 120)

    output$heatmap_ui <- renderUI({
      dims <- plot_dims()
      h_px <- paste0(round(dims$dev_h_px), "px")

      # Allow the drawing canvas to use overflow beyond the immediate box if needed
      div(style = "overflow: visible;",
          plotOutput(
            outputId = ns("heatmap_plot"),
            height   = h_px,
            width    = "100%"
          )
      )
    })


    output$plot_notes <- renderUI({
      htmltools::HTML("<small>NA values are grey (not expressed / absent in a dataset).</small>")
    })

    # ---------------- downloads ----------------
    output$dl_svg <- downloadHandler(
      filename = function() paste0("heatmap_", Sys.Date(), ".svg"),
      content = function(file) {
        dims <- plot_dims()
        w_in <- (session$clientData[[paste0("output_", ns("heatmap_plot"), "_width")]] %||% 1000) / 96
        h_in <- (dims$dev_h_px %||% 800) / 96  # use dev_h_px, not dims$height
        svglite::svglite(file, width = w_in, height = h_in)
        on.exit(grDevices::dev.off(), add = TRUE)
        # Re-draw with the same legend placement & padding
        ComplexHeatmap::draw(.last_ht(),
                             heatmap_legend_side="right",
                             annotation_legend_side="right",
                             padding = grid::unit(c(6,10,16,6), "mm"))
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
        ComplexHeatmap::draw(.last_ht(),
                             heatmap_legend_side="right",
                             annotation_legend_side="right",
                             padding = grid::unit(c(6,10,16,6), "mm"))
      }
    )



  })
}

