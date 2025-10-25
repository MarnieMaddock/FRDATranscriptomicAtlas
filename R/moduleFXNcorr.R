# FXN correlations — Sidebar UI
FXNcorrSidebarUI <- function(id) {
  ns <- NS(id)
  tagList(
    h4("FXN correlations"),
    textInput(
      ns("gene"),
      label = NULL,
      placeholder = "Type a gene symbol (e.g., PABPN1)",
      value = "PABPN1"
    ),
    checkboxInput(ns("color_by_dataset"), "Color points by dataset", value = FALSE)
  )
}


# FXN correlations — Main UI
FXNcorrMainUI <- function(id) {
  ns <- NS(id)
  tagList(
    # Plot
    shinycssloaders::withSpinner(plotOutput(ns("scatter"), height = "420px"),type = 4,  color = "#005249"),
    tags$div(
      style = "margin-top: .5rem;",
      downloadButton(ns("dl_plot_svg"), "Download plot (SVG)")
    ),
    tags$hr(),
    # Selected gene stats
    h5("Selected gene statistics"),
    tableOutput(ns("sel_stats")),
    tags$hr(),
    # Full table + download
    h5("All correlations (FXN vs gene)"),
    DT::dataTableOutput(ns("corr_table")),
    tags$div(
      style = "margin-top: .5rem;",
      downloadButton(ns("dl_corr_csv"), "Download table (CSV)")
    )
  )
}

.pkg_or_proj_path <- function(..., package, must_exist = TRUE) {
  local <- file.path("inst", ...)
  if (file.exists(local)) return(normalizePath(local, winslash = "/"))
  pkg <- system.file(..., package = package, mustWork = FALSE)
  if (nzchar(pkg)) return(normalizePath(pkg, winslash = "/"))
  if (isTRUE(must_exist)) stop("Could not locate file: ", file.path("inst", ...))
  ""
}

FXNcorrServer <- function(
    id,
    pkg = tryCatch(utils::packageName(), error = function(e) "") %||% ""
) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    `%||%` <- function(a, b) if (is.null(a)) b else a

    # ---- make pkg safe ----
    pkg <- tryCatch(pkg, error = function(e) "")
    if (!length(pkg) || !is.character(pkg) || !nzchar(pkg)) pkg <- "FRDATranscriptomicAtlas"
    pkg <- pkg[[1L]]

    # ---- pretty labels ----
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
    pretty_label <- function(id) {
      # return id if NA/empty or not in map
      if (is.na(id) || !nzchar(id)) return(id)
      lab <- unname(pretty_map[id])  # "[" returns NA instead of error for unknown/NA keys
      if (is.null(lab) || is.na(lab) || !length(lab)) id else lab
    }

    # Locate and source the theme file (works in both modes)
    if (!exists("theme_Marnie", inherits = TRUE)) {
      tp <- system.file("R", "utils_graphTheme.R", package = pkg, mustWork = FALSE)
      if (!nzchar(tp)) tp <- file.path("R", "utils_graphTheme.R")
      if (file.exists(tp)) source(tp)
    }

    requireNamespace("readr", quietly = TRUE)
    requireNamespace("dplyr", quietly = TRUE)
    requireNamespace("ggplot2", quietly = TRUE)
    requireNamespace("DT", quietly = TRUE)
    requireNamespace("svglite", quietly = TRUE)
    requireNamespace("tibble", quietly = TRUE)
    requireNamespace("purrr", quietly = TRUE)
    requireNamespace("stringr", quietly = TRUE)

    # ── load correlation CSV ────────────────────────────────────────────────────
    corr_file <- .pkg_or_proj_path("extdata", "corr", "correlations_data_FXN_v_X.csv",
                                   package = pkg)

    corr_tbl <- suppressMessages(
      readr::read_csv(corr_file, show_col_types = FALSE)
    ) |>
      dplyr::select(-dplyr::matches("^\\.\\.\\.")) |>
      dplyr::mutate(
        gene_name   = as.character(.data$gene_name),
        cor         = as.numeric(.data$cor),
        p_value     = as.numeric(.data$p_value),
        adj_p_value = as.numeric(.data$adj_p_value)
      ) |>
      dplyr::arrange(dplyr::desc(.data$cor))



    # populate gene choices (default = top row, e.g., PABPN1)
    observe({
      updateSelectizeInput(session, "gene",
                           choices = corr_tbl$gene_name,
                           selected = corr_tbl$gene_name[[1]], server = TRUE
      )
    })

    sel_row <- reactive({
      g <- input$gene
      if (is.null(g) || !nzchar(g)) return(NULL)
      dplyr::filter(corr_tbl, .data$gene_name == g) |> dplyr::slice(1)
    })

    output$sel_stats <- renderTable({
      x <- sel_row()
      if (is.null(x) || nrow(x) == 0) return(NULL)
      x |>
        dplyr::mutate(
          cor = round(cor, 3),
          p_value = format(p_value, digits = 3, scientific = TRUE),
          adj_p_value = format(adj_p_value, digits = 3, scientific = TRUE)
        )
    })

    # ── load and merge TPM matrices ────────────────────────────────────────────
    tpm_dir <- system.file("extdata/tpm", package = pkg, mustWork = FALSE)
    if (!nzchar(tpm_dir)) tpm_dir <- file.path("inst", "extdata", "tpm")

    manifest <- reactive({
      files <- if (nzchar(tpm_dir) && dir.exists(tpm_dir)) {
        list.files(tpm_dir, pattern = "_gene_tpm\\.rds$", full.names = TRUE)
      } else character(0)
      if (!length(files)) {
        return(tibble::tibble(path = character(), dataset = character()))
      }
      tibble::tibble(
        path = files,
        dataset = sub("_gene_tpm\\.rds$", "", basename(files))
      )
    })

    merged <- reactive({
      man <- manifest()
      if (!nrow(man)) return(NULL)

      # helper: turn one RDS df -> matrix with rownames = gene_name, numeric sample cols
      build_tpm_mat <- function(path) {
        df <- readRDS(path)
        # Must have gene_name + at least one sample col
        stopifnot("gene_name" %in% colnames(df))
        samp_cols <- setdiff(colnames(df), c("gene_id", "gene_name"))
        m <- as.matrix(df[, samp_cols, drop = FALSE])
        storage.mode(m) <- "numeric"
        rownames(m) <- df$gene_name
        m
      }

      mats <- lapply(seq_len(nrow(man)), function(i) {
        ds <- man$dataset[[i]]                # e.g., "Erwin"
        m  <- build_tpm_mat(man$path[[i]])    # rows = gene_name, cols = samples
        # If columns already start with "Erwin_", don't add any prefix
        if (!all(startsWith(colnames(m), paste0(ds, "_")))) {
          colnames(m) <- paste0(ds, ":", colnames(m))
        }
        m
      })

      gene_expr <- do.call(cbind, mats)

      # FXN vector by gene_name
      if ("FXN" %in% rownames(gene_expr)) {
        fxn_expr <- setNames(as.numeric(gene_expr["FXN", , drop = TRUE]), colnames(gene_expr))
      } else {
        fxn_expr <- setNames(rep(NA_real_, ncol(gene_expr)), colnames(gene_expr))
      }

      # sample -> dataset mapping
      sample_ids <- colnames(gene_expr)

      # robust: if no “_” or “:”, keep full sample name; never return NA
      dataset_id <- ifelse(grepl("[_:]", sample_ids),
                           sub("^([^:_]+)[_:].*$", "\\1", sample_ids),
                           sample_ids)

      sample_info <- tibble::tibble(
        sample_id   = sample_ids,
        dataset_id  = dataset_id,
        dataset_lab = vapply(dataset_id, pretty_label, FUN.VALUE = character(1))
      )

      list(gene_expr = gene_expr, fxn_expr = fxn_expr, sample_info = sample_info)
    })


    # ── plotting helper so we can reuse for download ───────────────────────────
    # Build the scatterplot for a selected gene (per-sample FXN vs gene TPM)
    # - Depends on `merged()` (reactive: list(gene_expr, fxn_expr, sample_info))
    # - Depends on `sel_row()` (reactive: one-row tibble with cor/p/FDR for the gene)
    # - Set color_by_ds = TRUE to color points by dataset (legend uses pretty labels)
    make_plot <- function(gene, color_by_ds = FALSE) {
      m <- merged()
      if (is.null(m)) {
        return(
          ggplot2::ggplot(data.frame(x = 0, y = 0), ggplot2::aes(x, y)) +
            ggplot2::annotate("text", x = 0, y = 0,
                              label = "No TPM matrices found in inst/extdata/tpm/",
                              size = 4.2) +
            ggplot2::theme_void()
        )
      }

      ge <- m$gene_expr     # matrix: rows = gene_name, cols = samples
      fx <- m$fxn_expr      # named numeric vector: names = samples
      si <- m$sample_info   # tibble: sample_id, dataset_id, dataset_lab

      # Check gene presence
      if (!gene %in% rownames(ge)) {
        return(
          ggplot2::ggplot(data.frame(x = 0, y = 0), ggplot2::aes(x, y)) +
            ggplot2::annotate("text", x = 0, y = 0,
                              label = sprintf("Gene '%s' not found in TPMs (row names are gene_name).", gene),
                              size = 4.2) +
            ggplot2::theme_void()
        )
      }

      # Check FXN availability
      if (all(!is.finite(fx))) {
        return(
          ggplot2::ggplot(data.frame(x = 0, y = 0), ggplot2::aes(x, y)) +
            ggplot2::annotate("text", x = 0, y = 0,
                              label = "FXN not found in TPMs (gene_name == 'FXN').",
                              size = 4.2) +
            ggplot2::theme_void()
        )
      }

      # Assemble per-sample data frame
      gx  <- as.numeric(ge[gene, , drop = TRUE])
      sam <- colnames(ge)
      ok  <- is.finite(gx) & is.finite(fx[sam])

      if (!any(ok)) {
        return(
          ggplot2::ggplot(data.frame(x = 0, y = 0), ggplot2::aes(x, y)) +
            ggplot2::annotate("text", x = 0, y = 0,
                              label = "No overlapping finite TPM values for FXN and selected gene.",
                              size = 4.2) +
            ggplot2::theme_void()
        )
      }

      df <- dplyr::tibble(
        FXN       = fx[sam][ok],
        GENE      = gx[ok],
        sample_id = sam[ok]
      ) |>
        dplyr::left_join(si, by = "sample_id")

      # Stats annotation from correlation table (keeps values consistent with table)
      sr <- sel_row()
      ann <- if (!is.null(sr) && nrow(sr)) {
        sprintf("r = %.3f\np = %s\nFDR = %s",
                sr$cor,
                format(sr$p_value, digits = 2, scientific = TRUE),
                format(sr$adj_p_value, digits = 2, scientific = TRUE))
      } else ""

      # --- build base plot ---
      p <- ggplot2::ggplot(df, ggplot2::aes(x = FXN, y = GENE)) +
        { if (isTRUE(color_by_ds))
          ggplot2::geom_point(ggplot2::aes(color = dataset_lab), alpha = 0.85, size = 2.8)
          else
            ggplot2::geom_point(alpha = 0.85, size = 2.8)
        } +
        ggplot2::geom_smooth(method = "lm", se = TRUE, color = "black") +
        ggplot2::labs(
          x = bquote(italic("FXN") ~ "expression (TPM)"),
          y = bquote(italic(.(gene)) ~ "expression (TPM)"),
          title = bquote(italic("FXN") ~ "vs" ~ italic(.(gene)) ~ "")
        ) +
        ggplot2::annotate("text", x = Inf, y = -Inf, label = ann,
                          hjust = 1.05, vjust = -0.6, size = 3.8) +
        theme_Marnie

      # ---- add your custom palette when coloring by dataset ----
      if (isTRUE(color_by_ds)) {
        colours <- c("#00359c", "#00b7c7", "#648fff", "#785ef0",
                     "#dc267f", "#fe6100", "#ffb000", "#FFEF00", "#72FF5F")

        p <- p +
          ggplot2::scale_color_manual(values = colours) +
          ggplot2::guides(color = ggplot2::guide_legend(title = "Dataset"))
      }

      p
    }


    # ── render plot ────────────────────────────────────────────────────────────
    output$scatter <- renderPlot({
      g <- input$gene
      if (is.null(g) || !nzchar(g)) return(NULL)
      suppressMessages(
        make_plot(gene = g, color_by_ds = isTRUE(input$color_by_dataset))
      )
    }, height = 420)

    # ── full table ────────────────────────────────────────────────────────────
    output$corr_table <- DT::renderDataTable({
      DT::datatable(
        corr_tbl,
        rownames = FALSE,
        options = list(pageLength = 10, scrollX = TRUE),
        filter = "top"
      )
    })

    # ── downloads ─────────────────────────────────────────────────────────────
    output$dl_plot_svg <- downloadHandler(
      filename = function() paste0("FXN_correlation_", input$gene %||% "plot", ".svg"),
      content  = function(file) {
        svglite::svglite(file, width = 10, height = 4.8, bg = "white")
        on.exit(grDevices::dev.off(), add = TRUE)
        print(make_plot(gene = input$gene %||% corr_tbl$gene_name[[1]],
                        color_by_ds = isTRUE(input$color_by_dataset)))
      }
    )

    output$dl_corr_csv <- downloadHandler(
      filename = function() "correlations_data_FXN_v_X.csv",
      content  = function(file) readr::write_csv(corr_tbl, file)
    )

    invisible(list(
      correlations = reactive(corr_tbl),
      merged_tpm   = merged,
      selected     = sel_row
    ))
  })
}
