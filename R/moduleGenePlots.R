#' @importFrom shiny HTML
#' @noRd
genePlotsSidebarUI <- function(id) {
  ns <- NS(id)
  tagList(
    h4("Gene plots"),
    checkboxGroupInput(
      ns("gp_datasets"),
      label   = "Datasets",
      choices = character(0)   # filled by server
    ),
    helpText("You may select multiple datasets only within the same study as TPMs are not comparable across studies."),
    uiOutput(ns("gp_datasets_note")),
    tags$br(),
    textInput(
      ns("gp_gene"),
      "Gene",
      value = "FXN",   # pre-filled
      placeholder = "Type a gene symbol to search..."
    ),
    tags$br(),
    checkboxInput(ns("gp_logy"), "Log10 Y-axis", value = FALSE),
    tags$hr(),
    strong("Download"),
    br(),
    downloadButton(ns("dl_points"), "Replicates (CSV)"),
    downloadButton(ns("dl_summary"), "Summary (CSV)"),
    downloadButton(ns("dl_plot"), "Plot (SVG)")
  )
}

#' Gene Plots - main UI
#' @noRd
genePlotsMainUI <- function(id) {
  ns <- NS(id)
  tagList(
    # if you don't have shinycssloaders installed yet, use plotOutput directly
    if (requireNamespace("shinycssloaders", quietly = TRUE)) {
      shinycssloaders::withSpinner(plotOutput(ns("gp_plot"), height = 620), type = 4,  color = "#005249")
    } else {
      plotOutput(ns("gp_plot"), height = 620)
    },
    br(),
    DT::dataTableOutput(ns("gp_table"))
  )
}

#' Gene Plots server - multi-dataset (checkboxes)
#' @noRd
genePlotsServer <- function(id, pkg = utils::packageName()) {
  moduleServer(id, function(input, output, session) {

    observeEvent(TRUE, {
      ensure_atlas_data("tpm_gene", package = pkg)
    }, once = TRUE)

    tpm_dir <- file.path(
      tools::R_user_dir(pkg, which = "cache"),
      "tpm"
    )

    # --- make pkg safe for both project + installed package modes ---
    pkg <- tryCatch(pkg, error = function(e) "")
    if (!length(pkg) || !is.character(pkg) || !nzchar(pkg)) pkg <- "FRDATranscriptomicAtlas"
    pkg <- pkg[[1L]]

    # ----- Pretty labels matching file base names -----
    pretty_map  <- c(
      "Chutake"      = "Chutake (Lymphoblastoid Cells)",
      "Erwin"        = "Erwin (Lymphoblastoid Cells)",
      "Indelicato"   = "Indelicato (Skeletal Muscle)",
      "Lai_iPSC"     = "Lai (iPSCs)",
      "Lai_CNS"      = "Lai (CNS neurons)",
      "Lai_PNS"      = "Lai (PNS neurons)",
      "Lees"         = "Lees (Cardiomyocytes)",
      "Li"           = "Li (Cardiomyocytes)",
      "Maddock_LMN"  = "Maddock (Lower Motor Neurons)",
      "Maddock_SN"   = "Maddock (Sensory Neurons)",
      "Maddock_NCC"  = "Maddock (Neural Crest Cells)",
      "Mishra"       = "Mishra (Neurons)",
      "Napierala"    = "Napierala (Fibroblasts)",
      "Vilema"       = "Vilema-Enriquez (Fibroblasts)",
      "Wang"         = "Wang (Fibroblasts)"
    )
    `%||%` <- function(a, b) if (is.null(a)) b else a
    pretty_label <- function(id) pretty_map[[id]] %||% id

    # Define dataset families based on prefix
    dataset_family <- function(ds) sub("_.*$", "", ds)

    observeEvent(input$gp_datasets, {
      req(input$gp_datasets)

      fams <- dataset_family(input$gp_datasets)

      # More than one family selected?
      if (length(unique(fams)) > 1) {

        # Drop the newly added dataset and keep only the original family
        original_family <- dataset_family(input$gp_datasets[1])
        valid <- input$gp_datasets[fams == original_family]

        showNotification(
          "Please select datasets from only one dataset family (e.g., only Lai, only Maddock).",
          type = "error",
          duration = 5
        )

        updateCheckboxGroupInput(
          session, "gp_datasets",
          selected = valid
        )
      }
    })

    # Locate and source the theme file (works in both modes)
    if (!exists("theme_Marnie", inherits = TRUE)) {
      tp <- system.file("R", "utils_graphTheme.R", package = pkg, mustWork = FALSE)
      if (!nzchar(tp)) tp <- file.path("R", "utils_graphTheme.R")
      if (file.exists(tp)) source(tp)
    }

    # ----- Locate TPM RDS files (package OR project) -----
    # tpm_dir <- system.file("extdata/tpm", package = pkg, mustWork = FALSE)
    # if (!nzchar(tpm_dir)) tpm_dir <- file.path("inst", "extdata", "tpm")

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

    # ----- Populate dataset checkboxes -----
    observe({
      m  <- manifest()
      ds <- sort(intersect(unique(m$dataset), names(pretty_map)))
      validate(need(length(ds), "No TPM RDS files found in extdata/tpm/."))

      choices_named <- stats::setNames(ds, unname(pretty_map[ds]))
      pre_sel <- if ("Maddock_LMN" %in% ds) "Maddock_LMN" else ds[1]

      updateCheckboxGroupInput(session, "gp_datasets",
                               choices  = choices_named,
                               selected = pre_sel)
    })

    output$gp_datasets_note <- renderUI({
      req(input$gp_datasets)
      labs <- vapply(input$gp_datasets, pretty_label, "", USE.NAMES = FALSE)
      htmltools::HTML(sprintf("<small>Selected: <b>%s</b></small>", paste(labs, collapse = ", ")))
    })

    # ----- Files -> load -> tidy -> plot/table (unchanged below this line) -----
    file_paths <- reactive({
      req(input$gp_datasets)
      m <- manifest()
      hits <- dplyr::filter(m, dataset %in% input$gp_datasets)
      validate(need(nrow(hits) >= 1, "No TPM files found for the selected datasets."))
      stats::setNames(hits$path, hits$dataset)
    })

    tpms_wide_list <- reactive({
      fps <- file_paths()
      out <- lapply(names(fps), function(ds) {
        x <- readRDS(fps[[ds]])
        if (!is.data.frame(x)) x <- as.data.frame(x)
        validate(need(all(c("gene_id", "gene_name") %in% names(x)),
                      sprintf("TPM RDS for '%s' must contain 'gene_id' and 'gene_name'.", ds)))
        x
      })
      names(out) <- names(fps)
      out
    })

    observeEvent(tpms_wide_list(), {
      wl <- tpms_wide_list()
      genes <- sort(unique(unlist(lapply(wl, function(x) x$gene_name))))

      current <- isolate(input$gp_gene)

      # keep user choice if it still exists in the new choices
      if (!is.null(current) && nzchar(current) && current %in% genes) {
        selected <- current
      } else {
        selected <- if ("FXN" %in% genes) "FXN" else genes[1]
      }

      updateSelectizeInput(
        session, "gp_gene",
        choices = genes,
        selected = selected,
        server = TRUE
      )
    }, ignoreInit = FALSE)


    tpms_long_all <- reactive({
      wl <- tpms_wide_list()
      if (!length(wl)) return(dplyr::tibble())
      parts <- lapply(names(wl), function(ds) {
        x <- wl[[ds]]
        lng <- tidyr::pivot_longer(
          x, cols = -c(gene_id, gene_name),
          names_to = "sample", values_to = "TPM"
        )
        lng$dataset <- ds
        lng
      })
      dplyr::bind_rows(parts)
    })

    base_of <- function(x) sub("_.*$", "", x)
    strip_dataset_prefix <- function(x, ds_id) {
      x2 <- sub(paste0("^", ds_id, "_"), "", x)
      if (identical(x2, x)) sub(paste0("^", base_of(ds_id), "_"), "", x) else x2
    }

    parse_conditions <- function(df_ds) {
      ds <- unique(df_ds$dataset)
      repnum <- suppressWarnings(as.integer(sub("(?i).*_rep([0-9]+)$", "\\1", df_ds$sample, perl = TRUE)))
      cond   <- sub("(?i)_rep[0-9]+$", "", df_ds$sample, perl = TRUE)
      repnum[is.na(repnum)] <- NA_integer_
      cond <- strip_dataset_prefix(cond, ds)

      if (identical(ds, "Maddock")) {
        celltype_order <- c("SN","LMN","NCC")
        tmp <- tidyr::extract(data.frame(condition = cond), "condition",
                              into = c("fa_num","ic","ctype"),
                              regex = "^FA(\\d+)(ic)?(SN|LMN|NCC)$", remove = FALSE)
        tmp$fa_num <- as.integer(tmp$fa_num)
        tmp$ic_first <- ifelse(is.na(tmp$ic), 1L, 0L)
        tmp$ctype_ord <- match(tmp$ctype, celltype_order)
        levs <- tmp |>
          dplyr::arrange(fa_num, ic_first, ctype_ord) |>
          dplyr::pull(condition) |>
          unique()
        cond <- factor(cond, levels = levs)
      } else if (identical(ds, "Lai")) {
        tissue_order <- c("iPSC","CNS","PNS")
        tmp <- tidyr::extract(data.frame(condition = cond), "condition",
                              into = c("status","tissue"),
                              regex = "^(FRDA2|FRDA|IC)_(CNS|iPSC|PNS)$", remove = FALSE)
        tmp$tissue_ord <- match(tmp$tissue, tissue_order)
        tmp$status_rank <- dplyr::case_when(
          tmp$status == "IC" ~ 1L, tmp$status == "FRDA" ~ 2L, tmp$status == "FRDA2" ~ 3L, TRUE ~ 9L)
        levs <- tmp |>
          dplyr::arrange(tissue_ord, status_rank) |>
          dplyr::pull(condition) |>
          unique()
        cond <- factor(cond, levels = levs)
      } else if (identical(ds, "Mishra")) {
        tmp <- tidyr::extract(data.frame(condition = cond), "condition",
                              into = c("status","code"),
                              regex = "^(FRDA|IC)_([A-Za-z0-9]+)$", remove = FALSE)
        unify_group <- function(code) if (code == "E35") "FF1" else code
        tmp$group <- vapply(tmp$code, unify_group, "", USE.NAMES = FALSE)
        groups <- unique(tmp$group)
        num_g <- groups[grepl("^\\d+$", groups)]
        non_g <- setdiff(groups, num_g)
        num_g <- num_g[order(as.integer(num_g))]
        preferred_non <- c("FF1","FF2")
        non_pref <- setdiff(non_g, preferred_non)
        non_g <- c(intersect(preferred_non, non_g), sort(non_pref))
        ordered_groups <- c(num_g, non_g)
        ic_code_for <- function(g) if (g == "FF1") "E35" else g
        frda_code_for <- function(g) g
        levs <- unlist(lapply(ordered_groups, function(g) {
          c(paste0("IC_", ic_code_for(g)), paste0("FRDA_", frda_code_for(g)))
        }))
        levs <- levs[levs %in% cond]
        cond <- factor(cond, levels = levs)
      } else {
        cond <- factor(cond, levels = sort(unique(cond)))
      }

      dplyr::mutate(df_ds, condition = cond, rep = repnum)
    }

    tpms_long_parsed <- reactive({
      all <- tpms_long_all()
      if (!nrow(all)) return(all)
      parsed <- lapply(split(all, all$dataset), parse_conditions)
      dplyr::bind_rows(parsed)
    })

    dat_points <- reactive({
      req(input$gp_datasets, input$gp_gene)
      gene_input <- trimws(input$gp_gene)
      df <- tpms_long_parsed() |>
        dplyr::filter(dataset %in% input$gp_datasets, gene_name == gene_input)
      if (nrow(df) == 0) {
        shiny::showNotification(sprintf("Warning - Gene '%s' was not found in the selected datasets.", gene_input),
                                type = "error", duration = 5)
      }
      df
    })

    dat_summary <- reactive({
      dat_points() |>
        dplyr::group_by(dataset, condition) |>
        dplyr::summarise(n = dplyr::n(), mean = mean(TPM, na.rm = TRUE),
                         sd = stats::sd(TPM, na.rm = TRUE), .groups = "drop")
    })

    build_plot <- function(points, summary, logy = FALSE, title_txt = "") {
      points$dataset  <- factor(points$dataset,  levels = input$gp_datasets)
      summary$dataset <- factor(summary$dataset, levels = input$gp_datasets)
      lab_ds <- ggplot2::labeller(dataset = ggplot2::as_labeller(
        function(v) vapply(v, pretty_label, "", USE.NAMES = FALSE)))
      p <- ggplot2::ggplot() +
        ggplot2::geom_crossbar(data = summary,
                               ggplot2::aes(x = condition, y = mean, ymin = mean, ymax = mean),
                               width = 0.6, linewidth = 1) +
        ggplot2::geom_errorbar(data = summary,
                               ggplot2::aes(x = condition, ymin = mean - sd, ymax = mean + sd),
                               width = 0.3, linewidth = 0.8) +
        ggplot2::geom_point(data = points,
                            ggplot2::aes(x = condition, y = TPM),
                            position = ggplot2::position_jitter(width = 0.15, height = 0, seed = 1),
                            size = 3.5, color = "#005249", na.rm = TRUE, show.legend = FALSE) +
        ggplot2::labs(title = title_txt,  x = NULL, y = "Transcripts Per Million") +
        theme_Marnie() +
        ggplot2::facet_wrap(dplyr::vars(dataset), scales = "free_x", nrow = 1, labeller = lab_ds)
      if (isTRUE(logy)) p <- p + ggplot2::scale_y_continuous(trans = "log10")
      p
    }

    output$gp_plot <- renderPlot({
      pts <- dat_points(); sms <- dat_summary()
      validate(need(nrow(pts) > 0, "No TPM values for this selection."))
      ds_lab <- if (length(input$gp_datasets) == 1) pretty_label(input$gp_datasets)
      else paste0(length(input$gp_datasets), " datasets")
      build_plot(pts, sms, input$gp_logy, paste(input$gp_gene, "-", ds_lab))
    })

    output$gp_table <- DT::renderDataTable({
      DT::datatable(
        dat_points() |>
          dplyr::arrange(dataset, condition, rep) |>
          dplyr::select(dataset, condition, rep, gene_id, gene_name, TPM),
        rownames = FALSE, options = list(pageLength = 5, scrollX = TRUE)
      )
    }, server = TRUE)

    safe_stem <- function(ds_vec) {
      if (!length(ds_vec)) return("none")
      if (length(ds_vec) <= 3) paste(ds_vec, collapse = "+") else
        paste0(ds_vec[1], "+", length(ds_vec) - 1, "more")
    }

    output$dl_points <- downloadHandler(
      filename = function() sprintf("TPM_%s_%s_replicates.csv",
                                    safe_stem(input$gp_datasets), input$gp_gene),
      content  = function(file) readr::write_csv(dat_points(), file)
    )
    output$dl_summary <- downloadHandler(
      filename = function() sprintf("TPM_%s_%s_summary.csv",
                                    safe_stem(input$gp_datasets), input$gp_gene),
      content  = function(file) readr::write_csv(dat_summary(), file)
    )
    output$dl_plot <- downloadHandler(
      filename = function() sprintf("TPM_%s_%s.svg",
                                    safe_stem(input$gp_datasets), input$gp_gene),
      content = function(file) {
        pts <- isolate(dat_points()); sms <- isolate(dat_summary())
        validate(need(nrow(pts) > 0, "No data available for this gene."))
        ds_lab <- isolate(if (length(input$gp_datasets) == 1)
          pretty_label(input$gp_datasets)
          else paste0(length(input$gp_datasets), " datasets"))
        grDevices::svg(file, width = 11, height = 7, onefile = TRUE)
        print(build_plot(pts, sms, isolate(input$gp_logy),
                         paste(isolate(input$gp_gene), "-", ds_lab)))
        grDevices::dev.off()
      }
    )
  })
}
