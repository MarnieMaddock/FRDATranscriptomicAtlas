# R/moduleGenePlots.R

#' Gene Plots – sidebar
#' @noRd
genePlotsSidebarUI <- function(id) {
  ns <- NS(id)
  tagList(
    h4("Gene plots"),
    selectizeInput(ns("gp_dataset"), "Dataset", choices = NULL, options = list(placeholder = "Select a dataset…")),
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

#' Gene Plots – main UI
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

#' Gene Plots server – one TPM RDS per dataset
#' @noRd
genePlotsServer <- function(id, pkg = utils::packageName()) {
  moduleServer(id, function(input, output, session) {

    # ----- Pretty labels matching file base names -----
    pretty_map <- c(
      "Erwin"      = "Erwin (Lymphoblastoid Cells)",
      "Indelicato" = "Indelicato (Skeletal Muscle)",
      "Lai"        = "Lai (iPSCs, CNS neurons, PNS neurons)",
      "Lees"       = "Lees (Cardiomyocytes)",
      "Maddock"    = "Maddock (Sensory Neurons, Lower Motor Neurons, Neural Crest Cells)",
      "Mishra"     = "Mishra (Neurons)",
      "Napierala"  = "Napierala (Fibroblasts)",
      "Vilema"     = "Vilema-Enriquez (Fibroblasts)"
    )
    `%||%` <- function(a, b) if (is.null(a)) b else a
    pretty_label <- function(id) pretty_map[[id]] %||% id

    # Locate and source the theme file
    theme_path <- system.file("R", "utils_graphTheme.R", package = pkg, mustWork = FALSE)
    if (file.exists(theme_path)) source(theme_path)

    # ----- Locate TPM RDS files -----
    tpm_dir <- system.file("extdata/tpm", package = pkg, mustWork = FALSE)
    manifest <- reactive({
      files <- if (nzchar(tpm_dir)) list.files(tpm_dir, pattern = "_gene_tpm\\.rds$", full.names = TRUE) else character(0)
      if (!length(files)) return(tibble::tibble(path = character(), dataset = character()))
      # Extract base name before "_gene_tpm.rds"
      tibble::tibble(
        path = files,
        dataset = sub("_gene_tpm\\.rds$", "", basename(files))
      )
    })

    # ----- Populate dataset dropdown -----
    observe({
      m <- manifest()
      ds <- sort(unique(m$dataset))
      validate(need(length(ds), "No TPM RDS files found in extdata/tpm/."))

      # fallback labeler: use pretty_map if present, else the id itself
      label_for <- function(ids) {
        labs <- unname(pretty_map[ids])
        labs[is.na(labs)] <- ids[is.na(labs)]
        labs
      }

      choices_named <- stats::setNames(ds, label_for(ds))

      sel <- if ("Maddock" %in% ds) "Maddock" else ds[1]
      updateSelectizeInput(session, "gp_dataset",
                           choices  = choices_named,
                           selected = sel,
                           server   = TRUE)
    })

    # ----- Select file & load it -----
    file_sel <- reactive({
      req(input$gp_dataset)
      hit <- dplyr::filter(manifest(), dataset == input$gp_dataset)
      validate(need(nrow(hit) >= 1, sprintf("No TPM file for dataset '%s'.", input$gp_dataset)))
      hit$path[1]
    })

    tpms_wide <- reactive({
      x <- readRDS(file_sel())
      if (!is.data.frame(x)) x <- as.data.frame(x)
      validate(need(all(c("gene_id", "gene_name") %in% names(x)),
                    "TPM RDS must contain 'gene_id' and 'gene_name' columns."))
      x
    })

    # ----- Gene dropdown -----
    observeEvent(tpms_wide(), {
      genes <- sort(unique(tpms_wide()$gene_name))
      updateSelectizeInput(
        session, "gp_gene",
        choices  = genes,        # used for autocomplete, but user can type anything
        selected = if ("FXN" %in% genes) "FXN" else genes[1],
        server   = TRUE
      )
    })


    # ----- Long format -----
    tpms_long <- reactive({
      tidyr::pivot_longer(
        tpms_wide(),
        cols = -c(gene_id, gene_name),
        names_to = "sample",
        values_to = "TPM"
      )
    })

    # ----- Parse sample name into condition & replicate -----
    tpms_long_parsed <- reactive({
      tl <- tpms_long()
      ds <- input$gp_dataset
      req(ds)

      # Extract replicate number (if present)
      repnum <- suppressWarnings(as.integer(sub(".*_REP([0-9]+)$", "\\1", tl$sample)))
      repnum[is.na(repnum)] <- NA_integer_

      # 1) remove trailing _REPX
      cond <- sub("_REP[0-9]+$", "", tl$sample)

      # 2) if the sample starts with "<DATASET>_", strip that prefix
      #    e.g., "Mishra_FRDA_223" -> "FRDA_223"
      cond <- sub(paste0("^", ds, "_"), "", cond)

      # If you ever get dots from make.names, normalize dots->underscores (optional):
      # cond <- gsub("\\.", "_", cond, fixed = FALSE)

      dplyr::mutate(
        tl,
        dataset   = ds,
        condition = cond,
        rep       = repnum
      )
    })


    # ----- Filter for current selection -----
    dat_points <- reactive({
      req(input$gp_dataset, input$gp_gene)
      gene_input <- trimws(input$gp_gene)

      df <- tpms_long_parsed() |>
        dplyr::filter(gene_name == gene_input)

      if (nrow(df) == 0) {
        shiny::showNotification(
          sprintf("⚠️ Gene '%s' does not exist in this dataset.", gene_input),
          type = "error",
          duration = 5
        )
      }

      if (identical(input$gp_dataset, "Maddock")) {
        # order of cell types:
        celltype_order <- c("SN", "LMN", "NCC")

        # Parse: FA number, ic flag, and cell type suffix
        tmp <- df |>
          tidyr::extract(
            condition,
            into = c("fa_num", "ic", "ctype"),
            regex = "^FA(\\d+)(ic)?(SN|LMN|NCC)$",
            remove = FALSE
          ) |>
          dplyr::mutate(
            fa_num    = as.integer(fa_num),
            ic_first  = ifelse(is.na(ic), 1L, 0L),                  # ic first (0 before 1)
            ctype_ord = match(ctype, celltype_order)
          )

        levs <- tmp |>
          dplyr::arrange(fa_num, ic_first, ctype_ord) |>
          dplyr::pull(condition) |>
          unique()

        df$condition <- factor(df$condition, levels = levs)
      } else {
        # Default: alphabetical (or add other dataset-specific rules as needed)
        df$condition <- factor(df$condition, levels = sort(unique(df$condition)))
      }

      if (identical(input$gp_dataset, "Lai")) {
        # Tissue order (tweak if you prefer)
        tissue_order <- c("iPSC", "CNS", "PNS")

        tmp <- df |>
          tidyr::extract(
            condition,
            into = c("status","tissue"),
            regex = "^(FRDA2|FRDA|IC)_(CNS|iPSC|PNS)$",
            remove = FALSE
          ) |>
          dplyr::mutate(
            tissue_ord  = match(tissue, tissue_order),
            status_rank = dplyr::case_when(
              status == "IC"    ~ 1L,   # IC first
              status == "FRDA"  ~ 2L,
              status == "FRDA2" ~ 3L,
              TRUE              ~ 9L
            )
          )

        levs <- tmp |>
          dplyr::arrange(tissue_ord, status_rank) |>
          dplyr::pull(condition) |>
          unique()

        df$condition <- factor(df$condition, levels = levs)
      }

      if (identical(input$gp_dataset, "Mishra")) {
        tmp <- df |>
          tidyr::extract(
            condition,
            into = c("status","code"),
            regex = "^(FRDA|IC)_([A-Za-z0-9]+)$",
            remove = FALSE
          )

        # Unify FRDA/IC into comparison groups
        # E35 is the isogenic control for FRDA_FF1
        unify_group <- function(code) if (code == "E35") "FF1" else code
        tmp$group <- vapply(tmp$code, unify_group, "", USE.NAMES = FALSE)

        # Determine group order: numeric ascending, then preferred named cohorts
        groups <- unique(tmp$group)
        num_g  <- groups[grepl("^\\d+$", groups)]
        non_g  <- setdiff(groups, num_g)
        num_g  <- num_g[order(as.integer(num_g))]

        preferred_non <- c("FF1", "FF2")          # edit if you want a different order
        non_pref      <- setdiff(non_g, preferred_non)
        non_g         <- c(intersect(preferred_non, non_g), sort(non_pref))

        ordered_groups <- c(num_g, non_g)

        # Map group -> actual IC/FRDA codes present in data
        ic_code_for <- function(g) if (g == "FF1") "E35" else g
        frda_code_for <- function(g) g

        # Build final factor levels with IC before FRDA per group
        levs <- unlist(lapply(ordered_groups, function(g) {
          cand <- c(paste0("IC_",   ic_code_for(g)),
                    paste0("FRDA_", frda_code_for(g)))
          cand[cand %in% df$condition]
        }))

        df$condition <- factor(df$condition, levels = levs)
      }

      df
    })


    dat_summary <- reactive({
      dat_points() |>
        dplyr::group_by(condition) |>
        dplyr::summarise(
          n    = dplyr::n(),
          mean = mean(TPM, na.rm = TRUE),
          sd   = sd(TPM,   na.rm = TRUE),
          .groups = "drop"
        )
    })

    # ----- Plot -----
    build_plot <- function(points, summary, logy = FALSE, title_txt = "") {
      p <- ggplot2::ggplot() +
        ggplot2::geom_crossbar(
          data = summary,
          ggplot2::aes(x = condition, y = mean, ymin = mean, ymax = mean),
          width = 0.6, linewidth = 1
        ) +
        ggplot2::geom_errorbar(
          data = summary,
          ggplot2::aes(x = condition, ymin = mean - sd, ymax = mean + sd),
          width = 0.3, linewidth = 0.8
        ) +
        ggplot2::geom_point(
          data = points,
          ggplot2::aes(x = condition, y = TPM),
          position = ggplot2::position_jitter(width = 0.15, height = 0, seed = 1),
          size = 4, color = "#005249", na.rm = TRUE, show.legend = FALSE
        ) +
        ggplot2::labs(title = title_txt, x = " ", y = "Transcripts Per Million") +
        theme_Marnie
      if (isTRUE(logy)) p <- p + ggplot2::scale_y_continuous(trans = "log10")
      p
    }

    output$gp_plot <- renderPlot({
      pts <- dat_points(); sms <- dat_summary()
      validate(need(nrow(pts) > 0, "No TPM values for this selection."))
      build_plot(pts, sms, input$gp_logy,
                 paste(input$gp_gene, "—", pretty_label(input$gp_dataset)))
    })


    # ----- Table -----
    output$gp_table <- DT::renderDataTable({
      DT::datatable(
        dat_points() |>
          dplyr::arrange(condition, rep) |>
          dplyr::select(condition, rep, gene_id, gene_name, TPM),
        rownames = FALSE,
        options = list(pageLength = 5, scrollX = TRUE)
      )
    }, server = TRUE)

    # ----- Downloads -----
    output$dl_points <- downloadHandler(
      filename = function() sprintf("TPM_%s_%s_replicates.csv", input$gp_dataset, input$gp_gene),
      content  = function(file) readr::write_csv(dat_points(), file)
    )
    output$dl_summary <- downloadHandler(
      filename = function() sprintf("TPM_%s_%s_summary.csv", input$gp_dataset, input$gp_gene),
      content  = function(file) readr::write_csv(dat_summary(), file)
    )
    # ---- downloads ----
    output$dl_plot <- downloadHandler(
      filename = function() sprintf("TPM_%s_%s.svg", input$gp_dataset, input$gp_gene),
      content = function(file) {
        pts <- isolate(dat_points())
        sms <- isolate(dat_summary())
        validate(need(nrow(pts) > 0, "No data available for this gene."))

        # Open SVG device instead of PNG
        grDevices::svg(file, width = 9, height = 6, onefile = TRUE)
        print(build_plot(pts, sms, isolate(input$gp_logy),
                         paste(isolate(input$gp_gene), "—", pretty_label(isolate(input$gp_dataset)))))
        grDevices::dev.off()
      }
    )

  })
}
