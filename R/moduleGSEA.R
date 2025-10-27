#moduleGSEA.R
# R/modules/goGSEA_module.R
# -------------------------

# Sidebar controls
GSEASidebarUI <- function(id) {
  ns <- shiny::NS(id)
  tagList(
    selectInput(ns("dataset"), "Dataset", choices = NULL),
    radioButtons(ns("ont"), "Ontology",
                 choices = c("BP", "CC", "MF"), inline = TRUE),
    sliderInput(ns("ncat"), "Show top categories",
                min = 5, max = 50, value = 10, step = 1),
    tags$hr(),
    strong("Download"),
    fluidRow(
      column(6,
             downloadButton(ns("dl_plot_gr_png"), "GeneRatio PNG", class = "btn-sm"),
             br(), br(),
             downloadButton(ns("dl_plot_gr_svg"), "GeneRatio SVG", class = "btn-sm")
      ),
      column(6,
             downloadButton(ns("dl_plot_nes_png"), "NES PNG", class = "btn-sm"),
             br(), br(),
             downloadButton(ns("dl_plot_nes_svg"), "NES SVG", class = "btn-sm")
      )
    ),
    br(),
    downloadButton(ns("dl_table_csv"), "Results table (CSV)", class = "btn-sm")
  )
}

# Main content (plots + table)
GSEAMainUI <- function(id) {
  ns <- shiny::NS(id)
  tagList(
    tabsetPanel(
      id = ns("tabs"),
      tabPanel(
        "Dotplot (GeneRatio)",
        shinycssloaders::withSpinner(
          plotOutput(ns("plot_gr"), height = 460),
          type = 4
        )
      ),
      tabPanel(
        "Dotplot (NES)",
        shinycssloaders::withSpinner(
          plotOutput(ns("plot_nes"), height = 460),
          type = 4
        )
      ),
      tabPanel(
        "Results table",
        shinycssloaders::withSpinner(
          DT::DTOutput(ns("tbl")),
          type = 4
        )
      )
    )
  )
}

# R/modules/goGSEA_module.R  (replace just the server)
# R/modules/goGSEA_module.R
# -------------------------
# Server for GSEA (CSV/RDS-on-disk, on-demand)
GSEAServer <- function(id, base_dir = NULL, pkg = utils::packageName()) {
  moduleServer(id, function(input, output, session) {

    # ---- deps ----
    requireNamespace("DT", quietly = TRUE)
    requireNamespace("ggplot2", quietly = TRUE)
    requireNamespace("svglite", quietly = TRUE)
    requireNamespace("enrichplot", quietly = TRUE)

    `%||%` <- function(x,y) if (is.null(x) || (is.character(x) && !nzchar(x))) y else x

    # ---- locate extdata/GSEA_results no matter how the app is run ----
    if (is.null(base_dir) || !nzchar(base_dir)) {
      # installed package path
      base_dir <- system.file("extdata", "GSEA_results", package = pkg, mustWork = FALSE)
      # project (dev) fallback
      if (!nzchar(base_dir) || !dir.exists(base_dir)) {
        base_dir <- file.path("inst", "extdata", "GSEA_results")
      }
    }

    # ---- index: dataset folder -> {BP,CC,MF} -> file path ----
    build_index <- function(root) {
      if (!dir.exists(root)) return(list())
      ds_dirs <- list.dirs(root, recursive = FALSE, full.names = TRUE)
      idx <- lapply(ds_dirs, function(d) {
        # prefer .rds; fallback to .csv (if you kept CSVs)
        map <- setNames(vector("list", 3), c("BP","CC","MF"))
        for (ont in c("BP","CC","MF")) {
          p_rds <- file.path(d, sprintf("gsea_GO_%s.rds", ont))
          p_csv <- file.path(d, sprintf("gsea_GO_%s.csv", ont))
          if (file.exists(p_rds)) {
            map[[ont]] <- p_rds
          } else if (file.exists(p_csv)) {
            map[[ont]] <- p_csv
          } else {
            map[[ont]] <- ""
          }
        }
        map
      })
      names(idx) <- basename(ds_dirs)
      # keep datasets that have at least one ontology file
      idx[vapply(idx, function(x) any(nzchar(unlist(x))), logical(1))]
    }

    idx <- build_index(base_dir)
    datasets_vec <- names(idx)

    # ---- pretty labels for the dataset select ----
    # expects a global/internal object `pretty_map` (named character vector)
    base_key <- function(x) sub("_.*$", "", x)   # "Erwin_0.05_all_genes" -> "Erwin"
    labels <- vapply(datasets_vec, function(d) {
      key <- base_key(d)
      lab <- tryCatch(pretty_map[[key]], error = function(e) NULL)
      lab %||% key
    }, character(1))
    choices_named <- stats::setNames(datasets_vec, labels)

    # ---- reactive: are we safe to render? ----
    ready <- reactive({
      ds_ok  <- !is.null(input$dataset) && input$dataset %in% datasets_vec
      ont_ok <- !is.null(input$ont) && input$ont %in% c("BP","CC","MF")
      file_ok <- if (ds_ok && ont_ok) nzchar(idx[[input$dataset]][[input$ont]]) else FALSE
      ds_ok && ont_ok && file_ok
    })

    # ---- populate controls ----
    observe({
      validate(need(length(datasets_vec) > 0, "No GSEA files found under extdata/GSEA_results."))
      updateSelectInput(session, "dataset", choices = choices_named, selected = datasets_vec[[1]])
    })
    observeEvent(input$dataset, {
      updateRadioButtons(session, "ont", choices = c("BP","CC","MF"), selected = "BP")
    }, ignoreInit = TRUE)

    # ---- current file path & object ----
    r_path <- reactive({
      req(ready())
      idx[[input$dataset]][[input$ont]]
    })

    read_any <- function(p) {
      if (grepl("\\.rds$", p, ignore.case = TRUE)) {
        readRDS(p)
      } else {
        # CSV fallback support (returns a df, not a gseaResult). We convert.
        df <- utils::read.csv(p, check.names = FALSE)
        df
      }
    }

    r_obj <- reactive({
      req(ready())
      read_any(r_path())
    })

    # ---- title ----
    title_txt <- reactive({
      paste0(input$dataset %||% "", " - GO ", input$ont %||% "", " GSEA")
    })

    # ---- df accessor (works for gseaResult or data.frame fallback) ----
    r_df <- reactive({
      x <- r_obj()
      df <- if (inherits(x, "gseaResult")) {
        as.data.frame(x)
      } else {
        # CSV fallback already returns a data.frame
        x
      }
      # order nicely
      if ("p.adjust" %in% names(df)) {
        df <- df[order(df$p.adjust, -abs(df$NES %||% 0)), , drop = FALSE]
      }
      df
    })

    # ---- plot helpers ----
    make_dotplot <- function(metric = c("GeneRatio","NES")) {
      metric <- match.arg(metric)
      x <- r_obj()
      if (inherits(x, "gseaResult")) {
        enrichplot::dotplot(x, x = metric, showCategory = input$ncat %||% 10) +
          ggplot2::ggtitle(title_txt())
      } else {
        # data.frame fallback (if you kept only CSVs) — manual dotplot
        df <- r_df()
        n  <- max(1, as.integer(input$ncat %||% 10))
        d  <- head(df, n)
        if (!"GeneRatio" %in% names(d)) {
          # compute GeneRatio from core_enrichment/setSize if present
          if ("core_enrichment" %in% names(d) && "setSize" %in% names(d)) {
            core_n <- vapply(strsplit(as.character(d$core_enrichment), "/", fixed = TRUE),
                             function(v) sum(nzchar(v)), integer(1))
            d$GeneRatio <- core_n / as.numeric(d$setSize)
          }
        }
        aes_x <- if (metric == "NES") ggplot2::aes(x = NES) else ggplot2::aes(x = GeneRatio)
        ggplot2::ggplot(d, aes_x + ggplot2::aes(y = stats::reorder(Description, !!rlang::sym(names(d)[match(metric, c("GeneRatio","NES"))])))) +
          ggplot2::geom_point(ggplot2::aes(size = setSize, alpha = -log10(p.adjust))) +
          { if (metric == "NES") ggplot2::geom_vline(xintercept = 0, linetype = 2) else NULL } +
          ggplot2::labs(x = if (metric == "NES") "Normalized Enrichment Score (NES)" else "GeneRatio (core / set size)",
                        y = NULL, title = title_txt(), alpha = "-log10(adj.P)") +
          ggplot2::theme_minimal(base_size = 12)
      }
    }

    # ---- outputs ----
    output$plot_gr  <- renderPlot({ req(ready()); make_dotplot("GeneRatio") }, res = 96, height = 460)
    output$plot_nes <- renderPlot({ req(ready()); make_dotplot("NES")       }, res = 96, height = 460)

    output$tbl <- DT::renderDT({
      req(ready())
      DT::datatable(r_df(), rownames = FALSE, filter = "top",
                    options = list(pageLength = 10, scrollX = TRUE))
    })

    # ---- downloads ----
    output$dl_plot_gr_png <- downloadHandler(
      filename = function() sprintf("dotplot_GeneRatio_%s_%s.png", input$dataset, input$ont),
      content  = function(file) ggplot2::ggsave(file, make_dotplot("GeneRatio"), width = 7, height = 5, dpi = 300)
    )
    output$dl_plot_gr_svg <- downloadHandler(
      filename = function() sprintf("dotplot_GeneRatio_%s_%s.svg", input$dataset, input$ont),
      content  = function(file) { svglite::svglite(file, 7, 5); print(make_dotplot("GeneRatio")); dev.off() }
    )
    output$dl_plot_nes_png <- downloadHandler(
      filename = function() sprintf("dotplot_NES_%s_%s.png", input$dataset, input$ont),
      content  = function(file) ggplot2::ggsave(file, make_dotplot("NES"), width = 7, height = 5, dpi = 300)
    )
    output$dl_plot_nes_svg <- downloadHandler(
      filename = function() sprintf("dotplot_NES_%s_%s.svg", input$dataset, input$ont),
      content  = function(file) { svglite::svglite(file, 7, 5); print(make_dotplot("NES")); dev.off() }
    )
    output$dl_table_csv <- downloadHandler(
      filename = function() sprintf("gsea_table_%s_%s.csv", input$dataset, input$ont),
      content  = function(file) utils::write.csv(r_df(), file, row.names = FALSE)
    )
  })
}
