#' @importFrom rlang .data
# a %||% b  ->  if a is NULL use b
`%||%` <- function(a, b) if (is.null(a)) b else a

# Discover labels from flat files in inst/extdata/DTU
# e.g. Maddock_SN_FA1_switch_summary.rds -> "Maddock_SN_FA1"
discover_dtu_labels_flat <- function(data_dir) {
  fs <- list.files(data_dir, pattern = "_switch_summary\\.rds$", full.names = FALSE)
  unique(sub("_switch_summary\\.rds$", "", fs))
}

# Build choices for selectize: names(shown) = pretty; values = machine labels
dtu_choice_vector <- function(labels, pretty_map = NULL) {
  if (is.null(pretty_map)) return(stats::setNames(labels, labels))
  disp <- unname(pretty_map[labels])
  # fallback to raw label if not in pretty_map
  disp[is.na(disp)] <- labels[is.na(disp)]
  stats::setNames(labels, disp)
}

# Paths for a given label (flat layout)
dtu_paths_for_label_flat <- function(data_dir, label) {
  list(
    summary = file.path(data_dir, paste0(label, "_switch_summary.rds")),
    genes   = file.path(data_dir, paste0(label, "_gene_switches.rds")),
    iso     = file.path(data_dir, paste0(label, "_significant_isoforms.rds"))
  )
}


# UI
dtuResultsSidebarUI <- function(id) {
  ns <- NS(id)
  tagList(
    h4("Differential Transcript Usage (DTU) (Using DEXSeq)"),
    tags$hr(style="margin:.5rem 0 1rem;"),
    selectizeInput(
      ns("label"), "Dataset / Contrast",
      choices = NULL, multiple = FALSE,
      options = list(placeholder = "Select a dataset...")
    ),
    radioButtons(
      ns("direction"), "Direction", inline = TRUE,
      choices = c("All", "Higher in Control", "Higher in FRDA"),
      selected = "All"
    ),
    tags$hr(style = "margin:1rem 0;"),

    # ---- help text ----
    h5("About this analysis"),
    tags$p(
      "Differential Transcript Usage (DTU) identifies genes whose",
      strong("isoforms (transcript variants)"),
      "change in relative abundance between experimental conditions."
    ),
    tags$p(
      "Isoform switching occurs when the dominant transcript of a gene differs ",
      "between conditions, for example, one isoform is more highly used in ",
      strong("FRDA") , "and another in", strong("Control"), "."
    ),
    tags$p(
      "DTU complements traditional differential expression by revealing changes ",
      "in transcript composition even when overall gene expression is unchanged."
    ),

    tags$hr(style = "margin:1rem 0;"),

    h5("Key terms"),
    tags$ul(
      tags$li(
        strong("Isoforms tested:"), " Number of transcript isoforms examined for DTU within this dataset."
      ),
      tags$li(
        strong("Switching isoforms:"), " Isoforms showing a significant difference in relative usage ",
        "(q ≤ 0.05 and |dIF| ≥ 0.1)."
      ),
      tags$li(
        strong("Genes with switch:"), " Genes containing at least one significantly switching isoform."
      ),
      tags$li(
        strong("IF (Isoform Fraction):"), " The proportion of a gene’s total expression accounted for by each isoform in a condition."
      ),
      tags$li(
        strong("dIF:"), " The change in isoform fraction between conditions ",
        "(e.g., IF_FRDA − IF_Control). A |dIF| ≥ 0.1 typically indicates a biologically relevant shift."
      ),
      tags$li(
        strong("q-value:"), " Adjusted p-value controlling the false discovery rate (FDR) for each isoform’s switch event."
      )
    ),

    h5("Tables:"),
    tags$li(
      strong("Top table = which genes are affected")
    ),
    tags$li(
      strong("Bottom table = which transcripts changed and how")
    ),

    tags$hr(style = "margin:1rem 0;"),

    h5("Downloads"),
    tags$div(
      style="display:flex; flex-direction:column; gap:8px;",
      downloadButton(ns("dl_summary"), "Summary CSV"),
      downloadButton(ns("dl_genes"),   "Gene-level CSV"),
      downloadButton(ns("dl_iso"),     "Isoform-level CSV")
    )
  )
}


dtuResultsMainUI <- function(id) {
  ns <- NS(id)
  tagList(
    br(),
    fluidRow(
      column(4, div(class="summary-box",
                    style="background:#f9fafb;border:1px solid #ddd;border-radius:12px;padding:16px;text-align:center;",
                    h5("Isoforms tested", style="margin-bottom:6px;"),
                    h2(textOutput(ns("v_iso")),   style="margin:0;color:#005249;"))),
      column(4, div(class="summary-box",
                    style="background:#f9fafb;border:1px solid #ddd;border-radius:12px;padding:16px;text-align:center;",
                    h5("Switching isoforms", style="margin-bottom:6px;"),
                    h2(textOutput(ns("v_sw")),    style="margin:0;color:#005249;"))),
      column(4, div(class="summary-box",
                    style="background:#f9fafb;border:1px solid #ddd;border-radius:12px;padding:16px;text-align:center;",
                    h5("Genes with switch", style="margin-bottom:6px;"),
                    h2(textOutput(ns("v_genes")), style="margin:0;color:#005249;")))
    ),
    br(),
    h5("Top switching genes"),
    DT::dataTableOutput(ns("tbl_genes")),
    br(),
    h5(uiOutput(ns("iso_title"))),
    DT::dataTableOutput(ns("tbl_iso"))
  )
}


# Server
dtuResultsServer <- function(id, pkg = utils::packageName(), labels = NULL, pretty_map = NULL){
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    # ---- ensure DTU data are available (Zenodo cache) ----
    ensure_atlas_data(
      keys    = "dtu_results",
      package = pkg
    )

    cache_root <- tools::R_user_dir(pkg, which = "cache")
    data_dir   <- file.path(cache_root, "dtu")

    # --- choices ---
    if (is.null(labels) || !length(labels)) labels <- discover_dtu_labels_flat(data_dir)
    updateSelectizeInput(
      session, "label",
      choices = dtu_choice_vector(labels, pretty_map),
      server = TRUE
    )

    # --- readers ---
    sum_tbl <- reactive({
      req(input$label)
      p <- dtu_paths_for_label_flat(data_dir, input$label)$summary
      req(file.exists(p))
      readRDS(p)
    })

    genes_tbl <- reactive({
      req(input$label)
      p <- dtu_paths_for_label_flat(data_dir, input$label)$genes
      req(file.exists(p))
      readRDS(p) |>
        dplyr::mutate(gene = dplyr::coalesce(.data$gene_name, .data$gene_id))
    })

    iso_tbl_raw <- reactive({
      req(input$label)
      p <- dtu_paths_for_label_flat(data_dir, input$label)$iso
      req(file.exists(p))
      readRDS(p)
    })

    # --- headline values ---
    output$v_iso   <- renderText({ format(sum_tbl()$nrIsoforms[[1]], big.mark = ",") })
    output$v_sw    <- renderText({ format(sum_tbl()$nrSwitches[[1]], big.mark = ",") })
    output$v_genes <- renderText({ format(sum_tbl()$nrGenes[[1]],   big.mark = ",") })

    # --- tables ---
    # 1) genes table (no buttons, use built-in search box)
    genes_tbl_view <- reactive({
      genes_tbl() |>
        dplyr::arrange(.data$gene_switch_q_value) |>
        dplyr::mutate(gene_switch_q_value = signif(.data$gene_switch_q_value, 3))
    })

    output$tbl_genes <- DT::renderDataTable({
      gt <- genes_tbl_view() |>
        # remove internal or redundant columns if present
        dplyr::select(-dplyr::any_of(c("gene_ref", "gene_id", "gene")))

      DT::datatable(
        gt,
        rownames = FALSE,
        selection = "single",
        options = list(
          pageLength = 10,
          scrollX = TRUE,
          dom = "tip"
        )
      )
    })


    # selected gene drives the isoform table filter
    sel_gene <- reactive({
      idx <- input$tbl_genes_rows_selected
      gt  <- genes_tbl_view()
      if (length(idx) == 1) gt[idx, c("gene_id","gene_name"), drop = FALSE] else NULL
    })

    # 2) isoform table (only Direction + selected gene)
    iso_tbl_view <- reactive({
      tt <- iso_tbl_raw()
      if (input$direction != "All") {
        tt <- tt |> dplyr::filter(.data$direction == input$direction)
      }
      sg <- sel_gene()
      if (!is.null(sg) && nrow(sg) == 1) {
        tt <- tt |>
          dplyr::filter(.data$gene_id == sg$gene_id | .data$gene_name == sg$gene_name)
      }
      tt |>
        dplyr::mutate(
          isoform_switch_q_value = signif(.data$isoform_switch_q_value, 3),
          dIF = signif(.data$dIF, 3),
          IF_condition_1 = signif(.data$IF_condition_1, 3),
          IF_condition_2 = signif(.data$IF_condition_2, 3)
        ) |>
        dplyr::arrange(.data$isoform_switch_q_value, dplyr::desc(abs(.data$dIF)))
    })

    output$iso_title <- renderUI({
      sg <- sel_gene()
      if (is.null(sg)) HTML("Isoform switches")
      else HTML(sprintf("Isoform switches for <b>%s</b> (%s)", sg$gene_name, sg$gene_id))
    })

    output$tbl_iso <- DT::renderDataTable({
      DT::datatable(
        iso_tbl_view(),
        rownames = FALSE,
        options = list(pageLength = 10, scrollX = TRUE, dom = "tip")
      )
    })

    # --- downloads (no references to removed inputs) ---
    output$dl_summary <- downloadHandler(
      filename = function() paste0(input$label, "_switch_summary.csv"),
      content  = function(file) readr::write_csv(sum_tbl(), file)
    )
    output$dl_genes <- downloadHandler(
      filename = function() paste0(input$label, "_gene_switches.csv"),
      content  = function(file) readr::write_csv(genes_tbl_view(), file)
    )
    output$dl_iso <- downloadHandler(
      filename = function() paste0(input$label, "_significant_isoforms.csv"),
      content  = function(file) readr::write_csv(iso_tbl_view(), file)
    )
  })
}
