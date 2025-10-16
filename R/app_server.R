#' Application Server Logic
#'
#' Defines the server-side logic for the FRDA Transcriptomic Atlas app.
#' @keywords internal
#' @import shiny
#' @importFrom ggplot2 ggplot aes geom_point
app_server <- function(input, output, session) {

  addResourcePath("pkgwww", system.file("www", package = "FRDATranscriptomicAtlas"))

  # Build long-format DEGs from your packaged results (padj ≤ 0.05 by default)
  deg_long_df <- build_deg_long(pkg = utils::packageName(), p_thr = 0.05, lfc_min = 0, level = "genes")

  ## Discover labels from flat files:
  # e.g. Maddock_SN_FA1_switch_summary.csv -> "Maddock_SN_FA1"
  discover_dtu_labels_flat <- function(data_dir) {
    fs <- list.files(data_dir, pattern = "_switch_summary\\.csv$", full.names = FALSE)
    unique(sub("_switch_summary\\.csv$", "", fs))
  }

  # Build shiny choices: show pretty names, keep label as value
  # pretty_map is your mapping named by label
  dtu_choice_vector <- function(labels, pretty_map) {
    disp <- unname(pretty_map[labels])
    # fallback to raw label if not in pretty_map
    disp[is.na(disp)] <- labels[is.na(pretty_map[labels])]
    stats::setNames(labels, disp)  # names = shown in UI, values = labels
  }

  # Paths for a given label (flat layout)
  dtu_paths_for_label_flat <- function(data_dir, label) {
    list(
      summary = file.path(data_dir, paste0(label, "_switch_summary.csv")),
      genes   = file.path(data_dir, paste0(label, "_gene_switches.csv")),
      iso     = file.path(data_dir, paste0(label, "_significant_isoforms.csv"))
    )
  }



  # --- Call modules -----------------------------------------------------
  degTablesServer("deg_tables")
  genePlotsServer("gene_plots")
  degVennServer("deg_venn", pkg = utils::packageName())
  volcanoServer("volc", level = "genes", pkg = "FRDATranscriptomicAtlas", custom_loader = deseq_loader)
  forestPlotsServer("forest")
  dtuResultsServer(
    id = "dtu",
    data_dir  = system.file("extdata/DTU", package = utils::packageName()),
    labels    = NULL,          # auto-discover from filenames
    pretty_map = pretty_map    # your mapping vector from utils
  )
  tpmHeatmapServer("tpm_hm",
    pkg = utils::packageName()
  )

  # ---- In server.R or inside your app_server() ----
  observe({
    shinyjs::toggle(id = "sidebar", condition = !(input$tabselected == 8))
  })

  observeEvent(input$tabselected, {
    # tab values arrive as strings
    if (identical(input$tabselected, "8")) {
      shinyjs::addClass(id = "main_wrap", class = "fullwidth")
    } else {
      shinyjs::removeClass(id = "main_wrap", class = "fullwidth")
    }
  })


}
