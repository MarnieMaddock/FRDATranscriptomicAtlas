#' Application Server Logic
#'
#' Defines the server-side logic for the FRDA Transcriptomic Atlas app.
#' @keywords internal
#' @import shiny
#' @importFrom ggplot2 ggplot aes geom_point
app_server <- function(input, output, session) {

  # --- make `pkg` robust for both project + installed package modes ----
  pkg <- tryCatch(utils::packageName(), error = function(e) "")
  if (!length(pkg) || !is.character(pkg) || !nzchar(pkg)) pkg <- "FRDATranscriptomicAtlas"
  pkg <- pkg[[1L]]  # ensure length 1

  pkg <- utils::packageName()
  pkg_www <- system.file("www", package = pkg, mustWork = FALSE)
  if (nzchar(pkg_www) && dir.exists(pkg_www)) {
    shiny::addResourcePath("pkgwww", pkg_www)  # /pkgwww → <package>/inst/www
  }

  # ---- helper for image ----
  # get_switchplot_example_src <- function(pkg) {
  #   # prefer packaged www
  #   pkg_file <- system.file("www", "switchplot_example.svg", package = pkg, mustWork = FALSE)
  #   if (nzchar(pkg_file) && file.exists(pkg_file)) {
  #     return("pkgwww/switchplot_example.svg")   # web-accessible path
  #   }
  #
  #   # fallback: if you’re running locally from the project root
  #   if (file.exists("www/switchplot_example.svg")) {
  #     return("switchplot_example.svg")          # Shiny serves /www automatically
  #   }
  #
  #   NULL
  # }

  # ---- Build long-format DEGs (uses robust `pkg`) ----
  deg_long_df <- build_deg_long(pkg = pkg, p_thr = 0.05, lfc_min = 0, level = "genes")

  # --- Optional helpers you already have ---
  discover_dtu_labels_flat <- function(data_dir) {
    fs <- if (nzchar(data_dir) && dir.exists(data_dir))
      list.files(data_dir, pattern = "_switch_summary\\.csv$", full.names = FALSE)
    else character(0)
    if (!length(fs)) return(character(0))
    unique(sub("_switch_summary\\.csv$", "", fs))
  }

  dtu_choice_vector <- function(labels, pretty_map) {
    disp <- unname(pretty_map[labels]); disp[is.na(disp)] <- labels[is.na(pretty_map[labels])]
    stats::setNames(labels, disp)
  }

  dtu_paths_for_label_flat <- function(data_dir, label) {
    list(
      summary = file.path(data_dir, paste0(label, "_switch_summary.csv")),
      genes   = file.path(data_dir, paste0(label, "_gene_switches.csv")),
      iso     = file.path(data_dir, paste0(label, "_significant_isoforms.csv"))
    )
  }

  get_DTU_dir <- function() {
    p <- file.path("inst", "extdata", "DTU")
    if (dir.exists(p)) return(p)
    system.file("extdata", "DTU", package = pkg)
  }

  # ---- Consequences helpers (flat layout; same style as your DTU code) ----

  discover_conseq_labels_flat <- function(data_dir) {
    fs <- if (nzchar(data_dir) && dir.exists(data_dir))
      list.files(
        data_dir,
        pattern = "_(genes|ISOFORMS)_with_consequences\\.csv$",
        full.names = FALSE
      )
    else character(0)
    if (!length(fs)) return(character(0))
    labs <- sub("_(genes|ISOFORMS)_with_consequences\\.csv$", "", fs)
    sort(unique(labs))
  }

  # reuse your existing dtu_choice_vector, or fall back safely if pretty_map is NULL
  conseq_choice_vector <- function(labels, pretty_map = NULL) {
    if (is.null(pretty_map)) return(stats::setNames(labels, labels))
    disp <- unname(pretty_map[labels])
    disp[is.na(disp)] <- labels[is.na(disp)]
    stats::setNames(labels, disp)
  }

  conseq_paths_for_label_flat <- function(data_dir, label) {
    list(
      genes = file.path(data_dir, paste0(label, "_genes_with_consequences.csv")),
      iso   = file.path(data_dir, paste0(label, "_ISOFORMS_with_consequences.csv"))
    )
  }

  # same resolver you used before
  get_DTU_dir <- function() {
    p <- file.path("inst", "extdata", "DTU")
    if (dir.exists(p)) return(p)
    # 'pkg' should be defined in your app_server (your package name)
    system.file("extdata", "DTU", package = pkg)
  }

  # --- ensure pretty_map exists (fallback is harmless) ---
  if (!exists("pretty_map", inherits = TRUE)) {
    pretty_map <- setNames(character(0), character(0))
  }

  # --- Call modules (pass the same `pkg`) -------------------------------
  aboutServer("about")
  pcaServer("pca")
  genePlotsServer("gene_plots", pkg = pkg)
  degTablesServer("deg_tables", pkg = pkg)
  degVennServer("deg_venn", pkg = pkg)
  volcanoServer("volc", level = "genes", pkg = pkg)
  forestPlotsServer("forest", pkg = pkg)
  dtuResultsServer(
    id         = "dtu",
    data_dir   = get_DTU_dir(),
    labels     = NULL,
    pretty_map = pretty_map
  )
  dtuVennServer("dtuVenn")
  tpmHeatmapServer("tpm_hm", pkg = pkg)
  consequencesServer(
    id         = "dtu_func_cons",
    data_dir   = get_DTU_dir(),
    labels     = NULL,              # we auto-discover with discover_conseq_labels_flat()
    pretty_map = pretty_map         # optional: same pretty labels you use elsewhere
  )

  # --- sidebar behavior for SwitchPlots full-width tab ------------------
  observe({
    shinyjs::toggle(
      id = "sidebar",
      condition = !(input$tabselected %in% c(0, 8))
    )
  })

  observeEvent(input$tabselected, {
    if (input$tabselected %in% c(0, 8)) {
      shinyjs::addClass(id = "main_wrap", class = "fullwidth")
    } else {
      shinyjs::removeClass(id = "main_wrap", class = "fullwidth")
    }
  })

}
