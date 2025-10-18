# ===== UI =====
#' Consequences — Sidebar UI
#' @noRd
# Discover dataset labels from files like "<LABEL>_genes_with_consequences.csv"
# and "<LABEL>_ISOFORMS_with_consequences.csv"
discover_conseq_labels_flat <- function(data_dir) {
  fs <- list.files(
    data_dir,
    pattern = "_(genes|ISOFORMS)_with_consequences\\.csv$",
    full.names = FALSE
  )
  labs <- sub("_(genes|ISOFORMS)_with_consequences\\.csv$", "", fs)
  sort(unique(labs))
}

conseq_paths_for_label_flat <- function(data_dir, label) {
  list(
    genes = file.path(data_dir, paste0(label, "_genes_with_consequences.csv")),
    iso   = file.path(data_dir, paste0(label, "_ISOFORMS_with_consequences.csv"))
  )
}



consequencesSidebarUI <- function(id) {
  ns <- NS(id)
  tagList(
    h4("Consequences"),
    selectInput(
      ns("label"),
      "Dataset",
      choices = NULL,  # populated in server from filenames
      multiple = FALSE
    ),
    radioButtons(
      ns("level"),
      label = "Show:",
      choices = c("Genes" = "genes", "Isoforms (transcripts)" = "isoforms"),
      selected = "genes",
      inline = TRUE
    ),
    tags$hr(),
    tags$strong("What are consequence results?"),
    helpText(HTML(
      paste0(
        "<p><b>DTU</b> tests whether the <i>relative usage</i> of isoforms differs between groups. ",
        "The <b>consequence</b> layer annotates the likely biological effect of those usage shifts, ",
        "such as coding status or domain changes.</p>"
      )
    )),
    tags$details(
      tags$summary("Key: consequence definitions (predictions)"),
      tags$p(
        HTML("<em>Important:</em> These are <strong>computational predictions</strong> derived from sequence and annotation tools (e.g., CPAT, PfamScan, SignalP, IsoformSwitchAnalyzeR’s NMD rule-set).
          They should be <strong>experimentally validated.</strong>")
      ),
      tags$ul(
        tags$li(HTML("<b>IR (Intron Retention)</b>: Evidence that an intron is retained/skipped between conditions; can alter the reading frame, localization, or mRNA stability.")),
        tags$li(HTML("<b>codingPotential</b>: switch alters likelihood a transcript encodes a protein (ORF gained/lost). Suggests a switch between coding vs non-coding output.")),
        tags$li(HTML("<b>ORF_seq_similarity</b>: ORF sequence is altered enough to affect similarity to the canonical protein (potentially impacting function/antibody binding).")),
        tags$li(HTML("<b>NMD Status</b>: Presence of a premature termination codon (PTC) or exon–exon junction context that predicts nonsense-mediated decay (reduced mRNA/protein).")),
        tags$li(HTML("<b>Domains Identified</b>: Gain/loss of annotated protein domains (e.g., Pfam) due to exon usage changes; may impact interaction or catalytic sites.")),
        tags$li(HTML("<b>PTC</b> (Premature Termination Codon): early stop; may trigger NMD and reduce protein output.")),
        tags$li(HTML("<b>Signal Peptide Identified</b>: Gain/loss of an N-terminal signal peptide (e.g., SignalP) suggesting altered secretion/ER targeting and subcellular routing."))
      )
    ),

    tags$hr(),
    tags$p(
      HTML(
        "<b>DTU vs Consequences:</b> DTU tells you <i>which genes/isoforms switch</i> (dIF, q-value). ",
        "Consequences explain <i>why it might matter</i> (coding status, domains, PTC, IR)."
      )
    ),
    tags$details(
      tags$summary("Key: statistical metrics"),
      tags$ul(
        tags$li(HTML("<b>dIF (delta Isoform Fraction):</b> Change in isoform usage between FRDA and Control. Positive = higher in FRDA, negative = higher in Control.")),
        tags$li(HTML("<b>max_abs_dIF:</b> The largest absolute dIF among a gene’s isoforms — summarises overall switch magnitude.")),
        tags$li(HTML("<b>q_use:</b> Adjusted p-value testing whether an isoform’s relative usage differs between groups.")),
        tags$li(HTML("<b>min_q:</b> The smallest q_use for any isoform of a gene — strongest evidence for switching at the gene level."))
      )
    )

  )
}

# Main panel: table (5 rows/page) + downloads
consequencesMainUI <- function(id) {
  ns <- NS(id)
  tagList(
    h4(textOutput(ns("title_text")), style = "margin-top: 0;"),
    DT::DTOutput(ns("tbl")),
    br(),
    fluidRow(
      column(
        width = 6,
        downloadButton(ns("dl_filtered"), "Download filtered (current view)")
      ),
      column(
        width = 6,
        downloadButton(ns("dl_full"), "Download full file")
      )
    )
  )
}
#' Consequences — Server (flat layout, same style as DTU module)
#' @param data_dir folder containing <LABEL>_genes_with_consequences.csv and <LABEL>_ISOFORMS_with_consequences.csv
#' @param labels optional vector of labels; if NULL they are auto-discovered from data_dir
#' @param pretty_map optional named vector: names = raw labels, values = display names
consequencesServer <- function(id, data_dir, labels = NULL, pretty_map = NULL) {
  moduleServer(id, function(input, output, session) {

    # --- populate dataset choices (labels) ---
    observe({
      labs <- labels
      if (is.null(labs) || !length(labs)) {
        labs <- discover_conseq_labels_flat(data_dir)
      }
      validate(need(length(labs) > 0, paste0(
        "No consequence files found in: ", normalizePath(data_dir, winslash = "/"),
        "\nExpected files like: <LABEL>_genes_with_consequences.csv and <LABEL>_ISOFORMS_with_consequences.csv"
      )))
      updateSelectInput(session, "label", choices = dtu_choice_vector(labs, pretty_map))

    })

    # --- title reflects dataset + level ---
    output$title_text <- renderText({
      req(input$label, input$level)
      lvl  <- if (identical(input$level, "genes")) "Gene" else "Isoform"
      disp <- if (!is.null(pretty_map) && !is.na(pretty_map[input$label])) pretty_map[input$label] else input$label
      sprintf("%s-level consequence results — %s", lvl, disp)
    })

    # --- data loader for the selected label + level ---
    dat <- reactive({
      req(input$label, input$level)
      paths <- conseq_paths_for_label_flat(data_dir, input$label)
      f <- if (identical(input$level, "genes")) paths$genes else paths$iso
      validate(
        need(file.exists(f),
             paste0("Missing file for ", input$label, " (", input$level, "): ", basename(f),
                    "\nSearched in: ", normalizePath(dirname(f), winslash = "/"),
                    "\nNote: filenames are case-sensitive on Linux."))
      )
      readr::read_csv(f, show_col_types = FALSE)
    })

    # --- table (5 rows/page) ---
    output$tbl <- DT::renderDT({
      DT::datatable(
        dat(),
        rownames = FALSE,
        filter   = "top",
        options  = list(
          pageLength = 5,
          lengthMenu = c(5, 10, 25, 50),
          scrollX    = TRUE,
          dom        = "ltip"
        )
      )
    })

    # --- downloads ---
    output$dl_filtered <- downloadHandler(
      filename = function() paste0("consequences_", input$label, "_", input$level, "_filtered_", Sys.Date(), ".csv"),
      content  = function(file) {
        idx <- input$tbl_rows_all  # respects current search/filter
        df  <- dat()
        readr::write_csv(if (!is.null(idx)) df[idx, , drop = FALSE] else df, file)
      }
    )

    output$dl_full <- downloadHandler(
      filename = function() paste0("consequences_", input$label, "_", input$level, "_full_", Sys.Date(), ".csv"),
      content  = function(file) {
        paths <- conseq_paths_for_label_flat(data_dir, input$label)
        f <- if (identical(input$level, "genes")) paths$genes else paths$iso
        file.copy(from = f, to = file, overwrite = TRUE)
      }
    )
  })
}
