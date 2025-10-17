#' @importFrom shiny HTML
#' @importFrom shiny tags
# =========================
# Switchplot help + rationale + OneDrive link + example image
# =========================
switchplotsHelpUI <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    tags$style(HTML("
      .sp-help h4 { margin-top: 0 }
      .sp-kbd { background:#f5f5f5; border:1px solid #e3e3e3; border-radius:6px; padding:2px 6px; font-family: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, 'Liberation Mono', monospace; }
      .sp-img { max-width:100%; height:auto; border:1px solid #eee; border-radius:12px; box-shadow:0 1px 4px rgba(0,0,0,0.06); }
      .sp-callout { background:#f8fbfa; border-left:4px solid #2a776d; padding:10px 12px; border-radius:8px; }
      .sp-grid { display:grid; grid-template-columns: 1fr; gap:14px; }
      @media (min-width: 900px) {
        .sp-grid { grid-template-columns: 0.6fr 1.4fr; align-items:start; }
      }
      .sp-li { margin-bottom:6px; }
      .sp-img {
            max-width: 1200px;
            width: 100%;
            height: auto;
            display: block;
            margin: 0 auto;
          }
    ")),

    # ---- Introductory well panel ----
    fluidRow(
      column(
        width = 12,
        shiny::wellPanel(
          shiny::HTML("
  <p>
    Switchplots have been generated for all significant differential transcript usage (DTU) events across every dataset.
    However, including thousands of plots directly within this web application would considerably increase load times and require approximately 2.5&nbsp;GB of additional storage.
    To maintain optimal performance, the complete collection of plots is hosted externally and can be accessed below.
  </p>
  <p style='margin-top:10px;'>
    <b>Switchplots:</b> To browse or download the full collection, visit the
    <a href='https://1drv.ms/f/s!abc123yourlink' target='_blank'
       style='color:#005249; text-decoration:underline; font-weight:bold;'>
       OneDrive folder <i class='fa fa-external-link-alt'></i>
    </a>

    Each dataset has its own subfolder containing the corresponding SVG plots for all significant isoform switches.
  </p>
"),
          HTML("
  <div style='text-align:center; margin-top:15px;'>
    <a href='https://1drv.ms/f/s!abc123yourlink' target='_blank'
       style='background:#005249; color:white; padding:10px 20px; border-radius:6px;
              text-decoration:none; font-weight:600;'>
       🔗 Open OneDrive Folder
    </a>
  </div>
")


        )
      )
    ),

    # ---- Detailed rationale and example ----
    shiny::wellPanel(class = "sp-help",
                     h4("Understanding Isoform Switchplots"),
                     div(class = "sp-grid",

                         # Left column: explanation
                         div(
                           div(class = "sp-callout",
                               tags$p(tags$b("Why this analysis is valuable")),
                               tags$ul(
                                 tags$li(class="sp-li", "Highlights transcript-level regulation that is hidden when analysing only gene-level expression."),
                                 tags$li(class="sp-li", "Reveals functional consequences of isoform switching, such as domain gain/loss, altered localisation signals, or nonsense-mediated decay (NMD) activation."),
                                 tags$li(class="sp-li", "Quantifies changes in isoform usage through ",
                                         tags$span(class="sp-kbd","ΔIF (dIF)"),
                                         ", providing a measure of the magnitude of the switch between conditions."),
                                 tags$li(class="sp-li", "Facilitates mechanistic interpretation by integrating transcript architecture with predicted protein features.")
                               )
                           ),

                           # ---- Replace your existing tags$details(...) block with this ----
                           tags$details(open = TRUE,
                                        tags$summary(tags$b("How to interpret the plot")),
                                        tags$ol(
                                          tags$li(
                                            tags$b("Isoform structure and protein features (top panel):"),
                                            " Each horizontal track represents a transcript isoform. The top and bottom tracks show isoforms that increase or decrease in usage between conditions. ",
                                            "Coloured blocks represent predicted protein domains or sequence features, matching the legends on the right (e.g., Signal Peptide, Cadherin domains, Transmembrane helix). ",
                                            "Arrows connect exons to indicate introns and alternative splicing patterns. ",
                                            "Labels such as 'Coding' or 'NMD' denote whether the isoform is protein-coding or predicted to undergo nonsense-mediated decay."
                                          ),
                                          tags$li(
                                            tags$b("Topology and domain legends (right-hand side):"),
                                            " These legends describe the spatial orientation (e.g., intracellular, extracellular, signal peptide, transmembrane helix) and specific Pfam or InterPro domains present in the protein. ",
                                            "Comparing the colours between isoforms reveals which functional domains are gained or lost."
                                          ),
                                          tags$li(
                                            tags$b("Gene Expression (bottom left):"),
                                            " Shows overall gene-level expression per condition. ",
                                            "This reflects total RNA abundance regardless of isoform composition, allowing you to see whether total gene expression changes alongside isoform switching."
                                          ),
                                          tags$li(
                                            tags$b("Isoform Expression (bottom middle):"),
                                            " Displays expression values for each transcript separately (e.g., TPM or normalized counts). ",
                                            "Helps determine which isoforms are strongly expressed and whether both are biologically relevant. "
                                          ),
                                          tags$li(
                                            tags$b("Isoform Usage (bottom right):"),
                                            " Shows the proportion of each isoform relative to total gene expression (Isoform Fraction, IF). ",
                                            "Bars and asterisks indicate ΔIF (difference in IF) and its statistical significance; larger |ΔIF| and more asterisks denote stronger, significant switching."
                                          ),
                                          tags$li(
                                            tags$b("Overall interpretation:"),
                                            " Identify which isoform becomes more prevalent in each condition and whether this corresponds to a change in predicted protein structure or localisation. ",
                                            "A biologically meaningful switch occurs when an isoform with different domains, NMD status, or topology replaces another, suggesting altered function without necessarily changing total gene expression."
                                          ),
                                          tags$li(
                                            tags$b("Additional notes:"),
                                            " Isoform switchplots are generated only for transcripts with significant differential usage (adjusted p-value ≤ 0.05) between the specified conditions. P values are not displayed on the gene expression or isoform expression graphs."
                                        )
                                        )
                           ),


                           div(
                             tags$p(
                               HTML("To explore all switchplots, access the&nbsp;<a href='https://1drv.ms/f/s!abc123yourlink' target='_blank'>OneDrive folder</a> (≈2.5&nbsp;GB), where each dataset’s results are organised into separate subfolders.")
                             )
                           )
                         ),

                         # Right column: example image
                         div(
                           tags$figure(
                             tags$img(
                               src = "pkgwww/switchplot_example.svg",  # note the 'pkgwww/' prefix
                               class = "sp-img",
                               alt   = "Example switchplot"
                             ),


                             tags$figcaption(
                                 HTML(
                                   "<b>Figure:</b> Example isoform switchplot for <i>PCDHGA10</i> (Control vs&nbsp;FRDA).
      The top panel shows the structural organisation of each transcript isoform, with coloured blocks indicating predicted protein domains and topology
      (extracellular, intracellular, signal peptide, transmembrane helix).
      In FRDA, isoform <b>ENST00000398810</b> is upregulated, while <b>ENST00000612503</b> decreases in usage, resulting in a significant isoform switch.
      The FRDA-enriched isoform retains multiple Cadherin domains and a transmembrane region, suggesting altered membrane localisation or adhesion potential.
      Lower panels show total gene expression, isoform-specific expression, and isoform fractions (IF) with significance asterisks denoting adjusted p-values."
                                 )

                             )
                           )
                         )
                     )
    )
  )
}
