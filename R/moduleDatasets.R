# R/modules/datasets_module.R
# ---------------------------
# R/modules/datasets_module.R
# ---------------------------

datasetsUI <- function(id) {
  ns <- shiny::NS(id)

  # helper to build external links safely
  alink <- function(label, href) {
    if (is.na(href) || !nzchar(href)) return("&mdash;")
    as.character(tags$a(label, href = href, target = "_blank", rel = "noopener"))
  }
       ena_link <- function(prjna) ifelse(
        nzchar(prjna),
         paste0("https://www.ebi.ac.uk/ena/browser/view/", prjna),
         NA_character_
      )

       geo_link <- function(gse) ifelse(
         nzchar(gse),
         paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", gse),
         NA_character_
       )

  sra_link <- function(prjna) ifelse(
    nzchar(prjna),
    paste0("https://www.ncbi.nlm.nih.gov/Traces/study/?acc=", prjna),
    NA_character_
  )

      df <- tibble::tribble(
        ~first_author_year, ~prjna,        ~gse,         ~sample_source,                                                                                         ~sample_size,                                                                                                   ~pub_title,                                                                                                                                      ~pub_url,
        "Chutake YK (2014)", "Data obtained from authors", "",           "Primary lymphoblastoid cell lines",                                                   "FRDA (N = 2, n = 1); CTRL (N = 2, n = 1) &#8225",                                                                  "Altered nucleosome positioning at the transcription start site and deficient transcriptional initiation in friedreich ataxia", "https://pubmed.ncbi.nlm.nih.gov/24737321/",
        "Erwin GS (2017)",  "PRJNA388394", "GSE99400",   "Primary lymphoblastoid cell lines",                                                                     "FRDA (N = 1, n = 4); CTRL (N = 1, n = 4) *",                                                                  "Synthetic transcription elongation factors license transcription across repressive chromatin",                                                    "https://pmc.ncbi.nlm.nih.gov/articles/PMC6037176/",
        "Indelicato E (2023)","PRJNA940990","GSE226646", "Primary skeletal muscle (gastrocnemius biopsy)",                                                        "FRDA (N = 7, n = 1); CTRL (N = 6, n = 1) &#8224;",                                                                   "Skeletal muscle transcriptomics dissects the pathogenesis of FRDA",                                                                               "https://academic.oup.com/hmg/article/32/13/2241/7110888",
        # Keep this Lai row on ONE LINE
        "Lai J (2018)",     "PRJNA495860", "",           "iPSCs; iPSC-derived PNS & CNS neurons",                                                                  "FRDA iPSCs (N = 1, n = 3); IC iPSCs (N = 1, n = 3); FRDA PNS (N = 1, n = 3); IC PNS (N = 1, n = 6); FRDA CNS (N = 2, n = 3); IC CNS (N = 1, n = 3) &#167;", "Transcriptional profiling of isogenic Friedreich ataxia neurons and effect of an HDAC inhibitor on disease signatures", "https://pmc.ncbi.nlm.nih.gov/articles/PMC6369281/",
        "Lees JG (2025*)",  "PRJNA1307399","GSE305638",  "iPSC-derived cardiomyocytes",                                                                            "FRDA (N = 3, n = 3); IC (N = 3, n = 3) &#167;",                                                                      "Frataxin deficiency drives cardiac dysfunction and transcriptional dysregulation in Friedreich ataxia iPSC model",                              "https://www.biorxiv.org/content/10.1101/2025.08.20.671405v1",
        "Li J (2019)", "Data obtained from authors", "", "iPSC-derived cardiomyocytes",                                                                 "FRDA (N = 1, n = 2); CTRL (N = 2, n = 1); IC (N = 1, n = 2);",                                                                  "Excision of the expanded GAA repeats corrects cardiomyopathy phenotypes of iPSC-derived Friedreich's ataxia cardiomyocytes", "https://pubmed.ncbi.nlm.nih.gov/31446150/",
        "Maddock ML (*)",   "PRJNA1425480",     "GSE319887", "iPSC-derived neural crest cells; sensory neurons; lower motor neurons",    "NCC: FRDA (N = 2, n = 4&ndash;5); IC (N = 2, n = 4); SN: FRDA (N = 2, n = 4&ndash;5); IC (N = 2, n = 4&ndash;5); LMN: FRDA (N = 1, n = 4); IC (N = 1, n = 4) &#167;",   "Friedreich ataxia transcriptomic dysregulation and identification of cell type-specific biomarkers: A systematic review and meta-analysis",  "https://doi.org/10.64898/2026.03.18.712785",
        "Mishra P (2024)",  "PRJNA1025300","GSE244886",  "iPSC-derived neurons",                                                                                   "FRDA (N = 4, n = 2&ndash;4); IC (N = 4, n = 2&ndash;4) &#182;",                                                                  "Gene editing improves endoplasmic reticulum-mitochondrial contacts and unfolded protein response in Friedreich's ataxia iPSC-derived neurons",   "https://pmc.ncbi.nlm.nih.gov/articles/PMC10899513/",
        "Napierala JS (2017)","PRJNA412241","GSE104288", "Primary fibroblasts",                                                                                    "FRDA (N = 18, n = 1); CTRL (N = 17, n = 1) &#8225;",                                                                 "Comprehensive analysis of gene expression patterns in Friedreich's ataxia fibroblasts by RNA sequencing reveals altered levels of protein synthesis factors and solute carriers", "https://pmc.ncbi.nlm.nih.gov/articles/PMC5719256/",
        "Vilema-Enriquez G (2020)","PRJNA606059","GSE145115","Primary fibroblasts",                                                                                "FRDA (N = 1, n = 3); CTRL (N = 1, n = 3) &#8225;",                                                                   "Inhibition of the SUV4-20 H1 histone methyltransferase increases frataxin expression in Friedreich's ataxia patient cells",                      "https://pmc.ncbi.nlm.nih.gov/articles/PMC7939392/",
        "Wang F (2022)", "PRJNA846268", "GSE205526", "Primary fibroblasts",                                                                                       "FRDA (N = 1, n = 3); CTRL (N = 1, n = 3) &#8756;",                                       "G-rich motifs within phosphorothioate-based antisense oligonucleotides (ASOs) drive activation of FXN expression through indirect effects",                        "https://pmc.ncbi.nlm.nih.gov/articles/PMC9825156/"
      )

      # ---- build linkified columns ----
      df_disp <- df |>
        dplyr::mutate(
          `Study (Year)`      = first_author_year,
          `Publication`       = purrr::pmap_chr(list(pub_title, pub_url), ~ alink(..1, ..2)),
          `GEO`               = purrr::map2_chr(gse, gse, ~ alink(.y, geo_link(.x))),
          `ENA (PRJNA)`       = purrr::map2_chr(prjna, prjna, ~ alink(.y, ena_link(.x))),
          `SRA (runs)`        = purrr::map2_chr(prjna, prjna, ~ alink("Run selector", sra_link(.x))),
          `Sampling & design` = sample_source,
          `Sample size`       = sample_size
        ) |>
        dplyr::select(`Study (Year)`, `ENA (PRJNA)`, GEO, Publication, `Sampling & design`, `Sample size`)

  tagList(
    h3("Datasets included"),

    DT::datatable(
      df_disp,
      escape = FALSE,
      rownames = FALSE,
      options = list(
        dom = "Bfrtip",
        pageLength = 10,
        buttons = c("copy", "excel"),
        scrollX = TRUE
      ),
      extensions = "Buttons"
    ),
    br(),
    br(),

    div(
      class = "text-muted",
      style = "margin-bottom:10px; font-size:0.9em;",
      shiny::HTML("
        <strong>Group definitions:</strong><br>
        * = FRDA = patient-derived; CTRL = healthy sibling<br>
        &#8224; = FRDA = patient-derived; CTRL = matched control<br>
        &#8756; = FRDA = patient-derived; CTRL = wild-type control<br>
        &#8225; = FRDA = patient-derived; CTRL = unaffected control<br>
        &#167; = FRDA = patient-derived; IC = isogenic control<br>
        &#182; = FRDA = patient-derived; IC = isogenic control (GAA excision, or GAA expansion replaced by healthy allele (E35))<br>
        <br>
        <strong>Sample size:</strong> N = biological replicates (different genetic backgrounds); n = technical replicates.
      ")
    ),
    br(),
    br()


  )
}

