utils::globalVariables(c(
  # --- identifiers ---
  "gene_id", "gene_name", "transcript_id", "isoform_switch_q_value",
  "ensembl_gene_id", "external_gene_name", "symbol", "Gene",

  # --- differential expression ---
  "padj", "pvalue", "qvalue", "log2FoldChange", "log2FC", "negLog10P",
  "status", "direction", "outline_class", "size_scaled",

  # --- DTU ---
  "dIF", "yi", "sei", "sd",

  # --- PCA ---
  "PC1", "PC2", "label",

  # --- TPM / expression ---
  "TPM", "dataset", "study", "Study", "Group",
  "sample_id", "sample_size", "sample_source",

  # --- heatmaps / ordering ---
  "grp_rank", "ds_rank", "status_rank",
  "tissue_ord", "ctype_ord", "gene_lab",

  # --- datasets table columns (non-syntactic names) ---
  "gse",
  "Study (Year)",
  "ENA (PRJNA)",
  "Sampling & design",
  "Sample size",

  # --- forest plots ---
  "gene",

  # --- genePlots ---
  "fa_num",
  "ic_first",
  "condition",

  # --- tables / metadata ---
  "Description", "Publication", "Sample", "Sum",
  "ENA", "GEO", "PRJNA", "prjna",
  "first_author_year", "pub_title", "pub_url",

  # --- volcano / hover ---
  "hover_txt", "display_label", "row_id", "P",

  # --- general ---
  "feature", "level", "path", "n", "n_studies", "score",
  "p_str", "desc", "up", "down",

  # --- app-level objects ---
  "PRETTY_MAP_LOCAL", "check_atlas_updates", "cran_install_if_missing"
))
