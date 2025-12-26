#' Loader for DESeq2 result objects saved as .rds
#' @param dataset key like "Erwin", "Mishra_FF1" (must match file stem)
#' @param level   "genes" or "transcripts"
#' @param pkg     package name (for system.file lookup)
#' @return data.frame with columns: gene, log2FC, padj (if present), pvalue (if present)
#' @keywords internal
#' @noRd
deseq_loader <- function(dataset, level = "genes", pkg = utils::packageName()) {

  stopifnot(level %in% c("genes", "transcripts"))

  # Folder inside the package
  dir_path <- system.file("extdata", "deg", level, package = pkg)
  if (!nzchar(dir_path) || !dir.exists(dir_path)) {
    stop("Data folder not found: ", dir_path)
  }

  # Files are like: DESEQ2_res_<DATASET>_0.05_all_<level>.rds
  pat <- paste0("^DESEQ2_res_", dataset, "_.*_all_", level, "\\.rds$")
  files <- list.files(dir_path, pattern = pat, full.names = TRUE)
  if (length(files) == 0L) {
    stop("No DESeq2 results found for dataset '", dataset, "' in: ", dir_path,
         "\nPattern: ", pat)
  }
  # If multiple match (different alphas), pick the first or the newest:
  f <- files[which.max(file.mtime(files))]

  obj <- readRDS(f)

  # Coerce to data.frame regardless of class
  df <- if (inherits(obj, "DESeqResults")) {
    as.data.frame(obj)
  } else if (is.data.frame(obj)) {
    obj
  } else {
    as.data.frame(obj)
  }

  # Try to find an ID column; fallback to rownames
  id_col <- NULL
  if (level == "genes") {
    for (cand in c("symbol", "gene", "gene_name", "gene_symbol", "Gene")) {
      if (cand %in% names(df)) { id_col <- cand; break }
    }
    id_label <- "gene"
  } else {
    for (cand in c("tx", "transcript", "transcript_id", "transcript_name", "target_id")) {
      if (cand %in% names(df)) { id_col <- cand; break }
    }
    id_label <- "gene"  # the volcano module expects column named 'gene'; reuse it for tx id
  }
  if (is.null(id_col)) {
    df[[id_label]] <- rownames(df)
  } else {
    df[[id_label]] <- df[[id_col]]
  }

  # Standardise effect / p columns
  if (!"log2FoldChange" %in% names(df)) {
    stop("Expected 'log2FoldChange' in DESeq2 results for ", dataset, " (", level, ")")
  }
  out <- data.frame(
    gene    = df[[id_label]],
    log2FC  = as.numeric(df$log2FoldChange),
    padj    = if ("padj"   %in% names(df)) as.numeric(df$padj)   else NA_real_,
    pvalue  = if ("pvalue" %in% names(df)) as.numeric(df$pvalue) else NA_real_,
    stringsAsFactors = FALSE
  )
  # Keep unique rows, drop complete NAs
  out <- out[!is.na(out$log2FC) & ( !is.na(out$padj) | !is.na(out$pvalue) ), ]
  rownames(out) <- NULL
  out
}
