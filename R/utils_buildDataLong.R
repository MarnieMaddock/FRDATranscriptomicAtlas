# Build long-format DEG table (dataset, gene, direction)
build_deg_long <- function(pkg = utils::packageName(),
                           p_thr = 0.05, lfc_min = 0, level = "genes") {
  stopifnot(level %in% c("genes", "transcripts"))

  dir_deg <- system.file(file.path("extdata", "deg", level), package = pkg, mustWork = FALSE)
  files <- if (nzchar(dir_deg)) {
    list.files(dir_deg, pattern = "^DESEQ2_res_.+_0\\.[0-9]+_all_(genes|transcripts)\\.rds$", full.names = TRUE)
  } else character(0)

  tx2_path <- system.file("extdata/maps/tx2gene.tsv", package = pkg, mustWork = FALSE)
  tx2 <- if (nzchar(tx2_path) && file.exists(tx2_path)) readr::read_tsv(tx2_path, col_types = "ccc") else NULL
  gene_map <- if (!is.null(tx2)) dplyr::distinct(tx2, gene_id, gene_name) else NULL

  if (!length(files)) return(tibble::tibble(dataset = character(), gene = character(), direction = character()))

  dfs <- lapply(files, function(fp) {
    m <- stringr::str_match(basename(fp),
                            "^DESEQ2_res_(.+)_(0\\.[0-9]+)_all_(genes|transcripts)\\.rds$")
    dataset <- m[2]; lev <- m[4]
    if (!identical(lev, level)) return(NULL)

    x <- readRDS(fp)
    if (!is.data.frame(x)) x <- as.data.frame(x)

    # normalize column names
    if (!"log2FoldChange" %in% names(x)) {
      if ("log2FC" %in% names(x))    x$log2FoldChange <- x$log2FC
      else if ("beta" %in% names(x)) x$log2FoldChange <- x$beta
    }
    if (!"padj" %in% names(x) && "qvalue" %in% names(x)) x$padj <- x$qvalue

    # IDs BEFORE filtering
    id_all <- if (identical(level, "genes")) {
      if ("ensembl_gene_id" %in% names(x)) x$ensembl_gene_id else rownames(x)
    } else {
      if ("transcript_id" %in% names(x)) x$transcript_id else rownames(x)
    }

    lfc  <- x$log2FoldChange
    keep <- !is.na(x$padj) & x$padj <= p_thr &
      !is.na(lfc)   & abs(lfc) >= lfc_min

    # Optional: drop exact 0 so direction is unambiguous when lfc_min == 0
    keep <- keep & (lfc != 0)

    if (!any(keep)) return(NULL)

    # subset EVERYTHING consistently
    x_f   <- x[keep, , drop = FALSE]
    id    <- id_all[keep]
    lfc_f <- lfc[keep]

    # map labels
    gene_lbl <- if (identical(level, "genes")) {
      if (!is.null(gene_map)) {
        sym <- gene_map$gene_name[match(id, gene_map$gene_id)]
        dplyr::coalesce(sym, id)
      } else id
    } else {
      id
    }

    direction <- dplyr::if_else(lfc_f > 0, "Up", "Down")

    tibble::tibble(
      dataset   = dataset,
      gene      = gene_lbl,
      direction = direction
    )
  })

  dplyr::bind_rows(dfs)
}
