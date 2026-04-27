get_deg_data <- function(
    dataset,
    level,
    file_path = NULL,
    data_mode = "local",
    padj_max = NULL,
    lfc_min = 0,
    direction = "both"
) {

  if (identical(data_mode, "local")) {
    validate(
      need(!is.null(file_path) && file.exists(file_path), "No results file found.")
    )

    x <- readRDS(file_path)
    if (!is.data.frame(x)) x <- as.data.frame(x)

    # Harmonise columns
    if (!"log2FoldChange" %in% names(x)) {
      if ("log2FC" %in% names(x)) {
        x$log2FoldChange <- x$log2FC
      } else if ("beta" %in% names(x)) {
        x$log2FoldChange <- x$beta
      }
    }

    if (!"padj" %in% names(x)) {
      if ("qvalue" %in% names(x)) {
        x$padj <- x$qvalue
      } else if ("adj.P.Val" %in% names(x)) {
        x$padj <- x$adj.P.Val
      }
    }

    # Apply p-value filter
    if (!is.null(padj_max) && !is.na(padj_max) && "padj" %in% names(x)) {
      x <- x[!is.na(x$padj) & x$padj <= padj_max, , drop = FALSE]
    }

    # Apply logFC / direction filter
    if ("log2FoldChange" %in% names(x)) {
      if (identical(direction, "up")) {
        x <- x[!is.na(x$log2FoldChange) & x$log2FoldChange >= lfc_min, , drop = FALSE]
      } else if (identical(direction, "down")) {
        x <- x[!is.na(x$log2FoldChange) & x$log2FoldChange <= -lfc_min, , drop = FALSE]
      } else {
        x <- x[!is.na(x$log2FoldChange) & abs(x$log2FoldChange) >= lfc_min, , drop = FALSE]
      }
    }

    return(x)
  }

  if (identical(data_mode, "cloud")) {
    return(get_deg_data_cloud_rds_cached(
      dataset_id = dataset,
      feature_level = level,
      padj_max = padj_max,
      lfc_min = lfc_min,
      direction = direction
    ))
  }

  stop("Unknown data_mode: ", data_mode)
}

#local
get_deg_manifest_local <- function(pkg = "FRDATranscriptomicAtlas") {

  ensure_atlas_data(
    keys = c("deg_genes", "deg_transcripts"),
    package = pkg
  )

  cache_root <- tools::R_user_dir(pkg, which = "cache")

  deg_dir_genes <- file.path(cache_root, "genes")
  deg_dir_transcripts <- file.path(cache_root, "txs")

  files <- c(
    if (dir.exists(deg_dir_genes)) {
      list.files(deg_dir_genes, full.names = TRUE)
    } else {
      character(0)
    },
    if (dir.exists(deg_dir_transcripts)) {
      list.files(deg_dir_transcripts, full.names = TRUE)
    } else {
      character(0)
    }
  )

  if (!length(files)) {
    return(tibble::tibble())
  }

  rx <- "^.*/DESEQ2_res_(.+)_(0\\.[0-9]+)_all_(genes|transcripts)\\.rds$"

  tibble::tibble(path = files) |>
    tidyr::extract(
      path,
      into = c("dataset", "p_str", "level"),
      regex = rx,
      remove = FALSE
    ) |>
    dplyr::mutate(p = suppressWarnings(as.numeric(p_str)))
}

get_deg_manifest <- function(pkg = "FRDATranscriptomicAtlas", data_mode = "local") {
  if (identical(data_mode, "local")) {
    return(get_deg_manifest_local(pkg = pkg))
  }

  if (identical(data_mode, "cloud")) {
    return(get_deg_manifest_cloud(pkg = pkg))
  }

  stop("Unknown data_mode: ", data_mode)
}

#cloud
get_deg_manifest_cloud <- function(pkg = "FRDATranscriptomicAtlas") {

  readr::read_csv(
    "https://frda-transcriptomic-atlas-835050295613-ap-southeast-2-an.s3.ap-southeast-2.amazonaws.com/metadata/deg_manifest.csv",
    show_col_types = FALSE
  ) |>
    dplyr::mutate(
      path = NA_character_,
      p_str = dplyr::case_when(
        threshold == "padj_0.10" ~ "0.10",
        threshold == "padj_0.05" ~ "0.05",
        threshold == "padj_0.01" ~ "0.01",
        threshold == "padj_0.001" ~ "0.001",
        threshold == "all" ~ NA_character_,
        TRUE ~ NA_character_
      ),
      p = suppressWarnings(as.numeric(p_str))
    ) |>
    dplyr::select(path, dataset, p_str, level, p, threshold)
}

get_deg_data_cloud <- function(
    dataset_id,
    feature_level,
    padj_max = NULL,
    lfc_min = 0,
    direction = "both"
) {
  threshold <- padj_max_to_threshold(padj_max)

  path <- paste0(
    "s3://frda-transcriptomic-atlas-835050295613-ap-southeast-2-an/deg_results/",
    "level=", feature_level,
    "/dataset=", dataset_id,
    "/threshold=", threshold
  )

  ds <- arrow::open_dataset(path)
  q <- ds

  if (!is.null(padj_max) && !is.na(padj_max) && threshold == "all") {
    q <- q |>
      dplyr::filter(!is.na(.data$padj), .data$padj <= !!padj_max)
  }

  # Safety filter. This is useful if padj_max is not one of the prebuilt thresholds.
  if (!is.null(padj_max) && !is.na(padj_max) && threshold == "all") {
    q <- q |>
      dplyr::filter(!is.na(.data$padj), .data$padj <= !!padj_max)
  }

  if (identical(direction, "up")) {
    q <- q |>
      dplyr::filter(.data$log2FoldChange >= !!lfc_min)
  } else if (identical(direction, "down")) {
    q <- q |>
      dplyr::filter(.data$log2FoldChange <= -!!lfc_min)
  } else {
    q <- q |>
      dplyr::filter(abs(.data$log2FoldChange) >= !!lfc_min)
  }

  x <- q |>
    dplyr::select(-dplyr::any_of(c("dataset", "level", "threshold", "p_str", "p"))) |>
    dplyr::collect()

  if (feature_level == "genes" && "transcript_id" %in% names(x)) {
    x <- dplyr::select(x, -transcript_id)
  }

  if (feature_level == "transcripts" && "ensembl_gene_id" %in% names(x)) {
    x <- dplyr::select(x, -ensembl_gene_id)
  }

  return(x)
}


padj_max_to_threshold <- function(padj_max) {
  if (is.null(padj_max) || is.na(padj_max)) {
    return("all")
  }

  dplyr::case_when(
    padj_max == 0.10 ~ "padj_0.10",
    padj_max == 0.05 ~ "padj_0.05",
    padj_max == 0.01 ~ "padj_0.01",
    padj_max == 0.001 ~ "padj_0.001",
    TRUE ~ "all"
  )
}

get_deg_data_cloud_rds <- function(
    dataset_id,
    feature_level,
    padj_max = NULL,
    lfc_min = 0,
    direction = "both"
) {
  threshold <- padj_max_to_threshold(padj_max)

  url <- paste0(
    "https://frda-transcriptomic-atlas-835050295613-ap-southeast-2-an.s3.ap-southeast-2.amazonaws.com/",
    "deg_results_rds/",
    feature_level, "/",
    dataset_id, "/",
    threshold, ".rds"
  )

  tf <- tempfile(fileext = ".rds")
  on.exit(unlink(tf), add = TRUE)

  utils::download.file(
    url,
    destfile = tf,
    mode = "wb",
    quiet = TRUE
  )

  x <- readRDS(tf)

  if ("log2FoldChange" %in% names(x)) {
    if (identical(direction, "up")) {
      x <- x[x$log2FoldChange >= lfc_min, , drop = FALSE]
    } else if (identical(direction, "down")) {
      x <- x[x$log2FoldChange <= -lfc_min, , drop = FALSE]
    } else {
      x <- x[abs(x$log2FoldChange) >= lfc_min, , drop = FALSE]
    }
  }

  x
}

deg_cache <- cachem::cache_mem(max_size = 500 * 1024^2)

get_deg_data_cloud_rds_cached <- memoise::memoise(
  get_deg_data_cloud_rds,
  cache = deg_cache
)


#biomarker server
biomarker_cache <- cachem::cache_mem(max_size = 200 * 1024^2)

get_biomarker_baseline_cloud <- function() {
  url <- paste0(
    "https://frda-transcriptomic-atlas-835050295613-ap-southeast-2-an.s3.ap-southeast-2.amazonaws.com/",
    "biomarker/baseline_long.rds"
  )

  tf <- tempfile(fileext = ".rds")
  on.exit(unlink(tf), add = TRUE)

  utils::download.file(url, tf, mode = "wb", quiet = TRUE)
  readRDS(tf)
}

get_biomarker_baseline_cloud_cached <- memoise::memoise(
  get_biomarker_baseline_cloud,
  cache = biomarker_cache
)


#PCA -----------------------------------------------------------
pca_cache <- cachem::cache_mem(max_size = 500 * 1024^2)

get_pca_manifest_cloud <- function() {
  readr::read_csv(
    "https://frda-transcriptomic-atlas-835050295613-ap-southeast-2-an.s3.ap-southeast-2.amazonaws.com/metadata/pca_manifest.csv",
    show_col_types = FALSE
  )
}

get_pca_data_cloud <- function(filename) {

  url <- paste0(
    "https://frda-transcriptomic-atlas-835050295613-ap-southeast-2-an.s3.ap-southeast-2.amazonaws.com/",
    "pca_input/",
    filename
  )

  tf <- tempfile(fileext = ".rds")
  on.exit(unlink(tf), add = TRUE)

  utils::download.file(url, tf, mode = "wb", quiet = TRUE)
  readRDS(tf)
}

get_pca_data_cloud_cached <- memoise::memoise(
  get_pca_data_cloud,
  cache = pca_cache
)


# Gene plots
get_tpm_gene_manifest_cloud <- function() {
  readr::read_csv(
    "https://frda-transcriptomic-atlas-835050295613-ap-southeast-2-an.s3.ap-southeast-2.amazonaws.com/metadata/tpm_gene_manifest.csv",
    show_col_types = FALSE
  )
}

get_tpm_gene_cloud <- function(filename) {

  url <- paste0(
    "https://frda-transcriptomic-atlas-835050295613-ap-southeast-2-an.s3.ap-southeast-2.amazonaws.com/",
    "tpm_gene/",
    filename
  )

  tmp <- tempfile(fileext = ".rds")

  download.file(
    url,
    destfile = tmp,
    mode = "wb",
    quiet = TRUE
  )

  readRDS(tmp)
}

get_tpm_gene_cloud_cached <- memoise::memoise(
  get_tpm_gene_cloud,
  cache = cachem::cache_disk(
    file.path(
      tools::R_user_dir("FRDATranscriptomicAtlas", "cache"),
      "cloud_tpm_gene"
    )
  )
)


get_tpm_transcript_manifest_cloud <- function() {
  readr::read_csv(
    "https://frda-transcriptomic-atlas-835050295613-ap-southeast-2-an.s3.ap-southeast-2.amazonaws.com/metadata/tpm_transcript_manifest.csv",
    show_col_types = FALSE
  )
}

get_tpm_transcript_cloud <- function(filename) {
  url <- paste0(
    "https://frda-transcriptomic-atlas-835050295613-ap-southeast-2-an.s3.ap-southeast-2.amazonaws.com/",
    "tpm_transcript/",
    filename
  )

  tmp <- tempfile(fileext = ".rds")

  download.file(url, tmp, mode = "wb", quiet = TRUE)

  readRDS(tmp)
}

get_tpm_transcript_cloud_cached <- memoise::memoise(
  get_tpm_transcript_cloud,
  cache = cachem::cache_disk(
    file.path(
      tools::R_user_dir("FRDATranscriptomicAtlas", "cache"),
      "cloud_tpm_transcript"
    )
  )
)

get_vsd_manifest_cloud <- function() {
  readr::read_csv(
    "https://frda-transcriptomic-atlas-835050295613-ap-southeast-2-an.s3.ap-southeast-2.amazonaws.com/metadata/vsd_manifest.csv",
    show_col_types = FALSE
  )
}

get_vsd_cloud <- function(filename) {
  url <- paste0(
    "https://frda-transcriptomic-atlas-835050295613-ap-southeast-2-an.s3.ap-southeast-2.amazonaws.com/",
    "vsd/",
    filename
  )

  tmp <- tempfile(fileext = ".rds")

  download.file(url, tmp, mode = "wb", quiet = TRUE)

  readRDS(tmp)
}

get_vsd_cloud_cached <- memoise::memoise(
  get_vsd_cloud,
  cache = cachem::cache_disk(
    file.path(
      tools::R_user_dir("FRDATranscriptomicAtlas", "cache"),
      "cloud_vsd"
    )
  )
)


# GSEA -----------------------------------------------------------
# GSEA -----------------------------------------------------------

gsea_cache <- cachem::cache_mem(max_size = 500 * 1024^2)

get_gsea_manifest_cloud <- function() {
  readr::read_csv(
    "https://frda-transcriptomic-atlas-835050295613-ap-southeast-2-an.s3.ap-southeast-2.amazonaws.com/metadata/gsea_manifest.csv",
    show_col_types = FALSE
  )
}

get_gsea_data_cloud <- function(filename) {

  url <- paste0(
    "https://frda-transcriptomic-atlas-835050295613-ap-southeast-2-an.s3.ap-southeast-2.amazonaws.com/",
    "GSEA_results/",
    filename
  )

  tf <- tempfile(fileext = ".rds")
  on.exit(unlink(tf), add = TRUE)

  utils::download.file(
    url,
    destfile = tf,
    mode = "wb",
    quiet = TRUE
  )

  readRDS(tf)
}

get_gsea_data_cloud_cached <- memoise::memoise(
  get_gsea_data_cloud,
  cache = gsea_cache
)
