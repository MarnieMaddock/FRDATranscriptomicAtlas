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

    return(x)
  }

  if (identical(data_mode, "cloud")) {
    return(get_deg_data_cloud(
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
  ds <- arrow::open_dataset("s3://frda-transcriptomic-atlas-835050295613-ap-southeast-2-an/deg_results")

  ds |>
    dplyr::distinct(dataset, p_str, level, p) |>
    dplyr::collect() |>
    dplyr::mutate(path = NA_character_) |>
    dplyr::select(path, dataset, p_str, level, p)
}



get_deg_data_cloud <- function(
    dataset_id,
    feature_level,
    padj_max = NULL,
    lfc_min = 0,
    direction = "both"
) {
  ds <- arrow::open_dataset(
    "s3://frda-transcriptomic-atlas-835050295613-ap-southeast-2-an/deg_results",
    partitioning = c("level", "dataset")
  )

  q <- ds |>
    dplyr::filter(
      .data$dataset == !!dataset_id,
      .data$level == !!feature_level
    )

  if (!is.null(padj_max) && !is.na(padj_max)) {
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
    dplyr::select(-dataset, -level, -p_str, -p) |>
    dplyr::collect()

  if (feature_level == "genes" && "transcript_id" %in% names(x)) {
    x <- dplyr::select(x, -transcript_id)
  }

  if (feature_level == "transcripts" && "ensembl_gene_id" %in% names(x)) {
    x <- dplyr::select(x, -ensembl_gene_id)
  }

  return(x)
}
