#' Cross-platform robust download using libcurl
#' @keywords internal
download_with_retry <- function(url, destfile, retries = 3) {

  options(timeout = max(6000, getOption("timeout")))

  has_libcurl <- isTRUE(capabilities("libcurl"))
  method <- if (has_libcurl) "libcurl" else "auto"

  for (i in seq_len(retries)) {

    message(sprintf("[Atlas] Attempt %d: %s", i, basename(url)))

    ok <- tryCatch({
      utils::download.file(
        url,
        destfile,
        mode   = "wb",
        method = method,
        quiet  = FALSE   # console progress
      )
      TRUE
    }, error = function(e) FALSE)

    if (ok && file.exists(destfile) && file.size(destfile) > 0) {
      return(invisible(TRUE))
    }

    Sys.sleep(2)
  }

  stop(
    "Failed to download file after ", retries, " attempts:\n", url,
    call. = FALSE
  )
}



#' Ensure FRDA atlas data are available locally (with progress bar + global banner)
#' @keywords internal
ensure_atlas_data <- function(
    keys,
    package = "FRDATranscriptomicAtlas"
) {

  manifest_path <- system.file(
    "extdata", "atlas_data_manifest.csv",
    package = package
  )

  if (!nzchar(manifest_path)) {
    stop("atlas_data_manifest.csv not found in installed package.")
  }

  manifest <- read.csv(manifest_path, stringsAsFactors = FALSE)

  rows <- manifest[manifest$key %in% keys, , drop = FALSE]
  if (!nrow(rows)) {
    stop("No matching data keys in manifest: ", paste(keys, collapse = ", "))
  }

  cache_root <- tools::R_user_dir(package, which = "cache")
  dir.create(cache_root, recursive = TRUE, showWarnings = FALSE)

  # Determine which datasets need downloading
  needs_download <- function(p, pattern = "\\.rds$") {
    if (!dir.exists(p)) return(TRUE)
    files <- list.files(p, pattern = pattern, recursive = TRUE)
    length(files) == 0
  }
  need_download <- vapply(
    rows$local_dir,
    function(d) {
      p <- file.path(cache_root, d)
      needs_download(p)
    },
    logical(1)
  )


  rows <- rows[need_download, , drop = FALSE]

  if (!nrow(rows)) {
    return(invisible(TRUE))
  }

  downloaded_any <- FALSE

  # ---- SHOW DOWNLOAD BANNER --------------------------------------------
  if (shiny::isRunning()) {
    shiny::showNotification(
      "Downloading atlas data… this may take several minutes.",
      type     = "warning",
      duration = NULL,
      id       = "atlas-download"
    )
  }

  shiny::withProgress(
    message = "Preparing atlas data",
    value   = 0,
    {

      n <- nrow(rows)

      for (i in seq_len(n)) {

        row <- rows[i, ]

        message(sprintf(
          "[Atlas] Downloading %s (%s)",
          row$key,
          row$description
        ))

        shiny::incProgress(
          1 / n,
          detail = paste("Downloading:", row$description)
        )

        zip_path <- file.path(cache_root, basename(row$url))
        download_with_retry(row$url, zip_path)

        utils::unzip(zip_path, exdir = cache_root)
        unlink(zip_path)

        target_dir <- file.path(cache_root, row$local_dir)
        if (!dir.exists(target_dir)) {
          stop("Expected data directory not created: ", target_dir)
        }

        message(sprintf(
          "[Atlas] Finished %s → %s",
          row$key,
          row$local_dir
        ))

        downloaded_any <- TRUE
      }
    }
  )

  # ---- REMOVE DOWNLOAD BANNER ------------------------------------------
  if (shiny::isRunning()) {
    shiny::removeNotification("atlas-download")
  }

  # ---- READY NOTIFICATION ----------------------------------------------
  if (downloaded_any && shiny::isRunning()) {
    shiny::showNotification(
      "Atlas data downloaded and ready.",
      type     = "message",
      duration = 5
    )
  }

  invisible(TRUE)
}
