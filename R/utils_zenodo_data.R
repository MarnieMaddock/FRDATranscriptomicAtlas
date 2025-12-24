#' Cross-platform robust download using libcurl
#' @keywords internal
download_with_retry <- function(url, destfile, retries = 3, quiet = TRUE) {

  options(timeout = max(600, getOption("timeout")))

  # Prefer libcurl everywhere
  has_libcurl <- isTRUE(capabilities("libcurl"))


  method <- if (has_libcurl) "libcurl" else "auto"

  for (i in seq_len(retries)) {
    ok <- tryCatch({
      utils::download.file(
        url,
        destfile,
        mode   = "wb",
        method = method,
        quiet  = quiet
      )
      TRUE
    }, error = function(e) FALSE, warning = function(w) FALSE)

    if (isTRUE(ok) && file.exists(destfile) && file.size(destfile) > 0) {
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

  # Determine which datasets actually need downloading
  need_download <- vapply(
    rows$local_dir,
    function(d) {
      p <- file.path(cache_root, d)
      !(dir.exists(p) && length(list.files(p)))
    },
    logical(1)
  )

  rows <- rows[need_download, , drop = FALSE]

  # Nothing to do → no banner, no progress
  if (!nrow(rows)) {
    return(invisible(TRUE))
  }

  downloaded_any <- FALSE

  shiny::withProgress(
    message = "Preparing atlas data",
    value = 0,
    {
      n <- nrow(rows)

      for (i in seq_len(n)) {
        row <- rows[i, ]

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

        downloaded_any <- TRUE
      }
    }
  )

  # ---- GLOBAL DATA READY BANNER -----------------------------------------
  if (downloaded_any) {
    shiny::showNotification(
      "Atlas data downloaded and ready.",
      type     = "message",
      duration = NULL,     # persistent
      id       = "atlas-data-ready"
    )
  }

  invisible(TRUE)
}
