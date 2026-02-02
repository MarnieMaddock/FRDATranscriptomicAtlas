#' Cross-platform robust download using libcurl
#' @keywords internal
#' @noRd
#' Cross-platform robust download using libcurl + sha256 integrity check
#' @keywords internal
#' @noRd
download_with_retry <- function(url, destfile, sha256 = NA_character_, retries = 3) {

  options(timeout = max(6000, getOption("timeout")))

  has_libcurl <- isTRUE(capabilities("libcurl"))
  method <- if (has_libcurl) "libcurl" else "auto"

  if (!requireNamespace("openssl", quietly = TRUE)) {
    stop("Package 'openssl' is required for checksum verification.", call. = FALSE)
  }

  sha_ok <- function(path, expected) {
    if (!is.character(expected) || !nzchar(expected) || is.na(expected)) return(TRUE)
    got <- as.character(openssl::sha256(file(path)))
    identical(tolower(got), tolower(expected))
  }

  for (i in seq_len(retries)) {

    message(sprintf("[Atlas] Attempt %d: %s", i, basename(url)))

    ok <- tryCatch({
      utils::download.file(
        url,
        destfile,
        mode   = "wb",
        method = method,
        quiet  = FALSE
      )
      TRUE
    }, error = function(e) FALSE)

    # basic file existence check
    if (!ok || !file.exists(destfile) || file.size(destfile) == 0) {
      Sys.sleep(2)
      next
    }

    # checksum check (critical)
    if (!sha_ok(destfile, sha256)) {
      message("[Atlas] SHA256 mismatch (download likely incomplete). Retrying...")
      unlink(destfile)
      Sys.sleep(2)
      next
    }

    return(invisible(TRUE))
  }

  stop("Failed to download file after ", retries, " attempts:\n", url, call. = FALSE)
}



#' Ensure FRDA atlas data are available locally (with progress bar + global banner)
#' @keywords internal
#' @noRd
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

  manifest <- utils::read.csv(manifest_path, stringsAsFactors = FALSE)

  rows <- manifest[manifest$key %in% keys, , drop = FALSE]
  if (!nrow(rows)) {
    stop("No matching data keys in manifest: ", paste(keys, collapse = ", "))
  }

  cache_root <- tools::R_user_dir(package, which = "cache")
  dir.create(cache_root, recursive = TRUE, showWarnings = FALSE)

  # Determine which datasets need downloading
  needs_download <- function(target_dir) {
    if (!dir.exists(target_dir)) return(TRUE)
    !file.exists(file.path(target_dir, ".complete"))
  }
  need_download <- vapply(
    rows$local_dir,
    function(d) needs_download(file.path(cache_root, d)),
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
      "Downloading atlas data... this may take several minutes.",
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
        # download + checksum validate
        download_with_retry(
          url      = row$url,
          destfile = zip_path,
          sha256   = if ("sha256" %in% names(row)) row$sha256 else NA_character_
        )


        utils::unzip(zip_path, exdir = cache_root)
        unlink(zip_path)

        target_dir <- file.path(cache_root, row$local_dir)
        if (!dir.exists(target_dir)) {
          stop("Expected data directory not created: ", target_dir)
        }

        # mark as complete ONLY after successful unzip + directory exists
        file.create(file.path(target_dir, ".complete"))

        message(sprintf(
          "[Atlas] Finished %s -> %s",
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
