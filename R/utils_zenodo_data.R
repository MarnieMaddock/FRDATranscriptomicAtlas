#' Cross-platform robust download using libcurl
#' @keywords internal
#' @noRd
download_with_retry <- function(url, destfile, md5 = NA_character_, retries = 5) {

  options(timeout = max(6000, getOption("timeout")))

  has_libcurl <- isTRUE(capabilities("libcurl"))
  method <- if (has_libcurl) "libcurl" else "auto"

  if (!requireNamespace("openssl", quietly = TRUE)) {
    stop("Package 'openssl' is required for checksum verification.", call. = FALSE)
  }

  md5_ok <- function(path, expected) {
    if (!is.character(expected) || !nzchar(expected) || is.na(expected)) return(TRUE)
    got <- tools::md5sum(path)
    identical(tolower(unname(got)), tolower(expected))
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
    if (!md5_ok(destfile, md5)) {
      message("[Atlas] MD5 mismatch (download likely incomplete). Retrying...")
      unlink(destfile)
      Sys.sleep(2)
      next
    }

    return(invisible(TRUE))
  }

  stop("Failed to download file after ", retries, " attempts:\n", url, call. = FALSE)
}


#Ensure FRDA atlas data are available locally (with progress bar + global banner)
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
  rows0 <- rows  # <--- keep original requested rows
  cache_root <- tools::R_user_dir(package, which = "cache")
  dir.create(cache_root, recursive = TRUE, showWarnings = FALSE)

  # Determine which datasets need downloading
  needs_download <- function(target_dir, key) {
    if (!dir.exists(target_dir)) return(TRUE)
    !file.exists(file.path(target_dir, paste0(".", key, ".complete")))
  }

  need_download <- vapply(
    seq_len(nrow(rows)),
    function(i) needs_download(file.path(cache_root, rows$local_dir[i]), rows$key[i]),
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
          md5   = if ("md5" %in% names(row)) row$md5 else NA_character_
        )


        utils::unzip(zip_path, exdir = cache_root)
        unlink(zip_path)

        target_dir <- file.path(cache_root, row$local_dir)
        if (!dir.exists(target_dir)) {
          stop("Expected data directory not created: ", target_dir)
        }

        # mark as complete ONLY after successful unzip + directory exists
        file.create(file.path(target_dir, paste0(".", row$key, ".complete")))

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
  all_complete <- all(vapply(
    seq_len(nrow(rows0)),
    function(i) {
      td <- file.path(cache_root, rows0$local_dir[i])
      file.exists(file.path(td, paste0(".", rows0$key[i], ".complete")))
    },
    logical(1)
  ))

  if (all_complete && shiny::isRunning()) {
    shiny::showNotification(
      "Atlas data downloaded and ready.",
      type     = "message",
      duration = 5
    )
  }

  invisible(TRUE)
}


check_atlas_updates <- function(package = "FRDATranscriptomicAtlas") {

  local_path <- system.file("extdata", "atlas_data_manifest.csv", package = package)
  if (!nzchar(local_path)) return(NULL)

  local <- utils::read.csv(local_path, stringsAsFactors = FALSE)

  remote_url <- "https://raw.githubusercontent.com/MarnieMaddock/FRDATranscriptomicAtlas/main/inst/extdata/atlas_data_manifest.csv"

  remote <- tryCatch(
    utils::read.csv(remote_url, stringsAsFactors = FALSE),
    error = function(e) NULL
  )
  if (is.null(remote)) return(NULL)

  # Require key column
  if (!("key" %in% names(local)) || !("key" %in% names(remote))) return(NULL)

  # Use version+md5 if present; otherwise fall back to key-only new dataset detection
  keep_cols <- intersect(c("key", "version", "md5", "description"), names(remote))
  remote2 <- remote[, keep_cols, drop = FALSE]

  keep_cols_local <- intersect(c("key", "version", "md5"), names(local))
  local2 <- local[, keep_cols_local, drop = FALSE]

  m <- merge(local2, remote2, by = "key", all = TRUE, suffixes = c("_local", "_remote"))

  is_new_key <- is.na(m$version_local) & is.na(m$md5_local)  # key missing locally

  version_changed <- ("version_local" %in% names(m) && "version_remote" %in% names(m)) &&
    !is.na(m$version_remote) && !is.na(m$version_local) && (m$version_remote != m$version_local)

  sha_changed <- ("md5_local" %in% names(m) && "md5_remote" %in% names(m)) &&
    !is.na(m$md5_remote) && !is.na(m$md5_local) && (tolower(m$md5_remote) != tolower(m$md5_local))

  # vectorise properly
  version_changed_vec <- if ("version_local" %in% names(m) && "version_remote" %in% names(m)) {
    !is.na(m$version_remote) & !is.na(m$version_local) & (m$version_remote != m$version_local)
  } else rep(FALSE, nrow(m))

  sha_changed_vec <- if ("md5_local" %in% names(m) && "md5_remote" %in% names(m)) {
    !is.na(m$md5_remote) & !is.na(m$md5_local) & (tolower(m$md5_remote) != tolower(m$md5_local))
  } else rep(FALSE, nrow(m))

  updates <- m[is_new_key | version_changed_vec | sha_changed_vec, , drop = FALSE]
  if (!nrow(updates)) return(NULL)

  updates
}
