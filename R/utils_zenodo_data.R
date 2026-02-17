#' Cross-platform robust download using libcurl
#' @keywords internal
#' @noRd
download_with_retry <- function(url, destfile, md5 = NA_character_, retries = 5) {

  options(timeout = max(6000, getOption("timeout")))

  has_libcurl <- isTRUE(capabilities("libcurl"))
  method <- if (has_libcurl) "libcurl" else "auto"

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

  # Require expected columns (prevents mysterious failures)
  req_cols <- c("key", "description", "url", "local_dir", "md5", "version")
  miss_cols <- setdiff(req_cols, names(manifest))
  if (length(miss_cols)) {
    stop("Manifest missing columns: ", paste(miss_cols, collapse = ", "), call. = FALSE)
  }

  rows <- manifest[manifest$key %in% keys, , drop = FALSE]
  if (!nrow(rows)) {
    stop("No matching data keys in manifest: ", paste(keys, collapse = ", "), call. = FALSE)
  }
  rows0 <- rows  # keep original requested rows

  cache_root <- tools::R_user_dir(package, which = "cache")
  dir.create(cache_root, recursive = TRUE, showWarnings = FALSE)

  # Determine which datasets need downloading (version/md5-aware marker)
  needs_download <- function(target_dir, key, md5, version) {

    marker <- file.path(target_dir, paste0(".", key, ".complete"))

    # No directory → must download
    if (!dir.exists(target_dir)) return(TRUE)

    # No marker → must download
    if (!file.exists(marker)) return(TRUE)

    # Read marker
    txt <- tryCatch(readLines(marker, warn = FALSE), error = function(e) character())

    stored_version <- sub("^version=", "", txt[grepl("^version=", txt)])
    stored_md5     <- sub("^md5=", "", txt[grepl("^md5=", txt)])

    # Malformed marker => refresh
    if (!length(stored_version) || !length(stored_md5)) return(TRUE)

    # If version differs → update
    if (!is.na(version) && nzchar(version) && stored_version != version) return(TRUE)

    # If md5 differs → update
    if (!is.na(md5) && nzchar(md5) && tolower(stored_md5) != tolower(md5)) return(TRUE)

    FALSE
  }

  need_download <- vapply(
    seq_len(nrow(rows)),
    function(i) needs_download(
      file.path(cache_root, rows$local_dir[i]),
      rows$key[i],
      rows$md5[i],
      rows$version[i]
    ),
    logical(1)
  )

  rows <- rows[need_download, , drop = FALSE]

  if (!nrow(rows)) {
    return(invisible(TRUE))
  }

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
          md5      = row$md5
        )

        target_dir <- file.path(cache_root, row$local_dir)

        # --- Option 1: full replacement on update/install ---
        if (dir.exists(target_dir)) {
          unlink(target_dir, recursive = TRUE, force = TRUE)
        }
        dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)

        # unzip into target_dir (robust across zip layouts + OS)
        utils::unzip(zip_path, exdir = target_dir)
        unlink(zip_path)

        # --- Robust flattening: handle target_dir/<same_name>/... nesting ---
        inner_same <- file.path(target_dir, basename(target_dir))
        if (dir.exists(inner_same)) {
          inner_items <- list.files(inner_same, full.names = TRUE, all.files = FALSE, no.. = TRUE)
          # (all.files = FALSE avoids .DS_Store / hidden stuff)

          ok <- file.rename(inner_items, file.path(target_dir, basename(inner_items)))
          if (any(!ok)) {
            warning("Some files could not be moved out of nested directory: ", inner_same)
          }


          unlink(inner_same, recursive = TRUE, force = TRUE)
        }

        # Fallback: if there is a single subdir and no top-level data files, flatten it
        subdirs <- list.dirs(target_dir, recursive = FALSE, full.names = TRUE)
        top_files <- list.files(target_dir, recursive = FALSE, all.files = FALSE)

        if (length(subdirs) == 1 && length(top_files) == 0) {
          nested <- subdirs[1]
          moved <- file.path(nested, list.files(nested, full.names = FALSE))
          file.rename(moved, file.path(target_dir, basename(moved)))
          unlink(nested, recursive = TRUE, force = TRUE)
        }

        # Validate that something was extracted
        if (!length(list.files(target_dir, recursive = TRUE, all.files = FALSE))) {
          stop("Unzip completed but no files were extracted into: ", target_dir, call. = FALSE)
        }

        # mark as complete ONLY after successful unzip
        marker <- file.path(target_dir, paste0(".", row$key, ".complete"))
        writeLines(
          c(
            paste0("version=", row$version),
            paste0("md5=", row$md5)
          ),
          marker
        )

        message(sprintf(
          "[Atlas] Finished %s -> %s",
          row$key,
          row$local_dir
        ))
      }
    }
  )

  # ---- REMOVE DOWNLOAD BANNER ------------------------------------------
  if (shiny::isRunning()) {
    shiny::removeNotification("atlas-download")
  }

  # ---- READY NOTIFICATION (currentness-aware) --------------------------
  all_complete <- all(!vapply(
    seq_len(nrow(rows0)),
    function(i) needs_download(
      file.path(cache_root, rows0$local_dir[i]),
      rows0$key[i],
      rows0$md5[i],
      rows0$version[i]
    ),
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


#' Check whether atlas data manifest indicates updates (GitHub vs installed)
#' @keywords internal
#' @noRd
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

  # Keep only useful columns (if present)
  remote_keep <- intersect(c("key", "version", "md5", "description"), names(remote))
  local_keep  <- intersect(c("key", "version", "md5"), names(local))

  remote2 <- remote[, remote_keep, drop = FALSE]
  local2  <- local[,  local_keep,  drop = FALSE]

  m <- merge(local2, remote2, by = "key", all = TRUE, suffixes = c("_local", "_remote"))

  # "new key" = not present locally at all (no version & no md5 locally)
  is_new_key <- is.na(m$version_local) & is.na(m$md5_local)

  version_changed_vec <- if (all(c("version_local", "version_remote") %in% names(m))) {
    !is.na(m$version_remote) & !is.na(m$version_local) & (m$version_remote != m$version_local)
  } else rep(FALSE, nrow(m))

  md5_changed_vec <- if (all(c("md5_local", "md5_remote") %in% names(m))) {
    !is.na(m$md5_remote) & !is.na(m$md5_local) &
      (tolower(m$md5_remote) != tolower(m$md5_local))
  } else rep(FALSE, nrow(m))

  updates <- m[is_new_key | version_changed_vec | md5_changed_vec, , drop = FALSE]
  if (!nrow(updates)) return(NULL)

  updates
}
