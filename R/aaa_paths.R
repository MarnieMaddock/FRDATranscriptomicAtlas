# R/aaa_paths.R
pkg_or_proj_path <- function(...,
                             package = "FRDATranscriptomicAtlas",
                             must_exist = TRUE) {
  local <- file.path("inst", ...)
  if (file.exists(local)) return(local)

  pkg <- system.file(..., package = package)
  if (nzchar(pkg)) return(pkg)

  if (isTRUE(must_exist)) {
    stop("Could not locate file: ", file.path("inst", ...))
  } else {
    ""
  }
}
