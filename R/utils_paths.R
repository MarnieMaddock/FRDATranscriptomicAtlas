# put this in a small utils file (e.g., R/utils_paths.R) and source it
pkg_or_proj_path <- function(..., package = "FRDATranscriptomicAtlas") {
  local <- file.path("inst", ...)
  if (file.exists(local)) return(local)
  pkg <- system.file(..., package = package)
  if (nzchar(pkg)) return(pkg)
  stop("Could not locate file: ", file.path("inst", ...))
}
