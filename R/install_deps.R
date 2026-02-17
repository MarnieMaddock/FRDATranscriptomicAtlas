#' Install optional dependencies for FRDATranscriptomicAtlas
#'
#' Installs Bioconductor packages needed for heatmap / SummarizedExperiment
#' functionality (e.g., ComplexHeatmap, SummarizedExperiment).
#'
#' @param update Logical; whether to allow updating already-installed packages.
#'   Default FALSE for reproducibility / minimal disruption.
#' @param ask Logical; passed to BiocManager::install(). Default FALSE.
#' @return Invisibly returns TRUE when complete.
#' @export
install_deps <- function(update = FALSE, ask = FALSE, extras = FALSE) {

  # CRAN helper
  cran_install_if_missing <- function(pkgs, quiet_fail = FALSE) {
    pkgs <- unique(pkgs)
    missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
    if (length(missing)) {
      ok <- tryCatch(
        { utils::install.packages(missing, quiet = TRUE); TRUE },
        error = function(e) FALSE
      )
      if (!ok && !quiet_fail) {
        stop("Failed to install CRAN packages: ", paste(missing, collapse = ", "))
      }
    }
    invisible(TRUE)
  }

  # Ensure BiocManager exists (CRAN)
  cran_install_if_missing("BiocManager")

  # Bioconductor deps you need for the app features
  bioc_pkgs <- c("SummarizedExperiment", "ComplexHeatmap", "DESeq2", "DOSE", "enrichplot")

  missing_bioc <- bioc_pkgs[!vapply(bioc_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_bioc)) {
    BiocManager::install(missing_bioc, ask = ask, update = update)
  }

    # best-effort: magick often fails due to system ImageMagick requirements
    cran_install_if_missing("magick", quiet_fail = TRUE)

    if (!requireNamespace("magick", quietly = TRUE)) {
      message(
        "Optional dependency 'magick' could not be installed.\n",
        "This is usually due to missing system requirements (ImageMagick).\n",
        "FRDATranscriptomicAtlas will still run; only magick-dependent features (if any) are unavailable."
      )
    }

  invisible(TRUE)
}
