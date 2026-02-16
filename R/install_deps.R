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
install_deps <- function(update = FALSE, ask = FALSE) {

  # CRAN helper
  cran_install_if_missing <- function(pkgs) {
    pkgs <- unique(pkgs)
    missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
    if (length(missing)) {
      utils::install.packages(missing)
    }
    invisible(TRUE)
  }

  # Ensure BiocManager exists (CRAN)
  cran_install_if_missing("BiocManager")

  # Bioconductor deps you need for the app features
  bioc_pkgs <- c("SummarizedExperiment", "ComplexHeatmap", "DESeq2")

  missing_bioc <- bioc_pkgs[!vapply(bioc_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_bioc)) {
    BiocManager::install(missing_bioc, ask = ask, update = update)
  }

  invisible(TRUE)
}
