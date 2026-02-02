suppressPackageStartupMessages({
  library(shiny)
})

# Source your package-like R/ files
invisible(lapply(list.files("R", pattern = "\\.R$", full.names = TRUE), source, local = TRUE))

# Run
run_app()
