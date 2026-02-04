suppressPackageStartupMessages({
  library(shiny)
})


# Source R/ files
invisible(lapply(list.files("R", pattern = "\\.R$", full.names = TRUE), source, local = TRUE))

# Run
run_app()
