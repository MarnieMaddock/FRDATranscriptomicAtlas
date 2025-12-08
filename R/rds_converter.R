library(readr)
library(dplyr)


convert_csv_to_rds <- function(in_dir) {
  csv_files <- list.files(in_dir, pattern = "\\.csv$", full.names = TRUE)

  for (f in csv_files) {
    message("Processing: ", f)
    dat <- readr::read_csv(f, show_col_types = FALSE)


    # RDS path
    rds_path <- sub("\\.csv$", ".rds", f)

    # Compressed RDS (xz is usually the best compression)
    saveRDS(dat, rds_path, compress = "xz")
  }
}

ind_dir <- "inst/extdata/DTU"
convert_csv_to_rds(ind_dir)

#remove csv files
csv_files <- list.files(ind_dir, pattern = "\\.csv$", full.names = TRUE)
file.remove(csv_files)

# recompress_rds <- function(in_dir) {
#   rds_files <- list.files(in_dir, pattern = "\\.rds$", full.names = TRUE)
#
#   for (f in rds_files) {
#     message("Recompressing: ", f)
#     obj <- readRDS(f)
#     saveRDS(obj, f, compress = "xz")  # overwrites with stronger compression
#   }
# }
# recompress_rds(ind_dir)
