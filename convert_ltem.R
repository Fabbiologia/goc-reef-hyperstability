# convert_ltem.R
# Convert the uploaded LTEM RDS to parquet so the analysis
# can read it instantly.
#
# Run once:
#   Rscript convert_ltem.R
#
# Output: data/ltem.parquet

if (!requireNamespace("arrow", quietly = TRUE)) {
  install.packages("arrow",
                   repos = "https://cloud.r-project.org")
}
suppressPackageStartupMessages({
  library(arrow)
})

# Stable, repo-relative location for the raw LTEM RDS (preferred).
# Place the raw file at data/raw/ltem_historic_updated_2025-12-08.RDS.
candidates <- c(
  "data/raw/ltem_historic_updated_2025-12-08.RDS",
  "data/ltem_historic_updated_2025-12-08.RDS"
)
src <- candidates[file.exists(candidates)][1]
if (is.na(src)) {
  stop("Cannot find the raw LTEM RDS. Place it at data/raw/ltem_historic_updated_2025-12-08.RDS ",
       "(the derived data/ltem.parquet is the analysis input; this script only rebuilds it).")
}
message("Reading: ", src)

x <- readRDS(src)
message("Class: ", paste(class(x), collapse = "/"))
if (is.data.frame(x)) {
  message("Dim:  ", paste(dim(x), collapse = " x "))
  message("Cols: ", paste(names(x), collapse = ", "))
  cat("\nHead:\n"); print(head(x, 3))
  dir.create("data", showWarnings = FALSE)
  write_parquet(x, "data/ltem.parquet")
  message("\nWrote data/ltem.parquet")
} else {
  message("Top-level object is not a data.frame; structure:")
  str(x, max.level = 2)
}
