# 00_run_all.R
# -----------------------------------------------------------
# Reproducibility orchestrator. Runs the full numbered pipeline
# in order, each step as an isolated Rscript process (the steps
# each setwd("..") internally, so they must not share a session).
#
# Run from the scripts/ directory:
#     Rscript 00_run_all.R
#
# A fresh checkout with the raw inputs in ../../data should
# rebuild every intermediate CSV in ../data and every figure in
# ../figures. This is the acceptance test for reproducibility.
# -----------------------------------------------------------

# Resolve this script's own directory so steps run from scripts/.
args <- commandArgs(FALSE)
fa <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(fa)) dirname(normalizePath(sub("^--file=", "", fa[1]))) else getwd()
setwd(script_dir)

# Step 01a ingests the raw CONAPESCA government files (~5 GB of xlsx and
# csv, ~20 min). Its output is cached, so it is skipped unless missing.
gz <- "../../data/conapesca_avisos_pacifico_2000_2026.csv.gz"
if (!file.exists(gz)) {
  message("\n========== 01a_ingest_conapesca_raw.py ==========")
  if (system2("python3", "01a_ingest_conapesca_raw.py") != 0)
    stop("Raw CONAPESCA ingest FAILED")
} else {
  message("Raw CONAPESCA ingest already cached at ", gz)
}

steps <- c(
  "01_data_preparation.R",     # artisanal landings + SST/MHW copy
  "02_climate_analysis.R",     # SST trend + warm-season anomaly
  "03_ltem_analysis.R",        # reef warming betas + Fig-2 trajectories
  "03b_productivity.R",        # production + turnover on the same panel (Fig 2c,d)
  "03c_pathways.R",            # energy pathways + subsidy interaction (Fig 3)
  "03d_buffer.R",              # buffer over time + regime rates
  "03e_buffer_climate.R",      # buffer strength Phi and the climate test (Fig 3d)
  "03f_buffer_nonlinear.R",    # smooth alternative + model-free arbiter
  "04_fishery_analysis.R",     # CPUE two-mode + lag scans
  "05_variation_partitioning.R",
  "06_economic_value.R",       # 5-state + La Paz/Loreto ex-vessel value
  "07_decoupling.R",           # survey vs landings: hyperstability test
  "08_gap_analysis.R",         # climate vs fishing vs combined; Figure S14
  "09_figures.R"               # all figures + in_text_statistics.csv (runs last)
)

t0 <- Sys.time()
for (s in steps) {
  message("\n========== ", s, " ==========")
  code <- system2("Rscript",
                  c("-e", shQuote(sprintf("setwd('%s'); source('%s')", script_dir, s))))
  if (code != 0) stop("Pipeline FAILED at ", s, " (exit ", code, ")")
}
message(sprintf("\nPipeline complete in %.1f min. Intermediates in ../data, figures in ../figures.",
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
