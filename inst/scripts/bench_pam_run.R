#!/usr/bin/env Rscript
# Single MCMC run for P(any missing) benchmark. Called by bench_pam_driver.R.
# Usage: Rscript bench_pam_run.R <data.rds> <out_dir>
# Writes: <out_dir>/time.txt with elapsed seconds.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: Rscript bench_pam_run.R <data.rds> <out_dir>")
data_path <- args[1]
out_dir <- args[2]

if (!file.exists(data_path)) stop("Data file not found: ", data_path)
if (!dir.exists(out_dir)) stop("Output dir not found: ", out_dir)

pkg_root <- Sys.getenv("BENCH_PAM_PKG_ROOT", getwd())
devtools::load_all(pkg_root, quiet = TRUE)
data <- readRDS(data_path)

burnin  <- as.integer(Sys.getenv("BENCH_PAM_BURNIN",  "500"))
samples <- as.integer(Sys.getenv("BENCH_PAM_SAMPLES", "500"))

elapsed <- system.time(
  invisible(moire::run_mcmc(data, burnin = burnin, samples_per_chain = samples, verbose = FALSE))
)[["elapsed"]]

writeLines(as.character(elapsed), file.path(out_dir, "time.txt"))
