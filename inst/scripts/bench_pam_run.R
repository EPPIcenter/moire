#!/usr/bin/env Rscript
# Single MCMC run for P(any missing) benchmark. Called by bench_pam_driver.R
# with MOIRE_USE_LEGACY_PAM unset (Gray) or set (legacy).
# Usage: Rscript bench_pam_run.R <data.rds> <gray|legacy> <out_dir>
# Writes: <out_dir>/time_<gray|legacy>.txt with elapsed seconds.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("Usage: Rscript bench_pam_run.R <data.rds> <gray|legacy> <out_dir>")
data_path <- args[1]
path_type <- match.arg(args[2], c("gray", "legacy"))
out_dir <- args[3]

if (!file.exists(data_path)) stop("Data file not found: ", data_path)
if (!dir.exists(out_dir)) stop("Output dir not found: ", out_dir)

# Assume run from package root (e.g. by bench_pam_driver.R)
pkg_root <- Sys.getenv("BENCH_PAM_PKG_ROOT", getwd())
devtools::load_all(pkg_root, quiet = TRUE)
data <- readRDS(data_path)

burnin   <- as.integer(Sys.getenv("BENCH_PAM_BURNIN",  "500"))
samples  <- as.integer(Sys.getenv("BENCH_PAM_SAMPLES", "500"))

elapsed <- system.time(
  invisible(moire::run_mcmc(data, burnin = burnin, samples_per_chain = samples, verbose = FALSE))
)[["elapsed"]]

out_file <- file.path(out_dir, paste0("time_", path_type, ".txt"))
writeLines(as.character(elapsed), out_file)
