#!/usr/bin/env Rscript
# Wall-clock MCMC benchmark for P(any missing) path (Gray-code inclusion-exclusion).
#
# Usage (from package root):
#   Rscript inst/scripts/bench_pam_driver.R [preset]
# Preset: small (default), medium (40 x 50), vignette. Same as profile_mcmc.R for small/vignette.
#
# For full profiler breakdown and committed baselines, prefer:
#   Rscript inst/scripts/bench_mcmc_baseline.R [preset] [--save|--compare]
#
# Optional env: BENCH_PAM_BURNIN, BENCH_PAM_SAMPLES (passed to child runs).

args <- commandArgs(trailingOnly = TRUE)
preset <- if (length(args) >= 1) args[1] else "small"

devtools::load_all(".", quiet = TRUE)

if (preset == "vignette") {
  message("Data: 100 samples x 100 loci (vignette-like)")
  num_samples <- 100L
  allele_counts <- c(rep(5L, 50L), rep(10L, 50L))
} else if (preset == "medium") {
  message("Data: 40 samples x 50 loci (medium)")
  num_samples <- 40L
  allele_counts <- c(rep(5L, 25L), rep(10L, 25L))
} else {
  message("Data: 20 samples x 10 loci (small)")
  num_samples <- 20L
  allele_counts <- c(rep(5L, 5L), rep(10L, 5L))
}

set.seed(42)
locus_freq_alphas <- lapply(allele_counts, function(a) rep(1, a))
data <- moire::simulate_data(
  mean_coi = 3,
  num_samples = num_samples,
  epsilon_pos = 0.01,
  epsilon_neg = 0.1,
  locus_freq_alphas = locus_freq_alphas,
  internal_relatedness_alpha = 0.1,
  internal_relatedness_beta = 1
)

out_dir <- tempdir()
data_path <- file.path(out_dir, "bench_pam_data.rds")
saveRDS(data, data_path)

script_dir <- "inst/scripts"
run_script <- file.path(script_dir, "bench_pam_run.R")
if (!file.exists(run_script)) run_script <- file.path(getwd(), script_dir, "bench_pam_run.R")
if (!file.exists(run_script)) stop("Cannot find bench_pam_run.R")

burnin  <- Sys.getenv("BENCH_PAM_BURNIN",  "500")
samples <- Sys.getenv("BENCH_PAM_SAMPLES", "500")
pkg_root <- getwd()

message("MCMC: burnin=", burnin, " samples_per_chain=", samples)
message("Running MCMC ...")
system2(
  R.home("bin/Rscript"),
  c(run_script, data_path, out_dir),
  env = c(
    paste0("BENCH_PAM_PKG_ROOT=", pkg_root),
    paste0("BENCH_PAM_BURNIN=", burnin),
    paste0("BENCH_PAM_SAMPLES=", samples)
  ),
  stdout = NULL,
  stderr = NULL
)

elapsed <- as.numeric(readLines(file.path(out_dir, "time.txt")))
message("")
message("=== P(any missing) full-MCMC benchmark ===")
message("Elapsed: ", round(elapsed, 2), " s")
message("")
message("Note: P(any missing) is only part of MCMC (transmission process). Observation process,")
message("update_p, update_samples, etc. dominate overall wall-clock.")
