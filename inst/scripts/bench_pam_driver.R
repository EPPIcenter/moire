#!/usr/bin/env Rscript
# Compare full-MCMC wall-clock time: Gray (default) vs legacy P(any missing).
# Runs the same MCMC twice in separate R processes (env must be set before R starts).
#
# Usage (from package root):
#   Rscript inst/scripts/bench_pam_driver.R [preset]
# Preset: small (default), medium (40 x 50), vignette. Same as profile_mcmc.R for small/vignette.
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
message("Running Gray (default) ...")
system2(
  R.home("bin/Rscript"),
  c(run_script, data_path, "gray", out_dir),
  env = c(
    paste0("BENCH_PAM_PKG_ROOT=", pkg_root),
    paste0("BENCH_PAM_BURNIN=", burnin),
    paste0("BENCH_PAM_SAMPLES=", samples)
  ),
  stdout = NULL,
  stderr = NULL
)

message("Running legacy (MOIRE_USE_LEGACY_PAM=1) ...")
system2(
  R.home("bin/Rscript"),
  c(run_script, data_path, "legacy", out_dir),
  env = c(
    "MOIRE_USE_LEGACY_PAM=1",
    paste0("BENCH_PAM_PKG_ROOT=", pkg_root),
    paste0("BENCH_PAM_BURNIN=", burnin),
    paste0("BENCH_PAM_SAMPLES=", samples)
  ),
  stdout = NULL,
  stderr = NULL
)

t_gray   <- as.numeric(readLines(file.path(out_dir, "time_gray.txt")))
t_legacy <- as.numeric(readLines(file.path(out_dir, "time_legacy.txt")))

message("")
message("=== P(any missing) full-MCMC comparison ===")
message("Gray (default):  ", round(t_gray, 2), " s")
message("Legacy (combo):  ", round(t_legacy, 2), " s")
message("Ratio (legacy/Gray): ", round(t_legacy / t_gray, 2))
if (t_legacy > t_gray) {
  message("Gray is ", round(100 * (t_legacy - t_gray) / t_legacy, 1), "% faster for this run.")
} else {
  message("Legacy is ", round(100 * (t_gray - t_legacy) / t_gray, 1), "% faster for this run.")
}
message("")
message("Note: P(any missing) is only part of MCMC (transmission process). Observation process,")
message("update_p, update_samples, etc. dominate, so a 2x per-call speedup yields ~10% overall.")
message("Build with MOIRE_ENABLE_PROFILER_REGISTRY to see the exact fraction in profiler stats.")
