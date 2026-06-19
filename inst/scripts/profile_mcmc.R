#!/usr/bin/env Rscript
# Profile MCMC C++ hotspots using simulated data.
#
# Prerequisite: Build the package with the profiler registry enabled so that
# moire_profiler_stats() returns data. In src/Makevars add:
#   PKG_CXXFLAGS += -DMOIRE_ENABLE_PROFILER_REGISTRY
# Then: devtools::install() or devtools::load_all()
#
# Usage:
#   Rscript inst/scripts/profile_mcmc.R              # default: vignette preset
#   Rscript inst/scripts/profile_mcmc.R minimal      # minimal long-form data (fast)
#   Rscript inst/scripts/profile_mcmc.R small        # simulated: 20 samples, 10 loci (default)
#   Rscript inst/scripts/profile_mcmc.R vignette     # simulated: 100 samples, 100 loci (like vignette)
#
# For repeatable baselines (wall-clock + profiler CSV), use:
#   Rscript inst/scripts/bench_mcmc_baseline.R small --save
#   Rscript inst/scripts/bench_mcmc_baseline.R vignette --pt --save
#   Rscript inst/scripts/bench_mcmc_baseline.R small --compare
# See inst/benchmarks/README.md for parallel modes and env vars.
#
# Env vars (optional): PROFILE_BURNIN, PROFILE_SAMPLES, PROFILE_SEED
#   e.g. PROFILE_SAMPLES=500 Rscript inst/scripts/profile_mcmc.R small

script_dir <- dirname(normalizePath(sub("--file=(.*)", "\\1", commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))][1], perl = TRUE)))
source(file.path(script_dir, "bench_helpers.R"))

args <- commandArgs(trailingOnly = TRUE)
preset <- if (length(args) >= 1) args[1] else bench_default_preset()

if (!requireNamespace("moire", quietly = TRUE)) {
  stop("Package 'moire' must be installed/loaded (e.g. devtools::load_all()).")
}

burnin   <- as.integer(Sys.getenv("PROFILE_BURNIN",  "200"))
samples  <- as.integer(Sys.getenv("PROFILE_SAMPLES", "200"))
seed     <- as.integer(Sys.getenv("PROFILE_SEED",   "42"))

set.seed(seed)

loaded <- bench_load_data(preset)
data <- loaded$data
parallel <- bench_parallel_config()

message("MCMC: burnin=", burnin, " samples_per_chain=", samples, " (verbose=FALSE)")
message("Parallel: ", bench_format_parallel_config(parallel))

result <- bench_run_once(data, burnin, samples, verbose = FALSE, parallel = parallel)

message("\nWall-clock time:")
print(c(elapsed = result$elapsed_s))

stats <- result$stats
if (nrow(stats) > 0) {
  message("\nC++ profiler stats (top by total_ms):")
  print(head(stats, 20))
  message("\nFraction of total C++ time per key (approx):")
  print(bench_profiler_fractions(stats)[, c("key", "calls", "total_ms", "avg_ms", "pct")])
} else {
  message("Profiler returned no data. Rebuild with PKG_CXXFLAGS += -DMOIRE_ENABLE_PROFILER_REGISTRY")
}
