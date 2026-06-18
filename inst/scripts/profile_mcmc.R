#!/usr/bin/env Rscript
# Profile MCMC C++ hotspots using simulated data.
#
# Prerequisite: Build the package with the profiler registry enabled so that
# moire_profiler_stats() returns data. In src/Makevars add:
#   PKG_CXXFLAGS += -DMOIRE_ENABLE_PROFILER_REGISTRY
# Then: devtools::install() or devtools::load_all()
#
# Usage:
#   Rscript inst/scripts/profile_mcmc.R              # default: simulated data, small preset
#   Rscript inst/scripts/profile_mcmc.R minimal      # minimal long-form data (fast)
#   Rscript inst/scripts/profile_mcmc.R small        # simulated: 20 samples, 10 loci (default)
#   Rscript inst/scripts/profile_mcmc.R vignette     # simulated: 100 samples, 100 loci (like vignette)
#
# Env vars (optional): PROFILE_BURNIN, PROFILE_SAMPLES, PROFILE_SEED
#   e.g. PROFILE_SAMPLES=500 Rscript inst/scripts/profile_mcmc.R small

args <- commandArgs(trailingOnly = TRUE)
preset <- if (length(args) >= 1) args[1] else "small"

if (!requireNamespace("moire", quietly = TRUE)) {
  stop("Package 'moire' must be installed/loaded (e.g. devtools::load_all()).")
}

burnin   <- as.integer(Sys.getenv("PROFILE_BURNIN",  "200"))
samples  <- as.integer(Sys.getenv("PROFILE_SAMPLES", "200"))
seed     <- as.integer(Sys.getenv("PROFILE_SEED",   "42"))

set.seed(seed)

if (preset == "minimal") {
  message("Data: minimal (long-form, 4 samples x 3 loci)")
  n_samples <- 4L
  n_loci   <- 3L
  sample_ids <- rep(seq_len(n_samples), each = n_loci * 2)
  loci       <- rep(rep(seq_len(n_loci), each = 2), times = n_samples)
  alleles    <- rep(c(1L, 2L), times = n_samples * n_loci)
  df <- data.frame(sample_id = sample_ids, locus = loci, allele = alleles)
  data <- moire::load_long_form_data(df)
} else {
  if (preset == "vignette") {
    message("Data: simulated (vignette-like), 100 samples x 100 loci")
    num_samples <- 100L
    allele_counts <- c(rep(5L, 50L), rep(10L, 50L))
  } else {
    message("Data: simulated (small), 20 samples x 10 loci")
    num_samples <- 20L
    allele_counts <- c(rep(5L, 5L), rep(10L, 5L))
  }
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
}

message("MCMC: burnin=", burnin, " samples_per_chain=", samples, " (verbose=FALSE)")

tryCatch(moire::moire_prof_reset(), error = function(e) {
  message("moire_prof_reset not found; build with MOIRE_ENABLE_PROFILER_REGISTRY")
})

elapsed <- system.time(
  invisible(moire::run_mcmc(data, burnin = burnin, samples_per_chain = samples, verbose = FALSE))
)

message("\nWall-clock time:")
print(elapsed)

tryCatch({
  stats <- moire::moire_prof_stats()
  if (nrow(stats) > 0) {
    message("\nC++ profiler stats (top by total_ms):")
    print(head(stats, 20))
    message("\nFraction of total C++ time per key (approx):")
    total_ms <- sum(stats$total_ms)
    if (total_ms > 0) {
      frac <- stats
      frac$pct <- round(100 * frac$total_ms / total_ms, 1)
      print(frac[order(-frac$total_ms), c("key", "calls", "total_ms", "avg_ms", "pct")])
    }
  } else {
    message("Profiler returned no data. Rebuild with PKG_CXXFLAGS += -DMOIRE_ENABLE_PROFILER_REGISTRY")
  }
}, error = function(e) {
  message("moire_prof_stats not available; rebuild with MOIRE_ENABLE_PROFILER_REGISTRY")
})
