#!/usr/bin/env Rscript
# End-to-end MCMC baseline benchmark: wall-clock + C++ profiler registry.
#
# Prerequisite: package built with MOIRE_ENABLE_PROFILER_REGISTRY (default in src/Makevars).
#   devtools::load_all() or devtools::install()
#
# Usage:
#   Rscript inst/scripts/bench_mcmc_baseline.R [preset] [--save|--compare]
#   Rscript inst/scripts/bench_mcmc_baseline.R vignette --pt [--save|--compare]
#   Rscript inst/scripts/bench_mcmc_baseline.R small --serial
#
# Presets: minimal, small, medium, vignette (default via BENCH_PRESET), all
#
# Parallel modes (see inst/benchmarks/README.md):
#   default / --single   one chain, inner TBB parallelism (num_threads)
#   --pt                 parallel tempering (BENCH_PT_CHAINS, default 20)
#   --serial             one chain, num_threads = 1 (algorithmic baseline)
#
# Env vars:
#   BENCH_PARALLEL_MODE   single | pt | serial (overridden by CLI flags)
#   BENCH_NUM_THREADS     pin TBB thread budget (default: physical cores - 1)
#   BENCH_PT_CHAINS       PT replica count when mode=pt (default: 20)
#   BENCH_BURNIN          burn-in iterations (default: 200)
#   BENCH_SAMPLES         post-burnin samples per chain (default: 200)
#   BENCH_SEED            RNG seed (default: 42)
#   BENCH_REPS            repetitions per preset (default: 3)
#   BENCH_BASELINE_FILE   baseline CSV path (default: inst/benchmarks/mcmc_baseline.csv)
#   BENCH_REGRESSION_THRESHOLD  max allowed wall-clock regression % (default: 5)

script_dir <- dirname(normalizePath(sub("--file=(.*)", "\\1", commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))][1], perl = TRUE)))
helpers <- file.path(script_dir, "bench_helpers.R")
if (!file.exists(helpers)) helpers <- file.path(getwd(), "inst/scripts/bench_helpers.R")
source(helpers)

args <- commandArgs(trailingOnly = TRUE)
mode <- "run"
parallel_mode <- NULL
presets <- character(0)
for (a in args) {
  if (a == "--save") {
    mode <- "save"
  } else if (a == "--compare") {
    mode <- "compare"
  } else if (a == "--pt") {
    parallel_mode <- "pt"
  } else if (a == "--serial") {
    parallel_mode <- "serial"
  } else if (a == "--single") {
    parallel_mode <- "single"
  } else if (!startsWith(a, "--")) {
    presets <- c(presets, a)
  }
}
if (length(presets) == 0) presets <- bench_default_preset()
if ("all" %in% presets) presets <- names(bench_preset_specs())

parallel <- bench_parallel_config(parallel_mode)

if (!requireNamespace("moire", quietly = TRUE)) {
  stop("Package 'moire' must be installed/loaded (e.g. devtools::load_all()).")
}

burnin <- as.integer(Sys.getenv("BENCH_BURNIN", "200"))
samples <- as.integer(Sys.getenv("BENCH_SAMPLES", "200"))
seed <- as.integer(Sys.getenv("BENCH_SEED", "42"))
reps <- as.integer(Sys.getenv("BENCH_REPS", "3"))
threshold <- as.numeric(Sys.getenv("BENCH_REGRESSION_THRESHOLD", "5"))
baseline_file <- Sys.getenv(
  "BENCH_BASELINE_FILE",
  file.path(getwd(), "inst/benchmarks/mcmc_baseline.csv")
)

set.seed(seed)
git_sha <- bench_git_sha()

message("MCMC baseline benchmark")
message("  mode: ", mode)
message("  presets: ", paste(presets, collapse = ", "))
message("  parallel: ", bench_format_parallel_config(parallel))
message("  burnin=", burnin, " samples_per_chain=", samples, " seed=", seed, " reps=", reps)
message("  git: ", git_sha)

all_rows <- list()

for (preset in presets) {
  loaded <- bench_load_data(preset)
  data <- loaded$data
  effective_preset <- bench_effective_preset(preset, parallel)
  elapsed <- numeric(reps)
  last_stats <- data.frame()

  for (rep in seq_len(reps)) {
    message("Running ", effective_preset, " (rep ", rep, "/", reps, ") ...")
    result <- bench_run_once(
      data,
      burnin,
      samples,
      verbose = FALSE,
      parallel = parallel
    )
    elapsed[rep] <- result$elapsed_s
    last_stats <- result$stats
    message("  elapsed: ", round(result$elapsed_s, 3), " s")
  }

  message(
    "Preset ", effective_preset, ": mean=", round(mean(elapsed), 3),
    " s, sd=", round(stats::sd(elapsed), 3), " s"
  )

  if (nrow(last_stats) > 0) {
    message("Top profiler keys (last rep):")
    print(head(bench_profiler_fractions(last_stats), 10))
  } else {
    message("Profiler returned no data. Rebuild with MOIRE_ENABLE_PROFILER_REGISTRY.")
  }

  all_rows[[length(all_rows) + 1L]] <- bench_summary_rows(
    effective_preset, burnin, samples, seed, reps, elapsed, git_sha, parallel
  )
  all_rows[[length(all_rows) + 1L]] <- bench_profiler_rows(effective_preset, last_stats)
}

current <- do.call(rbind, all_rows)

if (mode == "save") {
  if (file.exists(baseline_file) && parallel$mode == "single") {
    existing <- bench_read_baseline(baseline_file)
    keep <- existing[!existing$preset %in% current$preset, , drop = FALSE]
    current <- rbind(keep, current)
  } else if (file.exists(baseline_file) && parallel$mode != "single") {
    existing <- bench_read_baseline(baseline_file)
    keep <- existing[!existing$preset %in% current$preset, , drop = FALSE]
    current <- rbind(keep, current)
  }
  bench_write_csv(current, baseline_file)
  message("\nBaseline saved to: ", normalizePath(baseline_file, mustWork = FALSE))
} else if (mode == "compare") {
  baseline <- bench_read_baseline(baseline_file)
  effective_presets <- vapply(
    presets,
    function(p) bench_effective_preset(p, parallel),
    character(1)
  )
  comp <- bench_compare(
    current,
    baseline,
    threshold_pct = threshold,
    presets = effective_presets
  )
  bench_print_comparison(comp)
  if (bench_has_regression(comp)) {
    message("\nRegression detected (>", threshold, "% slower than baseline).")
    quit(status = 1)
  }
  message("\nNo wall-clock regression beyond threshold.")
} else {
  message("\nRun with --save to write ", baseline_file)
  message("Run with --compare to check against the saved baseline.")
}
