#!/usr/bin/env Rscript
# Wall-clock MCMC timing without profiler registry overhead.
#
# Usage:
#   MOIRE_DISABLE_PROFILER_REGISTRY=1 Rscript inst/scripts/bench_wall_clock.R [preset]
#   Rscript inst/scripts/bench_wall_clock.R small --compare-baseline
#   Rscript inst/scripts/bench_wall_clock.R vignette --pt
#
# Env: same as bench_mcmc_baseline.R (BENCH_* including BENCH_PARALLEL_MODE)

script_dir <- dirname(normalizePath(sub("--file=(.*)", "\\1", commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))][1], perl = TRUE)))
source(file.path(script_dir, "bench_helpers.R"))

args <- commandArgs(trailingOnly = TRUE)
compare_baseline <- "--compare-baseline" %in% args
parallel_mode <- NULL
if ("--pt" %in% args) parallel_mode <- "pt"
if ("--serial" %in% args) parallel_mode <- "serial"
if ("--single" %in% args) parallel_mode <- "single"
preset <- setdiff(args, c("--compare-baseline", "--pt", "--serial", "--single"))[1]
preset <- if (is.na(preset) || !nzchar(preset)) bench_default_preset()
parallel <- bench_parallel_config(parallel_mode)
effective_preset <- bench_effective_preset(preset, parallel)

if (!requireNamespace("moire", quietly = TRUE)) {
  stop("Package 'moire' must be installed/loaded (e.g. devtools::load_all()).")
}

burnin <- as.integer(Sys.getenv("BENCH_BURNIN", "200"))
samples <- as.integer(Sys.getenv("BENCH_SAMPLES", "200"))
seed <- as.integer(Sys.getenv("BENCH_SEED", "42"))
reps <- as.integer(Sys.getenv("BENCH_REPS", "5"))
baseline_file <- Sys.getenv(
  "BENCH_BASELINE_FILE",
  file.path(getwd(), "inst/benchmarks/mcmc_baseline.csv")
)

profiler_disabled <- identical(Sys.getenv("MOIRE_DISABLE_PROFILER_REGISTRY"), "1")
mode_label <- if (profiler_disabled) "no_profiler" else "profiler_on"

set.seed(seed)
loaded <- bench_load_data(preset)
data <- loaded$data

message("Wall-clock benchmark (profiler registry: ", if (profiler_disabled) "OFF" else "ON", ")")
message("  preset=", effective_preset, " parallel: ", bench_format_parallel_config(parallel))
message("  burnin=", burnin, " samples=", samples, " seed=", seed, " reps=", reps)

elapsed <- numeric(reps)
for (rep in seq_len(reps)) {
  result <- bench_run_once(
    data,
    burnin = burnin,
    samples = samples,
    verbose = FALSE,
    parallel = parallel
  )
  elapsed[rep] <- result$elapsed_s
  message("  rep ", rep, "/", reps, ": ", round(elapsed[rep], 3), " s")
}

mean_s <- mean(elapsed)
sd_s <- stats::sd(elapsed)
message("\nResult: mean=", round(mean_s, 3), " s, sd=", round(sd_s, 3), " s (",
        mode_label, ")")

if (compare_baseline) {
  baseline <- bench_read_baseline(baseline_file)
  base_row <- baseline[baseline$section == "summary" &
                         baseline$preset == effective_preset &
                         baseline$metric == "wall_elapsed_mean_s", ]
  if (nrow(base_row) == 0) {
    message("No baseline wall time for preset '", effective_preset, "' in ", baseline_file)
  } else {
    base_s <- base_row$value[1]
    pct <- 100 * (mean_s - base_s) / base_s
    message("\nBaseline (profiler on, from CSV): ", round(base_s, 3), " s")
    message("Delta: ", sprintf("%+.1f%%", pct),
            if (pct < 0) " (faster)" else if (pct > 0) " (slower)" else "")
    base_mode <- bench_meta_value(baseline, effective_preset, "parallel_mode")
    if (is.na(base_mode)) base_mode <- "single"
    if (!identical(parallel$mode, base_mode)) {
      message("Note: parallel mode ", parallel$mode, " vs baseline ", base_mode)
    }
    base_threads <- bench_meta_value(baseline, effective_preset, "num_threads", "value")
    if (!is.na(base_threads) && as.character(parallel$num_threads) != base_threads) {
      message("Note: num_threads ", parallel$num_threads, " vs baseline ", base_threads)
    }
  }
}

invisible(list(
  preset = effective_preset,
  mode = mode_label,
  parallel = parallel,
  mean_s = mean_s,
  sd_s = sd_s,
  elapsed = elapsed
))
