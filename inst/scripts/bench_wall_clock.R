#!/usr/bin/env Rscript
# Wall-clock MCMC timing without profiler registry overhead.
#
# Usage:
#   MOIRE_DISABLE_PROFILER_REGISTRY=1 Rscript inst/scripts/bench_wall_clock.R [preset]
#   Rscript inst/scripts/bench_wall_clock.R small --compare-baseline
#
# Env: BENCH_BURNIN, BENCH_SAMPLES, BENCH_SEED, BENCH_REPS (same as bench_mcmc_baseline.R)
#      BENCH_BASELINE_FILE

script_dir <- dirname(normalizePath(sub("--file=(.*)", "\\1", commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))][1], perl = TRUE)))
source(file.path(script_dir, "bench_helpers.R"))

args <- commandArgs(trailingOnly = TRUE)
compare_baseline <- "--compare-baseline" %in% args
preset <- setdiff(args, "--compare-baseline")[1]
preset <- if (is.na(preset) || !nzchar(preset)) bench_default_preset()

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
message("  preset=", preset, " burnin=", burnin, " samples=", samples,
        " seed=", seed, " reps=", reps)

elapsed <- numeric(reps)
for (rep in seq_len(reps)) {
  timing <- system.time(
    invisible(moire::run_mcmc(
      data,
      burnin = burnin,
      samples_per_chain = samples,
      verbose = FALSE
    ))
  )
  elapsed[rep] <- timing[["elapsed"]]
  message("  rep ", rep, "/", reps, ": ", round(elapsed[rep], 3), " s")
}

mean_s <- mean(elapsed)
sd_s <- stats::sd(elapsed)
message("\nResult: mean=", round(mean_s, 3), " s, sd=", round(sd_s, 3), " s (",
        mode_label, ")")

if (compare_baseline) {
  baseline <- bench_read_baseline(baseline_file)
  base_row <- baseline[baseline$section == "summary" &
                         baseline$preset == preset &
                         baseline$metric == "wall_elapsed_mean_s", ]
  if (nrow(base_row) == 0) {
    message("No baseline wall time for preset '", preset, "' in ", baseline_file)
  } else {
    base_s <- base_row$value[1]
    pct <- 100 * (mean_s - base_s) / base_s
    message("\nBaseline (profiler on, from CSV): ", round(base_s, 3), " s")
    message("Delta: ", sprintf("%+.1f%%", pct),
            if (pct < 0) " (faster)" else if (pct > 0) " (slower)" else "")
    message("Note: baseline CSV was recorded with profiler registry enabled.")
  }
}

invisible(list(
  preset = preset,
  mode = mode_label,
  mean_s = mean_s,
  sd_s = sd_s,
  elapsed = elapsed
))
