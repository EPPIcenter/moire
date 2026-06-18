# Shared helpers for MCMC wall-clock and profiler benchmarks.

bench_git_sha <- function() {
  out <- tryCatch(
    system2("git", c("rev-parse", "--short", "HEAD"), stdout = TRUE, stderr = FALSE),
    error = function(e) character(0)
  )
  if (length(out) == 0 || !nzchar(out[[1]])) "unknown" else out[[1]]
}

bench_preset_specs <- function() {
  list(
    minimal = list(
      label = "minimal (4 samples x 3 loci, long-form)",
      n_samples = 4L,
      n_loci = 3L,
      kind = "long_form"
    ),
    small = list(
      label = "small (20 samples x 10 loci)",
      n_samples = 20L,
      n_loci = 10L,
      kind = "simulated"
    ),
    medium = list(
      label = "medium (40 samples x 50 loci)",
      n_samples = 40L,
      n_loci = 50L,
      kind = "simulated"
    ),
    vignette = list(
      label = "vignette (100 samples x 100 loci)",
      n_samples = 100L,
      n_loci = 100L,
      kind = "simulated"
    )
  )
}

#' Default MCMC benchmark preset (override with BENCH_PRESET env).
#' Use vignette for wall-clock and cache-sensitive work; small for quick smoke tests.
bench_default_preset <- function() {
  preset <- Sys.getenv("BENCH_PRESET", "vignette")
  specs <- bench_preset_specs()
  if (!preset %in% names(specs)) {
    stop("Unknown BENCH_PRESET: ", preset, ". Choose one of: ", paste(names(specs), collapse = ", "))
  }
  preset
}

bench_load_data <- function(preset) {
  specs <- bench_preset_specs()
  if (!preset %in% names(specs)) {
    stop("Unknown preset: ", preset, ". Choose one of: ", paste(names(specs), collapse = ", "))
  }
  spec <- specs[[preset]]
  message("Data: ", spec$label)

  if (spec$kind == "long_form") {
    n_samples <- spec$n_samples
    n_loci <- spec$n_loci
    sample_ids <- rep(seq_len(n_samples), each = n_loci * 2)
    loci <- rep(rep(seq_len(n_loci), each = 2), times = n_samples)
    alleles <- rep(c(1L, 2L), times = n_samples * n_loci)
    df <- data.frame(sample_id = sample_ids, locus = loci, allele = alleles)
    data <- moire::load_long_form_data(df)
  } else {
    allele_counts <- if (spec$n_loci <= 10L) {
      c(rep(5L, spec$n_loci %/% 2L), rep(10L, spec$n_loci - spec$n_loci %/% 2L))
    } else if (spec$n_loci <= 50L) {
      c(rep(5L, spec$n_loci %/% 2L), rep(10L, spec$n_loci - spec$n_loci %/% 2L))
    } else {
      c(rep(5L, spec$n_loci %/% 2L), rep(10L, spec$n_loci - spec$n_loci %/% 2L))
    }
    locus_freq_alphas <- lapply(allele_counts, function(a) rep(1, a))
    data <- moire::simulate_data(
      mean_coi = 3,
      num_samples = spec$n_samples,
      epsilon_pos = 0.01,
      epsilon_neg = 0.1,
      locus_freq_alphas = locus_freq_alphas,
      internal_relatedness_alpha = 0.1,
      internal_relatedness_beta = 1
    )
  }

  list(data = data, spec = spec)
}

bench_profiler_fractions <- function(stats) {
  if (nrow(stats) == 0) return(stats)
  total_ms <- sum(stats$total_ms)
  if (total_ms <= 0) {
    stats$pct <- 0
    return(stats)
  }
  stats$pct <- round(100 * stats$total_ms / total_ms, 1)
  stats[order(-stats$total_ms), ]
}

bench_run_once <- function(data, burnin, samples, verbose = FALSE) {
  moire::moire_prof_reset()
  timing <- system.time(
    invisible(moire::run_mcmc(data, burnin = burnin, samples_per_chain = samples, verbose = verbose))
  )
  stats <- tryCatch(
    moire::moire_prof_stats(),
    error = function(e) data.frame()
  )
  list(elapsed_s = timing[["elapsed"]], stats = stats)
}

bench_summary_rows <- function(preset, burnin, samples, seed, reps, elapsed, git_sha) {
  rbind(
    data.frame(
      section = "summary",
      preset = preset,
      metric = c("wall_elapsed_mean_s", "wall_elapsed_sd_s", "burnin", "samples_per_chain", "seed", "reps"),
      value = c(mean(elapsed), stats::sd(elapsed), burnin, samples, seed, reps),
      unit = c("s", "s", "iter", "iter", "seed", "count"),
      calls = rep(NA_integer_, 6L),
      stringsAsFactors = FALSE
    ),
    data.frame(
      section = "meta",
      preset = preset,
      metric = "git_sha",
      value = NA_real_,
      unit = git_sha,
      calls = NA_integer_,
      stringsAsFactors = FALSE
    )
  )
}

bench_profiler_rows <- function(preset, stats) {
  if (nrow(stats) == 0) return(data.frame())
  stats <- bench_profiler_fractions(stats)
  data.frame(
    section = "profiler",
    preset = preset,
    metric = stats$key,
    value = stats$total_ms,
    unit = "ms_total",
    calls = as.integer(stats$calls),
    stringsAsFactors = FALSE
  )
}

bench_write_csv <- function(rows, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.csv(rows, path, row.names = FALSE)
  invisible(path)
}

bench_read_baseline <- function(path) {
  if (!file.exists(path)) {
    stop("Baseline file not found: ", path, "\nRun with --save to create one.")
  }
  df <- read.csv(path, stringsAsFactors = FALSE)
  required <- c("section", "preset", "metric", "value", "unit", "calls")
  missing <- setdiff(required, names(df))
  if (length(missing) > 0) {
    stop("Baseline file missing columns: ", paste(missing, collapse = ", "))
  }
  df
}

bench_compare <- function(current, baseline, threshold_pct = 5, presets = NULL) {
  if (!is.null(presets)) {
    baseline <- baseline[baseline$preset %in% presets, , drop = FALSE]
    current <- current[current$preset %in% presets, , drop = FALSE]
  }
  cur_summary <- current[current$section == "summary" & current$metric == "wall_elapsed_mean_s", , drop = FALSE]
  base_summary <- baseline[baseline$section == "summary" & baseline$metric == "wall_elapsed_mean_s", , drop = FALSE]

  merged <- merge(
    cur_summary[, c("preset", "value")],
    base_summary[, c("preset", "value")],
    by = "preset",
    suffixes = c("_current", "_baseline"),
    all = TRUE
  )
  names(merged) <- c("preset", "current_s", "baseline_s")
  merged$pct_change <- ifelse(
    is.finite(merged$baseline_s) & merged$baseline_s > 0,
    100 * (merged$current_s - merged$baseline_s) / merged$baseline_s,
    NA_real_
  )
  merged$status <- ifelse(
    is.na(merged$pct_change),
    "missing",
    ifelse(merged$pct_change > threshold_pct, "regression",
           ifelse(merged$pct_change < -threshold_pct, "improvement", "stable"))
  )

  cur_prof <- current[current$section == "profiler", , drop = FALSE]
  base_prof <- baseline[baseline$section == "profiler", , drop = FALSE]
  prof_merged <- merge(
    cur_prof[, c("preset", "metric", "value")],
    base_prof[, c("preset", "metric", "value")],
    by = c("preset", "metric"),
    suffixes = c("_current", "_baseline"),
    all = FALSE
  )
  names(prof_merged) <- c("preset", "metric", "current_ms", "baseline_ms")
  prof_merged$pct_change <- ifelse(
    prof_merged$baseline_ms > 0,
    100 * (prof_merged$current_ms - prof_merged$baseline_ms) / prof_merged$baseline_ms,
    NA_real_
  )
  prof_merged <- prof_merged[order(-abs(prof_merged$pct_change)), ]

  list(summary = merged, profiler = prof_merged, threshold_pct = threshold_pct)
}

bench_print_comparison <- function(comp) {
  message("\n=== Wall-clock comparison (threshold: ", comp$threshold_pct, "%) ===")
  print(comp$summary[, c("preset", "current_s", "baseline_s", "pct_change", "status")], row.names = FALSE)

  if (nrow(comp$profiler) > 0) {
    message("\n=== Profiler hotspots with largest absolute change ===")
    print(head(comp$profiler, 15), row.names = FALSE)
  }
}

bench_has_regression <- function(comp) {
  any(comp$summary$status == "regression", na.rm = TRUE)
}
