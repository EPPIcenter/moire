#!/usr/bin/env Rscript
# Summarize MCMC profiler hotspots for optimization targeting.
#
# Usage: Rscript inst/scripts/analyze_hotspots.R [preset]
# Env: PROFILE_BURNIN, PROFILE_SAMPLES, PROFILE_SEED (same as profile_mcmc.R)

script_dir <- dirname(normalizePath(sub("--file=(.*)", "\\1", commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))][1], perl = TRUE)))
source(file.path(script_dir, "bench_helpers.R"))

args <- commandArgs(trailingOnly = TRUE)
preset <- if (length(args) >= 1) args[1] else bench_default_preset()

if (!requireNamespace("moire", quietly = TRUE)) {
  stop("Package 'moire' must be installed/loaded.")
}

burnin <- as.integer(Sys.getenv("PROFILE_BURNIN", "200"))
samples <- as.integer(Sys.getenv("PROFILE_SAMPLES", "200"))
seed <- as.integer(Sys.getenv("PROFILE_SEED", "42"))

set.seed(seed)
loaded <- bench_load_data(preset)
result <- bench_run_once(loaded$data, burnin, samples, verbose = FALSE)

message("\nWall-clock: ", round(result$elapsed_s, 3), " s")
stats <- bench_profiler_fractions(result$stats)
if (nrow(stats) == 0) stop("No profiler data.")

focus <- c(
  "Chain::update_p",
  "Chain::update_p::build_logit",
  "Chain::update_p::expit",
  "Chain::update_p::recalc_transmission",
  "Chain::update_p::pam_groups",
  "Chain::recalculate_transmission_for_sample",
  "Chain::update_p::calc_post",
  "Chain::update_p::reject_restore",
  "Chain::update_p::accept_save",
  "Chain::calculate_transmission_likelihood",
  "Chain::calc_transmission_process::all",
  "Chain::calc_transmission_process::loop",
  "Chain::calc_transmission_process::p_update",
  "Chain::calc_transmission_process::p_inc",
  "Chain::calc_transmission_process::p_full",
  "Chain::calc_transmission_process::r_update",
  "Chain::calc_transmission_process::r_inc",
  "Chain::calculate_transmission_likelihood_after_p_change",
  "Chain::calculate_transmission_likelihood_after_r_change",
  "Chain::update_r::recalc_transmission",
  "Chain::pam_vec_cache_hit",
  "Chain::pam_vec_cache_miss",
  "Chain::pam_vec_fast_low_k",
  "Chain::pam_vec_fast_low_k_inline",
  "Chain::pam_vec_gray_code",
  "Chain::calc_new_likelihood",
  "Chain::calc_new_likelihood::transmission_sum",
  "Chain::calc_new_likelihood::obs_sum"
)

sub <- stats[stats$key %in% focus, c("key", "calls", "total_ms", "avg_ms", "pct")]
sub <- sub[order(-sub$total_ms), ]
message("\nHotspot breakdown:")
print(sub, row.names = FALSE)

update_p_ms <- sum(sub$total_ms[sub$key == "Chain::update_p" | grepl("^Chain::update_p::", sub$key)])
tx_ms <- sum(sub$total_ms[grepl("transmission", sub$key, ignore.case = TRUE)])
message("\nApprox update_p subtree: ", round(update_p_ms, 1), " ms (",
        round(100 * update_p_ms / sum(stats$total_ms), 1), "% of profiled C++)")
message("Approx transmission subtree: ", round(tx_ms, 1), " ms (",
        round(100 * tx_ms / sum(stats$total_ms), 1), "% of profiled C++)")
