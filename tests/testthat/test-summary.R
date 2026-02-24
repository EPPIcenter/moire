# Tests for summarize_coi, summarize_epsilon_neg, summarize_epsilon_pos (merge_chains TRUE/FALSE)

# Minimal barcode data for 3 samples, 2 loci, 2 alleles (for naive_coi / offset_naive_coi in summarize_coi)
mock_data_data <- list(
  list(c(1L, 0L), c(0L, 1L), c(1L, 1L)),   # locus 1: 3 samples
  list(c(1L, 0L), c(1L, 1L), c(0L, 1L))     # locus 2: 3 samples
)

mock_mcmc_results <- function(n_samples = 3L, n_chains = 2L) {
  list(
    args = list(
      data = list(
        sample_ids = paste0("S", seq_len(n_samples)),
        loci = c("L1", "L2"),
        data = mock_data_data
      )
    ),
    chains = lapply(seq_len(n_chains), function(chain_idx) {
      list(
        coi = lapply(seq_len(n_samples), function(s) runif(10L, 1, 3)),
        eps_neg = lapply(seq_len(n_samples), function(s) runif(10L, 0.01, 0.1)),
        eps_pos = lapply(seq_len(n_samples), function(s) runif(10L, 0.01, 0.1))
      )
    })
  )
}

test_that("summarize_coi returns correct structure with merge_chains TRUE", {
  mcmc <- mock_mcmc_results(3L, 2L)
  out <- moire::summarize_coi(mcmc, merge_chains = TRUE)
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 3L)
  expect_true(all(c("sample_id", "post_coi_lower", "post_coi_med", "post_coi_upper", "post_coi_mean",
                    "naive_coi", "offset_naive_coi", "prob_polyclonal") %in% names(out)))
  expect_equal(out$sample_id, mcmc$args$data$sample_ids)
  expect_true(all(out$post_coi_lower <= out$post_coi_med))
  expect_true(all(out$post_coi_med <= out$post_coi_upper))
})

test_that("summarize_coi returns correct structure with merge_chains FALSE", {
  mcmc <- mock_mcmc_results(3L, 2L)
  out <- moire::summarize_coi(mcmc, merge_chains = FALSE)
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 6L)  # 3 samples * 2 chains
  expect_true("chain" %in% names(out))
  expect_equal(sort(unique(out$chain)), c(1L, 2L))
})

test_that("summarize_epsilon_neg returns correct structure with merge_chains TRUE", {
  mcmc <- mock_mcmc_results(3L, 2L)
  out <- moire::summarize_epsilon_neg(mcmc, merge_chains = TRUE)
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 3L)
  expect_true(all(c("sample_id", "post_eps_neg_lower", "post_eps_neg_med", "post_eps_neg_upper", "post_eps_neg_mean") %in% names(out)))
  expect_equal(out$sample_id, mcmc$args$data$sample_ids)
  expect_true(all(out$post_eps_neg_lower <= out$post_eps_neg_med))
  expect_true(all(out$post_eps_neg_med <= out$post_eps_neg_upper))
})

test_that("summarize_epsilon_neg returns correct structure with merge_chains FALSE", {
  mcmc <- mock_mcmc_results(3L, 2L)
  out <- moire::summarize_epsilon_neg(mcmc, merge_chains = FALSE)
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 6L)
  expect_true("chain" %in% names(out))
})

test_that("summarize_epsilon_pos returns correct structure with merge_chains TRUE", {
  mcmc <- mock_mcmc_results(3L, 2L)
  out <- moire::summarize_epsilon_pos(mcmc, merge_chains = TRUE)
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 3L)
  expect_true(all(c("sample_id", "post_eps_pos_lower", "post_eps_pos_med", "post_eps_pos_upper", "post_eps_pos_mean") %in% names(out)))
  expect_true(all(out$post_eps_pos_lower <= out$post_eps_pos_med))
  expect_true(all(out$post_eps_pos_med <= out$post_eps_pos_upper))
})

test_that("summarize_epsilon_pos returns correct structure with merge_chains FALSE", {
  mcmc <- mock_mcmc_results(3L, 2L)
  out <- moire::summarize_epsilon_pos(mcmc, merge_chains = FALSE)
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 6L)
  expect_true("chain" %in% names(out))
})
