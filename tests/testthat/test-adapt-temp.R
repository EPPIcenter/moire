pt_data <- function() {
  set.seed(7)
  moire::simulate_data(
    mean_coi = 3,
    num_samples = 12,
    epsilon_pos = 0.01,
    epsilon_neg = 0.1,
    locus_freq_alphas = rep(list(rep(1, 5)), 6),
    missingness = 0
  )
}

run_pt <- function(adapt_temp, pt_chains = 8) {
  data <- pt_data()
  moire::run_mcmc(
    data,
    is_missing = data$is_missing,
    verbose = FALSE,
    burnin = 200,
    samples_per_chain = 20,
    pt_chains = pt_chains,
    adapt_temp = adapt_temp,
    pre_adapt_steps = 25,
    temp_adapt_steps = 25,
    seed = 3L
  )
}

test_that("adapted temperature ladder is finite, monotone, and keeps its ends", {
  res <- run_pt(adapt_temp = TRUE)
  ladder <- res$chains[[1]]$temp_gradient
  expect_length(ladder, 8)
  expect_true(all(is.finite(ladder)))
  expect_equal(ladder[1], 1)
  expect_equal(ladder[length(ladder)], 0)
  expect_true(all(diff(ladder) < 0))
  expect_false(isTRUE(all.equal(ladder, seq(1, 0, length.out = 8))))
})

test_that("temperature ladder is untouched when adapt_temp = FALSE", {
  res <- run_pt(adapt_temp = FALSE)
  expect_equal(
    res$chains[[1]]$temp_gradient, seq(1, 0, length.out = 8),
    tolerance = 1e-6
  )
})

test_that("adaptation works with a user supplied temperature vector", {
  temps <- c(1, 0.7, 0.4, 0.2, 0.05)
  res <- run_pt(adapt_temp = TRUE, pt_chains = temps)
  ladder <- res$chains[[1]]$temp_gradient
  expect_length(ladder, length(temps))
  expect_equal(ladder[1], 1)
  expect_equal(ladder[length(ladder)], 0.05)
  expect_true(all(diff(ladder) < 0))
})
