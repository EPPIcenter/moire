tiny_data <- function() {
  set.seed(1)
  moire::simulate_data(
    mean_coi = 2,
    num_samples = 3,
    epsilon_pos = 0.01,
    epsilon_neg = 0.1,
    locus_freq_alphas = list(rep(1, 3), rep(1, 3)),
    missingness = 0
  )
}

run_tiny <- function(...) {
  data <- tiny_data()
  args <- list(
    data = data,
    is_missing = data$is_missing,
    verbose = FALSE,
    burnin = 5,
    samples_per_chain = 5,
    adapt_temp = FALSE,
    max_coi = 8
  )
  dots <- list(...)
  args[names(dots)] <- dots
  do.call(moire::run_mcmc, args)
}

test_that("NULL seed leaves room for per-chain offsets", {
  res <- run_tiny(seed = NULL, num_chains = 3, num_cores = 1)
  expect_false(anyNA(res$chain_seeds))
  expect_equal(res$chain_seeds, res$seed + c(0L, 1000L, 2000L))
  expect_lte(max(res$chain_seeds), .Machine$integer.max)
})

test_that("user seed that would overflow chain offsets is rejected", {
  data <- tiny_data()
  expect_error(
    moire::run_mcmc(
      data,
      is_missing = data$is_missing,
      verbose = FALSE,
      burnin = 1,
      samples_per_chain = 1,
      seed = .Machine$integer.max,
      num_chains = 2
    ),
    "too large for num_chains"
  )
})

test_that("resolved seed is stored without leaking internals into args", {
  res <- run_tiny(seed = NULL)
  expect_equal(res$args$seed, res$seed)
  expect_false("run_completed" %in% names(res$args))
  expect_false(".moire_completed" %in% names(res$args))
  expect_false("chain_seeds" %in% names(res$args))
  expect_named(res$args, names(formals(moire::run_mcmc)), ignore.order = TRUE)
})

test_that("the same seed replays a serial chain", {
  a <- run_tiny(seed = 42L)
  b <- run_tiny(seed = 42L)
  expect_equal(a$chains[[1]]$llik_sample, b$chains[[1]]$llik_sample)
  expect_equal(a$seed, 42L)
})

test_that("PT multi-chain replay does not depend on num_cores", {
  a <- run_tiny(
    seed = 13L, num_chains = 2, num_cores = 1, pt_chains = 3, burnin = 15,
    samples_per_chain = 10
  )
  b <- run_tiny(
    seed = 13L, num_chains = 2, num_cores = 2, pt_chains = 3, burnin = 15,
    samples_per_chain = 10
  )
  expect_equal(a$chains[[1]]$llik_sample, b$chains[[1]]$llik_sample)
  expect_equal(a$chains[[2]]$llik_sample, b$chains[[2]]$llik_sample)
  expect_equal(a$chains[[1]]$swap_store, b$chains[[1]]$swap_store)
})
