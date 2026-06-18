skip_on_cran()

test_that("run_mcmc runs with count barcodes and observation_model = counts", {
  data <- list(
    sample_ids = c("s1", "s2"),
    loci = c("L1"),
    aggregate = "count",
    data = list(list(c(8L, 10L), c(5L, 0L))),
    is_missing = matrix(FALSE, nrow = 1, ncol = 2)
  )
  initial_allele_frequencies <- list(list(c(0.5, 0.5)))

  expect_error(
    res <- moire::run_mcmc(
      data,
      observation_model = "counts",
      burnin = 1L,
      samples_per_chain = 1L,
      verbose = FALSE,
      num_populations = 1L,
      initial_allele_frequencies = initial_allele_frequencies
    ),
    NA
  )
  expect_true(is.list(res))
  expect_length(res$chains, 1L)
})

test_that("run_mcmc errors when count barcodes are used with binary observation model", {
  data <- list(
    sample_ids = "s1",
    loci = "L1",
    aggregate = "count",
    data = list(list(c(8L, 10L))),
    is_missing = matrix(FALSE, nrow = 1, ncol = 1)
  )

  expect_error(
    moire::run_mcmc(
      data,
      observation_model = "binary",
      burnin = 1L,
      samples_per_chain = 1L,
      verbose = FALSE,
      num_populations = 1L
    ),
    "Count barcodes"
  )
})
