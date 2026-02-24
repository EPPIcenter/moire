# MCMC tests (validates C++ MultiVector/RaggedMultiVector path via run_mcmc)
# Skip on CRAN to avoid timeouts; MCMC run can be slow.
skip_on_cran()

test_that("run_mcmc runs without error with minimal data (multivector refactor compatibility)", {
  # Minimal data: 3 samples, 2 loci, 2 alleles each
  data <- list(
    sample_ids = c("s1", "s2", "s3"),
    loci = c("L1", "L2"),
    data = list(
      list(c(1L, 0L), c(0L, 1L), c(1L, 1L)),
      list(c(1L, 0L), c(1L, 1L), c(0L, 1L))
    )
  )
  # Provide initial allele frequencies so we do not depend on clustering init
  initial_allele_frequencies <- list(
    list(c(0.5, 0.5), c(0.5, 0.5))  # 1 population, 2 loci, 2 alleles each
  )
  expect_error(
    res <- moire::run_mcmc(
      data,
      burnin = 2L,
      samples_per_chain = 2L,
      verbose = FALSE,
      num_populations = 1L,
      initial_allele_frequencies = initial_allele_frequencies
    ),
    NA
  )
  expect_true(is.list(res))
  expect_true("chains" %in% names(res))
  expect_length(res$chains, 1L)
  expect_true("args" %in% names(res))
})
