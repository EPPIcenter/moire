test_that("initialization failure returns classed condition with diagnostics", {
  # Tiny synthetic dataset: one locus with many alleles tends to stress init.
  # Force failure with max_initialization_tries = 0 after the constructor draw
  # by using tries = 0 when the first draw is already non-finite — not reliable.
  # Instead, exercise the R formatter / condition shape directly, and a live
  # run_mcmc path that should succeed on a simple dataset.
  diagnostics <- list(
    max_initialization_tries = 10L,
    chains_attempted = 1L,
    total_ill_conditioned = 10L,
    unknown_genotyping_count = 0L,
    locus_counts = list(locus = 2L, n = 9L),
    sample_counts = list(sample = 1L, n = 9L),
    prior_nonfinite_counts = list(term = character(), n = integer()),
    examples = list(sample = 1L, locus = 2L),
    ill_conditioned_per_chain = 10L,
    classification = "consistent_failure_cause",
    concentration = "locus",
    dominant_locus = 2L,
    dominant_sample = NA_integer_,
    dominant_share = 0.9,
    secondary_sample = 1L
  )

  err <- tryCatch(
    stop_initialization_failure(
      diagnostics,
      sample_ids = c("S1", "S2"),
      loci = c("L1", "L2")
    ),
    error = function(e) e
  )

  expect_s3_class(err, "moire_initialization_failure")
  expect_true(!is.null(err$diagnostics))
  expect_equal(err$diagnostics$classification, "consistent_failure_cause")
  expect_match(conditionMessage(err), "Consistent failure cause: locus 2 \\(L2\\)")
  expect_match(conditionMessage(err), "Guidance: inspect this locus")
  expect_match(conditionMessage(err), "Secondary concentration on sample 1 \\(S1\\)")
})

test_that("hard starting set guidance is emitted", {
  diagnostics <- list(
    max_initialization_tries = 5L,
    chains_attempted = 1L,
    total_ill_conditioned = 5L,
    unknown_genotyping_count = 0L,
    locus_counts = list(locus = c(1L, 2L), n = c(2L, 3L)),
    sample_counts = list(sample = c(1L, 2L), n = c(2L, 3L)),
    prior_nonfinite_counts = list(term = character(), n = integer()),
    examples = list(sample = integer(), locus = integer()),
    ill_conditioned_per_chain = 5L,
    classification = "hard_starting_set",
    concentration = "",
    dominant_locus = 2L,
    dominant_sample = NA_integer_,
    dominant_share = 0.6,
    secondary_sample = NA_integer_
  )

  msg <- format_initialization_failure_message(diagnostics)
  expect_match(msg, "hard starting set")
  expect_match(msg, "increase max_initialization_tries")
})
