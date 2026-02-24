test_that("find_optimal_permutation finds correct permutation", {
  # Create two assignment matrices with known permutation
  # Step 1: pop1 has high values, pop2 has low values
  # Step 2: labels are swapped (pop1 in step2 = pop2 in step1)
  assignments_t <- matrix(c(
    0.9, 0.1,  # sample 1: mostly pop1
    0.8, 0.2,  # sample 2: mostly pop1
    0.1, 0.9   # sample 3: mostly pop2
  ), nrow = 3, ncol = 2, byrow = TRUE)
  
  # Step 2: labels swapped (what was pop1 is now pop2, and vice versa)
  assignments_t1 <- matrix(c(
    0.1, 0.9,  # sample 1: mostly pop2 (but should be pop1)
    0.2, 0.8,  # sample 2: mostly pop2 (but should be pop1)
    0.9, 0.1   # sample 3: mostly pop1 (but should be pop2)
  ), nrow = 3, ncol = 2, byrow = TRUE)
  
  permutation <- moire:::find_optimal_permutation(assignments_t, assignments_t1)
  
  # Should find that pop1 in step1 maps to pop2 in step2, and pop2 maps to pop1
  # So permutation[1] = 2 and permutation[2] = 1
  expect_equal(permutation[1], 2)
  expect_equal(permutation[2], 1)
})

test_that("apply_permutation correctly reorders columns", {
  assignments <- matrix(c(
    0.1, 0.9,  # sample 1
    0.2, 0.8,  # sample 2
    0.9, 0.1   # sample 3
  ), nrow = 3, ncol = 2, byrow = TRUE)
  
  # Permutation: pop1 -> pop2, pop2 -> pop1
  permutation <- c(2, 1)
  
  result <- moire:::apply_permutation(assignments, permutation)
  
  # After permutation, column 1 should have values from original column 2
  expect_equal(result[, 1], assignments[, 2])
  expect_equal(result[, 2], assignments[, 1])
})

test_that("correct_label_switching with iterative method works", {
  # Create population assignments with label switching
  # Step 1: sample1 -> pop1, sample2 -> pop2
  # Step 2: labels swapped - sample1 -> pop2, sample2 -> pop1
  # Step 3: labels swapped again - sample1 -> pop1, sample2 -> pop2
  
  # Step 1: log probabilities
  step1 <- list(
    list(log(0.9), log(0.1)),  # sample 1: mostly pop1
    list(log(0.1), log(0.9))   # sample 2: mostly pop2
  )
  
  # Step 2: labels swapped (log probabilities)
  step2 <- list(
    list(log(0.1), log(0.9)),  # sample 1: mostly pop2 (swapped)
    list(log(0.9), log(0.1))   # sample 2: mostly pop1 (swapped)
  )
  
  # Step 3: labels swapped back (log probabilities)
  step3 <- list(
    list(log(0.9), log(0.1)),  # sample 1: mostly pop1 (swapped back)
    list(log(0.1), log(0.9))   # sample 2: mostly pop2 (swapped back)
  )
  
  population_assignments <- list(step1, step2, step3)
  
  corrected <- moire::correct_label_switching(population_assignments, method = "iterative")
  
  # After correction, all steps should have consistent labels
  # Sample 1 should always have highest probability in pop1
  # Sample 2 should always have highest probability in pop2
  
  for (step_idx in seq_along(corrected)) {
    step <- corrected[[step_idx]]
    # Sample 1 should have pop1 > pop2
    expect_gt(exp(step[[1]][[1]]), exp(step[[1]][[2]]))
    # Sample 2 should have pop2 > pop1
    expect_gt(exp(step[[2]][[2]]), exp(step[[2]][[1]]))
  }
})

test_that("correct_label_switching with reference method works", {
  # Create population assignments with label switching
  step1 <- list(
    list(log(0.9), log(0.1)),  # sample 1: mostly pop1
    list(log(0.1), log(0.9))   # sample 2: mostly pop2
  )
  
  step2 <- list(
    list(log(0.1), log(0.9)),  # sample 1: mostly pop2 (swapped)
    list(log(0.9), log(0.1))   # sample 2: mostly pop1 (swapped)
  )
  
  step3 <- list(
    list(log(0.1), log(0.9)),  # sample 1: mostly pop2 (still swapped)
    list(log(0.9), log(0.1))   # sample 2: mostly pop1 (still swapped)
  )
  
  population_assignments <- list(step1, step2, step3)
  
  corrected <- moire::correct_label_switching(population_assignments, method = "reference")
  
  # After correction, all steps should align to step1 (reference)
  # Sample 1 should always have highest probability in pop1
  # Sample 2 should always have highest probability in pop2
  
  for (step_idx in seq_along(corrected)) {
    step <- corrected[[step_idx]]
    # Sample 1 should have pop1 > pop2
    expect_gt(exp(step[[1]][[1]]), exp(step[[1]][[2]]))
    # Sample 2 should have pop2 > pop1
    expect_gt(exp(step[[2]][[2]]), exp(step[[2]][[1]]))
  }
})

test_that("correct_label_switching handles edge cases", {
  # Empty assignments
  expect_equal(moire::correct_label_switching(list()), list())
  
  # Single step (no switching possible)
  single_step <- list(
    list(list(log(0.9), log(0.1)), list(log(0.1), log(0.9)))
  )
  result <- moire::correct_label_switching(single_step)
  expect_equal(length(result), 1)
  expect_equal(length(result[[1]]), 2)
  
  # Single population (no switching possible)
  single_pop <- list(
    list(list(log(1.0)), list(log(1.0))),
    list(list(log(1.0)), list(log(1.0)))
  )
  result <- moire::correct_label_switching(single_pop)
  expect_equal(length(result), 2)
})

test_that("correct_label_switching preserves structure", {
  step1 <- list(
    list(log(0.9), log(0.1)),
    list(log(0.1), log(0.9))
  )
  step2 <- list(
    list(log(0.1), log(0.9)),
    list(log(0.9), log(0.1))
  )
  
  population_assignments <- list(step1, step2)
  corrected <- moire::correct_label_switching(population_assignments)
  
  # Structure should be preserved
  expect_equal(length(corrected), length(population_assignments))
  expect_equal(length(corrected[[1]]), length(population_assignments[[1]]))
  expect_equal(length(corrected[[1]][[1]]), length(population_assignments[[1]][[1]]))
})

test_that("summarize_population_assignments with label switching correction works", {
  # Create mock MCMC results with label switching
  mock_mcmc_results <- list(
    args = list(
      data = list(
        sample_ids = c("S1", "S2"),
        loci = c("L1")
      )
    ),
    chains = list(
      list(
        # Step 1: S1 -> pop1, S2 -> pop2
        # Step 2: labels swapped - S1 -> pop2, S2 -> pop1
        population_assignment = list(
          # Sample 1
          list(
            c(log(0.9), log(0.1)),  # step 1: mostly pop1
            c(log(0.1), log(0.9))   # step 2: mostly pop2 (swapped)
          ),
          # Sample 2
          list(
            c(log(0.1), log(0.9)),  # step 1: mostly pop2
            c(log(0.9), log(0.1))   # step 2: mostly pop1 (swapped)
          )
        )
      )
    )
  )
  
  # Without correction, probabilities would be averaged incorrectly
  # With correction, S1 should have high probability for pop1, S2 for pop2
  result_corrected <- moire::summarize_population_assignments(
    mock_mcmc_results,
    merge_chains = TRUE,
    correct_label_switching = TRUE
  )
  
  # Check that correction worked: S1 should have higher prob for pop1 than pop2
  s1_pop1 <- result_corrected$most_likely_population[
    result_corrected$most_likely_population$sample_id == "S1" &
    result_corrected$most_likely_population$population == 1,
    "prob_most_likely"
  ]
  s1_pop2 <- result_corrected$most_likely_population[
    result_corrected$most_likely_population$sample_id == "S1" &
    result_corrected$most_likely_population$population == 2,
    "prob_most_likely"
  ]
  
  expect_gt(s1_pop1, s1_pop2)
  
  # S2 should have higher prob for pop2 than pop1
  s2_pop1 <- result_corrected$most_likely_population[
    result_corrected$most_likely_population$sample_id == "S2" &
    result_corrected$most_likely_population$population == 1,
    "prob_most_likely"
  ]
  s2_pop2 <- result_corrected$most_likely_population[
    result_corrected$most_likely_population$sample_id == "S2" &
    result_corrected$most_likely_population$population == 2,
    "prob_most_likely"
  ]
  
  expect_gt(s2_pop2, s2_pop1)
})

test_that("summarize_population_assignments respects correct_label_switching parameter", {
  # Create mock data
  mock_mcmc_results <- list(
    args = list(
      data = list(
        sample_ids = c("S1"),
        loci = c("L1")
      )
    ),
    chains = list(
      list(
        population_assignment = list(
          list(
            c(log(0.9), log(0.1)),
            c(log(0.1), log(0.9))
          )
        )
      )
    )
  )
  
  # Both should work without errors
  result_with <- moire::summarize_population_assignments(
    mock_mcmc_results,
    correct_label_switching = TRUE
  )
  result_without <- moire::summarize_population_assignments(
    mock_mcmc_results,
    correct_label_switching = FALSE
  )
  
  # Both should return same structure
  expect_equal(names(result_with), names(result_without))
  expect_equal(nrow(result_with$most_likely_population), nrow(result_without$most_likely_population))
})










