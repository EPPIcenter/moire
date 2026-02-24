test_that("summarize_population_assignments works with mock data", {
  # Create mock MCMC results structure
  mock_mcmc_results <- list(
    args = list(
      data = list(
        sample_ids = c("S1", "S2", "S3"),
        loci = c("L1", "L2")
      )
    ),
    chains = list(
      list(
        population_assignment = list(
          # Sample 1: mostly population 1, some population 2
          list(
            c(0.8, 0.2),  # step 1
            c(0.9, 0.1),  # step 2
            c(0.7, 0.3)   # step 3
          ),
          # Sample 2: mostly population 2, some population 1
          list(
            c(0.1, 0.9),  # step 1
            c(0.2, 0.8),  # step 2
            c(0.15, 0.85) # step 3
          ),
          # Sample 3: mixed assignments
          list(
            c(0.6, 0.4),  # step 1
            c(0.4, 0.6),  # step 2
            c(0.5, 0.5)   # step 3
          )
        )
      )
    )
  )
  
  # Test with merged chains
  result <- moire::summarize_population_assignments(mock_mcmc_results, merge_chains = TRUE)
  
  # Check that result is a list with two components
  expect_equal(length(result), 2)
  expect_equal(names(result), c("most_likely_population", "entropy_summary"))
  
  # Check most_likely_population structure
  most_likely <- result$most_likely_population
  expect_equal(ncol(most_likely), 3)
  expect_equal(colnames(most_likely), c("sample_id", "population", "prob_most_likely"))
  expect_equal(nrow(most_likely), 6)  # 3 samples × 2 populations
  
  # Check entropy_summary structure
  entropy <- result$entropy_summary
  expect_equal(ncol(entropy), 5)
  expect_equal(colnames(entropy), c("sample_id", "entropy_lower", "entropy_med", 
                                    "entropy_upper", "entropy_mean"))
  expect_equal(nrow(entropy), 3)  # 3 samples
  
  # Check that all sample IDs are present
  expect_equal(sort(unique(most_likely$sample_id)), c("S1", "S2", "S3"))
  expect_equal(sort(unique(entropy$sample_id)), c("S1", "S2", "S3"))
  
  # Check that all populations are present
  expect_equal(sort(unique(most_likely$population)), c(1, 2))
  
  # Check that probabilities are in valid range
  expect_true(all(most_likely$prob_most_likely >= 0 & most_likely$prob_most_likely <= 1))
  
  # Check that entropy values are non-negative
  expect_true(all(entropy$entropy_lower >= 0))
  expect_true(all(entropy$entropy_med >= 0))
  expect_true(all(entropy$entropy_upper >= 0))
  expect_true(all(entropy$entropy_mean >= 0))
  
  # Check that lower <= median <= upper for entropy
  expect_true(all(entropy$entropy_lower <= entropy$entropy_med))
  expect_true(all(entropy$entropy_med <= entropy$entropy_upper))
})

test_that("summarize_population_assignments works with separate chains", {
  # Create mock MCMC results with multiple chains
  mock_mcmc_results <- list(
    args = list(
      data = list(
        sample_ids = c("S1", "S2"),
        loci = c("L1", "L2")
      )
    ),
    chains = list(
      # Chain 1
      list(
        population_assignment = list(
          # Sample 1: mostly population 1
          list(c(0.8, 0.2), c(0.9, 0.1)),
          # Sample 2: mostly population 2
          list(c(0.1, 0.9), c(0.2, 0.8))
        )
      ),
      # Chain 2
      list(
        population_assignment = list(
          # Sample 1: mostly population 1
          list(c(0.7, 0.3), c(0.85, 0.15)),
          # Sample 2: mostly population 2
          list(c(0.15, 0.85), c(0.1, 0.9))
        )
      )
    )
  )
  
  # Test with separate chains
  result <- moire::summarize_population_assignments(mock_mcmc_results, merge_chains = FALSE)
  
  # Check that result is a list with two components
  expect_equal(length(result), 2)
  expect_equal(names(result), c("most_likely_population", "entropy_summary"))
  
  # Check most_likely_population structure
  most_likely <- result$most_likely_population
  expect_equal(ncol(most_likely), 4)
  expect_equal(colnames(most_likely), c("sample_id", "population", "prob_most_likely", "chain"))
  expect_equal(nrow(most_likely), 8)  # 2 samples × 2 populations × 2 chains
  
  # Check entropy_summary structure
  entropy <- result$entropy_summary
  expect_equal(ncol(entropy), 6)
  expect_equal(colnames(entropy), c("sample_id", "entropy_lower", "entropy_med", 
                                    "entropy_upper", "entropy_mean", "chain"))
  expect_equal(nrow(entropy), 4)  # 2 samples × 2 chains
  
  # Check that all sample IDs are present
  expect_equal(sort(unique(most_likely$sample_id)), c("S1", "S2"))
  expect_equal(sort(unique(entropy$sample_id)), c("S1", "S2"))
  
  # Check that all populations are present
  expect_equal(sort(unique(most_likely$population)), c(1, 2))
  
  # Check that all chains are present
  expect_equal(sort(unique(most_likely$chain)), c(1, 2))
  expect_equal(sort(unique(entropy$chain)), c(1, 2))
  
  # Check that probabilities are in valid range
  expect_true(all(most_likely$prob_most_likely >= 0 & most_likely$prob_most_likely <= 1))
  
  # Check that entropy values are non-negative
  expect_true(all(entropy$entropy_lower >= 0))
  expect_true(all(entropy$entropy_med >= 0))
  expect_true(all(entropy$entropy_upper >= 0))
  expect_true(all(entropy$entropy_mean >= 0))
})

test_that("summarize_population_assignments handles edge cases", {
  # Test with single population
  mock_mcmc_results_single_pop <- list(
    args = list(
      data = list(
        sample_ids = c("S1", "S2"),
        loci = c("L1")
      )
    ),
    chains = list(
      list(
        population_assignment = list(
          # Sample 1: always population 1
          list(c(1.0), c(1.0)),
          # Sample 2: always population 1
          list(c(1.0), c(1.0))
        )
      )
    )
  )
  
  result <- moire::summarize_population_assignments(mock_mcmc_results_single_pop, merge_chains = TRUE)
  
  # Check structure
  expect_equal(length(result), 2)
  expect_equal(names(result), c("most_likely_population", "entropy_summary"))
  
  # Check that we have results for all samples and the single population
  expect_equal(nrow(result$most_likely_population), 2)  # 2 samples × 1 population
  expect_equal(nrow(result$entropy_summary), 2)  # 2 samples
  
  # Check that all probabilities are 1.0 (since there's only one population)
  expect_true(all(result$most_likely_population$prob_most_likely == 1.0))
  
  # Check that entropy is 0 (since there's no uncertainty with single population)
  expect_true(all(result$entropy_summary$entropy_mean == 0))
})

test_that("summarize_population_assignments respects quantile parameters", {
  # Create mock data with known entropy distribution
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
          # Sample 1: varying entropy values
          list(
            c(0.1, 0.9),  # step 1 - low entropy (certain)
            c(0.2, 0.8),  # step 2 - low entropy
            c(0.3, 0.7),  # step 3 - medium entropy
            c(0.4, 0.6),  # step 4 - medium entropy
            c(0.5, 0.5)   # step 5 - high entropy (uncertain)
          )
        )
      )
    )
  )
  
  # Test with custom quantiles
  result <- moire::summarize_population_assignments(
    mock_mcmc_results, 
    lower_quantile = 0.1, 
    upper_quantile = 0.9, 
    merge_chains = TRUE
  )
  
  # Check that entropy quantiles are calculated correctly
  entropy_result <- result$entropy_summary
  expect_true(entropy_result$entropy_lower >= 0)
  expect_true(entropy_result$entropy_upper >= entropy_result$entropy_lower)
  expect_true(entropy_result$entropy_med >= entropy_result$entropy_lower)
  expect_true(entropy_result$entropy_upper >= entropy_result$entropy_med)
  
  # Check that most likely population probabilities are calculated correctly
  most_likely_result <- result$most_likely_population
  expect_true(all(most_likely_result$prob_most_likely >= 0))
  expect_true(all(most_likely_result$prob_most_likely <= 1))
})
