# Test for initial allele frequencies functionality
# ================================================

test_that("prepare_initial_allele_frequencies works correctly", {
  # Create simple test data using the correct structure
  data <- list(
    sample_ids = c("sample1", "sample2", "sample3", "sample4"),
    loci = c("locus1", "locus2"),
    data = list(
      # Locus 1: 2 alleles
      list(
        c(1, 0),  # sample1: allele 1 present
        c(0, 1),  # sample2: allele 2 present  
        c(1, 1),  # sample3: both alleles present
        c(1, 0)   # sample4: allele 1 present
      ),
      # Locus 2: 3 alleles
      list(
        c(1, 0, 0),  # sample1: allele 1 present
        c(0, 1, 0),  # sample2: allele 2 present
        c(0, 0, 1),  # sample3: allele 3 present
        c(1, 1, 0)   # sample4: alleles 1 and 2 present
      )
    )
  )
  
  # Test with 2 populations
  population_assignments <- c(1, 1, 2, 2)  # samples 1,2 in pop1; samples 3,4 in pop2
  num_populations <- 2
  
  result <- prepare_initial_allele_frequencies(
    data = data,
    population_assignments = population_assignments,
    num_populations = num_populations,
    pseudocount = 1
  )
  
  # Check structure
  expect_equal(length(result), num_populations)
  expect_equal(length(result[[1]]), 2)  # 2 loci
  expect_equal(length(result[[2]]), 2)  # 2 loci
  
  # Check that frequencies sum to 1
  for (pop_idx in 1:num_populations) {
    for (locus_idx in 1:2) {
      expect_equal(sum(result[[pop_idx]][[locus_idx]]), 1.0, tolerance = 1e-6)
    }
  }
  
  # Check specific values for population 1 (samples 1,2)
  # Locus 1: sample1 has [1,0], sample2 has [0,1], plus pseudocount [1,1]
  # Total: [1+0+1, 0+1+1] = [2,2], normalized: [2/4, 2/4] = [0.5, 0.5]
  expect_equal(result[[1]][[1]], c(0.5, 0.5), tolerance = 1e-6)
  
  # Locus 2: sample1 has [1,0,0], sample2 has [0,1,0], plus pseudocount [1,1,1]
  # Total: [1+0+1, 0+1+1, 0+0+1] = [2,2,1], normalized: [2/5, 2/5, 1/5] = [0.4, 0.4, 0.2]
  expect_equal(result[[1]][[2]], c(0.4, 0.4, 0.2), tolerance = 1e-6)
  
  # Check specific values for population 2 (samples 3,4)
  # Locus 1: sample3 has [1,1], sample4 has [1,0], plus pseudocount [1,1]
  # Total: [1+1+1, 1+0+1] = [3,2], normalized: [3/5, 2/5] = [0.6, 0.4]
  expect_equal(result[[2]][[1]], c(0.6, 0.4), tolerance = 1e-6)
  
  # Locus 2: sample3 has [0,0,1], sample4 has [1,1,0], plus pseudocount [1,1,1]
  # Total: [0+1+1, 0+1+1, 1+0+1] = [2,2,2], normalized: [2/6, 2/6, 2/6] = [1/3, 1/3, 1/3]
  expect_equal(result[[2]][[2]], c(1/3, 1/3, 1/3), tolerance = 1e-6)
})

test_that("prepare_initial_allele_frequencies handles edge cases", {
  # Test with empty population
  data <- list(
    sample_ids = c("sample1", "sample2"),
    loci = c("locus1"),
    data = list(
      list(
        c(1, 0),  # sample1
        c(0, 1)   # sample2
      )
    )
  )
  
  # Assign all samples to population 1, leaving population 2 empty
  population_assignments <- c(1, 1)
  num_populations <- 2
  
  result <- prepare_initial_allele_frequencies(
    data = data,
    population_assignments = population_assignments,
    num_populations = num_populations,
    pseudocount = 1
  )
  
  # Population 1 should have calculated frequencies
  expect_equal(sum(result[[1]][[1]]), 1.0, tolerance = 1e-6)
  
  # Population 2 should have uniform frequencies (fallback)
  expect_equal(result[[2]][[1]], c(0.5, 0.5), tolerance = 1e-6)
})

test_that("prepare_initial_allele_frequencies validates inputs", {
  data <- list(
    sample_ids = c("sample1", "sample2"),
    loci = c("locus1"),
    data = list(
      list(
        c(1, 0),
        c(0, 1)
      )
    )
  )
  
  # Test wrong length of population_assignments
  expect_error(
    prepare_initial_allele_frequencies(
      data = data,
      population_assignments = c(1, 2, 3),  # Wrong length
      num_populations = 2,
      pseudocount = 1
    ),
    "population_assignments must have length equal to number of samples"
  )
  
  # Test invalid population assignments
  expect_error(
    prepare_initial_allele_frequencies(
      data = data,
      population_assignments = c(1, 3),  # 3 > num_populations (2)
      num_populations = 2,
      pseudocount = 1
    ),
    "population_assignments must be integers from 1 to num_populations"
  )
  
  expect_error(
    prepare_initial_allele_frequencies(
      data = data,
      population_assignments = c(0, 1),  # 0 < 1
      num_populations = 2,
      pseudocount = 1
    ),
    "population_assignments must be integers from 1 to num_populations"
  )
})

