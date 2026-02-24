test_that("loading long form data works", {
  data <- data.frame(
    sample_id = c("S1", "S1", "S1", "S2", "S2", "S3", "S3", "S4", "S4", "S5"),
    locus = c("A", "A", "B", "A", "B", "A", "A", "B", "B", "A"),
    allele = c("1", "2", "1", "1", "2", "3", "1", "1", "2", "2")
  ) |> dplyr::slice_sample(prop = 1)

  res <- moire::load_long_form_data(data)

  expect_equal(res$sample_ids, c("S1", "S2", "S3", "S4", "S5"))
  expect_equal(res$loci, c("A", "B"))
  # Check that we have the expected number of samples and loci
  expect_equal(length(res$sample_ids), 5)
  expect_equal(length(res$loci), 2)
  # Check that data structure is correct
  expect_equal(length(res$data), 2)  # 2 loci
  expect_equal(length(res$data[[1]]), 5)  # 5 samples for locus A
  expect_equal(length(res$data[[2]]), 5)  # 5 samples for locus B
})

test_that("loading delimited data works", {
  data <- data.frame(
    sample_id = c("S1", "S2", "S3", "S4", "S5"),
    A = c("1;2", "1", "1;3", "2", "1;2;3"),
    B = c("1", "1", "2;3", "1;2", "2"),
    C = c("1;2", "1", "1", "2", "1;2")
  )

  res <- moire::load_delimited_data(data)

  expect_equal(res$sample_ids, c("S1", "S2", "S3", "S4", "S5"))
  expect_equal(res$loci, c("A", "B", "C"))
  # Check that we have the expected number of samples and loci
  expect_equal(length(res$sample_ids), 5)
  expect_equal(length(res$loci), 3)
  # Check that data structure is correct
  expect_equal(length(res$data), 3)  # 3 loci
  expect_equal(length(res$data[[1]]), 5)  # 5 samples for each locus
  expect_equal(length(res$data[[2]]), 5)
  expect_equal(length(res$data[[3]]), 5)
})


test_that("convert_to_long_form works with basic data", {
  # Create test data
  original_df <- data.frame(
    sample_id = c("S1", "S1", "S1", "S2", "S2", "S3", "S3"),
    locus = c("L1", "L1", "L2", "L1", "L2", "L1", "L2"),
    allele = c("A", "B", "A", "A", "B", "B", "A")
  )
  
  # Load data
  loaded_data <- moire::load_long_form_data(original_df)
  
  # Convert back to long form
  converted_df <- moire::convert_to_long_form(loaded_data)
  
  # Check structure
  expect_equal(ncol(converted_df), 3)
  expect_equal(colnames(converted_df), c("sample_id", "locus", "allele"))
  
  # Check that we get the expected number of observations
  expect_equal(nrow(converted_df), 7)
  
  # Check that all sample IDs are present
  expect_equal(sort(unique(converted_df$sample_id)), c("S1", "S2", "S3"))
  
  # Check that all loci are present
  expect_equal(sort(unique(converted_df$locus)), c("L1", "L2"))
  
  # Check that alleles are properly reconstructed
  expect_true(all(grepl("^Allele_", converted_df$allele)))
})


test_that("convert_to_long_form works with simulated data", {
  # Create locus frequency alphas for 2 loci with 2 alleles each
  locus_freq_alphas <- lapply(1:2, function(x) rep(1, 2))
  
  # Simulate data
  sim_data <- moire::simulate_data(
    mean_coi = 2,
    num_samples = 3,
    epsilon_pos = 0.1,
    epsilon_neg = 0.1,
    locus_freq_alphas = locus_freq_alphas,
    missingness = 0.05
  )
  
  # Convert to long form
  long_form <- moire::convert_to_long_form(sim_data)
  
  # Check structure
  expect_equal(ncol(long_form), 3)
  expect_equal(colnames(long_form), c("sample_id", "locus", "allele"))
  
  # Check that we have the expected samples and loci
  expect_equal(sort(unique(long_form$sample_id)), c("S1", "S2", "S3"))
  expect_equal(sort(unique(long_form$locus)), c("L1", "L2"))
  
  # Check that alleles are properly formatted
  expect_true(all(grepl("^Allele_", long_form$allele)))
  
  # Check that we have some data
  expect_true(nrow(long_form) > 0)
})

test_that("convert_to_long_form preserves data integrity", {
  # Create test data
  original_df <- data.frame(
    sample_id = c("S1", "S1", "S2", "S2"),
    locus = c("L1", "L1", "L1", "L1"),
    allele = c("A", "B", "A", "C")
  )
  
  # Load and convert
  loaded_data <- moire::load_long_form_data(original_df)
  converted_df <- moire::convert_to_long_form(loaded_data)
  
  # Check that the number of observations matches the original
  expect_equal(nrow(converted_df), nrow(original_df))
  
  # Check that all original sample IDs are preserved
  expect_equal(sort(unique(converted_df$sample_id)), sort(unique(original_df$sample_id)))
  
  # Check that all original loci are preserved
  expect_equal(sort(unique(converted_df$locus)), sort(unique(original_df$locus)))
  
  # Check that we have the right number of observations per sample-locus combination
  s1_l1_count <- sum(converted_df$sample_id == "S1" & converted_df$locus == "L1")
  s2_l1_count <- sum(converted_df$sample_id == "S2" & converted_df$locus == "L1")
  
  expect_equal(s1_l1_count, 2)  # S1 has alleles A and B at L1
  expect_equal(s2_l1_count, 2)  # S2 has alleles A and C at L1
})

test_that("convert_to_long_form works with delimited data", {
  # Create delimited data
  delimited_df <- data.frame(
    sample_id = c("S1", "S2", "S3"),
    L1 = c("A;B", "A", "B;C"),
    L2 = c("A", "B", "A;B")
  )
  
  # Load and convert
  loaded_data <- moire::load_delimited_data(delimited_df)
  converted_df <- moire::convert_to_long_form(loaded_data)
  
  # Check structure
  expect_equal(ncol(converted_df), 3)
  expect_equal(colnames(converted_df), c("sample_id", "locus", "allele"))
  
  # Check that we have the expected samples and loci
  expect_equal(sort(unique(converted_df$sample_id)), c("S1", "S2", "S3"))
  expect_equal(sort(unique(converted_df$locus)), c("L1", "L2"))
  
  # Check that alleles are properly formatted
  expect_true(all(grepl("^Allele_", converted_df$allele)))
})
