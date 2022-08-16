test_that("loading long form data works", {
  data <- data.frame(
    sample_id = c("S1", "S1", "S1", "S2", "S2", "S3", "S3"),
    locus = c("A", "A", "B", "A", "B", "A", "A"),
    allele = c("1", "2", "1", "1", "2", "3", "1")
  ) |> dplyr::slice_sample(prop = 1)

  res <- moire::load_long_form_data(data)

  expect_equal(res$sample_ids, c("S1", "S2", "S3"))
  expect_equal(res$loci, c("A", "B"))
  expect_equal(res$data[[1]][[1]], c(1, 1, 0), ignore_attr = TRUE)
  expect_equal(res$data[[1]][[2]], c(1, 0, 0), ignore_attr = TRUE)
  expect_equal(res$data[[1]][[3]], c(1, 0, 1), ignore_attr = TRUE)
  expect_equal(sum(res$is_missing), 1)
})

test_that("loading delimited data works", {
  data <- data.frame(
    sample_id = c("S1", "S2", "S3"),
    A = c("1;2", "1", "1;3"),
    B = c(NA, "1", "2;3")
  )

  res <- moire::load_delimited_data(data)

  expect_equal(res$sample_ids, c("S1", "S2", "S3"))
  expect_equal(res$loci, c("A", "B"))
  expect_equal(res$data[[1]][[1]], c(1, 1, 0), ignore_attr = TRUE)
  expect_equal(res$data[[1]][[2]], c(1, 0, 0), ignore_attr = TRUE)
  expect_equal(res$data[[1]][[3]], c(1, 0, 1), ignore_attr = TRUE)
  expect_equal(sum(res$is_missing), 1)
})
