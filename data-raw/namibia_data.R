namibia_data <- readxl::read_excel("data-raw/xls/elife-43510-supp1-v2.xlsx", skip = 1) |>
  dplyr::rename(sample_id = ID) |>
  tidyr::pivot_longer(cols = 6:31, names_to = "locus", values_to = "allele") |>
  tidyr::separate_rows(allele, sep = ";") |>
  tidyr::drop_na()

usethis::use_data(namibia_data, overwrite = TRUE, compress = "xz")
