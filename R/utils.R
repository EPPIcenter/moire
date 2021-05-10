#' Load long form data
#'
#' @details Long form data is a data frame with
#'  3 columns: `sample_id`, `locus`, `allele`. Returned data contains
#'  vectors `sample_ids` and `loci` that are ordered as the results
#'  will be ordered from running the MCMC algorithm.
#'
#' @export
#'
#' @param df data frame with 3 columns: `sample_id`, `locus`, `allele`.
#' Each row is a single observation of an allele at a particular
#' locus for a given sample.
#'
load_long_form_data <- function(df) {
  unique_alleles <- df %>%
    dplyr::group_by(locus) %>%
    dplyr::summarise(unique_alleles = list(sort(unique(allele))))

  num_loci <- df %>%
    dplyr::pull(locus) %>%
    dplyr::n_distinct()

  sample_locus_barcodes <- df %>%
    dplyr::group_by(sample_id, locus) %>%
    dplyr::summarize(alleles_grp = list(allele), .groups = "drop") %>%
    dplyr::mutate(missing = FALSE) %>%
    tidyr::complete(sample_id, locus,
      fill = list(alleles_grp = list(c(NULL)), missing = TRUE)
    )

  missing_vec <- sample_locus_barcodes %>%
    dplyr::arrange(sample_id, locus) %>%
    dplyr::pull(missing)

  sample_locus_barcodes <- sample_locus_barcodes %>%
    dplyr::left_join(unique_alleles, by = "locus") %>%
    dplyr::rowwise(sample_id, locus, missing) %>%
    dplyr::summarise(
      barcode = list(sapply(
        unique_alleles,
        function(x) {
          as.integer(x %in% alleles_grp)
        }
      )), .groups = "drop"
    ) %>%
    dplyr::group_by(locus) %>%
    dplyr::arrange(sample_id, by_group = TRUE) %>%
    dplyr::summarise(locus_barcodes = list(barcode))


  sample_ids <- df %>%
    dplyr::select(sample_id) %>%
    dplyr::arrange(sample_id) %>%
    dplyr::distinct() %>%
    dplyr::pull(sample_id)

  is_missing <- matrix(missing_vec, nrow = num_loci)

  return(list(
    sample_ids = sample_ids,
    data = sample_locus_barcodes$locus_barcodes,
    loci = sample_locus_barcodes$locus,
    is_missing = is_missing
  ))
}

#' Load delimited data
#'
#' @details Load `data.frame` with a `sample_id` column and the remaining
#'  columns are `loci`. Each cell contains a separator delimited string
#'  representing the observed alleles at that locus for that sample.
#'  Returned data contains vectors `sample_ids` and `loci` that are ordered
#'  as the results will be ordered from running the MCMC algorithm.
#'
#' @export
#'
#' @param data data.frame containing the described data
#' @param sep string used to separate alleles

load_delimited_data <- function(data, sep = ";") {
  df <- data %>%
    tidyr::pivot_longer(-sample_id,
      names_to = "locus",
      values_to = "allele"
    ) %>%
    tidyr::separate_rows(allele, sep = sep) %>%
    dplyr::filter(!is.na(allele))
  return(load_long_form_data(df))
}
