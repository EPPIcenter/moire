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
#' @param warn_uninformative boolean whether or not to print message when
#'  removing uninformative loci
#'
#' @importFrom rlang .data
load_long_form_data <- function(df, warn_uninformative=TRUE) {
  uninformative_loci <- df |>
    dplyr::group_by(.data$locus) |>
    dplyr::summarise(total_alleles = length(unique(.data$allele))) |>
    dplyr::filter(total_alleles == 1) |>
    dplyr::pull(.data$locus)

  if (length(uninformative_loci) > 0 ) {
    if (warn_uninformative) {
      message("Uninformative loci with only 1 allele included. Removing...")
      message(paste0("Loci: ", paste(uninformative_loci, collapse = ", ")))
    }
    df <- df |>
      dplyr::filter(!(.data$locus %in% uninformative_loci))
  }

  unique_alleles <- df |>
    dplyr::group_by(.data$locus) |>
    dplyr::summarise(unique_alleles = list(sort(unique(.data$allele))))

  uninformative_loci <- unique_alleles |>
    dplyr::mutate(total_alleles = length(.data$unique_alleles)) |>
    dplyr::filter(total_alleles == 1)



  num_loci <- df |>
    dplyr::pull(.data$locus) |>
    dplyr::n_distinct()

  sample_locus_barcodes <- df |>
    dplyr::group_by(.data$sample_id, .data$locus) |>
    dplyr::summarize(alleles_grp = list(.data$allele), .groups = "drop") |>
    dplyr::mutate(missing = FALSE) |>
    tidyr::complete(.data$sample_id, .data$locus,
      fill = list(alleles_grp = list(c(NULL)), missing = TRUE)
    )

  missing_vec <- sample_locus_barcodes |>
    dplyr::arrange(.data$sample_id, .data$locus) |>
    dplyr::pull(.data$missing)

  sample_locus_barcodes <- sample_locus_barcodes |>
    dplyr::left_join(unique_alleles, by = "locus") |>
    dplyr::rowwise(.data$sample_id, .data$locus, .data$missing) |>
    dplyr::summarise(
      barcode = list(sapply(
        unique_alleles,
        function(x) {
          as.integer(x %in% .data$alleles_grp)
        }
      )), .groups = "drop"
    ) |>
    dplyr::group_by(.data$locus) |>
    dplyr::arrange(.data$sample_id, by_group = TRUE) |>
    dplyr::summarise(locus_barcodes = list(.data$barcode))


  sample_ids <- df |>
    dplyr::select(.data$sample_id) |>
    dplyr::arrange(.data$sample_id) |>
    dplyr::distinct() |>
    dplyr::pull(.data$sample_id)

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
#' @param warn_uninformative boolean whether or not to print message when
#'  removing uninformative loci
#'
#' @importFrom rlang .data
load_delimited_data <- function(data, sep = ";", warn_uninformative = TRUE) {
  df <- data |>
    tidyr::pivot_longer(-.data$sample_id,
      names_to = "locus",
      values_to = "allele"
    ) |>
    tidyr::separate_rows(.data$allele, sep = sep) |>
    dplyr::filter(!is.na(.data$allele))
  return(load_long_form_data(df))
}
