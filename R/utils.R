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
load_long_form_data <- function(df, warn_uninformative = TRUE) {
  uninformative_loci <- df |>
    dplyr::ungroup() |>
    dplyr::group_by(.data$locus) |>
    dplyr::summarise(total_alleles = length(unique(.data$allele))) |>
    dplyr::filter(.data$total_alleles == 1) |>
    dplyr::pull(.data$locus)

  if (length(uninformative_loci) > 0) {
    if (warn_uninformative) {
      message("Uninformative loci with only 1 allele included. Removing...")
    }
    df <- df |>
      dplyr::filter(!(.data$locus %in% uninformative_loci))
  }

  unique_alleles <- df |>
    dplyr::group_by(.data$locus) |>
    dplyr::summarise(unique_alleles = list(sort(unique(.data$allele))))

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
    dplyr::rowwise("sample_id", "locus", "missing") |>
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
    dplyr::select("sample_id") |>
    dplyr::arrange(.data$sample_id) |>
    dplyr::distinct() |>
    dplyr::pull(.data$sample_id)

  is_missing <- matrix(missing_vec, nrow = num_loci)

  return(list(
    sample_ids = sample_ids,
    data = sample_locus_barcodes$locus_barcodes,
    loci = sample_locus_barcodes$locus,
    is_missing = is_missing,
    uninformative_loci = uninformative_loci
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
    tidyr::pivot_longer(-"sample_id",
      names_to = "locus",
      values_to = "allele"
    ) |>
    tidyr::separate_rows("allele", sep = sep) |>
    dplyr::filter(!is.na(.data$allele))
  return(load_long_form_data(df))
}


#' Plot chain swap acceptance rates
#'
#' @details Plot the swap acceptance rates for each chain.
#' The x-axis is the temperature, and the y-axis is the swap acceptance rate.
#' The dashed lines indicate the temperatures used for parallel tempering.
#'
#' @export
#'
#' @param mcmc_results list of results from `run_mcmc`
#'
#' @importFrom ggplot2 aes coord_cartesian geom_point geom_vline ggplot
#' @importFrom rlang .data
#'
#' @return list of ggplot objects
#'
plot_chain_swaps <- function(mcmc_results) {
  plots <- lapply(mcmc_results$chains, function(chain) {
    # swaps for a chain happen every 2 samples
    swaps_per_chain <- mcmc_results$args$samples_per_chain / 2
    swap_dist <- chain$swap_acceptances / swaps_per_chain
    temps <- temps <- mcmc_results$chains[[1]]$temp_gradient
    swap_idx <- (temps[1:length(temps) - 1] + temps[2:length(temps)]) / 2 # nolint: seq_linter.
    dat <- data.frame(swap_rate = swap_dist, temp = swap_idx)
    g <- ggplot(dat, aes(x = .data$temp, y = .data$swap_rate)) +
      geom_point() +
      geom_vline(data = data.frame(x = temps), aes(xintercept = .data$x), linetype = "dashed", alpha = 0.25) +
      coord_cartesian(ylim = c(0, 1))
    g
  })
  return(plots)
}
