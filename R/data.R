#' Simulated genotyping data
#'
#' A simulated dataset created using `simulate_data()`
#'
"simulated_data"

#' MCMC results from using the packaged simulated data and calling `run_mcmc()`
#'
"mcmc_results"

#' Genetic and epidemiological data from Namibia
#'
#' A dataset containing the genetic and epidemiological data from Namibia
#'
#' @format A data frame with 8 columns and 2585 rows:
#' \describe{
#'   \item{sample_id}{Sample ID}
#'   \item{HealthFacility}{Health facility}
#'   \item{HealthDistrict}{Health district}
#'   \item{Region}{Region}
#'   \item{Country}{Country}
#'   \item{locus}{Locus}
#'   \item{allele}{Allele}
#' }
#' @source \url{https://doi.org/10.7554/eLife.43510.018}
"namibia_data"

#' Allele frequencies for different regions
#'
#' A list of allele frequencies for different regions, estimated from the pf7k dataset.
#'
#' @format A list of lists, where each list element is a list of allele frequencies
#' for a specific region.
"regional_allele_frequencies"
