#' Simulated genotyping data
#'
#' A simulated dataset created using `simulate_data()`
#'
#' @format A list as returned by [simulate_data()], with 100 samples at 30
#'  loci (15 with 5 alleles, 15 with 10) and a mean COI of 5: observed `data`,
#'  `sample_ids`, `loci`, `is_missing`, and the simulated truth
#'  (`allele_freqs`, `sample_cois`, `sample_relatedness`, `true_genotypes`).
"simulated_data"

#' MCMC results from using the packaged simulated data and calling `run_mcmc()`
#'
#' @format A list as returned by [run_mcmc()] on [simulated_data]: a single
#'  chain with 1000 burnin and 1000 sampling iterations, using 80 parallel
#'  tempering replicas with adaptive temperatures. The `mcmc_demo` vignette
#'  shows the call.
"mcmc_results"

#' Genetic and epidemiological data from Namibia
#'
#' A dataset containing the genetic and epidemiological data from Namibia
#'
#' @format A data frame with 7 columns and 97214 rows:
#' \describe{
#'   \item{sample_id}{Sample ID}
#'   \item{HealthFacility}{Health facility}
#'   \item{HealthDistrict}{Health district}
#'   \item{Region}{Region}
#'   \item{Country}{Country}
#'   \item{locus}{Genetic locus}
#'   \item{allele}{Allele observed}
#' }
#' @source \doi{10.7554/eLife.43510.018}
"namibia_data"

#' Allele frequencies for different regions
#'
#' A list of allele frequencies for different regions, estimated from the pf7k dataset.
#'
#' @format A list of lists, where each list element is a list of allele frequencies
#' for a specific region.
"regional_allele_frequencies"
