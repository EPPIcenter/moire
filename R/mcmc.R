#' Sample from the target distribution using MCMC
#'
#' @export
#'
#' @param data Data to be used in MCMC, expected input is a list
#'  of lists of numeric vectors, where each list element is a
#'  collection of observations across samples at a single genetic locus.
#' @param sample_ids Vector containing the ordered unique identifier for samples
#' @param loci Vector containing the ordered unique identifer for loci
#' @param mean_coi Int prior mean complexity of infection for poisson prior
#' @param is_missing Boolean matrix indicating whether the observation
#'  should be treated as missing data and ignored. Number of rows equals the
#'  number of loci, number of columns equals the number samples. Alternatively,
#'  the user may pass in FALSE if no data should be considered missing.
#' @param thin Positive Integer. How often to sample from mcmc,
#'  1 means do not thin
#' @param burnin Positive Integer. Number of MCMC samples to
#'  discard as burnin
#' @param samples Positive Integer. Number of samples to take
#'  after burnin
#' @param importance_sampling_depth Positive Integer. Min number
#'  of samples to take during importance sampling. Larger values
#'  result in longer computation time, too small of values will
#'  result in poor mixing.
#' @param importance_sampling_scaling_factor Positive Ingeter.
#'  Total additional importance samples taken per increase in
#'  complexity of infection. Total number of
#'  samples = importance_sampling_depth + coi * importance_sampling_scaling_
#'  factor
#' @param verbose Logical indicating if progress is printed
#' @param eps_pos_0 0-1 Numeric. Initial eps_pos value
#' @param eps_pos_var 0-1 Numeric. Variance used in sampling eps_pos
#' @param eps_pos_alpha Positive Numeric. Alpha parameter in
#'  Beta distribution for eps_pos prior
#' @param eps_pos_beta Positive Numeric. Beta parameter in
#'  Beta distribution for eps_pos prior
#' @param eps_neg_0 0-1 Numeric. Initial eps_neg value
#' @param eps_neg_var 0-1 Numeric. Variance used in sampling eps_neg
#' @param eps_neg_alpha Positive Numeric. Alpha parameter in
#'  Beta distribution for eps_neg prior
#' @param eps_neg_beta Positive Numeric. Beta parameter in
#'  Beta distribution for eps_neg prior
#' @param max_eps_pos 0-1 Numeric. Maximum allowed value for eps_pos
#' @param max_eps_neg 0-1 Numeric. Maximum allowed value for eps_neg
#' @param allele_freq_var Positive Numeric. Variance used in sampling allele
#'  frequencies
run_mcmc <-
  function(data,
           sample_ids,
           loci,
           mean_coi,
           is_missing = FALSE,
           thin = 1,
           burnin = 1e4,
           samples = 1e4,
           importance_sampling_depth = 10,
           importance_sampling_scaling_factor = 10,
           verbose = TRUE,
           eps_pos_0 = .01,
           eps_pos_var = .001,
           eps_pos_alpha = 1,
           eps_pos_beta = 99,
           eps_neg_0 = .1,
           eps_neg_var = .005,
           eps_neg_alpha = 10,
           eps_neg_beta = 90,
           max_eps_pos = .5,
           max_eps_neg = .5,
           allele_freq_var = 1) {
    args <- as.list(environment())

    if (class(args$is_missing) == "logical" && args$is_missing == FALSE) {
      num_loci <- length(args$data)
      num_biological_samples <- length(args$data[[1]])
      args$is_missing <- matrix(
        FALSE,
        nrow = num_loci,
        ncol = num_biological_samples
      )
    }

    res <- run_mcmc_rcpp(args)
    res$args <- args
    res$total_samples <- args$samples / args$thin
    res
  }
