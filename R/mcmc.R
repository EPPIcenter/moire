#' Sample from the target distribution using MCMC
#'
#' @export
#'
#' @param data Data to be used in MCMC, expected input is a list
#'  of lists of numeric vectors, where each list element is a
#'  collection of observations across samples at a single genetic locus.
#' @param sample_ids Vector containing the ordered unique identifier for samples
#' @param loci Vector containing the ordered unique identifer for loci
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
#' @param complexity_limit Limit on the total number of computations before
#'  resorting to an importance sampling based approach of the integral. Total
#'  number of computations is approximately
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
#' @param mean_coi_prior_shape Positive Numeric. Shape parameter for gamma
#'  prior on mean complexity of infection.
#' @param mean_coi_prior_scale Positive Numeric. Scale parameter for gamma
#'  prior on mean complexity of infection.
#' @param mean_coi_var Positive Numeric. Variance used in sampling mean_coi
#' @param allele_freq_var Positive Numeric. Variance used in sampling allele
#'  frequencies
run_mcmc <-
  function(data,
           sample_ids,
           loci,
           is_missing = FALSE,
           thin = 1,
           burnin = 1e4,
           samples = 1e4,
           complexity_limit = 2050,
           importance_sampling_depth = 300,
           importance_sampling_scaling_factor = 100,
           verbose = TRUE,
           eps_pos_0 = .01,
           eps_pos_var = .001,
           eps_pos_alpha = 1,
           eps_pos_beta = 99,
           eps_neg_0 = .05,
           eps_neg_var = .005,
           eps_neg_alpha = 5,
           eps_neg_beta = 95,
           max_eps_pos = .5,
           max_eps_neg = .5,
           mean_coi_prior_shape = 1.5,
           mean_coi_prior_scale = .5,
           mean_coi_var = 1,
           allele_freq_var = .1) {
    args <- as.list(environment())

    ## if is_missing == FALSE, then generate a default FALSE matrix
    if (class(args$is_missing) == "logical" && args$is_missing == FALSE) {
      num_loci <- length(args$data)
      num_biological_samples <- length(args$data[[1]])
      args$is_missing <- matrix(
        FALSE,
        nrow = num_loci,
        ncol = num_biological_samples
      )
    }

    total_alleles <- lapply(data, function(x) {
      return(length(x[[1]]))
    })
    if (any(total_alleles < 2)) {
      stop("Loci with less than 2 alleles present, remove these loci")
    }

    res <- run_mcmc_rcpp(args)
    res$args <- args
    res$total_samples <- args$samples / args$thin
    res
  }
