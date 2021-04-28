# .run_prof_mcmc_rcpp <- function(args) {
#   .Call("start_profiler", "samples.log")
#
#   res <- .Call(`_moiR_run_mcmc`, args)
#
#   .Call("stop_profiler")
#
#   return(res)
# }

#' Sample from the target distribution using MCMC
#'
#' @export
#'
#' @param data Data to be used in MCMC, expected input is a list of lists of numeric vectors, where each list element is a collection of observations across samples at a single genetic locus.
#' @param thin Positive Integer. How often to sample from mcmc, 0 means do not thin
#' @param burnin Positive Integer. Number of MCMC samples to discard as burnin
#' @param samples Positive Integer. Number of samples to take after burnin
#' @param importance_sampling_depth Positive Integer. Min number of samples to take during importance sampling. Larger values result in longer computation time, too small of values will result in poor mixing.
#' @param importance_samping_scaling_factor Positive Ingeter. Total additional importance samples taken per increase in complexity of infection. Total number of samples = importance_sampling_depth + coi * importance_sampling_scaling_factor
#' @param verbose Logical indicating if progress is printed
#' @param max_coi Integer. Maximum allowed complexity of infection
#' @param eps_pos_0 0-1 Numeric. Initial eps_pos value
#' @param eps_pos_var 0-1 Numeric. Variance used in sampling eps_pos
#' @param eps_pos_alpha Positive Numeric. Alpha parameter in Beta distribution for eps_pos prior
#' @param eps_pos_beta Positive Numeric. Beta parameter in Beta distribution for eps_pos prior
#' @param eps_neg_0 0-1 Numeric. Initial eps_neg value
#' @param eps_neg_var 0-1 Numeric. Variance used in sampling eps_neg
#' @param eps_neg_alpha Positive Numeric. Alpha parameter in Beta distribution for eps_neg prior
#' @param eps_neg_beta Positive Numeric. Beta parameter in Beta distribution for eps_neg prior
#' @param max_eps_pos 0-1 Numeric. Maximum allowed value for eps_pos
#' @param max_eps_neg 0-1 Numeric. Maximum allowed value for eps_neg
run_mcmc <-
  function(data,
           is_missing = FALSE,
           thin = 0,
           burnin = 1e4,
           samples = 1e4,
           importance_sampling_depth = 10,
           importance_sampling_scaling_factor = 10,
           verbose = TRUE,
           max_coi = 30,
           eps_pos_0 = .01,
           eps_pos_var = .001,
           eps_pos_alpha = 1,
           eps_pos_beta = 99,
           eps_neg_0 = .1,
           eps_neg_var = .005,
           eps_neg_alpha = 10,
           eps_neg_beta = 90,
           max_eps_pos = .2,
           max_eps_neg = .2) {
    args <- as.list(environment())

    if (class(args$is_missing) == "logical" && args$is_missing == FALSE) {
      num_loci = length(args$data)
      num_biological_samples = length(args$data[[1]])
      args$is_missing = matrix(FALSE, nrow = num_loci, ncol = num_biological_samples)
    }

    res <- run_mcmc_rcpp(args)
    res
  }
