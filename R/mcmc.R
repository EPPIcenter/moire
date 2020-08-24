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
#' @param importance_sampling_depth Positive Integer. How many samples to take during importance sampling. Larger values result in longer computation time, too small of values will result in poor mixing.
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
           thin = 0,
           burnin = 1e4,
           samples = 1e4,
           importance_sampling_depth = 100,
           max_coi = 30,
           eps_pos_0 = .01,
           eps_pos_var = .25,
           eps_pos_alpha = 2,
           eps_pos_beta = 100,
           eps_neg_0 = .05,
           eps_neg_var = .25,
           eps_neg_alpha = 6,
           eps_neg_beta = 96,
           max_eps_pos = .2,
           max_eps_neg = .2) {
    args <- as.list(environment())
    # multiple chains not yet implemented
    args$num_chains = 1
    res <- run_mcmc_rcpp(args)
    res
  }
