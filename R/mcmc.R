.run_prof_mcmc_rcpp <- function(args) {
  .Call("start_profiler", "samples.log")

  res <- .Call(`_moiR_run_mcmc`, args)

  .Call("stop_profiler")

  return(res)
}

#' Run the MCMC
#' @export
#'
#' @param data Data to be used in MCMC
#' @param thin How often to sample from mcmc, 0 means do not thin
#' @param burnin Number of MCMC samples to discard as burnin
#' @param samples Number of samples to take after burnin
#' @param num_chains Number of concurrent chains to run (not implemented)
#' @param importance_sampling_depth How many samples to take during importance sampling
#' @param max_coi Maximum allowed complexity of infection
#' @param eps_pos_0 Initial eps_pos value
#' @param eps_pos_var Variance used in sampling eps_pos
#' @param eps_pos_alpha Alpha parameter in Beta distribution for eps_pos prior
#' @param eps_pos_beta Beta parameter in Beta distribution for eps_pos prior
#' @param eps_neg_0 Initial eps_neg value
#' @param eps_neg_var Variance used in sampling eps_neg
#' @param eps_neg_alpha Alpha parameter in Beta distribution for eps_neg prior
#' @param eps_neg_beta Beta parameter in Beta distribution for eps_neg prior
#' @param max_eps_pos Maximum allowed value for eps_pos
#' @param max_eps_neg Maximum allowed value for eps_neg
run_mcmc <- function(data, thin = 0, burnin = 1e4, samples = 1e4, num_chains = 1, importance_sampling_depth = 20, max_coi = 30, eps_pos_0 = .01, eps_pos_var=.25, eps_pos_alpha = 2, eps_pos_beta = 100, eps_neg_0 = .05, eps_neg_var=.25, eps_neg_alpha = 6, eps_neg_beta = 96, max_eps_pos = .2, max_eps_neg = .2) {
  args <- as.list(environment())
  res <- run_mcmc_rcpp(args)
  res
}
