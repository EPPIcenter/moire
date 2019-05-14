run_prof_mcmc_rcpp <- function(args) {
  .Call("start_profiler", "samples.log")

  res <- .Call(`_moiR_run_mcmc`, args)

  .Call("stop_profiler")

  return(res)
}

run_mcmc <- function(data, thin = 0, burnin = 1e4, samples = 1e4, num_chains = 1, importance_sampling_depth = 20, max_coi = 30, eps_pos_0 = .01, eps_pos_var=.25, eps_pos_alpha = 2, eps_pos_beta = 100, eps_neg_0 = .05, eps_neg_var=.25, eps_neg_alpha = 6, eps_neg_beta = 96, max_eps_pos = .2, max_eps_neg = .2, alpha = 200) {
  args <- as.list(environment())
  res <- run_mcmc_rcpp(args)
  res
}