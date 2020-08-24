## code to prepare `mcmc_results` dataset goes here

burnin = 1e3
num_samples = 1e3
mcmc_results <- run_mcmc(moire::simulated_data$data, burnin = burnin, samples = num_samples, eps_pos_0 = .05, eps_neg_0 = .1,
                eps_pos_var = .1, eps_neg_var = .1, max_coi = 20, importance_sampling_depth = 100,
                eps_pos_alpha = 1, eps_pos_beta = 1, eps_neg_alpha = 1, eps_neg_beta = 1)

usethis::use_data(mcmc_results, overwrite = TRUE)
