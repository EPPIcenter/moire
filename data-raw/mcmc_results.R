## ----simulate_data
set.seed(17325)

mean_moi <- 5
num_biological_samples <- 100
epsilon_pos <- .05
epsilon_neg <- .1

# Generate the number of alleles at each locus
allele_counts <- c(rep(5, 15), rep(10, 15), rep(25, 15))


# We'll use flat alpha vectors for our draws from the Dirichlet
locus_freq_alphas <- lapply(allele_counts, function(allele) rep(1, allele))

simulated_data <- moire::simulate_data(
  mean_moi, locus_freq_alphas,
  num_biological_samples,
  epsilon_pos, epsilon_neg
)


## ----run_mcmc
burnin <- 1e5
num_samples <- 1e4

mcmc_results <- moire::run_mcmc(
  simulated_data$data, simulated_data$sample_ids, simulated_data$loci,
  verbose = T, burnin = burnin, samples = num_samples, thin = 1,
  eps_pos_alpha = 1, eps_pos_beta = 99, complexity_limit = 1,
  eps_neg_alpha = 1, eps_neg_beta = 99, allele_freq_vars = 1,
  adapt_allele_freq_vars = TRUE
)

## ----save_results
usethis::use_data(mcmc_results, simulated_data, overwrite = TRUE)
