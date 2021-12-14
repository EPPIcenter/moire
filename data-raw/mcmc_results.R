## ----simulate_data
set.seed(17325)

mean_moi <- 15
num_biological_samples <- 100
epsilon_pos <- .1
epsilon_neg <- 1

# Generate the number of alleles at each locus
allele_counts <- c(rep(5, 15), rep(10, 15), rep(25, 15))
# allele_counts <- c(rep(2, 100))


# We'll use flat alpha vectors for our draws from the Dirichlet
locus_freq_alphas <- lapply(allele_counts, function(allele) rep(1, allele))

simulated_data <- moire::simulate_data(
  mean_moi, locus_freq_alphas,
  num_biological_samples,
  epsilon_pos, epsilon_neg
)


## ----run_mcmc
burnin <- 1e4
num_samples <- 1e3

mcmc_results <- moire::run_mcmc(
  simulated_data$data, simulated_data$sample_ids, simulated_data$loci,
  verbose = T, burnin = burnin, samples = num_samples, thin = 1, eps_pos_0 = 1, eps_neg_0 = 1,
  eps_pos_shape = .1, eps_pos_scale = 1, eps_pos_var = 1, eps_neg_var = 1,
  eps_neg_shape = .1, eps_neg_scale = 1, allele_freq_vars = 1,
  adapt_allele_freq_vars = TRUE
)

## ----save_results
usethis::use_data(mcmc_results, simulated_data, overwrite = TRUE)
