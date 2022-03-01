## ----load_settings
set.seed(17325)

mean_moi <- 4
num_biological_samples <- 100
epsilon_pos <- .05
epsilon_neg <- .1

# Generate the number of alleles at each locus
allele_counts <- c(rep(3, 15), rep(5, 15), rep(10, 15))


# We'll use flat alpha vectors for our draws from the Dirichlet
locus_freq_alphas <- lapply(allele_counts, function(allele) rep(1, allele))

## ----simulate_data
simulated_data <- moire::simulate_data(
  mean_moi,
  num_biological_samples,
  epsilon_pos, epsilon_neg,
  locus_freq_alphas = locus_freq_alphas
)


## ----run_mcmc
burnin <- 1e4
num_samples <- 1e2

mcmc_results <- moire::run_mcmc(
  simulated_data$data, simulated_data$sample_ids, simulated_data$loci,
  verbose = T, burnin = burnin, samples_per_chain = num_samples,
  adapt_allele_freq_vars = TRUE, num_chains = 10, num_cores = 10
)

## ----save_results
usethis::use_data(mcmc_results, simulated_data, overwrite = TRUE)
