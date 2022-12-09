## ----load_settings
set.seed(17325)

mean_moi <- 5
num_biological_samples <- 100
epsilon_pos <- .01
epsilon_neg <- .03

# Generate the number of alleles at each locus
allele_counts <- c(rep(5, 15), rep(10, 15), rep(25, 15), rep(50, 15))

# We'll use flat alpha vectors for our draws from the Dirichlet
locus_freq_alphas <- lapply(allele_counts, function(allele) rep(1, allele))

simulated_data <- moire::simulate_data(
  mean_moi,
  num_biological_samples,
  epsilon_pos, epsilon_neg,
  locus_freq_alphas = locus_freq_alphas,
  internal_relatedness_alpha = .1,
  internal_relatedness_beta = 1
)

## ----run_mcmc
burnin <- 1e4
num_samples <- 1e4
pt_chains <- seq(1, .5, length.out = 20)

mcmc_results <- moire::run_mcmc(
  simulated_data,
  verbose = T, burnin = burnin, samples_per_chain = num_samples,
  pt_chains = pt_chains, pt_num_threads = length(pt_chains), thin = 10
)

## ----save_results
usethis::use_data(mcmc_results, simulated_data, overwrite = TRUE, compress = "xz")
