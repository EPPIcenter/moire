## ----load_settings
set.seed(17325)

mean_moi <- 5
num_biological_samples <- 100
epsilon_pos <- .01
epsilon_neg <- .1
internal_relatedness_alpha <- .1
internal_relatedness_beta <- 1

# Generate 15 loci with 5 alleles each, and 15 loci with 10 alleles each
allele_counts <- c(rep(5, 15), rep(10, 15))

# We'll use flat alpha vectors for our draws from the Dirichlet
locus_freq_alphas <- lapply(allele_counts, function(allele) rep(1, allele))

simulated_data <- moire::simulate_data(
  mean_moi,
  num_biological_samples,
  epsilon_pos, epsilon_neg,
  locus_freq_alphas = locus_freq_alphas,
  internal_relatedness_alpha = internal_relatedness_alpha,
  internal_relatedness_beta = internal_relatedness_beta,
)

## ----run_mcmc
burnin <- 1e3
num_samples <- 1e3
pt_chains <- seq(1, 0, length.out = 80)
num_cores <- parallelly::availableCores() - 1 # number of threads to use for parallel processing

mcmc_results <- moire::run_mcmc(
  simulated_data,
  verbose = TRUE, burnin = burnin, samples_per_chain = num_samples,
  pt_chains = pt_chains, num_cores = num_cores,
  adapt_temp = TRUE
)

## ----save_results
usethis::use_data(mcmc_results, simulated_data, overwrite = TRUE, compress = "xz")
