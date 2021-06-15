## code to prepare `mcmc_results` dataset goes here

set.seed(17327)

mean_moi <- 4
num_biological_samples <- 100
epsilon_pos <- .01
epsilon_neg <- .03

# Generate the number of alleles at each locus
allele_counts <- c(rep(5, 10), rep(10, 10), rep(15, 10))

# We'll use flat alpha vectors for our draws from the Dirichlet
locus_freq_alphas <- lapply(allele_counts, function(allele) rep(1, allele))

simulated_data <- moire::simulate_data(
  mean_moi, locus_freq_alphas,
  num_biological_samples,
  epsilon_pos, epsilon_neg
)

burnin <- 1e3
num_samples <- 1e3

mcmc_results <- moire::run_mcmc(
  simulated_data$data, simulated_data$sample_ids, simulated_data$loci,
  verbose = T, burnin = burnin, samples = num_samples, thin = 1,
  eps_pos_alpha = 10, eps_pos_beta = 990,
  eps_neg_alpha = 30, eps_neg_beta = 970, allele_freq_var = .1
)

usethis::use_data(mcmc_results, simulated_data, overwrite = TRUE)
