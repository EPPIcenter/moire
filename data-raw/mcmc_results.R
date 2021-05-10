## code to prepare `mcmc_results` dataset goes here

set.seed(17327)

mean_moi <- 4
num_loci <- 25
num_biological_samples <- 40
epsilon_pos <- .03
epsilon_neg <- .1

# Generate the number of alleles at each locus
allele_counts = rep(15, num_loci)

# We'll use flat alpha vectors for our draws from the Dirichlet
locus_freq_alphas = lapply(allele_counts, function(allele) { rep(1, allele) })

simulated_data <- moire::simulate_data(mean_moi, locus_freq_alphas, num_biological_samples,
                                       epsilon_pos, epsilon_neg)

burnin = 5e3
num_samples = 1e4

mcmc_results <- moire::run_mcmc(simulated_data$data, verbose = TRUE,
                burnin = burnin, samples = num_samples,
                importance_sampling_scaling_factor = 10, thin = 10,
                eps_pos_alpha = 10, eps_pos_beta = 990, eps_neg_alpha = 100, eps_neg_beta = 900)

usethis::use_data(mcmc_results, simulated_data, overwrite = TRUE)
