## code to prepare `simulated_data` dataset goes here
mean_moi <- 4
num_loci <- 10
num_biological_samples <- 100
epsilon_pos <- .03
epsilon_neg <- .15

# Generate the number of alleles at each locus
allele_counts = rnbinom(num_loci, size = 15, mu = 10) + 2

# We'll use flat alpha vectors for our draws from the Dirichlet
locus_freq_alphas = lapply(allele_counts, function(allele) { rep(1, allele) })

simulated_data <- moire::simulate_data(mean_moi, locus_freq_alphas, num_biological_samples,
                                       epsilon_pos, epsilon_neg)

usethis::use_data(simulated_data, overwrite = TRUE)
