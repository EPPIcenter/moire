sim_locus_allele_frequencies <- function(alpha, num_loci) {
  dists <- MCMCpack::rdirichlet(num_loci, alpha)
  lapply(seq_len(num_loci), function(x) {
    dists[x,]
  })
}

sim_sample_moi <- function(num_samples, mean_moi) {
  qpois(runif(num_samples, dpois(0, mean_moi), 1), mean_moi)
}

sim_sample_genotype <- function(sample_mois, locus_allele_dist) {
  lapply(sample_mois, function(moi) {
    rmultinom(1, moi, locus_allele_dist)
  })
}

sim_observed_allele <- function(alleles, epsilon_pos, epsilon_neg) {
  sapply(alleles, function(allele) {
    if (allele) {
      rbinom(1, 1, prob = 1 - epsilon_pos ** allele)
    } else {
      rbinom(1, 1, epsilon_pos)
    }
  })
}

sim_observed_genotype <- function(true_genotypes, epsilon_pos, epsilon_neg) {
  lapply(true_genotypes, function(x) {
    sim_observed_allele(x, epsilon_pos, epsilon_neg)
  })
}

sim_data <- function(mean_moi, num_loci, locus_freq_alpha, num_samples, epsilon_pos, epsilon_neg) {
  # allele_freq_dists <- sim_locus_allele_frequencies(locus_freq_alpha, num_loci)
  allele_freq_dists <- c()
  for (alpha in locus_freq_alpha) {
    allele_freq_dists <- c(allele_freq_dists, sim_locus_allele_frequencies(alpha, num_loci))
  }

  sample_mois <- sim_sample_moi(num_samples, mean_moi)

  true_sample_genotypes <- lapply(allele_freq_dists, function(dist) {
    sim_sample_genotype(sample_mois, dist)
  })

  observed_sample_genotypes <- lapply(true_sample_genotypes, function(locus_genotypes) {
    sim_observed_genotype(locus_genotypes, epsilon_pos, epsilon_neg)
  })

  list(
    data =  observed_sample_genotypes,
    allele_freq_dists = allele_freq_dists,
    sample_mois = sample_mois,
    true_genotypes = true_sample_genotypes
  )
}
