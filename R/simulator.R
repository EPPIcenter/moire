.rdirichlet <- function(n, alpha) {
  len_alpha <- length(alpha)
  d <- matrix(rgamma(len_alpha * n, alpha), ncol = len_alpha, byrow = TRUE)
  m <- d %*% rep(1, len_alpha)
  d / as.vector(m)
}

.sim_locus_allele_frequencies <- function(alpha, num_loci) {
  dists <- rdirichlet(num_loci, alpha)
  lapply(seq_len(num_loci), function(x) {
    dists[x,]
  })
}

.sim_sample_moi <- function(num_samples, mean_moi) {
  qpois(runif(num_samples, dpois(0, mean_moi), 1), mean_moi)
}

.sim_sample_genotype <- function(sample_mois, locus_allele_dist) {
  lapply(sample_mois, function(moi) {
    rmultinom(1, moi, locus_allele_dist)
  })
}

.sim_observed_allele <- function(alleles, epsilon_pos, epsilon_neg) {
  sapply(alleles, function(allele) {
    if (allele) {
      rbinom(1, 1, prob = 1 - epsilon_pos ** allele)
    } else {
      rbinom(1, 1, epsilon_pos)
    }
  })
}

.sim_observed_genotype <- function(true_genotypes, epsilon_pos, epsilon_neg) {
  lapply(true_genotypes, function(x) {
    .sim_observed_allele(x, epsilon_pos, epsilon_neg)
  })
}

#' Simulate data that is compatible with run_mcmc
#'
#' @export
#'
#' @param mean_moi Mean multiplicity of infection drawn from a poisson
#' @param num_loci Number of genetic loci to simulate
#' @param locus_freq_alpha List of alphas to be used to simulate from a dirichlet distribution to generate allele frequencies. Generates num_loci distributions for each alpha
#' @param num_samples Total number of biological samples to simulate
#' @param epsilon_pos False positive rate
#' @param epsilon_neg False negative rate
#' @return Simulated data that is structured to go into the MCMC sampler
#'
sim_data <- function(mean_moi, num_loci, locus_freq_alpha, num_samples, epsilon_pos, epsilon_neg) {
  allele_freq_dists <- c()
  for (alpha in locus_freq_alpha) {
    allele_freq_dists <- c(allele_freq_dists, .sim_locus_allele_frequencies(alpha, num_loci))
  }

  sample_mois <- .sim_sample_moi(num_samples, mean_moi)

  true_sample_genotypes <- lapply(allele_freq_dists, function(dist) {
    .sim_sample_genotype(sample_mois, dist)
  })

  observed_sample_genotypes <- lapply(true_sample_genotypes, function(locus_genotypes) {
    .sim_observed_genotype(locus_genotypes, epsilon_pos, epsilon_neg)
  })

  list(
    data =  observed_sample_genotypes,
    allele_freq_dists = allele_freq_dists,
    sample_mois = sample_mois,
    true_genotypes = true_sample_genotypes
  )
}
