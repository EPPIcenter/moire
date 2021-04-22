#' Implementation of random sampling from a Dirichlet distribution
#'
#' @export
#'
#' @param n total number of draws
#' @param alpha vector controlling the concentration of simplex
rdirichlet <- function(n, alpha) {
  len_alpha <- length(alpha)
  d <- matrix(rgamma(len_alpha * n, alpha), ncol = len_alpha, byrow = TRUE)
  m <- d %*% rep(1, len_alpha)
  d / as.vector(m)
}

#' Simulate allele frequency vectors as a draw from a Dirichlet distribution
#'
#' @export
#'
#' @param alpha vector parameter controlling the Dirichlet distribution
#' @param num_loci total number of loci to draw
simulate_locus_allele_frequencies <- function(alpha, num_loci) {
  dists <- rdirichlet(num_loci, alpha)
  lapply(seq_len(num_loci), function(x) {
    dists[x, ]
  })
}

#' Simulate sample MOIs from a zero-truncated Poisson distribution
#'
#' @export
#'
#' @param num_samples the total number of biological samples to simulate
#' @param mean_moi mean multiplicity of infection
simulate_sample_moi <- function(num_samples, mean_moi) {
  qpois(runif(num_samples, dpois(0, mean_moi), 1), mean_moi)
}

#' Simulates sampling the genetics at a single locus given an allele frequency distribution and a vector of sample MOIs
#'
#' @param sample_mois Numeric vector indicating the multiplicity of infection for each biological sample
#' @param locus_allele_dist Allele frequencies -- simplex parameter of a multinomial distribution
simulate_sample_genotype <- function(sample_mois, locus_allele_dist) {
  lapply(sample_mois, function(moi) {
    rmultinom(1, moi, locus_allele_dist)
  })
}

#' Simulates the observation process. Takes a numeric value representing the number of strains contributing an allele and returns a binary vector indicating the presence or absence of the allele.
#'
#' @param alleles A numeric vector representing the number of strains contributing each allele
#' @param epsilon_pos false positive rate
#' @param epsilon_neg false negative rate
simulate_observed_allele <- function(alleles, epsilon_pos, epsilon_neg) {
  sapply(alleles, function(allele) {
    if (allele > 0) {
      rbinom(1, 1, prob = 1 - epsilon_neg)
    } else {
      rbinom(1, 1, epsilon_pos)
    }
  })
}


#' Simulate the observation process across a list of observation vectors
#'
#' @xport
#'
#' @param true_genotypes a list of numeric vectors that are input to sim_observed_allele
#' @param epsilon_pos false positive rate
#' @param epsilon_neg false negative rate
simulate_observed_genotype <- function(true_genotypes, epsilon_pos, epsilon_neg) {
  lapply(true_genotypes, function(x) {
    simulate_observed_allele(x, epsilon_pos, epsilon_neg)
  })
}

#' Simulate data that is compatible with run_mcmc and generated according to the assumed model
#'
#' @export
#'
#' @param mean_moi Mean multiplicity of infection drawn from a Poisson
#' @param locus_freq_alphas List of alpha vectors to be used to simulate from a Dirichlet distribution to generate allele frequencies.
#' @param num_samples Total number of biological samples to simulate
#' @param epsilon_pos False positive rate, between 0 and 1
#' @param epsilon_neg False negative rate, between 0 and 1
#' @return Simulated data that is structured to go into the MCMC sampler
#'
simulate_data <- function(mean_moi, locus_freq_alphas, num_samples, epsilon_pos, epsilon_neg) {
  allele_freq_dists <- c()
  for (alpha in locus_freq_alphas) {
    allele_freq_dists <- c(allele_freq_dists, simulate_locus_allele_frequencies(alpha, 1))
  }

  sample_mois <- simulate_sample_moi(num_samples, mean_moi)

  true_sample_genotypes <- lapply(allele_freq_dists, function(dist) {
    simulate_sample_genotype(sample_mois, dist)
  })

  observed_sample_genotypes <- lapply(true_sample_genotypes, function(locus_genotypes) {
    simulate_observed_genotype(locus_genotypes, epsilon_pos, epsilon_neg)
  })

  list(
    data = observed_sample_genotypes,
    allele_freq_dists = allele_freq_dists,
    sample_mois = sample_mois,
    true_genotypes = true_sample_genotypes,
    input = list(
      mean_moi = mean_moi,
      locus_freq_alphas = locus_freq_alphas,
      num_samples = num_samples,
      epsilon_pos = epsilon_pos,
      epsilon_neg = epsilon_neg
    )
  )
}
