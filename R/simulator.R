#' Dirichlet distribution
#'
#' @details Implementation of random sampling from a Dirichlet distribution
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


#' Simulate allele frequencies
#'
#' @details Simulate allele frequency vectors as a draw from a Dirichlet
#'  distribution
#'
#' @export
#'
#' @param alpha vector parameter controlling the Dirichlet distribution
#' @param num_loci total number of loci to draw
simulate_allele_frequencies <- function(alpha, num_loci) {
  dists <- rdirichlet(num_loci, alpha)
  sapply(seq_len(num_loci), function(x) {
    dists[x, ]
  })
}

#' Simulate sample COI
#'
#' @details Simulate sample COIs from a zero-truncated Poisson distribution
#'
#' @export
#'
#' @param num_samples the total number of biological samples to simulate
#' @param mean_coi mean multiplicity of infection
simulate_sample_coi <- function(num_samples, mean_coi) {
  qpois(runif(num_samples, dpois(0, mean_coi), 1), mean_coi)
}

#' Simulate sample genotype
#' @details Simulates sampling the genetics at a single locus given an allele
#'  frequency distribution and a vector of sample COIs
#'
#' @param sample_cois Numeric vector indicating the multiplicity of infection
#'  for each biological sample
#' @param locus_allele_dist Allele frequencies -- simplex parameter of a
#'  multinomial distribution
simulate_sample_genotype <- function(sample_cois, locus_allele_dist) {
  lapply(sample_cois, function(coi) {
    rmultinom(1, coi, locus_allele_dist)
  })
}

#' Simulates the observation process
#'
#' @details Takes a numeric value representing
#'  the number of strains contributing an allele and returns a binary vector
#'  indicating the presence or absence of the allele.
#'
#' @param alleles A numeric vector representing the number of strains
#'  contributing each allele
#' @param epsilon_pos expected number of false negatives
#' @param epsilon_neg expected number of false positives
simulate_observed_allele <- function(alleles, epsilon_pos, epsilon_neg) {
  positive_indices <- which(as.logical(alleles)) # True Positives
  negative_indices <- which(!as.logical(alleles)) # True Negatives

  # eps_pos_prob = epsilon_pos / length(negative_indices)
  # eps_neg_prob = epsilon_neg / length(positive_indices)
  eps_pos_prob = epsilon_pos / length(alleles)
  eps_neg_prob = epsilon_neg / length(alleles)


  alleles <- sapply(alleles, function(allele) {
    if (allele > 0) {
      rbinom(1, 1, prob = 1 - eps_neg_prob)
    } else {
      rbinom(1, 1, prob = eps_pos_prob)
      # allele
    }
  })

  # if (length(negative_indices) > 0) {
  #   for (idx in positive_indices) {
  #     if (rbinom(1, 1, prob = epsilon_pos) == 1) {
  #       fp_idx <- sample(negative_indices, 1)
  #       alleles[fp_idx] <- 1
  #     }
  #   }
  # }
  return(alleles)
}

#' Simulate observed genotypes
#'
#' @details Simulate the observation process across a list of observation
#'  vectors
#'
#' @export
#'
#' @param true_genotypes a list of numeric vectors that are input
#'  to sim_observed_allele
#' @param epsilon_pos expected number of false positives
#' @param epsilon_neg expected number of false negatives
simulate_observed_genotype <- function(true_genotypes,
                                       epsilon_pos,
                                       epsilon_neg) {
  lapply(true_genotypes, function(x) {
    simulate_observed_allele(x, epsilon_pos, epsilon_neg)
  })
}

#' Simulate data generated according to the assumed model
#'
#' @export
#'
#' @param mean_coi Mean multiplicity of infection drawn from a Poisson
#' @param locus_freq_alphas List of alpha vectors to be used to simulate
#'  from a Dirichlet distribution to generate allele frequencies.
#' @param num_samples Total number of biological samples to simulate
#' @param epsilon_pos False positive rate, expected number of false positives
#' @param epsilon_neg False negative rate, expected number of false negatives
#' @param allele_freqs List of allele frequencies to be used instead of
#'  simulating allele frequencies
#' @return Simulated data that is structured to go into the MCMC sampler
#'
simulate_data <- function(mean_coi,
                          num_samples,
                          epsilon_pos,
                          epsilon_neg,
                          locus_freq_alphas = NULL,
                          allele_freqs = NULL) {
  if(is.null(allele_freqs)) {
    allele_freqs <- list()
    for (i in 1:length(locus_freq_alphas)) {
      allele_freqs[[i]] <- simulate_allele_frequencies(locus_freq_alphas[[i]], 1)
    }
  }


  sample_cois <- simulate_sample_coi(num_samples, mean_coi)

  true_sample_genotypes <- lapply(allele_freqs, function(dist) {
    simulate_sample_genotype(sample_cois, dist)
  })

  observed_sample_genotypes <- lapply(
    true_sample_genotypes, function(locus_genotypes) {
      simulate_observed_genotype(locus_genotypes, epsilon_pos, epsilon_neg)
    }
  )

  list(
    data = observed_sample_genotypes,
    sample_ids = paste0("S", seq.int(1, num_samples)),
    loci = paste0("L", seq.int(1, length(allele_freqs))),
    allele_freqs = allele_freqs,
    sample_cois = sample_cois,
    true_genotypes = true_sample_genotypes,
    input = list(
      mean_coi = mean_coi,
      locus_freq_alphas = locus_freq_alphas,
      allele_freqs = allele_freqs,
      num_samples = num_samples,
      epsilon_pos = epsilon_pos,
      epsilon_neg = epsilon_neg
    )
  )
}
