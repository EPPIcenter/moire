#' Dirichlet distribution
#'
#' @details Implementation of random sampling from a Dirichlet distribution
#'
#' @importFrom stats rgamma
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
#' @importFrom stats qpois runif dpois
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
#'   frequency distribution and a vector of sample COIs
#'
#' @importFrom stats rmultinom rbinom quantile
#'
#' @export
#'
#' @param sample_cois Numeric vector indicating the multiplicity of infection
#'   for each biological sample
#' @param locus_allele_dist Allele frequencies -- simplex parameter of a
#'   multinomial distribution
#' @param internal_relatedness numeric 0-1 indicating the probability for a
#'   strain's allele to come from an existing lineage within host
#' @export
simulate_sample_genotype <- function(sample_cois, locus_allele_dist, internal_relatedness) {
  purrr::map2(sample_cois, internal_relatedness, function(coi, r) {
    related_strains <- rbinom(1, coi - 1, r)
    genotypes <- rmultinom(1, coi - related_strains, locus_allele_dist)
    t(genotypes)
  })
}

#' Simulates the observation process
#'
#' @details Takes a numeric value representing
#'  the number of strains contributing an allele and returns a binary vector
#'  indicating the presence or absence of the allele.
#'
#' @importFrom stats rbinom
#' @export
#'
#' @param alleles A numeric vector representing the number of strains
#'  contributing each allele
#' @param epsilon_pos expected number of false negatives
#' @param epsilon_neg expected number of false positives
#' @param missingness probability that the data is missing
simulate_observed_allele <- function(alleles, epsilon_pos, epsilon_neg, missingness) {
  # scale eps to the number of alleles so that given a fixed COI, there is a fixed
  # number of expected false positives or negatives across loci of varying
  # diversity.
  eps_pos_prob <- epsilon_pos / length(alleles)
  eps_neg_prob <- epsilon_neg / length(alleles)

  if (runif(1) < missingness) {
    obs_alleles <- rep(0, length(alleles))
  } else {
    obs_alleles <- sapply(alleles, function(allele) {
      if (allele > 0) {
        rbinom(1, 1, prob = 1 - eps_neg_prob)
      } else {
        rbinom(1, 1, prob = eps_pos_prob)
      }
    })
  }

  names(obs_alleles) <- colnames(alleles)

  return(obs_alleles)
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
#' @param missingness probability of data being missing
simulate_observed_genotype <- function(true_genotypes,
                                       epsilon_pos,
                                       epsilon_neg,
                                       missingness) {
  lapply(true_genotypes, function(x) {
    simulate_observed_allele(x, epsilon_pos, epsilon_neg, missingness)
  })
}

#' Simulate data generated according to the assumed model
#'
#' @export
#'
#' @importFrom stats rbeta
#'
#' @param mean_coi Mean multiplicity of infection drawn from a Poisson
#' @param locus_freq_alphas List of alpha vectors to be used to simulate
#'  from a Dirichlet distribution to generate allele frequencies.
#' @param num_samples Total number of biological samples to simulate
#' @param epsilon_pos False positive rate, expected number of false positives
#' @param epsilon_neg False negative rate, expected number of false negatives
#' @param sample_cois List of sample COIs to be used instead of simulating
#' @param allele_freqs List of allele frequencies to be used instead of
#'  simulating allele frequencies
#' @param internal_relatedness_alpha alpha parameter of beta distribution controlling
#'  the random relatedness draws for each sample
#' @param internal_relatedness_beta beta parameter of beta distribution controlling
#'  the random relatedness draws for each sample
#' @param internal_relatedness List of internal relatedness values to be used 
#'  instead of simulating
#' @param missingness probability of data being missing
#' @return Simulated data that is structured to go into the MCMC sampler
#'
simulate_data <- function(mean_coi = NULL,
                          num_samples,
                          epsilon_pos,
                          epsilon_neg,
                          sample_cois = NULL,
                          locus_freq_alphas = NULL,
                          allele_freqs = NULL,
                          internal_relatedness_alpha = 0,
                          internal_relatedness_beta = 1,
                          internal_relatedness = NULL,
                          missingness = 0) {
  if (is.null(allele_freqs)) {
    allele_freqs <- list()
    allele_freq_names <- paste0("L", 1:length(locus_freq_alphas))
    for (i in 1:length(locus_freq_alphas)) {
      total_alleles <- length(locus_freq_alphas[[i]])
      allele_names <- 1:total_alleles
      allele_freqs[[allele_freq_names[i]]] <- simulate_allele_frequencies(locus_freq_alphas[[i]], 1)[, 1]
      names(allele_freqs[[i]]) <- paste(allele_freq_names[i], allele_names, sep = "_")
    }
    names(allele_freqs) <- allele_freq_names
  }

  if (is.null(sample_cois)) {
    sample_cois <- simulate_sample_coi(num_samples, mean_coi)
  }

  if (is.null(internal_relatedness)) {
    internal_relatedness <- rbeta(num_samples, internal_relatedness_alpha, internal_relatedness_beta)
  }

  true_sample_genotypes <- lapply(allele_freqs, function(dist) {
    simulate_sample_genotype(sample_cois, dist, internal_relatedness)
  })

  observed_sample_genotypes <- lapply(
    true_sample_genotypes, function(locus_genotypes) {
      simulate_observed_genotype(locus_genotypes, epsilon_pos, epsilon_neg, missingness)
    }
  )

  suppressWarnings({
    is_missing <- !t(sapply(observed_sample_genotypes, function(loc) {
      sapply(loc, any)
    }))
  })

  list(
    data = observed_sample_genotypes,
    sample_ids = paste0("S", seq.int(1, num_samples)),
    loci = paste0("L", seq.int(1, length(allele_freqs))),
    is_missing = is_missing,
    allele_freqs = allele_freqs,
    sample_cois = sample_cois,
    sample_relatedness = internal_relatedness,
    true_genotypes = true_sample_genotypes,
    input = as.list(environment())
  )
}
