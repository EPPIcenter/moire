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
  allele_freq_dists <- sim_locus_allele_frequencies(locus_freq_alpha, num_loci)
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

run_prof_mcmc_rcpp <- function(args) {
  .Call("start_profiler", "samples.log")

  res <- .Call(`_moiR_run_mcmc`, args)

  .Call("stop_profiler")

  return(res)
}

run_mcmc <- function(data, thin = 0, burnin = 1e4, samples = 1e4, num_chains = 1, importance_sampling_depth = 20,
                     max_coi = 30, max_coi_delta = 5, eps_pos_0 = .01, eps_neg_0 = .05,
                     max_eps_pos = .2, max_eps_neg = .2, alpha = 200) {
  args <- as.list(environment())
  res <- run_mcmc_rcpp(args)
  res
}


mean_moi = 7
num_loci = 10
locus_freq_alpha = rep(1, 15)
num_samples = 50
epsilon_pos = .01
epsilon_neg = .10

res <- sim_data(mean_moi, num_loci, locus_freq_alpha, num_samples, epsilon_pos, epsilon_neg)
# ret <- run_mcmc(res$data)

Sys.getpid()
# run_mcmc(res)

