#' Calculate naive COI offset
#'
#' @details Estimates the complexity of infection using a naive approach
#'  that chooses the n'th highest number of observed alleles.
#'
#' @export
#'
#' @param data List of lists of numeric vectors, where each list
#' element is a collection of observations across samples at a
#' single genetic locus.
#' @param offset Numeric offset -- n'th highest number of observed alleles
calculate_naive_coi_offset <- function(data, offset) {
  num_alleles_by_locus <- lapply(data, function(locus) {
    lapply(locus, sum)
  })
  num_alleles_by_sample <- list()
  for (i in 1:length(num_alleles_by_locus)) {
    for (j in 1:length(num_alleles_by_locus[[i]])) {
      if (i == 1) {
        num_alleles_by_sample[[j]] <- list()
      }
      num_alleles_by_sample[[j]][[i]] <- num_alleles_by_locus[[i]][[j]]
    }
  }
  sapply(num_alleles_by_sample, function(sample) {
    max(1, sort(unlist(sample), decreasing = T)[offset])
  })
}

#' Calculate naive COI
#'
#' @details Estimates the complexity of infection using a naive approach
#' that chooses the highest number of observed alleles.
#'
#' @export
#'
#' @param data List of lists of numeric vectors, where each list
#' element is a collection of observations across samples at a
#' single genetic locus.
calculate_naive_coi <- function(data) {
  num_alleles_by_locus <- lapply(data, function(locus) {
    lapply(locus, sum)
  })
  num_alleles_by_sample <- list()
  for (i in 1:length(num_alleles_by_locus)) {
    for (j in 1:length(num_alleles_by_locus[[i]])) {
      if (i == 1) {
        num_alleles_by_sample[[j]] <- list()
      }
      num_alleles_by_sample[[j]][[i]] <- num_alleles_by_locus[[i]][[j]]
    }
  }
  sapply(num_alleles_by_sample, function(sample) {
    max(unlist(sample))
  })
}

#' Calculate naive allele frequencies
#'
#' @details Estimate naive allele frequencies from the empirical distribution
#'  of alleles
#'
#' @import purrr
#'
#' @export
#'
#' @param data List of lists of numeric vectors, where each list element
#' is a collection of observations across samples at a single genetic locus
calculate_naive_allele_frequencies <- function(data) {
  allele_freqs <- lapply(data, function(locus) {
    allele_counts <- purrr::reduce(locus, function(f1, f2) f1 + f2)
    return(allele_counts / sum(allele_counts))
  })
}

#' Calculate the expected heterozygosity from allele frequencies
#'
#' @export
#'
#' @param allele_freqs Simplex of allele frequencies
calculate_he <- function(allele_freqs) {
  return(1 - sum(allele_freqs**2))
}

#' Summarize COI
#'
#' @details Summarize complexity of infection results from MCMC. Returns
#'  a dataframe that contains summaries of the posterior
#'  distribution of COI for each biological sample, as well as naive
#'  estimates of COI.
#'
#'  @importFrom stats quantile
#'
#' @export
#'
#' @param mcmc_results Result of calling run_mcmc
#' @param lower_quantile The lower quantile of the posterior
#'  distribution to return
#' @param upper_quantile The upper quantile of the posterior
#'  distribution to return
#' @param naive_offset Offset used in calculate_naive_coi_offset
#' @param merge_chains boolean indicating that all chain results should be merged
summarize_coi <- function(mcmc_results, lower_quantile = .025,
                          upper_quantile = .975, naive_offset = 2, merge_chains = TRUE) {
  naive_coi <- calculate_naive_coi(mcmc_results$args$data)
  offset_naive_coi <- calculate_naive_coi_offset(mcmc_results$args$data, 2)

  if (merge_chains) {
    cois <- lapply(1:length(mcmc_results$args$sample_ids), function(x) c())
    for (idx in 1:length(mcmc_results$chains)) {
      chain <- mcmc_results$chains[[idx]]
      for (s in 1:length(chain$coi)) {
        cois[[s]] <- c(cois[[s]], chain$coi[[s]])
      }
      post_coi_lower <- sapply(cois, function(x) {
        quantile(x, lower_quantile)
      })
      post_coi_med <- sapply(cois, function(x) {
        quantile(x, .5)
      })
      post_coi_upper <- sapply(cois, function(x) {
        quantile(x, upper_quantile)
      })
      post_coi_mean <- sapply(cois, mean)
      return(data.frame(
        sample_id = mcmc_results$args$sample_ids,
        post_coi_lower, post_coi_med, post_coi_upper, post_coi_mean,
        naive_coi, offset_naive_coi
      ))
    }
  } else {
    chain_cois <- lapply(1:length(mcmc_results$chains), function(idx) {
      chain <- mcmc_results$chains[[idx]]
      cois <- chain$coi
      post_coi_lower <- sapply(cois, function(x) {
        quantile(x, lower_quantile)
      })
      post_coi_med <- sapply(cois, function(x) {
        quantile(x, .5)
      })
      post_coi_upper <- sapply(cois, function(x) {
        quantile(x, upper_quantile)
      })
      post_coi_mean <- sapply(cois, mean)
      return(data.frame(
        sample_id = mcmc_results$args$sample_ids,
        post_coi_lower, post_coi_med, post_coi_upper, post_coi_mean, chain = idx,
        naive_coi, offset_naive_coi
      ))
    })
    coi_data <- do.call(rbind, chain_cois)
    return(coi_data)
  }
}

#' Summarize epsilon_neg
#'
#' @details Summarize epsilon negative results from MCMC. Returns
#'  a dataframe that contains summaries of the posterior
#'  distribution of epsilon negative for each biological sample.
#'
#' @importFrom stats quantile
#'
#' @export
#'
#' @param mcmc_results Result of calling run_mcmc()
#' @param lower_quantile The lower quantile of the posterior
#'  distribution to return
#' @param upper_quantile The upper quantile of the posterior
#'  distribution to return
#' @param merge_chains boolean indicating that all chain results should be merged
summarize_epsilon_neg <- function(mcmc_results, lower_quantile = .025, upper_quantile = .975, merge_chains = TRUE) {
  if (merge_chains) {
    epsilon_neg <- lapply(1:length(mcmc_results$args$sample_ids), function(x) c())
    for (idx in 1:length(mcmc_results$chains)) {
      chain <- mcmc_results$chains[[idx]]
      for (s in 1:length(chain$eps_neg)) {
        epsilon_neg[[s]] <- c(epsilon_neg[[s]], chain$eps_neg[[s]])
      }
    }
    post_eps_neg_lower <- sapply(epsilon_neg, function(x) {
      quantile(x, lower_quantile)
    })
    post_eps_neg_med <- sapply(epsilon_neg, function(x) {
      quantile(x, .5)
    })
    post_eps_neg_upper <- sapply(epsilon_neg, function(x) {
      quantile(x, upper_quantile)
    })
    post_eps_neg_mean <- sapply(epsilon_neg, mean)

    return(data.frame(
      sample_id = mcmc_results$args$sample_ids,
      post_eps_neg_lower, post_eps_neg_med, post_eps_neg_upper, post_eps_neg_mean
    ))
  } else {
    chain_eps_neg <- lapply(1:length(mcmc_results$chains), function(idx) {
      epsilon_neg <- mcmc_results$chains[[idx]]$eps_neg
      post_eps_neg_lower <- sapply(epsilon_neg, function(x) {
        quantile(x, lower_quantile)
      })
      post_eps_neg_med <- sapply(epsilon_neg, function(x) {
        quantile(x, .5)
      })
      post_eps_neg_upper <- sapply(epsilon_neg, function(x) {
        quantile(x, upper_quantile)
      })
      post_eps_neg_mean <- sapply(epsilon_neg, mean)

      return(data.frame(
        sample_id = mcmc_results$args$sample_ids,
        post_eps_neg_lower, post_eps_neg_med, post_eps_neg_upper, post_eps_neg_mean,
        chain = idx
      ))
    })
    eps_neg_data <- do.call(rbind, chain_eps_neg)
    return(eps_neg_data)
  }
}

#' Summarize epsilon_pos
#'
#' @details Summarize epsilon positive results from MCMC. Returns
#'  a dataframe that contains summaries of the posterior
#'  distribution of epsilon positive for each biological sample.
#'
#' @importFrom stats quantile
#' @export
#'
#' @param mcmc_results Result of calling run_mcmc()
#' @param lower_quantile The lower quantile of the posterior
#'  distribution to return
#' @param upper_quantile The upper quantile of the posterior
#'  distribution to return
#' @param merge_chains boolean indicating that all chain results should be merged
summarize_epsilon_pos <- function(mcmc_results, lower_quantile = .025, upper_quantile = .975, merge_chains = TRUE) {
  if (merge_chains) {
    epsilon_pos <- lapply(1:length(mcmc_results$args$sample_ids), function(x) c())
    for (chain in mcmc_results$chains) {
      for (s in 1:length(chain$eps_pos)) {
        epsilon_pos[[s]] <- c(epsilon_pos[[s]], chain$eps_pos[[s]])
      }
    }
    post_eps_pos_lower <- sapply(epsilon_pos, function(x) {
      quantile(x, lower_quantile)
    })
    post_eps_pos_med <- sapply(epsilon_pos, function(x) {
      quantile(x, .5)
    })
    post_eps_pos_upper <- sapply(epsilon_pos, function(x) {
      quantile(x, upper_quantile)
    })
    post_eps_pos_mean <- sapply(epsilon_pos, mean)

    return(data.frame(
      sample_id = mcmc_results$args$sample_ids,
      post_eps_pos_lower, post_eps_pos_med, post_eps_pos_upper, post_eps_pos_mean
    ))
  } else {
    chain_eps_pos <- lapply(1:length(mcmc_results$chains), function(idx) {
      epsilon_pos <- mcmc_results$chains[[idx]]$eps_pos
      post_eps_pos_lower <- sapply(epsilon_pos, function(x) {
        quantile(x, lower_quantile)
      })
      post_eps_pos_med <- sapply(epsilon_pos, function(x) {
        quantile(x, .5)
      })
      post_eps_pos_upper <- sapply(epsilon_pos, function(x) {
        quantile(x, upper_quantile)
      })
      post_eps_pos_mean <- sapply(epsilon_pos, mean)

      return(data.frame(
        sample_id = mcmc_results$args$sample_ids,
        post_eps_pos_lower, post_eps_pos_med, post_eps_pos_upper, post_eps_pos_mean,
        chain = idx
      ))
    })
    eps_pos_data <- do.call(rbind, chain_eps_pos)
    return(eps_pos_data)
  }
}

#' Summarize Function of Allele Frequencies
#'
#' @details General function to summarize the posterior distribution of
#' functions of the sampled allele frequencies
#'
#' @importFrom stats quantile
#'
#' @export
#'
#' @param mcmc_results Result of calling run_mcmc()
#' @param fn Function that takes as input a simplex to apply to each
#'  allele frequency vector
#' @param lower_quantile The lower quantile of the posterior distribution
#'  to return
#' @param upper_quantile The upper quantile of the posterior distribution
#'  to return
#' @param merge_chains boolean indicating that all chain results should be merged
summarize_allele_freq_fn <- function(mcmc_results, fn,
                                     lower_quantile = .025,
                                     upper_quantile = .975, merge_chains = TRUE) {
  if (merge_chains) {
    post_allele_freqs <- lapply(1:length(mcmc_results$args$loci), function(x) c())

    for (chain in mcmc_results$chains) {
      for (l in 1:length(chain$allele_freqs)) {
        post_allele_freqs[[l]] <- c(post_allele_freqs[[l]], chain$allele_freqs[[l]])
      }
    }

    post_statistic <- lapply(post_allele_freqs, function(locus_posterior) {
      sapply(locus_posterior, function(allele_freq_sample) fn(allele_freq_sample))
    })

    res <- data.frame(
      locus = mcmc_results$args$loci,
      post_stat_lower = sapply(
        post_statistic,
        function(x) quantile(x, lower_quantile)
      ),
      post_stat_med = sapply(
        post_statistic,
        function(x) quantile(x, .5)
      ),
      post_stat_upper = sapply(
        post_statistic,
        function(x) quantile(x, upper_quantile)
      ),
      post_stat_mean = sapply(post_statistic, mean)
    )

    return(res)
  } else {
    res <- lapply(1:length(mcmc_results$chains), function(idx) {
      post_allele_freqs <- mcmc_results$chains[[idx]]$allele_freqs

      post_statistic <- lapply(post_allele_freqs, function(locus_posterior) {
        sapply(locus_posterior, function(allele_freq_sample) fn(allele_freq_sample))
      })

      chain_res <- data.frame(
        locus = mcmc_results$args$loci,
        post_stat_lower = sapply(
          post_statistic,
          function(x) quantile(x, lower_quantile)
        ),
        post_stat_med = sapply(
          post_statistic,
          function(x) quantile(x, .5)
        ),
        post_stat_upper = sapply(
          post_statistic,
          function(x) quantile(x, upper_quantile)
        ),
        post_stat_mean = sapply(post_statistic, mean),
        chain = idx
      )

      return(chain_res)
    })
    return(do.call(rbind, res))
  }
}

#' Summarize locus heterozygosity
#'
#' @details Summarize locus heterozygosity from the posterior distribution
#'  of sampled allele frequencies.
#'
#' @export
#'
#' @param mcmc_results Result of calling run_mcmc()
#' @param lower_quantile The lower quantile of the posterior distribution
#'  to return
#' @param upper_quantile The upper quantile of the posterior distribution
#' to return
#' @param merge_chains Merge the results of multiple chains into a single
#' summary
summarize_he <- function(mcmc_results,
                         lower_quantile = .025,
                         upper_quantile = .975, merge_chains = TRUE) {
  res <- summarize_allele_freq_fn(
    mcmc_results,
    fn = calculate_he,
    lower_quantile = lower_quantile,
    upper_quantile = upper_quantile,
    merge_chains = merge_chains
  )
  return(res)
}

#' Summarize allele frequencies
#'
#' @details Summarize individual allele frequencies from the posterior
#'  distribution of sampled allele frequencies
#'
#' @importFrom stats quantile
#' @import purrr
#'
#' @export
#'
#' @param mcmc_results Result of calling run_mcmc()
#' @param lower_quantile The lower quantile of the posterior distribution
#'  to return
#' @param upper_quantile The upper quantile of the posterior distribution
#' to return
#' @param merge_chains boolean indicating that all chain results should be merged
summarize_allele_freqs <- function(mcmc_results,
                                   lower_quantile = .025,
                                   upper_quantile = .975,
                                   merge_chains = TRUE) {
  locus_alleles <- do.call(
    rbind,
    purrr::imap(
      mcmc_results$args$data,
      ~ data.frame(locus = mcmc_results$args$loci[.y], allele = names(.x[[1]]))
    )
  )

  if (merge_chains) {
    allele_freq_matrices <- lapply(1:length(mcmc_results$args$loci), function(x) c())
    total_samples <- 0
    for (chain in mcmc_results$chains) {
      total_samples <- total_samples + length(chain$allele_freqs[[1]])
      for (l in 1:length(chain$allele_freqs)) {
        locus <- chain$allele_freqs[[l]]
        allele_freq_matrices[[l]] <- c(allele_freq_matrices[[l]], unlist(locus))
      }
    }
    allele_freq_matrices <- lapply(allele_freq_matrices, function(x) matrix(x, ncol = total_samples))
    res <- lapply(allele_freq_matrices, function(allele_freq_matrix) {
      post_allele_freqs_lower <- apply(
        allele_freq_matrix, 1, function(x) quantile(x, lower_quantile)
      )
      post_allele_freqs_med <- apply(
        allele_freq_matrix, 1, function(x) quantile(x, .5)
      )
      post_allele_freqs_upper <- apply(
        allele_freq_matrix, 1, function(x) quantile(x, upper_quantile)
      )
      post_allele_freqs_mean <- apply(
        allele_freq_matrix, 1, mean
      )

      return(data.frame(
        post_allele_freqs_lower = post_allele_freqs_lower,
        post_allele_freqs_med = post_allele_freqs_med,
        post_allele_freqs_upper = post_allele_freqs_upper,
        post_allele_freqs_mean = post_allele_freqs_mean
      ))
    })
    return(cbind(do.call(rbind, res), locus_alleles))
  } else {
    res <- list()
    for (idx in 1:length(mcmc_results$chains)) {
      chain <- mcmc_results$chains[[idx]]
      chain_res <- lapply(
        chain$allele_freqs,
        function(locus) {
          num_alleles <- length(locus[[1]])
          allele_freq_matrix <- matrix(unlist(locus), nrow = num_alleles)

          post_allele_freqs_lower <- apply(
            allele_freq_matrix, 1, function(x) quantile(x, lower_quantile)
          )
          post_allele_freqs_med <- apply(
            allele_freq_matrix, 1, function(x) quantile(x, .5)
          )
          post_allele_freqs_upper <- apply(
            allele_freq_matrix, 1, function(x) quantile(x, upper_quantile)
          )
          post_allele_freqs_mean <- apply(
            allele_freq_matrix, 1, mean
          )

          data.frame(
            post_allele_freqs_lower = post_allele_freqs_lower,
            post_allele_freqs_med = post_allele_freqs_med,
            post_allele_freqs_upper = post_allele_freqs_upper,
            post_allele_freqs_mean = post_allele_freqs_mean
          )
        }
      )
      chain_res <- cbind(do.call(rbind, chain_res), locus_alleles)
      chain_res$chain <- idx
      res[[idx]] <- chain_res
    }
    return(do.call(rbind, res))
  }
}



#' Summarize relatedness
#'
#' @details Summarize relatedness results from MCMC. Returns
#'  a dataframe that contains summaries of the posterior
#'  distribution of relatedness for each biological sample.
#'
#' @importFrom stats quantile
#' @export
#'
#' @param mcmc_results Result of calling run_mcmc()
#' @param lower_quantile The lower quantile of the posterior
#'  distribution to return
#' @param upper_quantile The upper quantile of the posterior
#'  distribution to return
#' @param merge_chains boolean indicating that all chain results should be merged
summarize_relatedness <- function(mcmc_results, lower_quantile = .025, upper_quantile = .975, merge_chains = TRUE) {
  if (merge_chains) {
    relatedness <- lapply(1:length(mcmc_results$args$sample_ids), function(x) c())
    for (chain in mcmc_results$chains) {
      for (s in 1:length(chain$relatedness)) {
        relatedness[[s]] <- c(relatedness[[s]], chain$relatedness[[s]])
      }
    }
    post_relatedness_lower <- sapply(relatedness, function(x) {
      quantile(x, lower_quantile)
    })
    post_relatedness_med <- sapply(relatedness, function(x) {
      quantile(x, .5)
    })
    post_relatedness_upper <- sapply(relatedness, function(x) {
      quantile(x, upper_quantile)
    })
    post_relatedness_mean <- sapply(relatedness, mean)

    return(data.frame(
      sample_id = mcmc_results$args$sample_ids,
      post_relatedness_lower, post_relatedness_med, post_relatedness_upper, post_relatedness_mean
    ))
  } else {
    chain_relatedness <- lapply(1:length(mcmc_results$chains), function(idx) {
      relatedness <- mcmc_results$chains[[idx]]$relatedness
      post_relatedness_lower <- sapply(relatedness, function(x) {
        quantile(x, lower_quantile)
      })
      post_relatedness_med <- sapply(relatedness, function(x) {
        quantile(x, .5)
      })
      post_relatedness_upper <- sapply(relatedness, function(x) {
        quantile(x, upper_quantile)
      })
      post_relatedness_mean <- sapply(relatedness, mean)

      return(data.frame(
        sample_id = mcmc_results$args$sample_ids,
        post_relatedness_lower, post_relatedness_med, post_relatedness_upper, post_relatedness_mean,
        chain = idx
      ))
    })
    relatedness_data <- do.call(rbind, chain_relatedness)
    return(relatedness_data)
  }
}

#' Summarize effective COI
#'
#' @details Summarize effective COI from MCMC. Returns
#'  a dataframe that contains summaries of the posterior
#'  distribution of effective COI for each biological sample.
#'
#' @importFrom stats quantile
#' @import purrr
#'
#' @export
#'
#' @param mcmc_results Result of calling run_mcmc()
#' @param lower_quantile The lower quantile of the posterior
#'  distribution to return
#' @param upper_quantile The upper quantile of the posterior
#'  distribution to return
#' @param merge_chains boolean indicating that all chain results should be merged
summarize_effective_coi <- function(mcmc_results, lower_quantile = .025, upper_quantile = .975, merge_chains = TRUE) {
  if (merge_chains) {
    effective_coi <- lapply(1:length(mcmc_results$args$sample_ids), function(x) c())
    for (chain in mcmc_results$chains) {
      for (s in 1:length(chain$relatedness)) {
        effective_coi[[s]] <- c(effective_coi[[s]], (1 - chain$relatedness[[s]]) * (chain$coi[[s]] - 1) + 1)
      }
    }
    post_effective_coi_lower <- sapply(effective_coi, function(x) {
      quantile(x, lower_quantile)
    })
    post_effective_coi_med <- sapply(effective_coi, function(x) {
      quantile(x, .5)
    })
    post_effective_coi_upper <- sapply(effective_coi, function(x) {
      quantile(x, upper_quantile)
    })
    post_effective_coi_mean <- sapply(effective_coi, mean)

    return(data.frame(
      sample_id = mcmc_results$args$sample_ids,
      post_effective_coi_lower, post_effective_coi_med, post_effective_coi_upper, post_effective_coi_mean
    ))
  } else {
    chain_effective_coi <- lapply(1:length(mcmc_results$chains), function(idx) {
      chain <- mcmc_results$chains[[idx]]
      r <- chain$relatedness
      coi <- chain$coi

      effective_coi <- purrr::map2(r, coi, ~ (1 - .x) * (.y - 1) + 1)

      # effective_coi <- (1 - mcmc_results$chains[[idx]]$relatedness) * mcmc_results$chains[[idx]]$coi
      post_effective_coi_lower <- sapply(effective_coi, function(x) {
        quantile(x, lower_quantile)
      })
      post_effective_coi_med <- sapply(effective_coi, function(x) {
        quantile(x, .5)
      })
      post_effective_coi_upper <- sapply(effective_coi, function(x) {
        quantile(x, upper_quantile)
      })
      post_effective_coi_mean <- sapply(effective_coi, mean)

      return(data.frame(
        sample_id = mcmc_results$args$sample_ids,
        post_effective_coi_lower, post_effective_coi_med, post_effective_coi_upper, post_effective_coi_mean,
        chain = idx
      ))
    })
    relatedness_data <- do.call(rbind, chain_effective_coi)
    return(relatedness_data)
  }
}
