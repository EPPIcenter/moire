#' Calculate the geometric median of the posterior distribution of allele
#' frequencies
#'
#' @details Returns the geometric median of the posterior distribution, defined
#' as the point minimizing the L2 distance from each sampled point.
#'
#' @import purrr
#' @importFrom stats dist
#'
#' @export
#'
#' @param mcmc_results Result of calling run_mcmc()
#'
#' @param merge_chains boolean indicating that all chain results should be merged
calculate_med_allele_freqs <- function(mcmc_results, merge_chains = TRUE) {
  if (merge_chains) {
    chains <- mcmc_results$chains
    post_af <- purrr::transpose(purrr::map(chains, ~ .x$allele_freqs))
    medians <- purrr::map(post_af, function(loc) {
      samples <- purrr::flatten(loc)
      mat <- matrix(unlist(samples), ncol = length(samples[[1]]))
      d <- dist(t(mat))
      samples[[which.min(rowSums(as.matrix(d)))]]
    })
    names(medians) <- mcmc_results$args$data$loci
  } else {
    chains <- mcmc_results$chains
    medians <- lapply(chains, function(chain) {
      post_af <- chain$allele_freqs
      res <- purrr::map(post_af, function(samples) {
        mat <- matrix(unlist(samples), ncol = length(samples[[1]]))
        d <- dist(t(mat))
        samples[[which.min(rowSums(as.matrix(d)))]]
      })
      names(res) <- mcmc_results$args$data$loci
      return(res)
    })
  }
  return(medians)
}


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
  for (i in seq_along(num_alleles_by_locus)) {
    for (j in seq_along(num_alleles_by_locus[[i]])) {
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
  for (i in seq_along(num_alleles_by_locus)) {
    for (j in seq_along(num_alleles_by_locus[[i]])) {
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
  return(allele_freqs)
}

# Internal helper: merge chains (or not), then apply quantile/mean to list of vectors per entity.
# chains: list of chain objects; get_vectors: function(chain) -> list of numeric vectors;
# ids: vector for id column; id_name, value_prefix: column naming; na.rm: for quantile/mean.
summarize_posterior_vectors <- function(chains, get_vectors, ids, id_name = "sample_id",
                                        value_prefix = "post_", lower_quantile = .025,
                                        upper_quantile = .975, merge_chains = TRUE, na.rm = FALSE) {
  if (merge_chains) {
    vectors <- lapply(seq_along(ids), function(x) c())
    for (chain in chains) {
      vecs <- get_vectors(chain)
      for (i in seq_along(vecs)) {
        vectors[[i]] <- c(vectors[[i]], vecs[[i]])
      }
    }
    lower_vec <- sapply(vectors, function(x) quantile(x, lower_quantile, na.rm = na.rm))
    med_vec <- sapply(vectors, function(x) quantile(x, .5, na.rm = na.rm))
    upper_vec <- sapply(vectors, function(x) quantile(x, upper_quantile, na.rm = na.rm))
    mean_vec <- sapply(vectors, function(x) mean(x, na.rm = na.rm))
    out <- data.frame(ids, lower_vec, med_vec, upper_vec, mean_vec, stringsAsFactors = FALSE)
    names(out) <- c(id_name, paste0(value_prefix, c("_lower", "_med", "_upper", "_mean")))
    return(out)
  } else {
    res <- lapply(seq_along(chains), function(idx) {
      vecs <- get_vectors(chains[[idx]])
      lower_vec <- sapply(vecs, function(x) quantile(x, lower_quantile, na.rm = na.rm))
      med_vec <- sapply(vecs, function(x) quantile(x, .5, na.rm = na.rm))
      upper_vec <- sapply(vecs, function(x) quantile(x, upper_quantile, na.rm = na.rm))
      mean_vec <- sapply(vecs, function(x) mean(x, na.rm = na.rm))
      out <- data.frame(ids, lower_vec, med_vec, upper_vec, mean_vec, chain = idx, stringsAsFactors = FALSE)
      names(out) <- c(id_name, paste0(value_prefix, c("_lower", "_med", "_upper", "_mean")), "chain")
      return(out)
    })
    return(do.call(rbind, res))
  }
}

#' Calculate the expected heterozygosity from allele frequencies
#'
#' @export
#'
#' @param allele_freqs Simplex of allele frequencies
calculate_he <- function(allele_freqs) {
  if (any(is.na(allele_freqs))) {
    allele_freqs <- replace(allele_freqs, which(is.na(allele_freqs)), 0)
    warning("NA values detected in allele frequency vector.This may indicate a problem with the MCMC chain or there are loci with no diversity. NA values will be replaced with 0.")
  }
  return(1 - sum(allele_freqs**2))
}

#' Calcuate the expected number of distinct alleles under a given multiplicity of infection (n)
#'
#' @export
#'
#' @param allele_freqs Simplex of allele frequencies
#' @param n multiplicity of infection
calculate_eda <- function(allele_freqs, n) {
  sum(1 - (1 - allele_freqs)^n)
}


#' Summarize COI
#'
#' @details Summarize complexity of infection results from MCMC. Returns
#'  a dataframe that contains summaries of the posterior
#'  distribution of COI for each biological sample, as well as naive
#'  estimates of COI.
#'
#' @importFrom stats quantile
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
  naive_coi <- calculate_naive_coi(mcmc_results$args$data$data)
  offset_naive_coi <- calculate_naive_coi_offset(mcmc_results$args$data$data, 2)
  base <- summarize_posterior_vectors(
    mcmc_results$chains,
    function(chain) chain$coi,
    mcmc_results$args$data$sample_ids,
    id_name = "sample_id",
    value_prefix = "post_coi",
    lower_quantile = lower_quantile,
    upper_quantile = upper_quantile,
    merge_chains = merge_chains
  )
  if (merge_chains) {
    cois <- lapply(seq_along(mcmc_results$args$data$sample_ids), function(x) c())
    for (chain in mcmc_results$chains) {
      for (s in seq_along(chain$coi)) {
        cois[[s]] <- c(cois[[s]], chain$coi[[s]])
      }
    }
    prob_polyclonal <- sapply(cois, function(x) mean(x > 1))
    return(cbind(base, naive_coi, offset_naive_coi, prob_polyclonal))
  } else {
    prob_polyclonal <- unlist(lapply(mcmc_results$chains, function(chain) {
      sapply(chain$coi, function(x) mean(x > 1))
    }))
    return(cbind(
      base,
      naive_coi = rep(naive_coi, length(mcmc_results$chains)),
      offset_naive_coi = rep(offset_naive_coi, length(mcmc_results$chains)),
      prob_polyclonal
    ))
  }
}

# Internal: unified epsilon summary parameterized by field name and column prefix.
summarize_epsilon_impl <- function(mcmc_results, field, value_prefix, lower_quantile = .025,
                                   upper_quantile = .975, merge_chains = TRUE) {
  summarize_posterior_vectors(
    mcmc_results$chains,
    function(chain) chain[[field]],
    mcmc_results$args$data$sample_ids,
    id_name = "sample_id",
    value_prefix = value_prefix,
    lower_quantile = lower_quantile,
    upper_quantile = upper_quantile,
    merge_chains = merge_chains
  )
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
  summarize_epsilon_impl(mcmc_results, "eps_neg", "post_eps_neg", lower_quantile, upper_quantile, merge_chains)
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
  summarize_epsilon_impl(mcmc_results, "eps_pos", "post_eps_pos", lower_quantile, upper_quantile, merge_chains)
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
    post_allele_freqs <- lapply(seq_along(mcmc_results$args$data$loci), function(x) c())

    for (chain in mcmc_results$chains) {
      for (l in seq_along(chain$allele_freqs)) {
        post_allele_freqs[[l]] <- c(post_allele_freqs[[l]], chain$allele_freqs[[l]])
      }
    }

    post_statistic <- lapply(post_allele_freqs, function(locus_posterior) {
      sapply(locus_posterior, function(allele_freq_sample) fn(allele_freq_sample))
    })

    res <- data.frame(
      locus = mcmc_results$args$data$loci,
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
    res <- lapply(seq_along(mcmc_results$chains), function(idx) {
      post_allele_freqs <- mcmc_results$chains[[idx]]$allele_freqs

      post_statistic <- lapply(post_allele_freqs, function(locus_posterior) {
        sapply(locus_posterior, function(allele_freq_sample) fn(allele_freq_sample))
      })

      chain_res <- data.frame(
        locus = mcmc_results$args$data$loci,
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

#' Summarize expected number of distinct alleles
#'
#' @details Summarize expected number of distinct alleles results from MCMC estimates of allele frequencies.
#' Returns a dataframe that contains summaries of the posterior distribution of expected number of distinct alleles
#'
#' @export
#'
#' @importFrom dplyr bind_rows
#' @importFrom dplyr arrange
#' @importFrom rlang .data
#'
#' @param mcmc_results Result of calling run_mcmc()
#' @param n values of multiplicity of infection to evaluate the expected number of distinct alleles
#' @param lower_quantile The lower quantile of the posterior distribution
#' to return
#' @param upper_quantile The upper quantile of the posterior distribution
#' to return
#' @param merge_chains boolean indicating that all chain results should be merged
summarize_eda <- function(mcmc_results, n, lower_quantile = .025, upper_quantile = .975, merge_chains = TRUE) {
  res <- lapply(n, function(n_) {
    out <- summarize_allele_freq_fn(
      mcmc_results,
      fn = function(x) calculate_eda(x, n_),
      lower_quantile = lower_quantile,
      upper_quantile = upper_quantile,
      merge_chains = merge_chains
    )
    out$n <- n_
    return(out)
  }) |>
    dplyr::bind_rows() |>
    dplyr::arrange(.data$locus, n)

  return(res)
}

names_or_idxs <- function(vec) {
  if (is.null(names(vec))) {
    return(sapply(seq(1, length(vec)), as.character))
  } else {
    return(names(vec))
  }
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
    purrr::map2(
      mcmc_results$args$data$data,
      seq_along(mcmc_results$args$data$data),
      ~ data.frame(locus = mcmc_results$args$data$loci[.y], allele = names_or_idxs(.x[[1]]))
    )
  )

  warning_flag <- FALSE

  if (merge_chains) {
    allele_freq_matrices <- lapply(seq_along(mcmc_results$args$data$loci), function(x) c())
    total_samples <- 0
    for (chain in mcmc_results$chains) {
      total_samples <- total_samples + length(chain$allele_freqs[[1]])
      for (l in seq_along(chain$allele_freqs)) {
        locus <- chain$allele_freqs[[l]]
        allele_freq_matrices[[l]] <- c(allele_freq_matrices[[l]], unlist(locus))
      }
    }
    allele_freq_matrices <- lapply(allele_freq_matrices, function(x) matrix(x, ncol = total_samples))
    res <- lapply(allele_freq_matrices, function(allele_freq_matrix) {
      if (any(is.na(allele_freq_matrix))) {
        if (!warning_flag) {
          warning_flag <- TRUE
          warning("NA values detected in allele frequency matrix. This may indicate a problem with the MCMC chain or there are loci with no diversity. NA values will be replaced with 0.")
        }
        allele_freq_matrix <- replace(allele_freq_matrix, which(is.na(allele_freq_matrix)), 0)
      }
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
    for (idx in seq_along(mcmc_results$chains)) {
      chain <- mcmc_results$chains[[idx]]
      chain_res <- lapply(
        chain$allele_freqs,
        function(locus) {
          num_alleles <- length(locus[[1]])
          allele_freq_matrix <- matrix(unlist(locus), nrow = num_alleles)

          if (any(is.na(allele_freq_matrix))) {
            if (!warning_flag) {
              warning_flag <- TRUE
              warning("NA values detected in allele frequency matrix. This may indicate a problem with the MCMC chain or there are loci with no diversity. NA values will be replaced with 0.")
            }
            allele_freq_matrix <- replace(allele_freq_matrix, which(is.na(allele_freq_matrix)), 0)
          }

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
get_relatedness_vectors <- function(chain) {
  lapply(seq_along(chain$relatedness), function(s) {
    rel <- chain$relatedness[[s]]
    mask <- chain$coi[[s]] > 1
    rel[!mask] <- NA
    rel
  })
}

summarize_relatedness <- function(mcmc_results, lower_quantile = .025, upper_quantile = .975, merge_chains = TRUE) {
  summarize_posterior_vectors(
    mcmc_results$chains,
    get_relatedness_vectors,
    mcmc_results$args$data$sample_ids,
    id_name = "sample_id",
    value_prefix = "post_relatedness",
    lower_quantile = lower_quantile,
    upper_quantile = upper_quantile,
    merge_chains = merge_chains,
    na.rm = TRUE
  )
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
get_effective_coi_vectors <- function(chain) {
  purrr::map2(chain$relatedness, chain$coi, ~ (1 - .x) * (.y - 1) + 1)
}

summarize_effective_coi <- function(mcmc_results, lower_quantile = .025, upper_quantile = .975, merge_chains = TRUE) {
  summarize_posterior_vectors(
    mcmc_results$chains,
    get_effective_coi_vectors,
    mcmc_results$args$data$sample_ids,
    id_name = "sample_id",
    value_prefix = "post_effective_coi",
    lower_quantile = lower_quantile,
    upper_quantile = upper_quantile,
    merge_chains = merge_chains
  )
}

#' Summarize population assignments
#'
#' @details Summarize population assignment probabilities from MCMC. Returns
#'  two dataframes: one with the distribution of most likely population assignments
#'  (probability that each population is the most likely) and another with entropy
#'  summaries (quantifying uncertainty around population assignment).
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
#' @param correct_label_switching boolean indicating whether to correct for label switching.
#'   If TRUE, population labels are aligned across MCMC steps using correlation-based matching.
#' @param label_switching_method Method for label switching correction: "iterative" (align each
#'   step to previous) or "reference" (align all steps to first step). Only used if
#'   correct_label_switching is TRUE.
#' @return A list with two dataframes:
#'   - `most_likely_population`: For each sample, the probability that each population
#'     is the most likely assignment
#'   - `entropy_summary`: For each sample, summary statistics of the entropy
#'     distribution (quantifying uncertainty)
summarize_population_assignments <- function(mcmc_results, lower_quantile = .025, upper_quantile = .975, merge_chains = TRUE, correct_label_switching = TRUE, label_switching_method = "iterative") {
  
  # Helper function to calculate entropy
  calculate_entropy <- function(log_probs) {
    entropy <- -sum(exp(log_probs) * log_probs)
    return(entropy)
  }
  
  # Helper function to transpose population assignments from sample-first to step-first format
  transpose_assignments <- function(chain_assignments) {
    num_samples <- length(chain_assignments)
    if (num_samples == 0) return(list())
    num_steps <- length(chain_assignments[[1]])
    if (num_steps == 0) return(list())
    
    # Transpose: from sample-first to step-first
    step_first <- lapply(seq_len(num_steps), function(step_idx) {
      lapply(seq_len(num_samples), function(sample_idx) {
        chain_assignments[[sample_idx]][[step_idx]]
      })
    })
    return(step_first)
  }
  
  # Helper function to transpose back from step-first to sample-first format
  transpose_assignments_back <- function(step_first_assignments) {
    num_steps <- length(step_first_assignments)
    if (num_steps == 0) return(list())
    num_samples <- length(step_first_assignments[[1]])
    if (num_samples == 0) return(list())
    
    # Transpose: from step-first to sample-first
    sample_first <- lapply(seq_len(num_samples), function(sample_idx) {
      lapply(seq_len(num_steps), function(step_idx) {
        step_first_assignments[[step_idx]][[sample_idx]]
      })
    })
    return(sample_first)
  }
  
  if (merge_chains) {
    # Get number of populations from the first chain
    num_populations <- length(mcmc_results$chains[[1]]$population_assignment[[1]][[1]])
    
    # Initialize storage for each sample
    most_likely_counts <- lapply(seq_along(mcmc_results$args$data$sample_ids), function(sample_idx) {
      rep(0, num_populations)
    })
    entropy_values <- lapply(seq_along(mcmc_results$args$data$sample_ids), function(sample_idx) {
      c()
    })
    total_steps <- length(mcmc_results$chains[[1]]$population_assignment[[1]])
    
    # Collect data across all chains
    for (chain in mcmc_results$chains) {
      # Apply label switching correction if requested
      chain_assignments <- chain$population_assignment
      if (correct_label_switching) {
        # Transpose to step-first format for correction
        step_first <- transpose_assignments(chain_assignments)
        # Apply correction
        corrected_step_first <- correct_label_switching(step_first, method = label_switching_method)
        # Transpose back to sample-first format
        chain_assignments <- transpose_assignments_back(corrected_step_first)
      }
      
      for (sample_idx in seq_along(chain_assignments)) {
        for (step_idx in seq_along(chain_assignments[[sample_idx]])) {
          step_log_probs <- chain_assignments[[sample_idx]][[step_idx]]
          
          # Aggregate probabilities across steps
          most_likely_counts[[sample_idx]] <- most_likely_counts[[sample_idx]] + (exp(step_log_probs) / total_steps)
          
          # Calculate entropy for this step
          entropy_val <- calculate_entropy(step_log_probs)
          entropy_values[[sample_idx]] <- c(entropy_values[[sample_idx]], entropy_val)
        }
      }
    }
    
    # Calculate most likely population probabilities
    most_likely_list <- list()
    for (sample_idx in seq_along(mcmc_results$args$data$sample_ids)) {
      sample_id <- mcmc_results$args$data$sample_ids[sample_idx]
      total_steps <- sum(most_likely_counts[[sample_idx]])
      
      for (pop_idx in seq_len(num_populations)) {
        prob_most_likely <- most_likely_counts[[sample_idx]][pop_idx] / total_steps
        
        most_likely_list[[length(most_likely_list) + 1]] <- data.frame(
          sample_id = sample_id,
          population = pop_idx,
          prob_most_likely = prob_most_likely
        )
      }
    }
    
    # Calculate entropy summaries
    entropy_list <- list()
    for (sample_idx in seq_along(mcmc_results$args$data$sample_ids)) {
      sample_id <- mcmc_results$args$data$sample_ids[sample_idx]
      entropy_vals <- entropy_values[[sample_idx]]
      
      # Filter out any remaining -Inf, Inf, or NaN values
      valid_entropy <- entropy_vals[is.finite(entropy_vals)]
      
      if (length(valid_entropy) == 0) {
        # If no valid entropy values, set all to 0
        entropy_lower <- entropy_med <- entropy_upper <- entropy_mean <- 0
      } else {
        entropy_lower <- quantile(valid_entropy, lower_quantile)
        entropy_med <- quantile(valid_entropy, .5)
        entropy_upper <- quantile(valid_entropy, upper_quantile)
        entropy_mean <- mean(valid_entropy)
      }
      
      entropy_list[[length(entropy_list) + 1]] <- data.frame(
        sample_id = sample_id,
        entropy_lower = entropy_lower,
        entropy_med = entropy_med,
        entropy_upper = entropy_upper,
        entropy_mean = entropy_mean
      )
    }
    
    return(list(
      most_likely_population = do.call(rbind, most_likely_list),
      entropy_summary = do.call(rbind, entropy_list)
    ))
    
  } else {
    # Handle separate chains
    chain_results_most_likely <- list()
    chain_results_entropy <- list()
    
    for (chain_idx in seq_along(mcmc_results$chains)) {
      chain <- mcmc_results$chains[[chain_idx]]
      num_populations <- length(chain$population_assignment[[1]][[1]])
      
      # Apply label switching correction if requested
      chain_assignments <- chain$population_assignment
      if (correct_label_switching) {
        # Transpose to step-first format for correction
        step_first <- transpose_assignments(chain_assignments)
        # Apply correction
        corrected_step_first <- correct_label_switching(step_first, method = label_switching_method)
        # Transpose back to sample-first format
        chain_assignments <- transpose_assignments_back(corrected_step_first)
      }
      
      # Initialize storage for this chain
      most_likely_counts <- lapply(seq_along(chain_assignments), function(sample_idx) {
        rep(0, num_populations)
      })
      entropy_values <- lapply(seq_along(chain_assignments), function(sample_idx) {
        c()
      })
      
      # Process this chain
      total_steps <- length(chain_assignments[[1]])
      for (sample_idx in seq_along(chain_assignments)) {
        for (step_idx in seq_along(chain_assignments[[sample_idx]])) {
          step_log_probs <- chain_assignments[[sample_idx]][[step_idx]]
          
          # Aggregate probabilities across steps
          most_likely_counts[[sample_idx]] <- most_likely_counts[[sample_idx]] + (exp(step_log_probs) / total_steps)
          
          # Calculate entropy for this step
          entropy_val <- calculate_entropy(step_log_probs)
          entropy_values[[sample_idx]] <- c(entropy_values[[sample_idx]], entropy_val)
        }
      }
      
      # Calculate results for this chain
      most_likely_list <- list()
      entropy_list <- list()
      
      for (sample_idx in seq_along(chain$population_assignment)) {
        sample_id <- mcmc_results$args$data$sample_ids[sample_idx]
        total_steps <- sum(most_likely_counts[[sample_idx]])
        
        # Most likely population probabilities
        for (pop_idx in seq_len(num_populations)) {
          prob_most_likely <- most_likely_counts[[sample_idx]][pop_idx] / total_steps
          
          most_likely_list[[length(most_likely_list) + 1]] <- data.frame(
            sample_id = sample_id,
            population = pop_idx,
            prob_most_likely = prob_most_likely,
            chain = chain_idx
          )
        }
        
        # Entropy summaries
        entropy_vals <- entropy_values[[sample_idx]]
        
        # Filter out any remaining -Inf, Inf, or NaN values
        valid_entropy <- entropy_vals[is.finite(entropy_vals)]
        
        if (length(valid_entropy) == 0) {
          # If no valid entropy values, set all to 0
          entropy_lower <- entropy_med <- entropy_upper <- entropy_mean <- 0
        } else {
          entropy_lower <- quantile(valid_entropy, lower_quantile)
          entropy_med <- quantile(valid_entropy, .5)
          entropy_upper <- quantile(valid_entropy, upper_quantile)
          entropy_mean <- mean(valid_entropy)
        }
        
        entropy_list[[length(entropy_list) + 1]] <- data.frame(
          sample_id = sample_id,
          entropy_lower = entropy_lower,
          entropy_med = entropy_med,
          entropy_upper = entropy_upper,
          entropy_mean = entropy_mean,
          chain = chain_idx
        )
      }
      
      chain_results_most_likely[[chain_idx]] <- do.call(rbind, most_likely_list)
      chain_results_entropy[[chain_idx]] <- do.call(rbind, entropy_list)
    }
    
    return(list(
      most_likely_population = do.call(rbind, chain_results_most_likely),
      entropy_summary = do.call(rbind, chain_results_entropy)
    ))
  }
}
