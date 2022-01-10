
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
#' @export
#'
#' @param mcmc_results Result of calling run_mcmc()
#' @param lower_quantile The lower quantile of the posterior
#'  distribution to return
#' @param upper_quantile The upper quantile of the posterior
#'  distribution to return
#' @param naive_offset Offset used in calculate_naive_coi_offset()
summarize_coi <- function(mcmc_results, lower_quantile = .025,
                          upper_quantile = .975, naive_offset = 2) {
  cois <- mcmc_results$coi
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

  naive_coi <- calculate_naive_coi(mcmc_results$args$data)
  offset_naive_coi <- calculate_naive_coi_offset(mcmc_results$args$data, 2)
  coi_data <- data.frame(
    sample_id = mcmc_results$args$sample_ids,
    post_coi_lower, post_coi_med, post_coi_upper,
    post_coi_mean, naive_coi, offset_naive_coi
  )
  return(coi_data)
}

#' Summarize Function of Allele Frequencies
#'
#' @details General function to summarize the posterior distribution of
#' functions of the sampled allele frequencies
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
summarize_allele_freq_fn <- function(mcmc_results, fn,
                                     lower_quantile = .025,
                                     upper_quantile = .975) {
  post_allele_freqs <- mcmc_results$allele_freqs
  post_statistic <- lapply(post_allele_freqs, function(locus_posterior) {
    sapply(locus_posterior, function(allele_freq_sample) fn(allele_freq_sample))
  })

  res <- data.frame(
    loci = mcmc_results$args$loci,
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
summarize_he <- function(mcmc_results,
                         lower_quantile = .025,
                         upper_quantile = .975) {
  res <- summarize_allele_freq_fn(
    mcmc_results,
    fn = calculate_he,
    lower_quantile = lower_quantile,
    upper_quantile = upper_quantile
  )
  return(res)
}


#' Summarize allele frequencies
#'
#' @details Summarize individual allele frequencies from the posterior
#'  distribution of sampled allele frequencies
#'
#' @export
#'
#' @param mcmc_results Result of calling run_mcmc()
#' @param lower_quantile The lower quantile of the posterior distribution
#'  to return
#' @param upper_quantile The upper quantile of the posterior distribution
#' to return
summarize_allele_freqs <- function(mcmc_results,
                                   lower_quantile = .025,
                                   upper_quantile = .975) {
  res <- lapply(
    mcmc_results$allele_freqs,
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
        post_allele_freqs_mean = post_allele_freqs_mean,
        num_alleles = num_alleles
      )
    }
  )
  return(do.call("rbind", res))
}
