#' Sample from the target distribution using MCMC
#'
#' @export
#'
#' @param data Data to be used in MCMC, expected input is a list
#'  of lists of numeric vectors, where each list element is a
#'  collection of observations across samples at a single genetic locus.
#' @param sample_ids Vector containing the ordered unique identifier for samples
#' @param loci Vector containing the ordered unique identifer for loci
#' @param is_missing Boolean matrix indicating whether the observation
#'  should be treated as missing data and ignored. Number of rows equals the
#'  number of loci, number of columns equals the number samples. Alternatively,
#'  the user may pass in FALSE if no data should be considered missing.
#' @param allow_relatedness Bool indicating whether or not to allow relatedness within
#' host
#' @param thin Positive Integer. How often to sample from mcmc,
#'  1 means do not thin
#' @param burnin Positive Integer. Number of MCMC samples to
#'  discard as burnin
#' @param samples_per_chain Positive Integer. Number of samples to take
#'  after burnin
#' @param verbose Logical indicating if progress is printed
#' @param allele_freq_threshold 0-1 Numeric. Lowest allowed value for an
#'  allele frequency.
#' @param eps_pos_0 Numeric. Initial eps_pos value
#' @param eps_pos_var Numeric. Variance used in sampling eps_pos
#' @param eps_pos_alpha Positive Numeric. Alpha parameter in
#'  Beta distribution for eps_pos prior
#' @param eps_pos_beta Positive Numeric. Beta parameter in
#'  Beta distribution for eps_pos prior
#' @param eps_neg_0 Numeric. Initial eps_neg value
#' @param eps_neg_var Numeric. Variance used in sampling eps_neg
#' @param eps_neg_alpha Positive Numeric. Alpha parameter in
#'  Beta distribution for eps_neg prior
#' @param eps_neg_beta Positive Numeric. Beta parameter in
#'  Beta distribution for eps_neg prior
#' @param max_eps_pos Numeric. Maximum allowed value for eps_pos
#' @param max_eps_neg Numeric. Maximum allowed value for eps_neg
#' @param max_coi Positive Numeric. Maximum allowed complexity of infection
#' @param allele_freq_vars Positive Numeric. Variance used in sampling allele
#'  frequencies
#' @param adapt_allele_freq_vars Logical indicating whether to adapt variance
#'  to achieve an acceptance rate of .27
#' @param num_chains Total number of chains to run, possibly simultaneously
#' @param num_cores Total cores to use to run chains
run_mcmc <-
  function(data,
           sample_ids,
           loci,
           is_missing = FALSE,
           allow_relatedness = TRUE,
           thin = 1,
           burnin = 1e4,
           samples_per_chain = 1e3,
           verbose = TRUE,
           allele_freq_threshold = 1e-12,
           eps_pos_0 = .1,
           eps_pos_var = 1,
           eps_pos_alpha = .1,
           eps_pos_beta = 9.9,
           eps_neg_0 = .1,
           eps_neg_var = 1,
           eps_neg_alpha = .1,
           eps_neg_beta = 9.9,
           max_eps_pos = 2,
           max_eps_neg = 2,
           max_coi = 20,
           allele_freq_vars = 1,
           adapt_allele_freq_vars = TRUE,
           num_chains = 1,
           num_cores = 1) {
    args <- as.list(environment())
    mcmc_args <- as.list(environment())

    ## if is_missing == FALSE, then generate a default FALSE matrix
    suppressWarnings({
      if (inherits(mcmc_args$is_missing, "logical") && mcmc_args$is_missing == FALSE) {
        num_loci <- length(mcmc_args$data)
        num_biological_samples <- length(mcmc_args$data[[1]])
        mcmc_args$is_missing <- matrix(
          FALSE,
          nrow = num_loci,
          ncol = num_biological_samples
        )
      }
    })

    if (length(mcmc_args$allele_freq_vars)) {
      mcmc_args$allele_freq_vars <- rep(mcmc_args$allele_freq_vars, length(loci))
    }

    total_alleles <- lapply(data, function(x) {
      return(length(x[[1]]))
    })

    if (any(total_alleles < 2)) {
      stop("Loci with less than 2 alleles present, remove these loci")
    }

    if (max_coi < 1) {
      stop("Max COI must be greater than 1")
    }

    res <- list()
    chains <- list()
    if (num_chains > 1) {
      chains <- parallel::mclapply(
        1:num_chains,
        function(i) {
          mcmc_args$chain_number <- i
          mcmc_args$simple_verbose <- (mcmc_args$num_chains > 1)
          mcmc_args$samples <- round(mcmc_args$samples_per_chain)
          chain <- run_mcmc_rcpp(mcmc_args)
          return(chain)
        },
        mc.cores = num_cores
      )
    } else {
      mcmc_args$chain_number <- 1
      mcmc_args$simple_verbose <- FALSE
      mcmc_args$samples <- mcmc_args$samples_per_chain
      chain <- run_mcmc_rcpp(mcmc_args)
      chains[[1]] <- chain
    }


    res$chains <- chains
    res$args <- args
    res
  }
