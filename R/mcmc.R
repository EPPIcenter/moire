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
#' @param thin Positive Integer. How often to sample from mcmc,
#'  1 means do not thin
#' @param burnin Positive Integer. Number of MCMC samples to
#'  discard as burnin
#' @param samples Positive Integer. Number of samples to take
#'  after burnin
#' @param verbose Logical indicating if progress is printed
#' @param allele_freq_threshold 0-1 Numeric. Lowest allowed value for an
#'  allele frequency.
#' @param eps_pos_0 Numeric. Initial eps_pos value
#' @param eps_pos_var Numeric. Variance used in sampling eps_pos
#' @param eps_pos_shape Positive Numeric. Shape parameter in
#'  Gamma distribution for eps_pos prior
#' @param eps_pos_scale Positive Numeric. Scale parameter in
#'  Gamma distribution for eps_pos prior
#' @param eps_neg_0 Numeric. Initial eps_neg value
#' @param eps_neg_var Numeric. Variance used in sampling eps_neg
#' @param eps_neg_shape Positive Numeric. Shape parameter in
#'  Gamma distribution for eps_neg prior
#' @param eps_neg_scale Positive Numeric. Scale parameter in
#'  Gamma distribution for eps_neg prior
#' @param max_eps_pos Numeric. Maximum allowed value for eps_pos
#' @param max_eps_neg Numeric. Maximum allowed value for eps_neg
#' @param max_coi Positive Numeric. Maximum allowed complexity of infection
#' @param allele_freq_vars Positive Numeric. Variance used in sampling allele
#'  frequencies
#' @param adapt_allele_freq_vars Logical indicating whether to adapt variance
#'  to achieve an acceptance rate of .27
run_mcmc <-
  function(data,
           sample_ids,
           loci,
           is_missing = FALSE,
           thin = 1,
           burnin = 1e4,
           samples = 1e4,
           verbose = TRUE,
           allele_freq_threshold = 1e-5,
           eps_pos_0 = .005,
           eps_pos_var = 1,
           eps_pos_alpha = .5,
           eps_pos_beta = 99.5,
           eps_neg_0 = .005,
           eps_neg_var = 1,
           eps_neg_alpha = .5,
           eps_neg_beta = 99.5,
           max_eps_pos = 2,
           max_eps_neg = 2,
           max_coi = 20,
           allele_freq_vars = 1,
           adapt_allele_freq_vars = TRUE) {
    args <- as.list(environment())

    ## if is_missing == FALSE, then generate a default FALSE matrix
    if (class(args$is_missing) == "logical" && args$is_missing == FALSE) {
      num_loci <- length(args$data)
      num_biological_samples <- length(args$data[[1]])
      args$is_missing <- matrix(
        FALSE,
        nrow = num_loci,
        ncol = num_biological_samples
      )
    }

    if (length(args$allele_freq_vars)) {
      args$allele_freq_vars <- rep(args$allele_freq_vars, length(loci))
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

    res <- run_mcmc_rcpp(args) # nolint
    res$args <- args
    res$total_samples <- args$samples / args$thin
    res
  }
