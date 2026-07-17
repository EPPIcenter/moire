#' Build and raise a classed Initialization failure condition
#'
#' @param diagnostics List from C++ Initialization diagnostics
#' @param sample_ids Character vector of sample IDs (optional names)
#' @param loci Character vector of locus IDs (optional names)
#' @noRd
stop_initialization_failure <- function(diagnostics, sample_ids = NULL, loci = NULL) {
  msg <- format_initialization_failure_message(diagnostics, sample_ids, loci)
  err <- structure(
    list(
      message = msg,
      call = NULL,
      diagnostics = diagnostics
    ),
    class = c("moire_initialization_failure", "error", "condition")
  )
  stop(err)
}

#' @noRd
format_location <- function(index, names = NULL) {
  if (is.null(index) || length(index) == 0 || is.na(index)) {
    return(NULL)
  }
  label <- as.character(index)
  if (!is.null(names) && index >= 1 && index <= length(names)) {
    label <- sprintf("%s (%s)", index, names[[index]])
  }
  label
}

#' @noRd
format_initialization_failure_message <- function(diagnostics,
                                                  sample_ids = NULL,
                                                  loci = NULL) {
  n_ill <- if (is.null(diagnostics$total_ill_conditioned)) 0L else diagnostics$total_ill_conditioned
  n_tries <- if (is.null(diagnostics$max_initialization_tries)) NA_integer_ else diagnostics$max_initialization_tries
  n_chains <- if (is.null(diagnostics$chains_attempted)) NA_integer_ else diagnostics$chains_attempted
  classification <- if (is.null(diagnostics$classification)) "hard_starting_set" else diagnostics$classification
  concentration <- if (is.null(diagnostics$concentration)) "" else diagnostics$concentration
  share <- if (is.null(diagnostics$dominant_share)) NA_real_ else diagnostics$dominant_share

  lines <- c(
    "Initialization failed: no finite genotyping log-likelihood after Initialization retries.",
    sprintf(
      "Attempted %s chain(s); recorded %s Ill-conditioned start(s) (max_initialization_tries = %s).",
      n_chains, n_ill, n_tries
    )
  )

  if (identical(classification, "consistent_failure_cause")) {
    if (identical(concentration, "locus")) {
      loc <- format_location(diagnostics$dominant_locus, loci)
      lines <- c(
        lines,
        sprintf(
          "Consistent failure cause: locus %s in %.0f%% of Ill-conditioned starts.",
          loc, 100 * share
        ),
        "Guidance: inspect this locus; high allelic diversity or extremely rare alleles often destabilize the genotyping likelihood. Consider dropping rare alleles or this locus."
      )
      sec <- format_location(diagnostics$secondary_sample, sample_ids)
      if (!is.null(sec)) {
        lines <- c(
          lines,
          sprintf("Secondary concentration on sample %s was also observed.", sec)
        )
      }
    } else if (identical(concentration, "sample")) {
      samp <- format_location(diagnostics$dominant_sample, sample_ids)
      lines <- c(
        lines,
        sprintf(
          "Consistent failure cause: sample %s in %.0f%% of Ill-conditioned starts.",
          samp, 100 * share
        ),
        "Guidance: inspect this sample's missingness or call patterns; consider excluding it for a sanity run."
      )
    }
  } else {
    lines <- c(
      lines,
      "Classification: hard starting set (no single locus or sample dominates failures).",
      "Guidance: increase max_initialization_tries; if it still fails, simplify the dataset (fewer loci/alleles) and retry."
    )
  }

  ex_s <- diagnostics$examples$sample
  ex_l <- diagnostics$examples$locus
  if (length(ex_s) > 0 && length(ex_l) > 0) {
    example_bits <- vapply(seq_along(ex_s), function(i) {
      sprintf(
        "sample %s, locus %s",
        format_location(ex_s[[i]], sample_ids),
        format_location(ex_l[[i]], loci)
      )
    }, character(1))
    lines <- c(lines, paste("Example Failure loci:", paste(example_bits, collapse = "; ")))
  }

  prior_terms <- diagnostics$prior_nonfinite_counts$term
  prior_n <- diagnostics$prior_nonfinite_counts$n
  if (length(prior_terms) > 0) {
    prior_bits <- paste0(prior_terms, "=", prior_n, collapse = ", ")
    lines <- c(
      lines,
      paste("Secondary non-finite prior terms observed:", prior_bits)
    )
  }

  paste(lines, collapse = "\n")
}
