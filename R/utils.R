#' Prepare initial allele frequencies for MCMC
#'
#' @details This function helps prepare initial allele frequencies for the MCMC
#' algorithm. It takes a list of population assignments and the genotyping data,
#' then calculates allele frequencies for each population. This allows users to
#' use custom clustering strategies (e.g., hierarchical clustering, k-means, etc.)
#' to initialize populations before running the MCMC.
#'
#' @export
#'
#' @param data The genotyping data object (as returned by `load_*_data` functions)
#' @param population_assignments A vector of population assignments for each sample.
#' Should be integers from 1 to `num_populations`, with length equal to the number
#' of samples in the data.
#' @param num_populations The number of populations to initialize
#' @param pseudocount Numeric value to add as pseudocount to avoid zero frequencies
#' @return A list of length `num_populations`, where each element is a list of
#' length `num_loci`, and each locus element is a numeric vector of allele frequencies
#' that sum to 1.
#'
#' @examples
#' # Example with simulated data
#' data <- simulate_data(
#'   num_samples = 100,
#'   mean_coi = 2,
#'   epsilon_pos = 0.01,
#'   epsilon_neg = 0.01,
#'   locus_freq_alphas = list(c(1, 1), c(1, 1, 1))
#' )
#' # Use custom clustering (e.g., hierarchical clustering)
#' # population_assignments <- cutree(hclust(dist_matrix), k = 3)
#' # initial_freqs <- prepare_initial_allele_frequencies(data, population_assignments, 3)
cat_progress_header <- function(title, data) {
  cat("=== ", title, " ===\n", sep = "")
  cat("Number of samples:", length(data$sample_ids), "\n")
  cat("Number of loci:", length(data$loci), "\n")
}

prepare_initial_allele_frequencies <- function(data, population_assignments, num_populations, pseudocount = 1) {
  cat_progress_header("Preparing Initial Allele Frequencies", data)
  cat("Number of populations:", num_populations, "\n")
  cat("Population assignments:", table(population_assignments), "\n")
  
  if (length(population_assignments) != length(data$sample_ids)) {
    stop("population_assignments must have length equal to number of samples")
  }
  
  if (max(population_assignments) > num_populations || min(population_assignments) < 1) {
    stop("population_assignments must be integers from 1 to num_populations")
  }
  
  num_loci <- length(data$loci)
  result <- vector("list", num_populations)
  
  for (pop_idx in 1:num_populations) {
    pop_samples <- which(population_assignments == pop_idx)
    cat("Population", pop_idx, "has", length(pop_samples), "samples\n")
    
    if (length(pop_samples) == 0) {
      warning(paste("No samples assigned to population", pop_idx, "- using uniform frequencies"))
      # Use uniform frequencies if no samples assigned to this population
      pop_frequencies <- vector("list", num_loci)
      for (locus_idx in 1:num_loci) {
        num_alleles <- length(data$data[[locus_idx]][[1]])
        pop_frequencies[[locus_idx]] <- rep(1/num_alleles, num_alleles)
      }
      cat("  -> Using uniform frequencies for empty population\n")
    } else {
      # Calculate allele frequencies for this population
      pop_frequencies <- vector("list", num_loci)
      
      for (locus_idx in 1:num_loci) {
        num_alleles <- length(data$data[[locus_idx]][[1]])
        allele_counts <- rep(pseudocount, num_alleles)  # Start with pseudocount
        
        # Sum up allele counts across samples in this population
        for (sample_idx in pop_samples) {
          if (sample_idx <= length(data$data[[locus_idx]])) {
            sample_alleles <- data$data[[locus_idx]][[sample_idx]]
            for (allele_idx in 1:num_alleles) {
              if (allele_idx <= length(sample_alleles)) {
                allele_counts[allele_idx] <- allele_counts[allele_idx] + sample_alleles[allele_idx]
              }
            }
          }
        }
        
        # Normalize to get frequencies
        total_count <- sum(allele_counts)
        if (total_count > 0) {
          pop_frequencies[[locus_idx]] <- allele_counts / total_count
        } else {
          pop_frequencies[[locus_idx]] <- rep(1/num_alleles, num_alleles)
        }
      }
      cat("  -> Calculated frequencies from", length(pop_samples), "samples\n")
    }
    
    result[[pop_idx]] <- pop_frequencies
  }
  
  # Log summary of calculated frequencies
  cat("\n=== Initial Allele Frequencies Summary ===\n")
  for (pop_idx in 1:num_populations) {
    cat("Population", pop_idx, ":\n")
    for (locus_idx in 1:min(3, num_loci)) {  # Show first 3 loci
      freqs <- result[[pop_idx]][[locus_idx]]
      cat("  Locus", locus_idx, ":", sprintf("%.3f", freqs), "(sum:", sprintf("%.3f", sum(freqs)), ")\n")
    }
    if (num_loci > 3) {
      cat("  ... and", num_loci - 3, "more loci\n")
    }
  }
  cat("==========================================\n\n")
  
  return(result)
}

# ---- Profiling helpers ----
#' Get C++ profiler stats
#' @return data.frame with key, calls, total_ms, avg_ms
#' @export
moire_prof_stats <- function() {
  if (!exists("moire_profiler_stats")) {
    stop("Profiler not available; rebuild and load package.")
  }
  df <- moire_profiler_stats()
  df[order(-df$total_ms), ]
}

#' Reset C++ profiler stats
#' @export
moire_prof_reset <- function() {
  if (!exists("moire_profiler_reset")) {
    stop("Profiler not available; rebuild and load package.")
  }
  invisible(moire_profiler_reset())
}

#' Load long form data
#'
#' @details Long form data is a data frame with
#'  3 columns: `sample_id`, `locus`, `allele`. Returned data contains
#'  vectors `sample_ids` and `loci` that are ordered as the results
#'  will be ordered from running the MCMC algorithm.
#'
#' @export
#'
#' @param df data frame with 3 columns: `sample_id`, `locus`, `allele`.
#' Each row is a single observation of an allele at a particular
#' locus for a given sample.
#' @param warn_uninformative boolean whether or not to print message when
#'  removing uninformative loci
#'
#' @importFrom rlang .data
load_long_form_data <- function(df, warn_uninformative = TRUE) {
  uninformative_loci <- df |>
    dplyr::ungroup() |>
    dplyr::group_by(.data$locus) |>
    dplyr::summarise(total_alleles = length(unique(.data$allele))) |>
    dplyr::filter(.data$total_alleles == 1) |>
    dplyr::pull(.data$locus)

  if (length(uninformative_loci) > 0) {
    if (warn_uninformative) {
      message("Uninformative loci with only 1 allele included. Removing...")
    }
    df <- df |>
      dplyr::filter(!(.data$locus %in% uninformative_loci))
  }

  unique_alleles <- df |>
    dplyr::group_by(.data$locus) |>
    dplyr::summarise(unique_alleles = list(sort(unique(.data$allele))))

  num_loci <- df |>
    dplyr::pull(.data$locus) |>
    dplyr::n_distinct()

  sample_locus_barcodes <- df |>
    dplyr::group_by(.data$sample_id, .data$locus) |>
    dplyr::summarize(alleles_grp = list(.data$allele), .groups = "drop") |>
    dplyr::mutate(missing = FALSE) |>
    tidyr::complete(.data$sample_id, .data$locus,
      fill = list(alleles_grp = list(c(NULL)), missing = TRUE)
    )

  missing_vec <- sample_locus_barcodes |>
    dplyr::arrange(.data$sample_id, .data$locus) |>
    dplyr::pull(.data$missing)

  sample_locus_barcodes <- sample_locus_barcodes |>
    dplyr::left_join(unique_alleles, by = "locus") |>
    dplyr::rowwise("sample_id", "locus", "missing") |>
    dplyr::summarise(
      barcode = list(sapply(
        unique_alleles,
        function(x) {
          as.integer(x %in% .data$alleles_grp)
        }
      )), .groups = "drop"
    ) |>
    dplyr::group_by(.data$locus) |>
    dplyr::arrange(.data$sample_id, by_group = TRUE) |>
    dplyr::summarise(locus_barcodes = list(.data$barcode))


  sample_ids <- df |>
    dplyr::select("sample_id") |>
    dplyr::arrange(.data$sample_id) |>
    dplyr::distinct() |>
    dplyr::pull(.data$sample_id)

  is_missing <- matrix(missing_vec, nrow = num_loci)

  return(list(
    sample_ids = sample_ids,
    data = sample_locus_barcodes$locus_barcodes,
    loci = sample_locus_barcodes$locus,
    is_missing = is_missing,
    uninformative_loci = uninformative_loci
  ))
}

#' Load delimited data
#'
#' @details Load `data.frame` with a `sample_id` column and the remaining
#'  columns are `loci`. Each cell contains a separator delimited string
#'  representing the observed alleles at that locus for that sample.
#'  Returned data contains vectors `sample_ids` and `loci` that are ordered
#'  as the results will be ordered from running the MCMC algorithm.
#'
#' @export
#'
#' @param data data.frame containing the described data
#' @param sep string used to separate alleles
#' @param warn_uninformative boolean whether or not to print message when
#'  removing uninformative loci
#'
#' @importFrom rlang .data
load_delimited_data <- function(data, sep = ";", warn_uninformative = TRUE) {
  df <- data |>
    tidyr::pivot_longer(-"sample_id",
      names_to = "locus",
      values_to = "allele"
    ) |>
    tidyr::separate_rows("allele", sep = sep) |>
    dplyr::filter(!is.na(.data$allele))
  return(load_long_form_data(df))
}


#' Convert loaded data back to long form
#'
#' @details This function converts the internal data structure back to long form
#' format (data frame with columns `sample_id`, `locus`, `allele`). This is useful
#' for exporting data, visualization, or further analysis in standard formats.
#' The function reconstructs the original allele observations from the barcode
#' format used internally by the MCMC algorithm.
#'
#' @export
#'
#' @param data The loaded data object (as returned by `load_long_form_data` or `load_delimited_data`)
#' @return A data frame with columns `sample_id`, `locus`, `allele` representing
#' the original long form data structure
#'
#' @examples
#' # Load some data
#' df <- data.frame(
#'   sample_id = c("S1", "S1", "S2", "S2"),
#'   locus = c("L1", "L1", "L1", "L1"),
#'   allele = c("A", "B", "A", "C")
#' )
#' loaded_data <- load_long_form_data(df)
#' 
#' # Convert back to long form
#' long_form <- convert_to_long_form(loaded_data)
#' print(long_form)
convert_to_long_form <- function(data) {
  cat_progress_header("Converting Data to Long Form", data)

  # Initialize result data frame
  result_rows <- list()
  row_count <- 0
  
  # Get unique alleles for each locus (reconstruct from barcode data)
  unique_alleles_per_locus <- list()
  for (locus_idx in seq_along(data$loci)) {
    # Find all unique alleles by looking at the barcode data
    all_alleles <- c()
    for (sample_idx in seq_along(data$sample_ids)) {
      if (sample_idx <= length(data$data[[locus_idx]])) {
        barcode <- data$data[[locus_idx]][[sample_idx]]
        # Convert barcode back to allele names
        # The barcode is a binary vector indicating presence/absence of alleles
        # We need to reconstruct the original allele names
        for (allele_idx in seq_along(barcode)) {
          if (barcode[allele_idx] > 0) {
            # Create allele name based on position
            allele_name <- paste0("Allele_", allele_idx)
            all_alleles <- c(all_alleles, allele_name)
          }
        }
      }
    }
    unique_alleles_per_locus[[locus_idx]] <- sort(unique(all_alleles))
  }
  
  # Process each locus and sample
  for (locus_idx in seq_along(data$loci)) {
    locus_name <- data$loci[locus_idx]
    cat("Processing locus:", locus_name, "\n")
    
    for (sample_idx in seq_along(data$sample_ids)) {
      sample_id <- data$sample_ids[sample_idx]
      
      # Check if this sample-locus combination is missing
      is_missing_sample_locus <- if (nrow(data$is_missing) > 0 && ncol(data$is_missing) > 0) {
        data$is_missing[locus_idx, sample_idx]
      } else {
        FALSE
      }
      
      if (is_missing_sample_locus) {
        # Skip missing data
        next
      }
      
      if (sample_idx <= length(data$data[[locus_idx]])) {
        barcode <- data$data[[locus_idx]][[sample_idx]]
        
        # Convert barcode to allele observations
        for (allele_idx in seq_along(barcode)) {
          if (barcode[allele_idx] > 0) {
            # Create allele name
            allele_name <- paste0("Allele_", allele_idx)
            
            # Add one row for each copy of the allele
            for (copy in 1:barcode[allele_idx]) {
              row_count <- row_count + 1
              result_rows[[row_count]] <- data.frame(
                sample_id = sample_id,
                locus = locus_name,
                allele = allele_name,
                stringsAsFactors = FALSE
              )
            }
          }
        }
      }
    }
  }
  
  # Combine all rows
  if (row_count > 0) {
    result_df <- do.call(rbind, result_rows)
  } else {
    result_df <- data.frame(
      sample_id = character(0),
      locus = character(0),
      allele = character(0),
      stringsAsFactors = FALSE
    )
  }
  
  cat("Generated", nrow(result_df), "observations\n")
  cat("==========================================\n\n")
  
  return(result_df)
}

#' Plot chain swap acceptance rates
#'
#' @details Plot the swap acceptance rates for each chain.
#' The x-axis is the temperature, and the y-axis is the swap acceptance rate.
#' The dashed lines indicate the temperatures used for parallel tempering.
#'
#' @export
#'
#' @param mcmc_results list of results from `run_mcmc`
#'
#' @importFrom ggplot2 aes coord_cartesian geom_point geom_vline ggplot
#' @importFrom rlang .data
#'
#' @return list of ggplot objects
#'
plot_chain_swaps <- function(mcmc_results) {
  plots <- lapply(mcmc_results$chains, function(chain) {
    # swaps for a chain happen every 2 samples
    swaps_per_chain <- mcmc_results$args$samples_per_chain / 2
    swap_dist <- chain$swap_acceptances / swaps_per_chain
    temps <- mcmc_results$chains[[1]]$temp_gradient
    swap_idx <- (temps[1:length(temps) - 1] + temps[2:length(temps)]) / 2 # nolint: seq_linter.
    dat <- data.frame(swap_rate = swap_dist, temp = swap_idx)
    g <- ggplot(dat, aes(x = .data$temp, y = .data$swap_rate)) +
      geom_point() +
      geom_vline(data = data.frame(x = temps), aes(xintercept = .data$x), linetype = "dashed", alpha = 0.25) +
      coord_cartesian(ylim = c(0, 1))
    g
  })
  return(plots)
}

# ---- Label Switching Correction ----

#' Find optimal permutation between two population assignment matrices
#'
#' @details Finds the permutation of population labels that maximizes the
#' correlation between two assignment matrices. Uses a greedy matching algorithm
#' that finds the best matching for each population by maximizing correlation.
#'
#' @param assignments_t Matrix of population assignments at step t.
#'   Rows are samples, columns are populations. Values should be probabilities.
#' @param assignments_t1 Matrix of population assignments at step t+1.
#'   Same structure as assignments_t.
#' @return A permutation vector where \code{perm[i]} gives the population index in
#'   \code{assignments_t1} that corresponds to population \code{i} in \code{assignments_t}.
#'
#' @importFrom stats cor
#' @keywords internal
find_optimal_permutation <- function(assignments_t, assignments_t1) {
  num_populations <- ncol(assignments_t)
  
  # Compute correlation matrix between populations
  # cor[i, j] = correlation between population i in step t and population j in step t+1
  cor_matrix <- matrix(0, nrow = num_populations, ncol = num_populations)
  for (i in seq_len(num_populations)) {
    for (j in seq_len(num_populations)) {
      cor_matrix[i, j] <- cor(assignments_t[, i], assignments_t1[, j], use = "pairwise.complete.obs")
    }
  }
  
  # Replace NA with 0 (occurs when there's no variation)
  cor_matrix[is.na(cor_matrix)] <- 0
  
  # Greedy matching: find best permutation
  # Start with highest correlations and match them
  permutation <- rep(NA_integer_, num_populations)
  used_target <- rep(FALSE, num_populations)
  
  # Sort correlations in descending order
  cor_list <- list()
  for (i in seq_len(num_populations)) {
    for (j in seq_len(num_populations)) {
      cor_list[[length(cor_list) + 1]] <- list(
        source = i,
        target = j,
        cor = cor_matrix[i, j]
      )
    }
  }
  cor_list <- cor_list[order(sapply(cor_list, function(x) x$cor), decreasing = TRUE)]
  
  # Greedily assign matches
  for (item in cor_list) {
    if (is.na(permutation[item$source]) && !used_target[item$target]) {
      permutation[item$source] <- item$target
      used_target[item$target] <- TRUE
    }
  }
  
  # Handle any unmatched populations (shouldn't happen, but just in case)
  unmatched_source <- which(is.na(permutation))
  unmatched_target <- which(!used_target)
  if (length(unmatched_source) > 0 && length(unmatched_target) > 0) {
    for (i in seq_along(unmatched_source)) {
      if (i <= length(unmatched_target)) {
        permutation[unmatched_source[i]] <- unmatched_target[i]
      }
    }
  }
  
  return(permutation)
}

#' Apply permutation to population assignments
#'
#' @details Reorders the columns of an assignment matrix according to a permutation.
#'
#' @param assignments Matrix of population assignments. Rows are samples,
#'   columns are populations.
#' @param permutation Vector where \code{permutation[i]} gives the new index for
#'   population \code{i}. If NULL or identity, returns assignments unchanged.
#' @return Matrix with permuted columns.
#'
#' @keywords internal
apply_permutation <- function(assignments, permutation) {
  if (is.null(permutation) || all(permutation == seq_along(permutation))) {
    return(assignments)
  }
  
  # Reorder columns according to permutation
  # permutation[i] = j means: population i in original maps to population j in permuted
  # So we need the inverse: which column in original corresponds to each column in permuted
  num_populations <- length(permutation)
  inverse_permutation <- integer(num_populations)
  for (i in seq_len(num_populations)) {
    inverse_permutation[permutation[i]] <- i
  }
  
  return(assignments[, inverse_permutation, drop = FALSE])
}

#' Correct label switching in population assignments
#'
#' @details Corrects label switching by aligning population labels across MCMC steps.
#' Uses correlation between population assignment vectors to identify the optimal
#' permutation of labels for each step.
#'
#' @param population_assignments List of population assignment vectors, one per MCMC step.
#'   Each element should be a list of vectors (one per sample), where each vector contains
#'   log probabilities for each population.
#' @param method Alignment method: "iterative" (align each step to previous) or
#'   "reference" (align all steps to first step).
#' @return List of corrected population assignment vectors with same structure as input.
#'
#' @export
correct_label_switching <- function(population_assignments, method = "iterative") {
  if (length(population_assignments) == 0) {
    return(population_assignments)
  }
  
  num_samples <- length(population_assignments[[1]])
  if (num_samples == 0) {
    return(population_assignments)
  }
  
  num_steps <- length(population_assignments)
  num_populations <- length(population_assignments[[1]][[1]])
  
  # Convert log probabilities to probabilities for correlation computation
  # Build matrices for each step: rows = samples, columns = populations
  step_matrices <- lapply(population_assignments, function(step_assignments) {
    mat <- matrix(0, nrow = num_samples, ncol = num_populations)
    for (sample_idx in seq_len(num_samples)) {
      log_probs <- step_assignments[[sample_idx]]
      log_probs <- as.numeric(unlist(log_probs))  # ensure numeric vector (handles list format)
      # Handle -Inf values by setting to very small number before exp
      log_probs_finite <- pmax(log_probs, -700)  # exp(-700) is very small but finite
      probs <- exp(log_probs_finite)
      # Normalize to ensure probabilities sum to 1 (handles numerical issues)
      probs <- probs / sum(probs)
      mat[sample_idx, ] <- probs
    }
    return(mat)
  })
  
  # Apply label switching correction (only when we have 2+ steps)
  if (num_steps == 1) {
    corrected_matrices <- step_matrices
  } else if (method == "iterative") {
    # Align each step to the previous step
    corrected_matrices <- list(step_matrices[[1]])  # First step unchanged
    
    for (step_idx in 2:num_steps) {
      # Find optimal permutation to align step t to step t-1
      permutation <- find_optimal_permutation(
        corrected_matrices[[step_idx - 1]],
        step_matrices[[step_idx]]
      )
      
      # Apply permutation
      corrected_matrices[[step_idx]] <- apply_permutation(
        step_matrices[[step_idx]],
        permutation
      )
    }
  } else if (method == "reference") {
    # Align all steps to the first step
    reference_matrix <- step_matrices[[1]]
    corrected_matrices <- list(reference_matrix)
    
    for (step_idx in 2:num_steps) {
      # Find optimal permutation to align step t to reference (first step)
      permutation <- find_optimal_permutation(
        reference_matrix,
        step_matrices[[step_idx]]
      )
      
      # Apply permutation
      corrected_matrices[[step_idx]] <- apply_permutation(
        step_matrices[[step_idx]],
        permutation
      )
    }
  } else {
    stop("method must be 'iterative' or 'reference'")
  }
  
  # Convert back to original format (log probabilities)
  corrected_assignments <- lapply(corrected_matrices, function(mat) {
    lapply(seq_len(nrow(mat)), function(sample_idx) {
      probs <- mat[sample_idx, ]
      # Avoid log(0) by adding small epsilon
      probs <- pmax(probs, .Machine$double.eps)
      log_probs <- log(probs)
      # Normalize log probabilities (subtract log-sum-exp)
      log_sum_exp <- log(sum(exp(log_probs)))
      log_probs <- log_probs - log_sum_exp
      return(log_probs)
    })
  })
  
  return(corrected_assignments)
}
