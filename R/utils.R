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
#' data <- simulate_data(num_samples = 100, num_loci = 10)
#' # Use custom clustering (e.g., hierarchical clustering)
#' # population_assignments <- cutree(hclust(dist_matrix), k = 3)
#' # initial_freqs <- prepare_initial_allele_frequencies(data, population_assignments, 3)
prepare_initial_allele_frequencies <- function(data, population_assignments, num_populations, pseudocount = 1) {
  cat("=== Preparing Initial Allele Frequencies ===\n")
  cat("Number of samples:", length(data$sample_ids), "\n")
  cat("Number of loci:", length(data$loci), "\n")
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
  cat("=== Converting Data to Long Form ===\n")
  cat("Number of samples:", length(data$sample_ids), "\n")
  cat("Number of loci:", length(data$loci), "\n")
  
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
    temps <- temps <- mcmc_results$chains[[1]]$temp_gradient
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
