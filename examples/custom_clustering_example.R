# Example: Multiple Populations with Correlated Allele Frequencies
# =================================================================

# This example demonstrates how to generate correlated allele frequency
# distributions for multiple populations and use them in the simulate_data function

devtools::load_all(".")
library(moire)
library(ggplot2)

# Function to generate correlated allele frequency distributions
# ============================================================

#' Generate correlated allele frequency distributions for multiple populations
#' 
#' @param num_loci Number of loci
#' @param num_alleles_per_locus Number of alleles per locus (can be vector)
#' @param num_populations Number of populations
#' @param correlation_strength Strength of correlation between populations (0-1)
#' @param base_alpha Base alpha parameter for Dirichlet distribution
#' @return List of allele frequency distributions for each population
generate_correlated_allele_frequencies <- function(num_loci, 
                                                  num_alleles_per_locus, 
                                                  num_populations, 
                                                  correlation_strength = 1,
                                                  base_alpha = .5) {
  
  # Ensure num_alleles_per_locus is a vector
  if (length(num_alleles_per_locus) == 1) {
    num_alleles_per_locus <- rep(num_alleles_per_locus, num_loci)
  }
  
  # Generate base allele frequencies for the first population
  base_frequencies <- list()
  for (locus_idx in 1:num_loci) {
    alpha_vec <- rep(base_alpha, num_alleles_per_locus[locus_idx])
    base_frequencies[[locus_idx]] <- simulate_allele_frequencies(alpha_vec, 1)[, 1]
  }
  
  # Generate correlated frequencies for other populations
  all_frequencies <- list()
  
  for (pop_idx in 1:num_populations) {
    pop_frequencies <- list()
    
    for (locus_idx in 1:num_loci) {
      num_alleles <- num_alleles_per_locus[locus_idx]
      
      if (pop_idx == 1) {
        # Use base frequencies for first population
        pop_frequencies[[locus_idx]] <- base_frequencies[[locus_idx]]
      } else {
        # Generate correlated frequencies
        # Method: Mix base frequencies with random noise
        base_freq <- base_frequencies[[locus_idx]]
        
        # # Generate random frequencies
        # alpha_vec <- rep(base_alpha, num_alleles)
        # random_freq <- simulate_allele_frequencies(alpha_vec, 1)[, 1]
        
        # # Mix base and random frequencies based on correlation strength
        # mixed_freq <- correlation_strength * base_freq + (1 - correlation_strength) * random_freq

        corr_freq <- simulate_allele_frequencies(base_freq * correlation_strength, 1)
        
        # Ensure frequencies sum to 1
        pop_frequencies[[locus_idx]] <- corr_freq 
      }
    }
    
    all_frequencies[[pop_idx]] <- pop_frequencies
  }
  
  return(all_frequencies)
}

# Function to simulate data for multiple populations
# ==================================================

#' Simulate data for multiple populations with correlated allele frequencies
#' 
#' @param num_populations Number of populations
#' @param samples_per_population Number of samples per population
#' @param num_loci Number of loci
#' @param num_alleles_per_locus Number of alleles per locus
#' @param correlation_strength Correlation strength between populations
#' @param mean_coi Mean multiplicity of infection
#' @param epsilon_pos False positive rate
#' @param epsilon_neg False negative rate
#' @param missingness Missing data rate
#' @return Combined data from all populations
simulate_multi_population_data <- function(num_populations = 3,
                                          samples_per_population = 50,
                                          num_loci = 20,
                                          num_alleles_per_locus = 10,
                                          correlation_strength = 1,
                                          mean_coi = 3,
                                          epsilon_pos = 0.01,
                                          epsilon_neg = 0.1,
                                          missingness = 0) {
  
  # Generate correlated allele frequencies
  cat("Generating correlated allele frequencies for", num_populations, "populations...\n")
  correlated_frequencies <- generate_correlated_allele_frequencies(
    num_loci = num_loci,
    num_alleles_per_locus = num_alleles_per_locus,
    num_populations = num_populations,
    correlation_strength = correlation_strength
  )
  
  # Simulate data for each population
  population_data <- list()
  
  for (pop_idx in 1:num_populations) {
    cat("Simulating data for population", pop_idx, "...\n")
    
    # Convert frequencies to the format expected by simulate_data
    allele_freqs <- list()
    for (locus_idx in 1:num_loci) {
      locus_name <- paste0("L", locus_idx)
      # Handle both single value and vector for num_alleles_per_locus
      if (length(num_alleles_per_locus) == 1) {
        num_alleles <- num_alleles_per_locus
      } else {
        num_alleles <- num_alleles_per_locus[locus_idx]
      }
      allele_names <- paste0(locus_name, "_", 1:num_alleles)
      allele_freqs[[locus_name]] <- correlated_frequencies[[pop_idx]][[locus_idx]]
      names(allele_freqs[[locus_name]]) <- allele_names
    }
    
    # Simulate data for this population
    pop_data <- simulate_data(
      mean_coi = mean_coi,
      num_samples = samples_per_population,
      epsilon_pos = epsilon_pos,
      epsilon_neg = epsilon_neg,
      allele_freqs = allele_freqs,
      internal_relatedness_alpha = 1,
      internal_relatedness_beta = 9,
      missingness = missingness
    )
    
    # Add population labels to sample IDs
    pop_data$sample_ids <- paste0("Pop", pop_idx, "_", pop_data$sample_ids)
    pop_data$population <- rep(pop_idx, samples_per_population)
    
    population_data[[pop_idx]] <- pop_data
  }
  
  # Combine all population data
  cat("Combining data from all populations...\n")
  combined_data <- combine_simulated_data(population_data[[1]], population_data[[2]])
  
  if (num_populations > 2) {
    for (pop_idx in 3:num_populations) {
      combined_data <- combine_simulated_data(combined_data, population_data[[pop_idx]])
    }
  }
  
  # Add population information
  combined_data$population_assignments <- unlist(lapply(population_data, function(x) x$population))
  combined_data$original_populations <- population_data
  
  return(combined_data)
}

# Generate the multi-population data
# ==================================
# Generate multi-population data with correlated allele frequencies
data <- simulate_multi_population_data(
  num_populations = 3,
  samples_per_population = 100,
  num_loci = 10,
  num_alleles_per_locus = 10,
  correlation_strength = 4,
  mean_coi = 3,
  epsilon_pos = 0.01,
  epsilon_neg = 0.1,
  missingness = 0
)
readr::write_rds(data, "examples/data.rds")
data <- readr::read_rds("examples/data.rds")

# Display summary of the generated data
cat("=== Multi-Population Data Summary ===\n")
cat("Total samples:", length(data$sample_ids), "\n")
cat("Number of loci:", length(data$loci), "\n")
cat("Population distribution:\n")
print(table(data$population_assignments))

# Visualize correlation between populations
# =========================================

# Function to calculate allele frequency correlation between populations
calculate_frequency_correlation <- function(pop1_freqs, pop2_freqs) {
  correlations <- numeric(length(pop1_freqs))
  
  for (locus_idx in seq_along(pop1_freqs)) {
    freq1 <- pop1_freqs[[locus_idx]]
    freq2 <- pop2_freqs[[locus_idx]]
    
    if (length(freq1) == length(freq2)) {
      correlations[locus_idx] <- cor(freq1, freq2, use = "complete.obs")
    } else {
      correlations[locus_idx] <- NA
    }
  }
  
  return(correlations)
}

# Calculate correlations between populations
cat("\n=== Allele Frequency Correlations ===\n")
for (pop1 in 1:2) {
  for (pop2 in (pop1+1):3) {
    corr_values <- calculate_frequency_correlation(
      data$original_populations[[pop1]]$allele_freqs,
      data$original_populations[[pop2]]$allele_freqs
    )
    mean_corr <- mean(corr_values, na.rm = TRUE)
    cat("Pop", pop1, "vs Pop", pop2, "correlation:", round(mean_corr, 3), "\n")
  }
}

# Clustering Analysis on Multi-Population Data
# =============================================

# Method 1: Hierarchical Clustering with Jaccard Distance
# --------------------------------------------------------

# Calculate pairwise Jaccard similarity matrix
calculate_jaccard_similarity <- function(data) {
  num_samples <- length(data$sample_ids)
  similarity_matrix <- matrix(0, nrow = num_samples, ncol = num_samples)
  
  for (i in 1:num_samples) {
    for (j in 1:num_samples) {
      if (i != j) {
        # Calculate Jaccard similarity across all loci
        total_shared <- 0
        total_possible <- 0
        for (locus_idx in seq_along(data$loci)) {
          alleles_i <- data$data[[locus_idx]][[i]]
          alleles_j <- data$data[[locus_idx]][[j]]
          
          # Convert to binary presence/absence
          present_i <- alleles_i > 0
          present_j <- alleles_j > 0

          
          # Jaccard similarity = intersection / union
          intersection <- sum(present_i & present_j)
          union_size <- sum(present_i | present_j)

          if (i == 1 && j == 2) {
            cat("alleles_i: ", alleles_i, "\n")
            cat("alleles_j: ", alleles_j, "\n")
            cat("present_i: ", present_i, "\n")
            cat("present_j: ", present_j, "\n")
            cat("intersection: ", intersection, "\n")
            cat("union_size: ", union_size, "\n")
          }
          if (union_size > 0) {
            total_shared <- total_shared + intersection
            total_possible <- total_possible + union_size
          }
        }
        similarity_matrix[i, j] <- total_shared / total_possible
      } else {
        similarity_matrix[i, i] <- 1.0
      }
    }
  }
  
  return(similarity_matrix)
}


# Calculate similarity matrix
similarity_matrix <- calculate_jaccard_similarity(data)

# Convert to distance matrix
distance_matrix <- 1 - similarity_matrix

# Perform hierarchical clustering
hc_result <- hclust(as.dist(distance_matrix), method = "ward.D2")

# reorder the similarity matrix
similarity_matrix_ord <- similarity_matrix[order(data$population_assignments), order(data$population_assignments)]
similarity_matrix_hclust <- similarity_matrix[order(hc_result$order), order(hc_result$order)]

# plot the similarity matrix
image(similarity_matrix_ord)
image(similarity_matrix_hclust)

# Plot the dendrogram with tips colored by true population assignments
original_population_assignments <- data$population_assignments
cols <- c("red", "blue", "green")

# Convert to dendrogram object to get the proper ordering
dend <- as.dendrogram(hc_result)

# Get the order of samples in the dendrogram
dend_order <- order.dendrogram(dend)

# Create colors in the correct order for the dendrogram
tip_colors <- cols[original_population_assignments[dend_order]]

# Create the dendrogram plot
plot(dend, main = "Hierarchical Clustering with Jaccard Distance", )

# Add colored points at the tips in the correct positions
n_samples <- length(original_population_assignments)
for (i in 1:n_samples) {
  points(i, 0, col = tip_colors[i], pch = 19, cex = 1.2)
}

# Add legend
legend("topright", 
       legend = paste("Population", 1:3), 
       col = cols, 
       pch = 19, 
       title = "True Populations")


# Cut tree to get 5 populations
num_populations <- 5
population_assignments <- cutree(hc_result, k = num_populations)

# Compare clustering results with true population assignments
cat("\n=== Clustering Performance Analysis ===\n")
true_assignments <- data$population_assignments
clustered_assignments <- population_assignments

# Create confusion matrix
confusion_matrix <- table(True = true_assignments, Clustered = clustered_assignments)
cat("Confusion Matrix (True vs Clustered):\n")
print(confusion_matrix)

# Calculate clustering accuracy
accuracy <- sum(diag(confusion_matrix)) / sum(confusion_matrix)
cat("Clustering Accuracy:", round(accuracy, 3), "\n")

# Prepare initial allele frequencies
initial_freqs <- prepare_initial_allele_frequencies(
  data = data,
  population_assignments = population_assignments,
  num_populations = num_populations,
  pseudocount = 1
)

# Run MCMC with custom initialization
cat("\n=== Running MCMC with Clustering Initialization ===\n")
devtools::load_all(".")
start_profiler("mcmc")
mcmc_result2 <- run_mcmc(
  data = data,
  num_populations = num_populations,
  initial_allele_frequencies = initial_freqs,
  burnin = 50,
  r_alpha = 5,
  r_beta = 5,
  samples_per_chain = 50,
  verbose = TRUE
)
stop_profiler()


readr::write_rds(mcmc_result2, "examples/multi_pop_mcmc_result.rds")
mcmc_result2 <- readr::read_rds("examples/multi_pop_mcmc_result.rds")

mcmc_result_single_pop <- run_mcmc(
  data = data,
  num_populations = 1,
  initial_allele_frequencies = initial_freqs,
  r_alpha = 5,
  r_beta = 5,
  burnin = 100,
  samples_per_chain = 1000,
  verbose = TRUE
)
readr::write_rds(mcmc_result_single_pop, "examples/single_pop_mcmc_result.rds")
mcmc_result_single_pop <- readr::read_rds("examples/single_pop_mcmc_result.rds")


ecoi_summary_single_pop <- moire::summarize_effective_coi(mcmc_result_single_pop) |>
  dplyr::select(sample_id, single_pop_post_ecoi_med = post_effective_coi_med, single_pop_post_ecoi_lower = post_effective_coi_lower, single_pop_post_ecoi_upper = post_effective_coi_upper)
ecoi_summary_multi_pop <- moire::summarize_effective_coi(mcmc_result2) |>
  dplyr::select(sample_id, multi_pop_post_ecoi_med = post_effective_coi_med, multi_pop_post_ecoi_lower = post_effective_coi_lower, multi_pop_post_ecoi_upper = post_effective_coi_upper)

ecoi_summary <- data.frame(
  true_ecoi = c(
    1 + (data$original_populations[[1]]$sample_cois - 1) * (1 - data$original_populations[[1]]$sample_relatedness),
    1 + (data$original_populations[[2]]$sample_cois - 1) * (1 - data$original_populations[[2]]$sample_relatedness),
    1 + (data$original_populations[[3]]$sample_cois - 1) * (1 - data$original_populations[[3]]$sample_relatedness)
  ),
  sample_id = data$sample_ids
) |>
  dplyr::left_join(ecoi_summary_single_pop, by = "sample_id") |>
  dplyr::left_join(ecoi_summary_multi_pop, by = "sample_id") |>
  tidyr::pivot_longer(cols = c(single_pop_post_ecoi_med, multi_pop_post_ecoi_med), names_to = "model", values_to = "ecoi") |>
  dplyr::mutate(model = ifelse(model == "single_pop_post_ecoi_med", "Single population", "Multi-population"))

ggplot(ecoi_summary) +
  geom_point(aes(x = true_ecoi, y = ecoi)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_smooth(aes(x = true_ecoi, y = ecoi), method = "lm", se = FALSE) +
  theme_minimal() +
  labs(x = "True effective COI", y = "Estimated effective COI") +
  facet_wrap(~model)

pop_summary <- moire::summarize_population_assignments(mcmc_result2)

inferred_population_assignments <- pop_summary$most_likely_population |>
  dplyr::group_by(sample_id) |>
  dplyr::summarize(inferred_population = population[which.max(prob_most_likely)]) |>
  dplyr::left_join(pop_summary$entropy_summary, by = "sample_id")

long_form_data <- moire::convert_to_long_form(data) |>
  dplyr::left_join(tibble::tibble(sample_id = data$sample_ids, population = data$population_assignments)) |>
  dplyr::left_join(inferred_population_assignments, by = "sample_id")


source("examples/parallel_ibd.r")
dsmp <- dcifer::formatDat(long_form_data, "sample_id", "locus", "allele")
coi <- dcifer::getCOI(dsmp, lrank = 2)
afreq <- dcifer::calcAfreq(dsmp, coi, tol = 1e-5)
combined_dcifer_res <- parallel_ibdDat(dsmp, coi, afreq) |>
  dplyr::rename(
    ibd_est_single_pop = estimate, 
    ibd_p_value_single_pop = p_value, 
    ibd_CI_lower_single_pop = CI_lower, 
    ibd_CI_upper_single_pop = CI_upper
  )

populations <- data$population_assignments |> unique()
true_ibd_res <- list()
for(pop in 1:length(populations)) {
  print(paste("Running DCIFER for population", pop))
  pop_long_form <- long_form_data |>
    dplyr::filter(population == pop)
  dsmp <- dcifer::formatDat(pop_long_form, "sample_id", "locus", "allele")
  coi <- dcifer::getCOI(dsmp, lrank = 2)
  afreq <- dcifer::calcAfreq(dsmp, coi, tol = 1e-5)
  dcifer_res <- parallel_ibdDat(dsmp, coi, afreq)
  dcifer_res$population <- pop
  true_ibd_res[[pop]] <- dcifer_res
}
true_ibd_res_df <- do.call(rbind, true_ibd_res) |>
  dplyr::rename(
    ibd_est_true_pop = estimate, 
    ibd_p_value_true_pop = p_value, 
    ibd_CI_lower_true_pop = CI_lower, 
    ibd_CI_upper_true_pop = CI_upper
  )

inferred_populations <- long_form_data$inferred_population |> unique()
inferred_ibd_res <- list()
for(pop in inferred_populations) {
  print(paste("Running DCIFER for population", pop))
  pop_long_form <- long_form_data |>
    dplyr::filter(inferred_population == pop)
  print(paste("Number of samples in population", pop, ":", pop_long_form$sample_id |> unique() |> length()))
  if (pop_long_form$sample_id |> unique() |> length() < 5) {
    next
  }
  dsmp <- dcifer::formatDat(pop_long_form, "sample_id", "locus", "allele")
  coi <- dcifer::getCOI(dsmp, lrank = 2)
  afreq <- dcifer::calcAfreq(dsmp, coi, tol = 1e-5)
  dcifer_res <- parallel_ibdDat(dsmp, coi, afreq)
  dcifer_res$inferred_population <- pop
  inferred_ibd_res[[pop]] <- dcifer_res
}
inferred_ibd_res_df <- do.call(rbind, inferred_ibd_res) |>
  dplyr::rename(
    ibd_est_inferred_pop = estimate,
    ibd_p_value_inferred_pop = p_value,
    ibd_CI_lower_inferred_pop = CI_lower,
    ibd_CI_upper_inferred_pop = CI_upper
  )


ibd_res_df <- true_ibd_res_df |>
  dplyr::left_join(combined_dcifer_res, by = c("x", "y")) |>
  dplyr::left_join(inferred_ibd_res_df, by = c("x", "y")) |>
  dplyr::mutate(
    inferred_sig_misspec = ibd_p_value_inferred_pop < 0.05 & ibd_p_value_true_pop > 0.05,
    single_pop_sig_misspec = ibd_p_value_single_pop < 0.05 & ibd_p_value_true_pop > 0.05
  )


ggplot(ibd_res_df |> dplyr::filter(!is.na(inferred_sig_misspec), !is.na(inferred_population))) +
  geom_point(aes(x = ibd_est_inferred_pop, y = ibd_est_true_pop, color = inferred_sig_misspec)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", alpha = 0.5) +
  facet_wrap(~inferred_population) +
  theme_minimal() +
  labs(x = "Inferred IBD", y = "True IBD", color = "Sig. Misspec")

ggplot(ibd_res_df) +
  geom_point(aes(x = ibd_est_single_pop, y = ibd_est_true_pop, color = single_pop_sig_misspec)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", alpha = 0.5) +
  facet_wrap(~population) +
  theme_minimal() +
  labs(x = "Inferred IBD", y = "True IBD", color = "Sig. Misspec")

combined_sample_ids <- c(combined_dcifer_res$x, combined_dcifer_res$y) |> unique(na.rm = TRUE) 
combined_ibd_matrix <- matrix(NA, nrow = length(combined_sample_ids), ncol = length(combined_sample_ids))
rownames(combined_ibd_matrix) <- combined_sample_ids
colnames(combined_ibd_matrix) <- combined_sample_ids
for(i in 1:nrow(combined_dcifer_res)) {
    x <- combined_dcifer_res$x[i]
    y <- combined_dcifer_res$y[i]
    combined_ibd_matrix[x, y] <- combined_dcifer_res$ibd_est_single_pop[i]
    combined_ibd_matrix[y, x] <- combined_dcifer_res$ibd_est_single_pop[i]
}

mds <- cmdscale(dist(combined_ibd_matrix), k = 2)
mds_df <- tibble::tibble(x = mds[, 1], y = mds[, 2], sample_id = rownames(combined_ibd_matrix)) |>
  dplyr::left_join(tibble::tibble(sample_id = data$sample_ids, population = data$population_assignments), by = "sample_id") |>
  dplyr::left_join(inferred_population_assignments, by = "sample_id")
ggplot(mds_df, aes(x = x, y = y, label = sample_id)) +
  geom_point(aes(color = as.factor(inferred_population), size = entropy_mean), alpha = 0.7) +
  theme_minimal() +
  labs(x = "MDS1", y = "MDS2", color = "Inferred Population", size = "Entropy")


ggplot(mds_df, aes(x = x, y = y, label = sample_id)) +
  geom_point(aes(color = as.factor(population)), alpha = 0.7) +
  theme_minimal() +
  labs(x = "MDS1", y = "MDS2", color = "True Population")

























# Method 2: K-means Clustering
# ----------------------------

# Convert data to matrix format for k-means
# (This is a simplified approach - you might want to use more sophisticated methods)
convert_to_matrix <- function(data) {
  num_samples <- length(data$sample_ids)
  num_loci <- length(data$loci)
  
  # Create a binary matrix: samples x (loci * alleles)
  matrix_data <- matrix(0, nrow = num_samples, ncol = 0)
  
  for (locus_idx in 1:num_loci) {
    num_alleles <- length(data$data[[locus_idx]][[1]])
    for (allele_idx in 1:num_alleles) {
      column <- numeric(num_samples)
      for (sample_idx in 1:num_samples) {
        if (sample_idx <= length(data$data[[locus_idx]])) {
          column[sample_idx] <- data$data[[locus_idx]][[sample_idx]][allele_idx] > 0
        }
      }
      matrix_data <- cbind(matrix_data, column)
    }
  }
  
  return(matrix_data)
}

# Convert data and perform k-means
matrix_data <- convert_to_matrix(data)
kmeans_result <- kmeans(matrix_data, centers = num_populations, nstart = 20)
population_assignments_kmeans <- kmeans_result$cluster

# Prepare initial allele frequencies for k-means result
initial_freqs_kmeans <- prepare_initial_allele_frequencies(
  data = data,
  population_assignments = population_assignments_kmeans,
  num_populations = num_populations,
  pseudocount = 1
)

# Method 3: Using External Clustering Libraries
# ---------------------------------------------

# Example using PAM (Partitioning Around Medoids) from cluster package
if (require(cluster, quietly = TRUE)) {
  pam_result <- pam(distance_matrix, k = num_populations)
  population_assignments_pam <- pam_result$clustering
  
  initial_freqs_pam <- prepare_initial_allele_frequencies(
    data = data,
    population_assignments = population_assignments_pam,
    num_populations = num_populations,
    pseudocount = 1
  )
}

# Method 4: Custom Distance Metrics
# ----------------------------------

# Example using a custom distance metric based on allele sharing
calculate_allele_sharing_distance <- function(data) {
  num_samples <- length(data$sample_ids)
  distance_matrix <- matrix(0, nrow = num_samples, ncol = num_samples)
  
  for (i in 1:num_samples) {
    for (j in 1:num_samples) {
      if (i != j) {
        total_shared <- 0
        total_possible <- 0
        
        for (locus_idx in seq_along(data$loci)) {
          alleles_i <- data$data[[locus_idx]][[i]]
          alleles_j <- data$data[[locus_idx]][[j]]
          
          # Count shared alleles (minimum of the two counts)
          shared <- sum(pmin(alleles_i, alleles_j))
          possible <- sum(pmax(alleles_i, alleles_j))
          
          total_shared <- total_shared + shared
          total_possible <- total_possible + possible
        }
        
        if (total_possible > 0) {
          distance_matrix[i, j] <- 1 - (total_shared / total_possible)
        } else {
          distance_matrix[i, j] <- 1
        }
      }
    }
  }
  
  return(distance_matrix)
}

# Use custom distance metric
custom_distance <- calculate_allele_sharing_distance(data)
hc_custom <- hclust(as.dist(custom_distance), method = "ward.D2")
population_assignments_custom <- cutree(hc_custom, k = num_populations)

initial_freqs_custom <- prepare_initial_allele_frequencies(
  data = data,
  population_assignments = population_assignments_custom,
  num_populations = num_populations,
  pseudocount = 1
)

# Compare different clustering approaches
cat("\n=== Clustering Method Comparison ===\n")
cat("True population distribution:", table(true_assignments), "\n")
cat("Hierarchical clustering (Jaccard):", table(population_assignments), "\n")
cat("K-means clustering:", table(population_assignments_kmeans), "\n")
if (exists("population_assignments_pam")) {
  cat("PAM clustering:", table(population_assignments_pam), "\n")
}
cat("Custom distance clustering:", table(population_assignments_custom), "\n")

# Final Summary and Analysis
# ==========================

cat("\n=== Final Analysis Summary ===\n")
cat("1. Generated", num_populations, "populations with", length(data$sample_ids), "total samples\n")
cat("2. Average correlation between populations:", round(mean_corr, 3), "\n")
cat("3. Clustering accuracy:", round(accuracy, 3), "\n")
cat("4. MCMC completed with both clustering and true population initializations\n")

# Optional: Create a simple visualization of population structure
if (require(ggplot2, quietly = TRUE)) {
  # Create a data frame for plotting
  plot_data <- data.frame(
    Sample = seq_along(data$sample_ids),
    True_Population = as.factor(true_assignments),
    Clustered_Population = as.factor(population_assignments),
    Sample_ID = data$sample_ids
  )
  
  # Plot true vs clustered assignments
  p1 <- ggplot(plot_data, aes(x = Sample, y = True_Population, color = True_Population)) +
    geom_point(size = 2) +
    labs(title = "True Population Assignments", y = "Population") +
    theme_minimal()
  
  p2 <- ggplot(plot_data, aes(x = Sample, y = Clustered_Population, color = Clustered_Population)) +
    geom_point(size = 2) +
    labs(title = "Clustered Population Assignments", y = "Population") +
    theme_minimal()
  
  # Print plots
  print(p1)
  print(p2)
}

cat("\n=== Example Complete ===\n")
cat("This example demonstrates:\n")
cat("- Generation of correlated allele frequency distributions\n")
cat("- Multi-population data simulation\n")
cat("- Various clustering methods for population structure inference\n")
cat("- MCMC analysis with different initialization strategies\n")
