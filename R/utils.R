### Convert alleles list strings to barcode datastructure, requires the set
### of all alleles observed a each locus
.convert_to_barcode <- function(alleles, loci, all_alleles_list) {
  alleles_lists <- str_split(alleles, ';')
  barcodes <- list()
  for (i in 1:length(alleles_lists)) {
    alleles_list <- unlist(alleles_lists[i])
    locus <- loci[i]
    all_alleles <- unlist(all_alleles_list[[locus]])
    barcodes[[i]] = rep(0, length(all_alleles))
    for (j in 1:length(alleles_list)) {
      allele <- alleles_list[j]
      if (!is.na(allele)) {
        barcodes[[i]] = barcodes[[i]] + as.integer(all_alleles == allele)
      }
    }
  }
  barcodes
}


#' Estimates the complexity of infection using a naive approach that chooses the n'th highest number of observed alleles.
#'
#' @export
#'
#' @param data List of lists of numeric vectors, where each list element is a collection of observations across samples at a single genetic locus.
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
  sapply(num_alleles_by_sample, function(sample) {sort(unlist(sample), decreasing = T)[offset]})
}


#' Estimates the complexity of infection using a naive approach that chooses the highest number of observed alleles.
#'
#' @export
#'
#' @param data List of lists of numeric vectors, where each list element is a collection of observations across samples at a single genetic locus.
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
  sapply(num_alleles_by_sample, function(sample) {max(unlist(sample))})
}


#' Calculate the expected heterozygosity from allele frequencies
#'
#' @export
#'
#' @param allele_freqs Simplex of allele frequencies
calculate_he <- function(allele_freqs) {
  return(1 - sum(allele_freqs**2))
}
