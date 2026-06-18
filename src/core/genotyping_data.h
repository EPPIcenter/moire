#pragma once

#ifndef DATA_H_
#define DATA_H_

#include "multivector.h"
#include "observation_model_kind.h"

#include <Rcpp.h>
#include <string>
#include <vector>
#include <span>


//------------------------------------------------
// class containing multi allelic genotyping data
class GenotypingData
{
   public:
    // indexed by sample, locus, allele
    RaggedMultiVector<int, 3> observed_alleles;
    // indexed by sample, locus
    MultiVector<int, 2> is_missing_;
    // indexed by locus
    std::vector<std::size_t> num_alleles;
    // indexed by sample
    std::vector<std::size_t> observed_coi;
    std::size_t num_samples;
    std::size_t num_loci;
    std::string aggregate;

    MultiVector<float, 2> jaccard_similarity_matrix;

    // constructors
    GenotypingData(){};
    GenotypingData(const Rcpp::List &args);

    std::span<int const> get_observed_alleles(std::size_t sample, std::size_t locus) const;
    bool is_missing(std::size_t sample, std::size_t locus) const;
    bool has_count_barcodes() const;
    void validate_for_observation_model(ObservationModelKind observation_model_kind) const;
};

#endif  // DATA_H_
