#pragma once

#ifndef DATA_H_
#define DATA_H_

#include <Rcpp.h>
#include <vector>
#include <span>

//------------------------------------------------
// class containing multi allelic genotyping data
class GenotypingData
{
   public:
    // Data are ordered by Locus, then Sample
    std::vector<std::vector<std::vector<int>>> observed_alleles;
    std::vector<std::vector<bool>> is_missing_;
    std::vector<std::size_t> num_alleles;
    std::vector<std::size_t> observed_coi;
    std::size_t num_samples;
    std::size_t num_loci;

    // constructors
    GenotypingData(){};
    GenotypingData(const Rcpp::List &args);

    std::span<int const> get_observed_alleles(int locus, int sample) const;
    bool is_missing(int locus, int sample) const;
};

#endif  // DATA_H_
