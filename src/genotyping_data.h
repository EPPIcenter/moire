#pragma once

#ifndef DATA_H_
#define DATA_H_

#include <Rcpp.h>
#include <vector>

//------------------------------------------------
// class containing multi allelic genotyping data
class GenotypingData
{
   public:
    // Data are ordered by Locus, then Sample
    static std::vector<std::vector<std::vector<int>>> observed_alleles;
    static std::vector<std::vector<bool>> is_missing_;
    static std::vector<int> num_alleles;
    static std::vector<int> observed_coi;
    static std::size_t num_samples;
    static std::size_t num_loci;
    static int max_alleles;

    // constructors
    GenotypingData(){};
    GenotypingData(const Rcpp::List &args);

    const std::vector<int> &get_observed_alleles(int locus, int sample) const;
    bool is_missing(int locus, int sample) const;
};

#endif  // DATA_H_
