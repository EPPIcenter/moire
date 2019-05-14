#pragma once

#ifndef DATA_H_
#define DATA_H_

#include <Rcpp.h>


//------------------------------------------------
// class containing multi allelic genotyping data
class GenotypingData {
public:
  // Data are ordered by Locus, then Sample
  static std::vector<std::vector<std::vector<int> > > observed_alleles;
  
  static std::vector<int> num_alleles;
  static std::vector<int> observed_coi;
  static size_t num_samples;
  static size_t num_loci;
  static int max_alleles;
  
  // constructors
  GenotypingData() {};
  GenotypingData(const Rcpp::List &args);

};

#endif // DATA_H_
