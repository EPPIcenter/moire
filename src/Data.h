#pragma once

#include <Rcpp.h>


//------------------------------------------------
// class containing multi-allelic data
class Data_multiallelic {

public:

  static std::vector<std::vector<std::vector<int8_t> > > data;
  static std::vector<int> observed_COI;
  static std::vector<int> alleles;
  static int n;
  static int L;

  // constructors
  Data_multiallelic() {};
  Data_multiallelic(const Rcpp::List &args);

};
