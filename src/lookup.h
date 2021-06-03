#pragma once

#ifndef LOOKUP_H_
#define LOOKUP_H_

#include "parameters.h"

#include <Rcpp.h>

class Lookup
{
   private:
    bool lgamma_initialized;
    bool sample_depth_initialized;
    int max_coi;
    int max_alleles;

   public:
    Lookup(){};
    Lookup(int max_coi, int max_alleles);
    std::vector<double> lookup_lgamma;
    std::vector<double> lookup_sampling_depth;

    void init_lgamma();
    void init_sampling_depth();

    long get_sampling_depth(int coi, int num_alleles);
};

#endif  // LOOKUP_H_
