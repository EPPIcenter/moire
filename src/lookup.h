#pragma once

#ifndef LOOKUP_H_
#define LOOKUP_H_

#include "parameters.h"

#include <Rcpp.h>
#include <map>
#include <utility>
#include <vector>

class Lookup
{
    using SamplingDepthEntry = std::pair<int, int>;

   private:
    bool lgamma_initialized;
    bool sample_depth_initialized;
    int max_alleles;

   public:
    Lookup(int max_alleles);
    std::vector<double> lookup_lgamma;
    std::map<SamplingDepthEntry, double> lookup_sampling_depth{};

    void init_lgamma();
    void init_sampling_depth();

    double get_sampling_depth(int coi, int num_alleles);
};

#endif  // LOOKUP_H_
