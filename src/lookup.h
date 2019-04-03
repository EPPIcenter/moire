#pragma once

#ifndef LOOKUP_H_
#define LOOKUP_H_

#include <Rcpp.h>
#include "parameters.h"

class Lookup {
    private:
    static bool lgamma_initialized;
    static bool sample_depth_initialized;
    
    public:
    Lookup() {};
    Lookup(int max_coi, int max_alleles);
    static std::vector<double> lookup_lgamma;
    static std::vector<std::vector<int > > lookup_sampling_depth;

    void init_lgamma(int max_coi);
    void init_sampling_depth(int max_coi, int max_alleles);
};


#endif // LOOKUP_H_