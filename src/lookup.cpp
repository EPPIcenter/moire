#include <algorithm>

#include "lookup.h"
#include "mcmc_utils.h"

std::vector<double> Lookup::lookup_lgamma;
std::vector<std::vector<int > > Lookup::lookup_sampling_depth;
bool Lookup::lgamma_initialized = false;

Lookup::Lookup(int max_coi, int max_alleles) {
    init_lgamma(max_coi, max_alleles);
};


void Lookup::init_lgamma(int max_coi, int max_alleles) {
    
    lookup_lgamma = std::vector<double>(max_coi + max_alleles + 10);
    
    for(size_t i = 0; i < lookup_lgamma.size(); i++) {
        lookup_lgamma[i] = std::lgamma(i);
    }

};
