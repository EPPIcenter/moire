
#include "lookup.h"
#include "mcmc_utils.h"

std::vector<double> Lookup::lookup_lgamma;
bool Lookup::lgamma_initialized = false;

Lookup::Lookup(int max_coi) {
    init_lgamma(max_coi);
};


void Lookup::init_lgamma(int max_coi) {
    if(!lgamma_initialized) {
        lookup_lgamma = std::vector<double>(max_coi + 2);
        
        for(size_t i = 0; i < lookup_lgamma.size(); i++) {
            lookup_lgamma[i] = std::lgamma(i);
        }

        lgamma_initialized = true;
    }
};