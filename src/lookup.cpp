
#include <algorithm>

#include "lookup.h"
#include "mcmc_utils.h"

std::vector<double> Lookup::lookup_lgamma;
std::vector<std::vector<int > > Lookup::lookup_sampling_depth;
bool Lookup::lgamma_initialized = false;
bool Lookup::sample_depth_initialized = false;

Lookup::Lookup(int max_coi, int max_alleles) {
    init_lgamma(max_coi);
    init_sampling_depth(max_coi, max_alleles);
};


void Lookup::init_lgamma(int max_coi) {
    if (!lgamma_initialized) {
        lookup_lgamma = std::vector<double>(max_coi + 2);
        
        for(size_t i = 0; i < lookup_lgamma.size(); i++) {
            lookup_lgamma[i] = std::lgamma(i);
        }

        lgamma_initialized = true;
    }
};

void Lookup::init_sampling_depth(int max_coi, int max_alleles) {
    if (!sample_depth_initialized) {
        for(int i = 0; i <= max_coi + 1; i++) {
            lookup_sampling_depth.push_back(std::vector<int>(max_alleles));
            for(int j = 0; j <= max_alleles + 1; j++) {
                if(i == 0 || j == 0) {
                    lookup_sampling_depth[i][j] = 0;
                } else {
                    lookup_sampling_depth[i][j] = std::min(1e8, tgamma(i + j) / (tgamma(j - 1) * tgamma(i)));
                    UtilFunctions::print(i, j, lookup_sampling_depth[i][j]);
                }
            }
            
        }
        
        sample_depth_initialized = true;
    }
}