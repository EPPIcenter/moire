#include "lookup.h"

#include "mcmc_utils.h"

#include <cmath>
#include <algorithm>

Lookup::Lookup(int max_alleles) : max_alleles(max_alleles)
{
    init_lgamma();
    init_sampling_depth();
};

void Lookup::init_lgamma()
{
    lookup_lgamma = std::vector<double>(max_alleles + 2);

    for (size_t i = 0; i < lookup_lgamma.size(); i++)
    {
        lookup_lgamma[i] = std::lgamma(i);
    }
    lgamma_initialized = true;
};

void Lookup::init_sampling_depth()
{
    for (int ii = 1; ii <= max_alleles; ii++)
    {
        for (int jj = 1; jj <= ii; jj++)
        {
            SamplingDepthEntry key{ii, jj};
            lookup_sampling_depth[key] = 0;
            for (int kk = jj; kk > 0; --kk)
            {
                lookup_sampling_depth[key] += std::exp(lookup_lgamma[ii + 1] -
                                              lookup_lgamma[kk + 1] -
                                              lookup_lgamma[(ii - kk) + 1]);
            }
            lookup_sampling_depth[key] = std::log(lookup_sampling_depth[key]);
            UtilFunctions::print("Sampling Depth:", ii, jj, lookup_sampling_depth[key]);
        }
    }
}

double Lookup::get_sampling_depth(int coi, int num_alleles)
{
    assert(num_alleles <= max_alleles);
    SamplingDepthEntry key{num_alleles, std::min(coi, num_alleles)};
    return lookup_sampling_depth[key];
}
