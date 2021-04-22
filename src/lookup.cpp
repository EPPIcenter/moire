#include "lookup.h"

#include "mcmc_utils.h"

#include <algorithm>

Lookup::Lookup(int max_coi, int max_alleles)
    : max_coi(max_coi), max_alleles(max_alleles)
{
    init_lgamma();
    init_sampling_depth();
};

void Lookup::init_lgamma()
{
    lookup_lgamma = std::vector<double>(max_coi + max_alleles + 10);

    for (size_t i = 0; i < lookup_lgamma.size(); i++)
    {
        lookup_lgamma[i] = std::lgamma(i);
    }
    lgamma_initialized = true;
};

void Lookup::init_sampling_depth()
{
    // 2D vector with rows of length max_coi, max_alleles total rows
    // max_coi is leq max_alleles
    // this is approximate due to rounding errors
    lookup_sampling_depth =
        std::vector<long>((max_coi + 1) * (max_alleles + 1), 0);

    for (int i = 1; i < max_alleles + 1; i++)
    {
        for (int j = 1; j < max_coi + 1; j++)
        {
            // total samples required is the number of samples for
            // i choose 1 + i choose 2 + ... + i choose j
            for (int k = 1; k <= j; k++)
            {
                lookup_sampling_depth[max_coi * i + j] +=
                    std::exp(lookup_lgamma[i + 1] - lookup_lgamma[k + 1] -
                             lookup_lgamma[(i - k) + 1]);
            }
        }
    }
    sample_depth_initialized = true;
}

long Lookup::get_sampling_depth(int coi, int num_alleles)
{
    assert(num_alleles <= max_alleles);
    return lookup_sampling_depth[max_coi * num_alleles + coi];
}
