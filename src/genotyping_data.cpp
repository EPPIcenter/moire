
#include "genotyping_data.h"

#include "mcmc_utils.h"

#include <span>

GenotypingData::GenotypingData(const Rcpp::List &args)
{
    observed_alleles = UtilFunctions::r_to_array_int(args["data"]);
    is_missing_ = UtilFunctions::r_to_mat_bool(args["is_missing"]);

    num_loci = observed_alleles.size();
    num_samples = observed_alleles[0].size();

    num_alleles = std::vector<std::size_t>(num_loci);
    observed_coi = std::vector<std::size_t>(num_samples, 0);

    for (size_t i = 0; i < num_loci; i++)
    {
        num_alleles[i] = observed_alleles[i][0].size();

        for (size_t j = 0; j < num_samples; j++)
        {
            std::size_t total_alleles = 0;
            for (size_t k = 0; k < observed_alleles[i][j].size(); k++)
            {
                total_alleles += observed_alleles[i][j][k];
            }

            if (total_alleles > observed_coi[j])
            {
                observed_coi[j] = total_alleles;
            }
        }
    }
}

std::span<int const> GenotypingData::get_observed_alleles(int locus,
                                                             int sample) const
{
    return observed_alleles.at(locus).at(sample);
}

bool GenotypingData::is_missing(int locus, int sample) const
{
    return is_missing_.at(locus).at(sample);
}
