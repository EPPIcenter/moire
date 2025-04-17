
#include "genotyping_data.h"

#include "mcmc_utils.h"

#include <span>

GenotypingData::GenotypingData(const Rcpp::List &args)
{
    // indexed by locus, sample, allele
    std::vector<std::vector<std::vector<int>>> observed_alleles_input = UtilFunctions::r_to_array_int(args["data"]);
    // indexed by locus, sample
    std::vector<std::vector<int>> is_missing_input = UtilFunctions::r_to_mat_int(args["is_missing"]);

    num_loci = observed_alleles_input.size();
    num_samples = observed_alleles_input[0].size();
    num_alleles = std::vector<std::size_t>(num_loci, 0);

    for (std::size_t locus = 0; locus < num_loci; ++locus)
    {
        num_alleles[locus] = observed_alleles_input[locus][0].size();
    }

    observed_alleles = RaggedMultiVector<int, 3>({num_samples, num_loci}, num_alleles);
    is_missing_ = MultiVector<int, 2>({num_samples, num_loci});


    for (std::size_t locus = 0; locus < num_loci; ++locus)
    {
        for (std::size_t sample = 0; sample < num_samples; ++sample)
        {
            for (std::size_t allele = 0; allele < num_alleles[locus]; ++allele)
            {
                observed_alleles.at({sample, locus, allele}) = observed_alleles_input[locus][sample][allele];
            }
        }
    }

    for (std::size_t locus = 0; locus < num_loci; ++locus)
    {
        for (std::size_t sample = 0; sample < num_samples; ++sample)
        {
            is_missing_.at({sample, locus}) = is_missing_input[locus][sample];
        }
    }

    observed_coi = std::vector<std::size_t>(num_samples, 0);

    for (size_t sample_idx = 0; sample_idx < num_samples; sample_idx++)
    {
        for (size_t locus_idx = 0; locus_idx < num_loci; locus_idx++)
        {
            std::size_t total_alleles = 0;
            for (size_t allele_idx = 0; allele_idx < num_alleles[locus_idx]; allele_idx++)
            {
                total_alleles += observed_alleles.at({sample_idx, locus_idx, allele_idx});
            }

            if (total_alleles > observed_coi[sample_idx])
            {
                observed_coi[sample_idx] = total_alleles;
            }
        }
    }
}

std::span<int const> GenotypingData::get_observed_alleles(std::size_t sample, std::size_t locus) const
{
    auto [start, end] = observed_alleles.inner_iterators({sample, locus});
    return std::span<int const>(start, end);
}

bool GenotypingData::is_missing(std::size_t sample, std::size_t locus) const
{
    return is_missing_.at({sample, locus});
}
