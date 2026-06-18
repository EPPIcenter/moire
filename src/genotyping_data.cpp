
#include "genotyping_data.h"

#include "mcmc_utils.h"

#include <span>
#include <stdexcept>

GenotypingData::GenotypingData(const Rcpp::List &args)
{
    if (args.containsElementNamed("aggregate") && !Rcpp::RObject(args["aggregate"]).isNULL())
    {
        aggregate = UtilFunctions::r_to_string(args["aggregate"]);
    }
    else
    {
        aggregate = "binary";
    }

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
            observed_alleles.inner_fill(
                {sample, locus},
                std::span<const int>(
                    observed_alleles_input[locus][sample].data(),
                    observed_alleles_input[locus][sample].size()));
        }
    }

    for (std::size_t locus = 0; locus < num_loci; ++locus)
    {
        for (std::size_t sample = 0; sample < num_samples; ++sample)
        {
            is_missing_.unchecked_at({sample, locus}) = is_missing_input[locus][sample];
        }
    }

    observed_coi = std::vector<std::size_t>(num_samples, 0);

    for (size_t sample_idx = 0; sample_idx < num_samples; sample_idx++)
    {
        for (size_t locus_idx = 0; locus_idx < num_loci; locus_idx++)
        {
            std::size_t present_alleles = 0;
            const auto [begin, end] = observed_alleles.inner_iterators({sample_idx, locus_idx});
            for (auto it = begin; it != end; ++it)
            {
                if (*it > 0)
                {
                    ++present_alleles;
                }
            }

            if (present_alleles > observed_coi[sample_idx])
            {
                observed_coi[sample_idx] = present_alleles;
            }
        }
    }

    jaccard_similarity_matrix = MultiVector<float, 2>({num_samples, num_samples});
    jaccard_similarity_matrix.fill(0.0);

    for (size_t i = 0; i < num_samples; ++i) {
        for (size_t j = i + 1; j < num_samples; ++j) {
            for (size_t locus_idx = 0; locus_idx < num_loci; ++locus_idx) {
                jaccard_similarity_matrix.unchecked_at({i, j}) += UtilFunctions::jaccard_similarity<float>(get_observed_alleles(i, locus_idx), get_observed_alleles(j, locus_idx));
            }
            jaccard_similarity_matrix.unchecked_at({i, j}) /= num_loci;
            jaccard_similarity_matrix.unchecked_at({j, i}) = jaccard_similarity_matrix.unchecked_at({i, j});
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
    return is_missing_.unchecked_at({sample, locus});
}

bool GenotypingData::has_count_barcodes() const
{
    for (std::size_t sample_idx = 0; sample_idx < num_samples; ++sample_idx)
    {
        for (std::size_t locus_idx = 0; locus_idx < num_loci; ++locus_idx)
        {
            const auto [begin, end] = observed_alleles.inner_iterators({sample_idx, locus_idx});
            for (auto it = begin; it != end; ++it)
            {
                if (*it > 1)
                {
                    return true;
                }
            }
        }
    }
    return false;
}

void GenotypingData::validate_for_observation_model(
    ObservationModelKind observation_model_kind) const
{
    if (observation_model_kind == ObservationModelKind::Binary && has_count_barcodes())
    {
        throw std::runtime_error(
            "Count barcodes (allele values > 1) require observation_model = \"counts\".");
    }

    if (observation_model_kind == ObservationModelKind::Counts && aggregate == "binary"
        && has_count_barcodes())
    {
        throw std::runtime_error(
            "Count barcodes require aggregate = \"count\" when loading data.");
    }
}
