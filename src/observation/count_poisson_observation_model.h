#pragma once

#include "observation_model.h"
#include "observation_model_math.h"

#include <algorithm>
#include <cmath>
#include <span>
#include <vector>

class CountPoissonObservationModel : public ObservationModel {
   public:
    CountPoissonObservationModel(float max_eps_neg, float max_eps_pos)
        : max_eps_neg_(max_eps_neg), max_eps_pos_(max_eps_pos) {}

    float log_likelihood(std::span<int const> latent_allele_indices,
                         std::span<int const> observed_barcode,
                         float epsilon_neg,
                         float epsilon_pos) const override
    {
        const int num_alleles = static_cast<int>(observed_barcode.size());
        if (num_alleles == 0)
        {
            return 0.0f;
        }

        const int total_depth = observation_model_math::total_reads(observed_barcode);
        const float allele_count = static_cast<float>(num_alleles);
        const float dropout_prob = std::min(1.0f, epsilon_neg * max_eps_neg_);
        const float noise_rate = epsilon_pos * max_eps_pos_ * (static_cast<float>(total_depth) / allele_count);
        const float noise_mean = observation_model_math::floor_mean(noise_rate / allele_count);

        const std::size_t support_alleles = observation_model_math::support_size(latent_allele_indices);
        const float signal_mean = support_alleles > 0
            ? observation_model_math::floor_mean(static_cast<float>(total_depth)
                                               / static_cast<float>(support_alleles))
            : observation_model_math::kPoissonMeanFloor;

        float log_likelihood_total = 0.0f;
        for (int allele_idx = 0; allele_idx < num_alleles; ++allele_idx)
        {
            const int observed_count = observed_barcode[static_cast<std::size_t>(allele_idx)];
            if (observation_model_math::allele_in_support(allele_idx, latent_allele_indices))
            {
                log_likelihood_total += observation_model_math::present_mixture_log_pmf(
                    observed_count, signal_mean, dropout_prob, noise_mean);
            }
            else
            {
                log_likelihood_total +=
                    observation_model_math::poisson_log_pmf(observed_count, noise_mean);
            }
        }

        return log_likelihood_total;
    }

    LatentGenotype sample_latent_genotype(Sampler &sampler,
                                          std::span<int const> observed_barcode,
                                          int coi,
                                          float epsilon_pos,
                                          float epsilon_neg) const override;

    LatentGenotype propose_latent_genotype_marginal(
        Sampler &sampler, std::span<int const> observed_barcode, float epsilon_pos,
        float epsilon_neg) const override;

    float latent_genotype_log_prob_marginal(
        std::span<int const> latent_allele_indices, std::span<int const> observed_barcode,
        float epsilon_pos, float epsilon_neg) const override;

   private:
    // Per-allele Bernoulli presence probabilities for the marginal proposal.
    void marginal_presence_probs(std::span<int const> observed_barcode, float epsilon_pos,
                                 float epsilon_neg, std::vector<float> &rho_out) const;

    float max_eps_neg_;
    float max_eps_pos_;
};
