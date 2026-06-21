#pragma once

#include "observation_model.h"
#include "observation_model_math.h"

#include <algorithm>
#include <cmath>
#include <span>
#include <vector>

class BinaryObservationModel : public ObservationModel {
   public:
    BinaryObservationModel(float max_eps_neg, float max_eps_pos)
        : max_eps_neg_(max_eps_neg), max_eps_pos_(max_eps_pos) {}

    float log_likelihood(std::span<int const> latent_allele_indices,
                         std::span<int const> observed_barcode,
                         float epsilon_neg,
                         float epsilon_pos) const override
    {
        const float norm_factor = 1.0f / static_cast<float>(observed_barcode.size());
        const float neg_rate = epsilon_neg * max_eps_neg_ * norm_factor;
        const float pos_rate = epsilon_pos * max_eps_pos_ * norm_factor;
        const float log_neg_correct = std::log(1.0f - neg_rate);
        const float log_neg_wrong = std::log(neg_rate);
        const float log_pos_correct = std::log(1.0f - pos_rate);
        const float log_pos_wrong = std::log(pos_rate);

        unsigned int fp = 0;
        unsigned int tp = 0;
        unsigned int fn = 0;
        unsigned int tn = 0;

        auto first_instance =
            std::find(latent_allele_indices.begin(), latent_allele_indices.end(), -1);
        std::span<int const> allele_index_vec =
            latent_allele_indices.subspan(0, first_instance - latent_allele_indices.begin());

        const unsigned int total_alleles = static_cast<unsigned int>(allele_index_vec.size());
        unsigned int vec_pointer = 0;
        int j = 0;
        for (const auto &e : observed_barcode)
        {
            const int next_allele_index =
                (vec_pointer < total_alleles) ? allele_index_vec[vec_pointer] : -1;
            const unsigned int mask = (j == next_allele_index) ? 1u : 0u;
            fp += (e & 1) & (1u - mask);
            tp += (e & 1) & mask;
            fn += static_cast<unsigned int>(!(e & 1)) & mask;
            tn += static_cast<unsigned int>(!(e & 1)) & (1u - mask);
            vec_pointer += mask;
            ++j;
        }

        return log_neg_correct * static_cast<float>(tp) + log_neg_wrong * static_cast<float>(fn)
             + log_pos_correct * static_cast<float>(tn) + log_pos_wrong * static_cast<float>(fp);
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
