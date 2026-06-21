#include "binary_observation_model.h"
#include "count_poisson_observation_model.h"
#include "observation_model_math.h"
#include "sampler.h"

#include <boost/math/special_functions/binomial.hpp>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <numeric>
#include <span>
#include <vector>

namespace {

// Shared sampler for the COI-independent marginal proposal: draw each allele
// present with probability rho_a, conditioned on a non-empty genotype, and
// return the genotype together with its (non-empty conditioned) log density that
// observation_model_math::marginal_proposal_log_prob reproduces exactly.
LatentGenotype marginal_sample_from_rho(Sampler &sampler, const std::vector<float> &rho)
{
    std::vector<int> present;
    constexpr int kMaxAttempts = 128;
    for (int attempt = 0; attempt < kMaxAttempts && present.empty(); ++attempt)
    {
        present.clear();
        for (std::size_t a = 0; a < rho.size(); ++a)
        {
            if (sampler.sample_unif() < rho[a])
            {
                present.push_back(static_cast<int>(a));
            }
        }
    }
    if (present.empty() && !rho.empty())
    {
        std::size_t best = 0;
        for (std::size_t a = 1; a < rho.size(); ++a)
        {
            if (rho[a] > rho[best]) best = a;
        }
        present.push_back(static_cast<int>(best));
    }
    // present is already ascending since alleles are visited in order.
    const float log_prob = observation_model_math::marginal_proposal_log_prob(
        std::span<const float>(rho.data(), rho.size()),
        std::span<int const>(present.data(), present.size()));
    return LatentGenotype{present, log_prob};
}

}  // namespace

namespace {

struct WeightedSelection {
    std::vector<int> selected;
    float log_prob = 0.0f;
};

WeightedSelection weighted_sample_k(Sampler &sampler, const std::vector<int> &pool,
                                    const std::vector<float> &log_weights, int k)
{
    WeightedSelection result;
    if (k <= 0)
    {
        return result;
    }

    std::vector<int> remaining_pool = pool;
    std::vector<float> remaining_log_weights = log_weights;

    for (int pick = 0; pick < k; ++pick)
    {
        float max_log_weight = -std::numeric_limits<float>::infinity();
        for (const float log_weight : remaining_log_weights)
        {
            max_log_weight = std::max(max_log_weight, log_weight);
        }

        float total_weight = 0.0f;
        std::vector<float> weights(remaining_log_weights.size());
        for (std::size_t idx = 0; idx < remaining_log_weights.size(); ++idx)
        {
            weights[idx] = std::exp(remaining_log_weights[idx] - max_log_weight);
            total_weight += weights[idx];
        }

        const float log_total_weight = max_log_weight + std::log(total_weight);
        const float draw = sampler.sample_unif() * total_weight;

        float cumulative_weight = 0.0f;
        std::size_t chosen_idx = 0;
        for (std::size_t idx = 0; idx < weights.size(); ++idx)
        {
            cumulative_weight += weights[idx];
            if (draw <= cumulative_weight)
            {
                chosen_idx = idx;
                break;
            }
        }

        result.selected.push_back(remaining_pool[chosen_idx]);
        result.log_prob += remaining_log_weights[chosen_idx] - log_total_weight;

        remaining_pool.erase(remaining_pool.begin() + static_cast<std::ptrdiff_t>(chosen_idx));
        remaining_log_weights.erase(
            remaining_log_weights.begin() + static_cast<std::ptrdiff_t>(chosen_idx));
    }

    return result;
}

}  // namespace

LatentGenotype BinaryObservationModel::sample_latent_genotype(
    Sampler &sampler, std::span<int const> observed_barcode, int coi, float epsilon_pos,
    float epsilon_neg) const
{
    const int total_alleles = static_cast<int>(observed_barcode.size());
    int total_obs_positives = 0;
    int total_obs_negatives = 0;
    std::vector<int> obs_positive_indices{};
    std::vector<int> obs_negative_indices{};

    for (int allele_idx = 0; allele_idx < total_alleles; ++allele_idx)
    {
        if (observed_barcode[static_cast<std::size_t>(allele_idx)] == 1)
        {
            total_obs_positives++;
            obs_positive_indices.push_back(allele_idx);
        }
        else
        {
            total_obs_negatives++;
            obs_negative_indices.push_back(allele_idx);
        }
    }

    const int min_false_negatives = total_obs_positives == 0;
    const int max_false_negatives =
        std::max(min_false_negatives, std::min(coi, total_obs_negatives));

    int total_false_negatives = min_false_negatives;
    for (int allele_idx = total_false_negatives; allele_idx < max_false_negatives; ++allele_idx)
    {
        total_false_negatives +=
            (sampler.sample_unif() < (epsilon_neg / static_cast<float>(total_alleles)));
    }
    const int total_true_negatives = total_obs_negatives - total_false_negatives;

    const float per_allele_fn_rate = epsilon_neg / static_cast<float>(total_alleles);
    const float log_prob_total_false_negatives =
        std::log(boost::math::binomial_coefficient<float>(
            total_obs_negatives, total_false_negatives - min_false_negatives))
        + static_cast<float>(total_false_negatives - min_false_negatives)
              * std::log(per_allele_fn_rate)
        + static_cast<float>(total_true_negatives) * std::log(1.0f - per_allele_fn_rate);

    const int min_false_positives =
        std::max(0, (total_obs_positives + total_false_negatives) - coi);
    const int max_false_positives =
        std::min(total_obs_positives, total_false_negatives / 2);

    int total_false_positives = min_false_positives;
    for (int allele_idx = total_false_positives; allele_idx < max_false_positives; ++allele_idx)
    {
        total_false_positives +=
            (sampler.sample_unif() < (epsilon_pos / static_cast<float>(total_alleles)));
    }
    const int total_true_positives = total_obs_positives - total_false_positives;

    const float per_allele_fp_rate = epsilon_pos / static_cast<float>(total_alleles);
    const float log_prob_total_false_positives =
        std::log(boost::math::binomial_coefficient<float>(
            total_obs_positives, total_false_positives - min_false_positives))
        + static_cast<float>(total_false_positives - min_false_positives)
              * std::log(per_allele_fp_rate)
        + static_cast<float>(total_true_positives) * std::log(1.0f - per_allele_fp_rate);

    sampler.shuffle_vec(obs_positive_indices);
    sampler.shuffle_vec(obs_negative_indices);

    std::vector<int> allele_index_vec{};
    allele_index_vec.insert(
        allele_index_vec.end(), obs_positive_indices.begin(),
        obs_positive_indices.begin() + total_true_positives);
    allele_index_vec.insert(
        allele_index_vec.end(), obs_negative_indices.begin(),
        obs_negative_indices.begin() + total_false_negatives);

    std::sort(allele_index_vec.begin(), allele_index_vec.end());

    const float log_prob_positive_indices =
        -std::log(boost::math::binomial_coefficient<float>(total_obs_positives,
                                                           total_true_positives));
    const float log_prob_negative_indices =
        -std::log(boost::math::binomial_coefficient<float>(total_obs_negatives,
                                                           total_false_negatives));

    const float log_prob = log_prob_positive_indices + log_prob_negative_indices
                         + log_prob_total_false_positives + log_prob_total_false_negatives;

    assert(allele_index_vec.size() > 0);
    assert(allele_index_vec.size() <= static_cast<std::size_t>(coi));

    return LatentGenotype{allele_index_vec, log_prob};
}

LatentGenotype CountPoissonObservationModel::sample_latent_genotype(
    Sampler &sampler, std::span<int const> observed_barcode, int coi, float epsilon_pos,
    float epsilon_neg) const
{
    const int total_alleles = static_cast<int>(observed_barcode.size());
    const int total_depth = observation_model_math::total_reads(observed_barcode);
    const float allele_count = static_cast<float>(total_alleles);
    const float dropout_prob = std::min(1.0f, epsilon_neg * max_eps_neg_);
    const float fp_rate = std::min(1.0f, epsilon_pos * max_eps_pos_ / allele_count);
    const float noise_rate =
        epsilon_pos * max_eps_pos_ * (static_cast<float>(total_depth) / allele_count);
    const float noise_mean = observation_model_math::floor_mean(noise_rate / allele_count);

    int total_obs_positives = 0;
    int total_obs_negatives = 0;
    std::vector<int> obs_positive_indices{};
    std::vector<int> obs_negative_indices{};

    for (int allele_idx = 0; allele_idx < total_alleles; ++allele_idx)
    {
        if (observed_barcode[static_cast<std::size_t>(allele_idx)] > 0)
        {
            total_obs_positives++;
            obs_positive_indices.push_back(allele_idx);
        }
        else
        {
            total_obs_negatives++;
            obs_negative_indices.push_back(allele_idx);
        }
    }

    const int min_false_negatives = total_obs_positives == 0 ? 1 : 0;
    const int max_false_negatives =
        std::max(min_false_negatives, std::min(coi, total_obs_negatives));

    int total_false_negatives = min_false_negatives;
    for (int allele_idx = total_false_negatives; allele_idx < max_false_negatives; ++allele_idx)
    {
        total_false_negatives += (sampler.sample_unif() < dropout_prob);
    }
    const int total_true_negatives = total_obs_negatives - total_false_negatives;

    const float log_prob_total_false_negatives =
        std::log(boost::math::binomial_coefficient<float>(
            total_obs_negatives, total_false_negatives - min_false_negatives))
        + static_cast<float>(total_false_negatives - min_false_negatives) * std::log(dropout_prob)
        + static_cast<float>(total_true_negatives) * std::log(1.0f - dropout_prob);

    const int min_false_positives =
        std::max(0, (total_obs_positives + total_false_negatives) - coi);
    const int max_false_positives =
        std::min(total_obs_positives, total_false_negatives / 2);

    int total_false_positives = min_false_positives;
    for (int allele_idx = total_false_positives; allele_idx < max_false_positives; ++allele_idx)
    {
        total_false_positives += (sampler.sample_unif() < fp_rate);
    }
    const int total_true_positives = total_obs_positives - total_false_positives;

    const float log_prob_total_false_positives =
        std::log(boost::math::binomial_coefficient<float>(
            total_obs_positives, total_false_positives - min_false_positives))
        + static_cast<float>(total_false_positives - min_false_positives) * std::log(fp_rate)
        + static_cast<float>(total_true_positives) * std::log(1.0f - fp_rate);

    const std::size_t support_alleles =
        static_cast<std::size_t>(total_true_positives + total_false_negatives);
    const float signal_mean = observation_model_math::floor_mean(
        static_cast<float>(total_depth) / static_cast<float>(support_alleles));

    std::vector<float> positive_log_weights(obs_positive_indices.size());
    for (std::size_t idx = 0; idx < obs_positive_indices.size(); ++idx)
    {
        const int allele_idx = obs_positive_indices[idx];
        positive_log_weights[idx] = observation_model_math::present_mixture_log_pmf(
            observed_barcode[static_cast<std::size_t>(allele_idx)], signal_mean, dropout_prob,
            noise_mean);
    }

    std::vector<float> negative_log_weights(obs_negative_indices.size());
    for (std::size_t idx = 0; idx < obs_negative_indices.size(); ++idx)
    {
        negative_log_weights[idx] = observation_model_math::present_mixture_log_pmf(
            0, signal_mean, dropout_prob, noise_mean);
    }

    const WeightedSelection positive_selection = weighted_sample_k(
        sampler, obs_positive_indices, positive_log_weights, total_true_positives);
    const WeightedSelection negative_selection = weighted_sample_k(
        sampler, obs_negative_indices, negative_log_weights, total_false_negatives);

    std::vector<int> allele_index_vec = positive_selection.selected;
    allele_index_vec.insert(allele_index_vec.end(), negative_selection.selected.begin(),
                            negative_selection.selected.end());
    std::sort(allele_index_vec.begin(), allele_index_vec.end());

    const float log_prob = positive_selection.log_prob + negative_selection.log_prob
                         + log_prob_total_false_positives + log_prob_total_false_negatives;

    assert(allele_index_vec.size() > 0);
    assert(allele_index_vec.size() <= static_cast<std::size_t>(coi));

    return LatentGenotype{allele_index_vec, log_prob};
}

// ---------------------------------------------------------------------------
// COI-independent marginal proposals.
// ---------------------------------------------------------------------------

void BinaryObservationModel::marginal_presence_probs(std::span<int const> observed_barcode,
                                                     float epsilon_pos, float epsilon_neg,
                                                     std::vector<float> &rho_out) const
{
    const std::size_t n = observed_barcode.size();
    const float norm_factor = n > 0 ? 1.0f / static_cast<float>(n) : 1.0f;
    const float neg_rate = epsilon_neg * max_eps_neg_ * norm_factor;
    const float pos_rate = epsilon_pos * max_eps_pos_ * norm_factor;
    rho_out.resize(n);
    for (std::size_t a = 0; a < n; ++a)
    {
        const bool observed_positive = (observed_barcode[a] & 1) != 0;
        rho_out[a] = observation_model_math::marginal_presence_prob(observed_positive, neg_rate,
                                                                    pos_rate);
    }
}

LatentGenotype BinaryObservationModel::propose_latent_genotype_marginal(
    Sampler &sampler, std::span<int const> observed_barcode, float epsilon_pos,
    float epsilon_neg) const
{
    std::vector<float> rho;
    marginal_presence_probs(observed_barcode, epsilon_pos, epsilon_neg, rho);
    return marginal_sample_from_rho(sampler, rho);
}

float BinaryObservationModel::latent_genotype_log_prob_marginal(
    std::span<int const> latent_allele_indices, std::span<int const> observed_barcode,
    float epsilon_pos, float epsilon_neg) const
{
    std::vector<float> rho;
    marginal_presence_probs(observed_barcode, epsilon_pos, epsilon_neg, rho);
    return observation_model_math::marginal_proposal_log_prob(
        std::span<const float>(rho.data(), rho.size()), latent_allele_indices);
}

void CountPoissonObservationModel::marginal_presence_probs(std::span<int const> observed_barcode,
                                                           float epsilon_pos, float epsilon_neg,
                                                           std::vector<float> &rho_out) const
{
    const std::size_t n = observed_barcode.size();
    const float allele_count = n > 0 ? static_cast<float>(n) : 1.0f;
    const float neg_rate = std::min(1.0f, epsilon_neg * max_eps_neg_);
    const float pos_rate = std::min(1.0f, epsilon_pos * max_eps_pos_ / allele_count);
    rho_out.resize(n);
    for (std::size_t a = 0; a < n; ++a)
    {
        const bool observed_positive = observed_barcode[a] > 0;
        rho_out[a] = observation_model_math::marginal_presence_prob(observed_positive, neg_rate,
                                                                    pos_rate);
    }
}

LatentGenotype CountPoissonObservationModel::propose_latent_genotype_marginal(
    Sampler &sampler, std::span<int const> observed_barcode, float epsilon_pos,
    float epsilon_neg) const
{
    std::vector<float> rho;
    marginal_presence_probs(observed_barcode, epsilon_pos, epsilon_neg, rho);
    return marginal_sample_from_rho(sampler, rho);
}

float CountPoissonObservationModel::latent_genotype_log_prob_marginal(
    std::span<int const> latent_allele_indices, std::span<int const> observed_barcode,
    float epsilon_pos, float epsilon_neg) const
{
    std::vector<float> rho;
    marginal_presence_probs(observed_barcode, epsilon_pos, epsilon_neg, rho);
    return observation_model_math::marginal_proposal_log_prob(
        std::span<const float>(rho.data(), rho.size()), latent_allele_indices);
}
