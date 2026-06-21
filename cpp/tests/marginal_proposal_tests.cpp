#include <catch2/catch_test_macros.hpp>

#include "observation_model_math.h"

#include <cmath>
#include <map>
#include <random>
#include <vector>

using observation_model_math::marginal_presence_prob;
using observation_model_math::marginal_proposal_log_prob;

namespace {

std::vector<int> mask_to_indices(int mask, int n)
{
    std::vector<int> idx;
    for (int a = 0; a < n; ++a) {
        if (mask & (1 << a)) idx.push_back(a);
    }
    return idx;
}

}  // namespace

// The marginal latent-genotype proposal must be a proper distribution over the
// non-empty genotypes (it is conditioned on >=1 strain present), so the density
// returned by marginal_proposal_log_prob must sum to 1 over all non-empty subsets.
TEST_CASE("marginal proposal density normalizes over non-empty genotypes")
{
    const std::vector<float> rho = {0.8f, 0.3f, 0.6f, 0.1f};
    const int n = static_cast<int>(rho.size());

    double total = 0.0;
    for (int mask = 1; mask < (1 << n); ++mask) {  // mask 0 == empty, excluded
        const auto idx = mask_to_indices(mask, n);
        total += std::exp(static_cast<double>(marginal_proposal_log_prob(
            std::span<const float>(rho.data(), rho.size()),
            std::span<int const>(idx.data(), idx.size()))));
    }
    REQUIRE(std::abs(total - 1.0) < 1e-5);
}

// The sampler (draw allele a present w.p. rho_a, reject empty draws) must realize
// exactly the density evaluator, so empirical genotype frequencies match the
// closed-form non-empty-conditioned density. This is the contract that makes the
// reverse-move MH ratio correct.
TEST_CASE("marginal proposal sampling matches its density")
{
    const std::vector<float> rho = {0.7f, 0.4f, 0.55f};
    const int n = static_cast<int>(rho.size());

    std::mt19937 rng(12345u);
    std::uniform_real_distribution<float> unif(0.0f, 1.0f);

    std::map<int, long> counts;
    long total = 0;
    for (int iter = 0; iter < 600000; ++iter) {
        std::vector<int> present;
        for (int attempt = 0; attempt < 128 && present.empty(); ++attempt) {
            present.clear();
            for (int a = 0; a < n; ++a) {
                if (unif(rng) < rho[static_cast<std::size_t>(a)]) present.push_back(a);
            }
        }
        if (present.empty()) continue;
        int mask = 0;
        for (int a : present) mask |= (1 << a);
        ++counts[mask];
        ++total;
    }

    for (int mask = 1; mask < (1 << n); ++mask) {
        const auto idx = mask_to_indices(mask, n);
        const double empirical = static_cast<double>(counts[mask]) / static_cast<double>(total);
        const double theoretical = std::exp(static_cast<double>(marginal_proposal_log_prob(
            std::span<const float>(rho.data(), rho.size()),
            std::span<int const>(idx.data(), idx.size()))));
        REQUIRE(std::abs(empirical - theoretical) < 0.005);
    }
}

// Presence probabilities are clamped into (0,1) and ordered sensibly: an observed
// allele is more likely present than an unobserved one under the same error rates.
TEST_CASE("marginal presence probability is bounded and ordered")
{
    const float neg_rate = 0.05f;
    const float pos_rate = 0.02f;
    const float rho_pos = marginal_presence_prob(true, neg_rate, pos_rate);
    const float rho_neg = marginal_presence_prob(false, neg_rate, pos_rate);
    REQUIRE(rho_pos > rho_neg);
    REQUIRE(rho_pos > 0.0f);
    REQUIRE(rho_pos < 1.0f);
    REQUIRE(rho_neg > 0.0f);
    REQUIRE(rho_neg < 1.0f);
}
