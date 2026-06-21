#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numeric>
#include <span>
#include <vector>

namespace observation_model_math {

inline constexpr float kPoissonMeanFloor = 1e-10f;

inline float floor_mean(float mean) { return std::max(mean, kPoissonMeanFloor); }

/// Memoized lgamma(n+1) for non-negative integer read counts. The observed
/// counts are fixed data, so lgamma is recomputed for the same handful of
/// integer values on every likelihood evaluation; a per-thread table removes
/// those repeated calls. Falls back to a direct lgamma for atypically large
/// counts to keep the table bounded.
inline float log_factorial(int n)
{
    constexpr int kMaxCached = 1 << 16;  // ~256 KB/thread upper bound
    if (n < 0 || n >= kMaxCached) {
        return static_cast<float>(std::lgamma(static_cast<double>(n) + 1.0));
    }
    static thread_local std::vector<float> cache;
    if (static_cast<std::size_t>(n) >= cache.size()) {
        const std::size_t old_size = cache.size();
        cache.resize(static_cast<std::size_t>(n) + 1);
        for (std::size_t i = old_size; i < cache.size(); ++i) {
            cache[i] = static_cast<float>(std::lgamma(static_cast<double>(i) + 1.0));
        }
    }
    return cache[static_cast<std::size_t>(n)];
}

inline float log_sum_exp(float a, float b)
{
    if (a == -std::numeric_limits<float>::infinity())
    {
        return b;
    }
    if (b == -std::numeric_limits<float>::infinity())
    {
        return a;
    }
    if (a > b)
    {
        return a + std::log1p(std::exp(b - a));
    }
    return b + std::log1p(std::exp(a - b));
}

inline float poisson_log_pmf(int count, float mean)
{
    const float mu = floor_mean(mean);
    return static_cast<float>(count) * std::log(mu) - mu - log_factorial(count);
}

inline float present_mixture_log_pmf(int count, float signal_mean, float dropout_prob,
                                     float noise_mean)
{
    const float mu = floor_mean(signal_mean);
    const float noise = floor_mean(noise_mean);
    const float pi = std::clamp(dropout_prob, 0.0f, 1.0f);
    const float log_one_minus_pi = std::log(std::max(1.0f - pi, kPoissonMeanFloor));

    const float log_signal = log_one_minus_pi + poisson_log_pmf(count, mu);
    const float log_dropout = std::log(std::max(pi, kPoissonMeanFloor)) + poisson_log_pmf(count, noise);
    return log_sum_exp(log_signal, log_dropout);
}

inline int total_reads(std::span<int const> observed_barcode)
{
    return std::accumulate(observed_barcode.begin(), observed_barcode.end(), 0);
}

inline std::size_t support_size(std::span<int const> latent_allele_indices)
{
    auto first_sentinel =
        std::find(latent_allele_indices.begin(), latent_allele_indices.end(), -1);
    return static_cast<std::size_t>(first_sentinel - latent_allele_indices.begin());
}

inline bool allele_in_support(int allele_index, std::span<int const> latent_allele_indices)
{
    const auto support_end =
        latent_allele_indices.begin()
        + static_cast<std::ptrdiff_t>(support_size(latent_allele_indices));
    return std::find(latent_allele_indices.begin(), support_end, allele_index) != support_end;
}

// ---------------------------------------------------------------------------
// Marginal (eCOI) latent-genotype proposal.
//
// The fully-collapsed eCOI sampler integrates COI out of the model, so its
// latent-genotype proposal must be COI-independent AND have a density that can
// be evaluated for an *arbitrary* genotype (needed for the reverse move). We use
// an independent per-allele Bernoulli: allele a is present with probability
// rho_a, the posterior presence probability under the per-allele error channel
// and a uniform present/absent prior. rho_a is clamped away from {0,1} so the
// chain stays ergodic and reverse-proposal densities are finite. The proposal is
// conditioned on a non-empty genotype (every observed locus carries >=1 strain).
// ---------------------------------------------------------------------------

inline constexpr float kMarginalRhoClamp = 1e-4f;

inline float marginal_presence_prob(bool observed_positive, float neg_rate, float pos_rate)
{
    const float nr = std::clamp(neg_rate, 0.0f, 1.0f);
    const float pr = std::clamp(pos_rate, 0.0f, 1.0f);
    float rho;
    if (observed_positive) {
        // P(present | obs=positive) propto P(obs=pos | present) = 1 - FN rate
        const float denom = (1.0f - nr) + pr;
        rho = denom > 0.0f ? (1.0f - nr) / denom : 0.5f;
    } else {
        // P(present | obs=negative) propto P(obs=neg | present) = FN rate
        const float denom = nr + (1.0f - pr);
        rho = denom > 0.0f ? nr / denom : 0.5f;
    }
    return std::clamp(rho, kMarginalRhoClamp, 1.0f - kMarginalRhoClamp);
}

// log P(empty genotype) under independent Bernoulli(rho_a).
inline float marginal_log_prob_empty(std::span<const float> rho)
{
    double s = 0.0;
    for (const float r : rho) {
        s += std::log1p(-static_cast<double>(r));
    }
    return static_cast<float>(s);
}

// log q(G | non-empty) for an arbitrary latent genotype G (sentinel -1 padded ok).
inline float marginal_proposal_log_prob(std::span<const float> rho,
                                        std::span<int const> latent_allele_indices)
{
    const std::size_t k = support_size(latent_allele_indices);
    std::vector<char> present(rho.size(), 0);
    for (std::size_t i = 0; i < k; ++i) {
        const int a = latent_allele_indices[i];
        if (a >= 0 && static_cast<std::size_t>(a) < rho.size()) {
            present[static_cast<std::size_t>(a)] = 1;
        }
    }
    double lp = 0.0;
    for (std::size_t a = 0; a < rho.size(); ++a) {
        lp += present[a] ? std::log(static_cast<double>(rho[a]))
                         : std::log1p(-static_cast<double>(rho[a]));
    }
    const double log_p_empty = static_cast<double>(marginal_log_prob_empty(rho));
    // log(1 - P(empty)); guard against log(0) when P(empty) -> 1 (rho all tiny).
    const double one_minus = -std::expm1(log_p_empty);  // = 1 - exp(log_p_empty)
    const double log_p_nonempty =
        one_minus > 0.0 ? std::log(one_minus) : std::log(kPoissonMeanFloor);
    return static_cast<float>(lp - log_p_nonempty);
}

}  // namespace observation_model_math
