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

}  // namespace observation_model_math
