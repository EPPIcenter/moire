#pragma once

#include "prob_any_missing_cache.h"
#include "pam_fast_paths.h"
#include "sampler.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <iterator>
#include <limits>
#include <span>
#include <vector>

namespace transmission_incremental {

inline float log_sum_exp(std::span<const float> values)
{
    if (values.empty()) {
        return -std::numeric_limits<float>::infinity();
    }
    float max_el = values[0];
    for (float v : values) {
        max_el = std::max(max_el, v);
    }
    if (max_el == -std::numeric_limits<float>::infinity()) {
        return max_el;
    }
    float sum = 0.f;
    for (float v : values) {
        sum += std::exp(v - max_el);
    }
    return max_el + std::log(sum);
}

inline float log_sum_exp(const std::vector<float>& values)
{
    return log_sum_exp(std::span<const float>(values));
}

inline bool build_constrained_q(std::span<const int> allele_indices,
                                std::span<const float> p_full,
                                std::vector<float>& q_out,
                                float& sum_out)
{
    q_out.clear();
    sum_out = 0.f;
    for (int e : allele_indices) {
        if (e < 0 || static_cast<std::size_t>(e) >= p_full.size()) {
            return false;
        }
        q_out.push_back(p_full[e]);
        sum_out += q_out.back();
    }
    if (sum_out <= 0.f || !std::isfinite(sum_out)) {
        return false;
    }
    for (float& v : q_out) {
        v /= sum_out;
    }
    return true;
}

inline std::uint32_t float_bits(float x) noexcept
{
    std::uint32_t bits = 0;
    std::memcpy(&bits, &x, sizeof(bits));
    return bits;
}

/// Reuse log Binomial(i, coi-1, r) weights across transmission evals with the same (coi, r).
inline std::span<const float> cached_log_binomial_weights(Sampler& sampler,
                                                          int coi,
                                                          float relatedness,
                                                          std::size_t num_terms)
{
    thread_local int cached_coi = -1;
    thread_local std::uint32_t cached_r_bits = 0;
    thread_local std::vector<float> log_w;

    const std::uint32_t r_bits = float_bits(relatedness);
    if (cached_coi != coi || cached_r_bits != r_bits || log_w.size() < num_terms) {
        cached_coi = coi;
        cached_r_bits = r_bits;
        const std::size_t reserve_n = std::max(num_terms, static_cast<std::size_t>(coi));
        log_w.resize(reserve_n);
        for (std::size_t i = 0; i < reserve_n; ++i) {
            log_w[i] = sampler.dbinom(static_cast<int>(i), coi - 1, relatedness);
        }
    }
    return std::span<const float>(log_w.data(), num_terms);
}

inline float transmission_log_no_relatedness(std::span<const float> log_one_minus_pam,
                                             int coi,
                                             float log_sum)
{
    constexpr float kNegInf = -std::numeric_limits<float>::infinity();
    const std::size_t pam_idx = static_cast<std::size_t>(coi - 1);
    if (pam_idx >= log_one_minus_pam.size()) {
        return kNegInf;
    }
    return log_one_minus_pam[pam_idx] + log_sum * static_cast<float>(coi);
}

inline float transmission_log_no_relatedness(const std::vector<double>& pam_vec,
                                             int coi,
                                             float log_sum)
{
    const std::size_t pam_idx = static_cast<std::size_t>(coi - 1);
    if (pam_idx >= pam_vec.size()) {
        return -std::numeric_limits<float>::infinity();
    }
    const float log_one_minus =
        static_cast<float>(std::log1p(-pam_vec[pam_idx]));
    return log_one_minus + log_sum * static_cast<float>(coi);
}

/// K=1 on latent support: P(any missing)=0 for n>=1, so log(1-PAM)=0.
inline float transmission_log_with_relatedness_k1(Sampler& sampler,
                                                  int coi,
                                                  float relatedness,
                                                  float log_sum)
{
    constexpr float kNegInf = -std::numeric_limits<float>::infinity();
    if (coi < 1) {
        return kNegInf;
    }

    const std::size_t loop_upper = static_cast<std::size_t>(coi - 1);
    const auto log_w =
        cached_log_binomial_weights(sampler, coi, relatedness, loop_upper + 1);

    constexpr std::size_t kMaxTerms = 128;
    float terms[kMaxTerms];
    const std::size_t n_terms = loop_upper + 1;
    if (n_terms > kMaxTerms) {
        return kNegInf;
    }

    for (std::size_t i = 0; i < n_terms; ++i) {
        terms[i] = log_w[i] + log_sum * static_cast<float>(coi - static_cast<int>(i));
    }
    return log_sum_exp(std::span<const float>(terms, n_terms));
}

inline float transmission_log_with_relatedness(std::span<const double> pam_vec,
                                               std::span<const float> log_one_minus_pam,
                                               Sampler& sampler,
                                               int coi,
                                               std::size_t total_alleles,
                                               float relatedness,
                                               float log_sum)
{
    constexpr float kNegInf = -std::numeric_limits<float>::infinity();

    if (static_cast<int>(coi) < static_cast<int>(total_alleles)) {
        return kNegInf;
    }
    if (coi > 100 || total_alleles > 30 || (coi > 50 && total_alleles > 20)) {
        return kNegInf;
    }

    if (total_alleles == 1 && pam_fast_paths::tx_opts_enabled()) {
        return transmission_log_with_relatedness_k1(sampler, coi, relatedness, log_sum);
    }

    const std::size_t pam_idx = static_cast<std::size_t>(coi - 1);
    if (pam_idx >= pam_vec.size() || pam_idx >= log_one_minus_pam.size()) {
        return kNegInf;
    }

    const float log_one_minus_pam_at_coi = log_one_minus_pam[pam_idx];
    const float log_sum_scaled = log_sum * static_cast<float>(coi);

    if (static_cast<int>(total_alleles) == coi && pam_fast_paths::tx_opts_enabled()) {
        const float pr = cached_log_binomial_weights(sampler, coi, relatedness, 1)[0];
        return pr + log_one_minus_pam_at_coi + log_sum_scaled;
    }

    const std::size_t loop_upper = static_cast<std::size_t>(coi - total_alleles);
    const auto log_w =
        cached_log_binomial_weights(sampler, coi, relatedness, loop_upper + 1);

    constexpr std::size_t kMaxTerms = 128;
    float terms[kMaxTerms];
    const std::size_t n_terms = loop_upper + 1;
    if (n_terms > kMaxTerms) {
        return kNegInf;
    }

    for (std::size_t i = 0; i < n_terms; ++i) {
        const std::size_t idx = static_cast<std::size_t>(coi - i - 1);
        if (idx >= log_one_minus_pam.size()) {
            terms[i] = kNegInf;
            continue;
        }
        terms[i] =
            log_w[i] + log_one_minus_pam[idx] +
            log_sum * static_cast<float>(coi - static_cast<int>(i));
    }

    return log_sum_exp(std::span<const float>(terms, n_terms));
}

inline float transmission_log_with_relatedness(const std::vector<double>& pam_vec,
                                               Sampler& sampler,
                                               int coi,
                                               std::size_t total_alleles,
                                               float relatedness,
                                               float log_sum)
{
    thread_local std::vector<float> log_scratch;
    log_scratch.resize(pam_vec.size());
    for (std::size_t i = 0; i < pam_vec.size(); ++i) {
        log_scratch[i] = static_cast<float>(std::log1p(-pam_vec[i]));
    }
    return transmission_log_with_relatedness(
        std::span<const double>(pam_vec.data(), pam_vec.size()),
        std::span<const float>(log_scratch.data(), log_scratch.size()),
        sampler, coi, total_alleles, relatedness, log_sum);
}

} // namespace transmission_incremental
