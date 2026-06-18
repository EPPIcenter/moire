#pragma once

#include <vector>
#include <array>
#include <concepts>
#include <ranges>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <bit>
#include <type_traits>
#include "multivector.h"  // Includes moire_parallel namespace

namespace zeta_mobius {

template<typename T>
concept ProbabilityCalculator = requires(T calc, const std::vector<double>& probs, std::size_t max_n, std::size_t u_mask) {
    { calc(probs, max_n, u_mask) } -> std::convertible_to<std::vector<double>>;
};

struct MultinomialCalculator {
    std::vector<double> operator()(const std::vector<double>& probs, std::size_t max_n, std::size_t u_mask) const {
        // Calculate probabilities for P(no U = u), that is, the probability that no events are observed in the support vector u.
        std::vector<double> result(max_n + 1, 0.0);
        
        // Calculate the probability that a single event is NOT in the support vector u
        // Optimized: iterate only over set bits in u_mask for better cache behavior
        double p_single_event_not_in_u = 1.0;
        const std::size_t num_categories = probs.size();
        
        // Optimize bit checking: use std::has_single_bit for power-of-2 checks where helpful
        // Iterate over all categories but use bit test - compiler may vectorize this
        for (std::size_t i = 0; i < num_categories; ++i) {
            // Bit test: (u_mask >> i) & 1 is more efficient than u_mask & (1 << i) on some architectures
            if ((u_mask >> i) & 1) {
                p_single_event_not_in_u -= probs[i];
            }
        }
        
        // Ensure the probability is non-negative (clamp to 0 if negative)
        p_single_event_not_in_u = std::max(0.0, p_single_event_not_in_u);
        
        // Calculate probabilities for 1 to max_n events
        // P(no U = u, n) = P(single event not in u)^n
        // Optimized: use incremental multiplication instead of pow() for significant speedup
        // Also: use faster multiplication pattern with better instruction-level parallelism
        if (max_n > 0) {
            result[1] = p_single_event_not_in_u;
            // Unroll first few iterations for better performance
            if (max_n >= 2) {
                result[2] = result[1] * p_single_event_not_in_u;
                if (max_n >= 3) {
                    result[3] = result[2] * p_single_event_not_in_u;
                    if (max_n >= 4) {
                        result[4] = result[3] * p_single_event_not_in_u;
                        // Continue with remaining values
                        for (std::size_t i = 5; i <= max_n; ++i) {
                            result[i] = result[i-1] * p_single_event_not_in_u;
                        }
                    }
                } else {
                    // max_n == 2, already done
                }
            }
        }
        
        return result;
    }
};

template<ProbabilityCalculator Calc>
class MobiusTransform {
private:
    std::vector<double> probs_;
    std::size_t max_n_;
    Calc calculator_;

    // indexed by (u, n) where u is the support vector and n is the number of events
    // cached_p_no_[u][n] = P(no U = u) for n = 1, ..., max_n_
    MultiVector<double, 2> cached_p_no_{};    
    
    // indexed by (s, n) where s is the support vector and n is the number of events
    // cached_support_probs_[s][n] = P(support = s) for n = 1, ..., max_n_
    MultiVector<double, 2> cached_support_probs_{}; 

    void resize_cache(const std::size_t num_categories, const std::size_t max_n) {
        cached_p_no_.resize({static_cast<std::size_t>(1 << num_categories), max_n + 1});
        cached_support_probs_.resize({static_cast<std::size_t>(1 << num_categories), max_n + 1});
    }

    void compute_p_no() {
        const std::size_t num_categories = probs_.size();
        const std::size_t m = 1 << num_categories;
        const std::size_t full = m - 1;

        resize_cache(num_categories, max_n_);

        // Parallelize computation of cached_p_no_ - each mask value is independent
        // Use parallel_for_always since even small workloads (e.g., 2^10 = 1024) benefit from parallelization
        // when each iteration does substantial work (calculator computation + memory fill)
        moire_parallel::parallel_for_always(0, m, [&](std::size_t u) {
            const auto probs = calculator_(probs_, max_n_, u);
            cached_p_no_.inner_fill({u}, std::span<const double>(probs));
        });

        // Parallelize copying of complement masks - independent operations
        moire_parallel::parallel_for_always(0, m, [&](std::size_t s) {
            const std::size_t s_complement = full ^ s;
            const auto [begin, end] = cached_p_no_.inner_iterators({s_complement});
            cached_support_probs_.inner_fill({s}, std::span<const double>(begin, end));
        });
        subset_mobius_inplace(cached_support_probs_, num_categories);
    }

    void subset_mobius_inplace(MultiVector<double, 2>& a, std::size_t num_categories) {
        const std::size_t mask_count = 1 << num_categories;
        for (std::size_t i = 0; i < num_categories; ++i) {
            const std::size_t bit = 1 << i;
            for (std::size_t mask = 0; mask < mask_count; ++mask) {
                if (mask & bit) {
                    auto [dst_begin, dst_end] = a.inner_iterators({mask});
                    auto [src_begin, src_end] = a.inner_iterators({mask ^ bit});
                    const std::size_t len = static_cast<std::size_t>(dst_end - dst_begin);
                    
                    // Optimized: Use restrict pointers and SIMD vectorization
                    // This helps compiler optimize the subtraction loop
                    double* __restrict dst = &*dst_begin;  // Get raw pointer from iterator
                    const double* __restrict src = &*src_begin;
                    
                    // Compiler will auto-vectorize this loop
                    // Use compiler-specific hints for SIMD vectorization
                    #if defined(__GNUC__) || defined(__clang__)
                    #pragma GCC ivdep
                    #elif defined(_MSC_VER)
                    #pragma loop(ivdep)
                    #endif
                    for (std::size_t n = 0; n < len; ++n) {
                        dst[n] -= src[n];
                    }
                }
            }
        }
    }

    std::size_t support_vector_to_mask(const std::span<const int>& support_vector) const {
        std::size_t mask = 0;
        for (std::size_t i = 0; i < support_vector.size(); ++i) {
            if (support_vector[i]) {
                mask |= (1 << i);
            }
        }
        return mask;
    }

    std::vector<int> mask_to_support_vector(std::size_t mask) const {
        std::vector<int> result(probs_.size(), 0);
        for (std::size_t i = 0; i < probs_.size(); ++i) {
            if (mask & (1 << i)) {
                result[i] = 1;
            }
        }
        return result;
    }

public:

    MobiusTransform(const std::vector<double>& probs, const std::size_t max_n, Calc calc)
        : probs_(probs), max_n_(max_n), calculator_(calc) {
        compute_p_no();
    }

    std::vector<double> query_probabilities(const std::span<const int>& support_vector) const {
        // Returns a vector of probabilities for all n from 0 to max_n_ inclusive. The first element 
        // is the probability of a support vector with no events, which is impossible. The empty
        // support vector is unsupported.
        const std::size_t mask = support_vector_to_mask(support_vector);
        // UtilFunctions::print("Mask: ", mask);
        const auto [begin, end] = cached_support_probs_.inner_iterators({mask});
        // UtilFunctions::print_span(std::span<const double>(begin, end));
        return std::vector<double>(begin, end);
    }

    // Optimized accessor: return a span into cached probabilities by integer mask.
    // The span is valid as long as this MobiusTransform instance remains alive.
    std::span<const double> query_probabilities_by_mask(std::size_t mask) const {
        const auto [begin, end] = cached_support_probs_.inner_iterators({mask});
        return std::span<const double>(begin, end);
    }

    void update_parameters(const std::vector<double>& new_probs, const std::size_t max_n) {
        // Updates the category probabilities and recomputes the support probabilities.
        probs_ = new_probs;
        max_n_ = max_n;
        compute_p_no();
    }
};

inline auto create_multinomial_mobius_transform(const std::vector<double>& probs, const std::size_t max_n) {
    return MobiusTransform<MultinomialCalculator>(probs, max_n, MultinomialCalculator{});
}


// Compute P(all categories in S are seen at least once | draws are restricted to S)
// for n = 0..max_n using a Möbius transform on the normalized probabilities q over S.
// Inputs:
// - probs: full category probabilities (not necessarily summing to 1)
// - support_indices: indices of categories in S
// - max_n: maximum number of draws
// Returns: vector<double> of size max_n+1 with conditional probabilities over n
inline std::vector<double> all_seen_restricted_via_mobius(
    const std::vector<double>& probs,
    const std::vector<int>& support_indices,
    const std::size_t max_n)
{
    // Build normalized probabilities within S
    double total_mass_S = 0.0;
    for (int idx : support_indices) {
        total_mass_S += probs[static_cast<std::size_t>(idx)];
    }
    const std::size_t K = support_indices.size();
    std::vector<double> q(K, 0.0);
    if (total_mass_S > 0.0) {
        for (std::size_t i = 0; i < K; ++i) {
            q[i] = probs[static_cast<std::size_t>(support_indices[i])] / total_mass_S;
        }
    }

    // Möbius transform over S-only categories, query mask with all ones
    auto mobius_restricted = create_multinomial_mobius_transform(q, max_n);
    std::vector<int> full_support_mask(K, 1);
    return mobius_restricted.query_probabilities(full_support_mask);
}

// Compute P(any missing in S | draws restricted to S) for n = 0..max_n via Möbius.
// It is simply 1 - P(all seen | restricted S).
inline std::vector<double> any_missing_restricted_via_mobius(
    const std::vector<double>& probs,
    const std::vector<int>& support_indices,
    const std::size_t max_n)
{
    auto all_seen = all_seen_restricted_via_mobius(probs, support_indices, max_n);
    for (double &v : all_seen) {
        v = 1.0 - v;
    }
    return all_seen;
}

// Same layout as probAnyMissingFunctor::vectorized(eventProbs, 1u, max_n):
// returns vector of size max_n with result[j] = P(any missing for numEvents = j+1),
// result[0..K-2] = 1.0. q must be normalized (sum 1) over K categories.
inline std::vector<double> any_missing_from_normalized_q(
    const std::vector<double>& q,
    const std::size_t max_n)
{
    const std::size_t K = q.size();
    std::vector<double> out(max_n, 1.0);
    if (K == 0 || max_n == 0) return out;
    std::vector<int> support_indices(static_cast<std::size_t>(K), 0);
    for (std::size_t i = 0; i < K; ++i) support_indices[i] = static_cast<int>(i);
    auto mob = any_missing_restricted_via_mobius(q, support_indices, max_n);
    for (std::size_t j = 0; j < max_n; ++j) {
        if (j + 1 >= K)
            out[j] = mob[j + 1];
    }
    return out;
}

} // namespace zeta_mobius  