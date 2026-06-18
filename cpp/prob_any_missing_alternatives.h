#pragma once

// Reference algorithms for benchmarks and cross-validation only.
// Production MCMC uses probAnyMissingFunctor::vectorized() (Gray-code).

#include "combination_indices_generator.h"
#include "prob_any_missing.h"
#include "zeta_mobius.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

namespace prob_any_missing_bench {

namespace detail {

inline bool q_equal(const std::vector<double>& a, const std::vector<double>& b) noexcept
{
    if (a.size() != b.size()) return false;
    for (std::size_t i = 0; i < a.size(); ++i)
        if (a[i] != b[i]) return false;
    return true;
}

} // namespace detail

struct ProbAnyMissingCache {
    static constexpr unsigned DEFAULT_MAX_N = 100u;

    void clear() noexcept
    {
        cached_q_.clear();
        cached_max_n_ = 0;
        cached_result_.clear();
    }

    const double* get(const std::vector<double>& q, unsigned coi,
                      unsigned* out_max_n) const noexcept
    {
        if (cached_q_.empty() || !detail::q_equal(q, cached_q_) || cached_max_n_ < coi)
            return nullptr;
        if (out_max_n) *out_max_n = cached_max_n_;
        return cached_result_.data();
    }

    void set(std::vector<double> q, unsigned max_n, std::vector<double> result)
    {
        cached_q_ = std::move(q);
        cached_max_n_ = max_n;
        cached_result_ = std::move(result);
    }

private:
    std::vector<double> cached_q_{};
    unsigned cached_max_n_{0};
    std::vector<double> cached_result_{};
};

struct ReferenceFunctor {
    CombinationIndicesGenerator c;

    std::vector<double> vectorized_combination(const std::vector<float>& eventProbs,
                                               unsigned int maxNumEvents)
    {
        const std::size_t totalEvents = eventProbs.size();
        std::vector<double> probVec(maxNumEvents, 0.0);
        if (maxNumEvents < totalEvents) {
            std::fill_n(probVec.begin(), maxNumEvents, 1.0);
            return probVec;
        }
        std::fill_n(probVec.begin(), totalEvents - 1, 1.0);

        int sign = -1;
        for (std::size_t i = 1; i <= totalEvents; ++i) {
            sign = -sign;
            c.reset(static_cast<int>(totalEvents), static_cast<int>(i));
            for (std::size_t k = 0; k < c.numCombinations; ++k) {
                float base = 1.0f;
                for (const auto j : c.curr)
                    base -= eventProbs[j];
                c.next();

                float r = static_cast<float>(sign);
                for (std::size_t j = 0; j < totalEvents - 1; ++j)
                    r *= base;
                for (std::size_t j = totalEvents - 1; j < maxNumEvents; ++j) {
                    r *= base;
                    probVec[j] += static_cast<double>(r);
                }
            }
        }
        return probVec;
    }

    std::vector<double> vectorized_mobius(const std::vector<float>& eventProbs,
                                          unsigned int maxNumEvents)
    {
        const std::size_t K = eventProbs.size();
        if (maxNumEvents < K) {
            return std::vector<double>(maxNumEvents, 1.0);
        }
        std::vector<double> q(K);
        for (std::size_t i = 0; i < K; ++i)
            q[i] = static_cast<double>(eventProbs[i]);
        return zeta_mobius::any_missing_from_normalized_q(q, static_cast<std::size_t>(maxNumEvents));
    }

    std::vector<double> vectorized_cached(const std::vector<float>& eventProbs,
                                          unsigned int maxNumEvents,
                                          ProbAnyMissingCache* cache)
    {
        const std::size_t K = eventProbs.size();
        if (maxNumEvents < K) {
            return std::vector<double>(maxNumEvents, 1.0);
        }
        std::vector<double> q(K);
        for (std::size_t i = 0; i < K; ++i)
            q[i] = static_cast<double>(eventProbs[i]);

        probAnyMissingFunctor production;

        if (cache) {
            unsigned cached_max_n = 0;
            const double* ptr = cache->get(q, maxNumEvents, &cached_max_n);
            if (ptr) {
                return std::vector<double>(ptr, ptr + maxNumEvents);
            }
            const unsigned store_max_n = (maxNumEvents > ProbAnyMissingCache::DEFAULT_MAX_N)
                                             ? maxNumEvents
                                             : ProbAnyMissingCache::DEFAULT_MAX_N;
            std::vector<double> result =
                zeta_mobius::any_missing_from_normalized_q(q, static_cast<std::size_t>(store_max_n));
            cache->set(std::move(q), store_max_n, result);
            if (store_max_n >= maxNumEvents)
                return std::vector<double>(result.begin(),
                                           result.begin() + static_cast<std::ptrdiff_t>(maxNumEvents));
            return std::vector<double>(result.begin(), result.end());
        }

        return production.vectorized(eventProbs, maxNumEvents);
    }
};

} // namespace prob_any_missing_bench
