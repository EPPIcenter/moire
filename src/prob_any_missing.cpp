#include "prob_any_missing.h"
#include "profiler.h"
#include "zeta_mobius.h"

#include <bit>
#include <cmath>
#include <cstdint>
#include <vector>

namespace {

bool q_equal(const std::vector<double>& a, const std::vector<double>& b) noexcept
{
    if (a.size() != b.size()) return false;
    for (std::size_t i = 0; i < a.size(); ++i)
        if (a[i] != b[i]) return false;
    return true;
}

} // namespace

void ProbAnyMissingCache::clear() noexcept
{
    cached_q_.clear();
    cached_max_n_ = 0;
    cached_result_.clear();
}

const double* ProbAnyMissingCache::get(const std::vector<double>& q, unsigned coi,
                                       unsigned* out_max_n) const noexcept
{
    if (cached_q_.empty() || !q_equal(q, cached_q_) || cached_max_n_ < coi)
        return nullptr;
    if (out_max_n) *out_max_n = cached_max_n_;
    return cached_result_.data();
}

void ProbAnyMissingCache::set(std::vector<double> q, unsigned max_n,
                              std::vector<double> result)
{
    cached_q_ = std::move(q);
    cached_max_n_ = max_n;
    cached_result_ = std::move(result);
}

// todo: try implementing this with simd intrinsics by precomputing the event
// probabilities
double probAnyMissingFunctor::operator()(const std::vector<float> &eventProbs,
                                        int numEvents)
{
    const int totalEvents = eventProbs.size();
    if (numEvents < totalEvents)
    {
        // early exit if impossible
        return 1.0;
    }

    double prob = 0.0;

    // Calculate via inclusion-exclusion principle
    int sign = -1;
    for (int i = 1; i <= totalEvents; ++i)
    {
        sign = -sign;
        c.reset(totalEvents, i);
        while (!c.completed)
        {
            double base = 1.0;

            for (const auto j : c.curr)
            {
                base -= eventProbs[j];
            }
            c.next();

            double r = sign;
            int multCounter = static_cast<signed>(numEvents);
            // squared exponentiation
            while (multCounter > 0)
            {
                if (multCounter & 1)
                {
                    r *= base;
                }
                base = (base * base);
                multCounter >>= 1;
            }
            prob += r;
        }
    }
    return prob;
}

std::vector<double> probAnyMissingFunctor::vectorized(
    const std::vector<float> &eventProbs, unsigned int numEvents)
{
    return vectorized(eventProbs, 1, numEvents);
}

std::vector<double> probAnyMissingFunctor::vectorized(
    const std::vector<float> &eventProbs, unsigned int minNumEvents,
    unsigned int maxNumEvents)
{
    const std::size_t totalEvents = eventProbs.size();

    std::vector<double> probVec(maxNumEvents - minNumEvents + 1, 0.0);
    
    if (maxNumEvents < totalEvents) {
        std::fill_n(probVec.begin(), maxNumEvents - minNumEvents + 1, 1.0);
        return probVec;
    }

    std::fill_n(probVec.begin(), totalEvents - 1, 1.0);

    // Fast path: minNumEvents == 1 → Gray-code iteration (O(1) running sum per subset)
    if (minNumEvents == 1) {
        const std::size_t n = totalEvents;
        if (n > 0 && n <= 31) {
            const std::uint32_t numMasks = (n == 31) ? 0x7FFFFFFFu : (static_cast<std::uint32_t>(1) << n) - 1;
            double running_sum = 0.0;
            std::uint32_t prev_mask = 0;

            for (std::uint32_t i = 1; i <= numMasks; ++i)
            {
                const std::uint32_t mask = i ^ (i >> 1);

                if (i == 1) {
                    running_sum = static_cast<double>(eventProbs[0]);
                } else {
                    const std::uint32_t changed = mask ^ prev_mask;
                    const int idx = static_cast<int>(std::countr_zero(changed));
                    if ((mask & changed) != 0) {
                        running_sum += static_cast<double>(eventProbs[static_cast<std::size_t>(idx)]);
                    } else {
                        running_sum -= static_cast<double>(eventProbs[static_cast<std::size_t>(idx)]);
                    }
                }
                prev_mask = mask;

                const double base = 1.0 - running_sum;
                const int sign = (std::popcount(mask) & 1) ? 1 : -1;

                double r = static_cast<double>(sign);
                for (std::size_t j = 0; j < totalEvents; ++j)
                    r *= base;
                for (std::size_t j = totalEvents - 1; j < maxNumEvents; ++j) {
                    probVec[j] += r;
                    r *= base;
                }
            }
            return probVec;
        }
    }

    // Fallback: minNumEvents > 1 or n > 31 → combination-based inclusion-exclusion
    int sign = -1;
    for (std::size_t i = minNumEvents; i <= totalEvents; ++i)
    {
        sign = -sign;
        c.reset(static_cast<int>(totalEvents), static_cast<int>(i));

        for (std::size_t k = 0; k < c.numCombinations; ++k)
        {
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

std::vector<double> probAnyMissingFunctor::vectorized_mobius(
    const std::vector<float> &eventProbs, unsigned int maxNumEvents)
{
    const std::size_t K = eventProbs.size();
    if (maxNumEvents < K) {
        std::vector<double> out(maxNumEvents, 1.0);
        return out;
    }
    std::vector<double> q(K);
    for (std::size_t i = 0; i < K; ++i)
        q[i] = static_cast<double>(eventProbs[i]);
    return zeta_mobius::any_missing_from_normalized_q(q, static_cast<std::size_t>(maxNumEvents));
}

std::vector<double> probAnyMissingFunctor::vectorized_cached(
    const std::vector<float> &eventProbs, unsigned int maxNumEvents,
    ProbAnyMissingCache* cache)
{
    const std::size_t K = eventProbs.size();
    if (maxNumEvents < K) {
        std::vector<double> out(maxNumEvents, 1.0);
        return out;
    }
    std::vector<double> q(K);
    for (std::size_t i = 0; i < K; ++i)
        q[i] = static_cast<double>(eventProbs[i]);

    if (cache) {
        unsigned cached_max_n = 0;
        const double* ptr = cache->get(q, maxNumEvents, &cached_max_n);
        if (ptr) {
            ProfileScope _("prob_any_missing::vectorized_cached::cache_hit");
            return std::vector<double>(ptr, ptr + maxNumEvents);
        }
        {
            ProfileScope _("prob_any_missing::vectorized_cached::mobius_miss");
            const unsigned store_max_n = (maxNumEvents > ProbAnyMissingCache::DEFAULT_MAX_N)
                                             ? maxNumEvents
                                             : ProbAnyMissingCache::DEFAULT_MAX_N;
            std::vector<double> result = zeta_mobius::any_missing_from_normalized_q(q, static_cast<std::size_t>(store_max_n));
            cache->set(std::move(q), store_max_n, result);
            if (store_max_n >= maxNumEvents)
                return std::vector<double>(result.begin(), result.begin() + static_cast<std::ptrdiff_t>(maxNumEvents));
            return std::vector<double>(result.begin(), result.end());
        }
    }

    {
        ProfileScope _("prob_any_missing::vectorized_cached::gray");
        return vectorized(eventProbs, 1u, maxNumEvents);
    }
}

std::vector<double> probAnyMissingFunctor::vectorized_combination(
    const std::vector<float> &eventProbs, unsigned int maxNumEvents)
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