#pragma once

#ifndef PROBANYMISSING_H
#define PROBANYMISSING_H

#include "combination_indices_generator.h"

#include <cstddef>
#include <vector>

/// Optional cache for P(any missing) when the same normalized q is reused
/// (e.g. same locus/population across samples). Use with vectorized_cached.
/// Not thread-safe; use one per chain or per thread.
struct ProbAnyMissingCache
{
    static constexpr unsigned DEFAULT_MAX_N = 100u;

    void clear() noexcept;

    /// Returns cached vector for this q and coi if we have a valid entry (same q, cached_max_n >= coi).
    /// Otherwise returns nullptr and *out_max_n is unchanged.
    const double* get(const std::vector<double>& q, unsigned coi, unsigned* out_max_n) const noexcept;

    /// Store result for this q so future get(q, coi) with coi <= max_n can reuse.
    void set(std::vector<double> q, unsigned max_n, std::vector<double> result);

private:
    std::vector<double> cached_q_{};
    unsigned cached_max_n_{0};
    std::vector<double> cached_result_{};
};

struct probAnyMissingFunctor
{
    probAnyMissingFunctor() = default;

    double operator()(const std::vector<float> &eventProbs, int numEvents);

    std::vector<double> vectorized(const std::vector<float> &eventProbs,
                                   unsigned int numEvents);

    std::vector<double> vectorized(const std::vector<float> &eventProbs,
                                   unsigned int minNumEvents,
                                   unsigned int maxNumEvents);

    /// When cache is non-null: use cache on hit; on miss compute via Möbius and store.
    /// When cache is null: same as vectorized() (Gray-code). Use when q is stable across calls.
    std::vector<double> vectorized_cached(const std::vector<float> &eventProbs,
                                         unsigned int maxNumEvents,
                                         ProbAnyMissingCache* cache);

    /// Same result as vectorized(eventProbs, 1u, maxNumEvents) but always uses
    /// the combination-based inclusion-exclusion (for benchmarking vs Gray-code).
    std::vector<double> vectorized_combination(const std::vector<float> &eventProbs,
                                               unsigned int maxNumEvents);

    /// Same result as vectorized(eventProbs, 1u, maxNumEvents) via Möbius transform
    /// (parallel over 2^K masks). For benchmarking and optional use in hot path.
    std::vector<double> vectorized_mobius(const std::vector<float> &eventProbs,
                                         unsigned int maxNumEvents);

    CombinationIndicesGenerator c;
};

#endif /* PROBANYMISSING_H */
