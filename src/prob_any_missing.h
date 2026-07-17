#pragma once

#ifndef PROBANYMISSING_H
#define PROBANYMISSING_H

#include "combination_indices_generator.h"

#include <vector>

struct probAnyMissingFunctor
{
    probAnyMissingFunctor() = default;

    float operator()(const std::vector<float> &eventProbs, int numEvents);

    std::vector<float> vectorized(const std::vector<float> &eventProbs,
                                  unsigned int numEvents);

    std::vector<float> vectorized(const std::vector<float> &eventProbs,
                                  unsigned int minNumEvents,
                                  unsigned int maxNumEvents);

    std::vector<float> baseVec{};
    CombinationIndicesGenerator c;

    // Scratch for EGF coverage path (avoids per-call allocation).
    std::vector<double> egf_a{};
    std::vector<double> egf_b{};
    std::vector<double> egf_p_pow{};
};

#endif /* PROBANYMISSING_H */
