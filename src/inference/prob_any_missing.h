#pragma once

#include "combination_indices_generator.h"

#include <cstddef>
#include <vector>

struct probAnyMissingFunctor {
    probAnyMissingFunctor() = default;

    double operator()(const std::vector<float>& eventProbs, int numEvents);

    std::vector<double> vectorized(const std::vector<float>& eventProbs,
                                   unsigned int numEvents);

    std::vector<double> vectorized(const std::vector<float>& eventProbs,
                                   unsigned int minNumEvents,
                                   unsigned int maxNumEvents);

    CombinationIndicesGenerator c;
};
