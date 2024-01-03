#pragma once

#ifndef PROBANYMISSING_H
#define PROBANYMISSING_H

#include "combination_indices_generator.h"

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
};

#endif /* PROBANYMISSING_H */
