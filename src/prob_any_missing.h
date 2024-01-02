#pragma once

#ifndef PROBANYMISSING_H
#define PROBANYMISSING_H

#include "combination_indices_generator.h"

struct probAnyMissingFunctor
{
    probAnyMissingFunctor() = default;

    double operator()(const std::vector<double> &eventProbs, int numEvents);

    std::vector<double> vectorized(const std::vector<double> &eventProbs,
                                   unsigned int numEvents);

    std::vector<double> vectorized(const std::vector<double> &eventProbs,
                                   unsigned int minNumEvents,
                                   unsigned int maxNumEvents);

    std::vector<double> baseVec{};
    CombinationIndicesGenerator c;
};

#endif /* PROBANYMISSING_H */
