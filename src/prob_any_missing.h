#pragma once

#ifndef PROBANYMISSING_H
#define PROBANYMISSING_H

#include "combination_indices_generator.h"

struct probAnyMissingFunctor
{
    probAnyMissingFunctor() = default;

    double operator()(const std::vector<double> &eventProbs, int numEvents);

    double eventCombo{};
    CombinationIndicesGenerator c;
};

#endif /* PROBANYMISSING_H */
