#include "prob_any_missing.h"

#include <cmath>
#include <vector>

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

    //      Calculate via inclusion-exclusion principle
    int sign = -1;
    for (std::size_t i = minNumEvents; i <= totalEvents; ++i)
    {
        sign = -sign;
        c.reset(totalEvents, i);

        for (std::size_t k = 0; k < c.numCombinations; ++k)
        {
            float base = 1.0;

            for (const auto j : c.curr)
            {
                base -= eventProbs[j];
            }
            c.next();

            float r = sign;
            for (std::size_t j = 0; j < totalEvents - 1; ++j)
            {
                r *= base;
            }
            for (std::size_t j = totalEvents - 1; j < maxNumEvents; ++j)
            {
                r *= base;
                probVec[j] += r;
            }
        }
    }
    return probVec;
}