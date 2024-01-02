#include "prob_any_missing.h"

#include <cmath>
#include <vector>

// todo: try implementing this with simd intrinsics by precomputing the event
// probabilities
double probAnyMissingFunctor::operator()(const std::vector<double> &eventProbs,
                                         int numEvents)
{
    const std::size_t totalEvents = eventProbs.size();
    if (numEvents < totalEvents)
    {
        // early exit if impossible
        return 1.0;
    }

    double prob = 0.0;

    //      Calculate via inclusion-exclusion principle
    int sign = -1;
    for (std::size_t i = 1; i <= totalEvents; ++i)
    {
        sign = -sign;
        c.reset(totalEvents, i);
        baseVec.resize(0);
        while (!c.completed)
        {
            double base = 1.0;

            for (const auto j : c.curr)
            {
                base -= eventProbs[j];
            }
            baseVec.push_back(base);
            c.next();
        }

        for (double base : baseVec)
        {
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
    const std::vector<double> &eventProbs, unsigned int numEvents)
{
    return vectorized(eventProbs, 1, numEvents);
}

std::vector<double> probAnyMissingFunctor::vectorized(
    const std::vector<double> &eventProbs, unsigned int minNumEvents,
    unsigned int maxNumEvents)
{
    const std::size_t totalEvents = eventProbs.size();

    std::vector<double> probVec(maxNumEvents - minNumEvents + 1, 0.0);
    std::fill_n(probVec.begin(), totalEvents - 1, 1.0);

    //      Calculate via inclusion-exclusion principle
    int sign = -1;
    for (std::size_t i = minNumEvents; i <= totalEvents; ++i)
    {
        sign = -sign;
        c.reset(totalEvents, i);
        baseVec.clear();
        baseVec.reserve(c.numCombinations);
        while (!c.completed)
        {
            double base = 1.0;

            for (const auto j : c.curr)
            {
                base -= eventProbs[j];
            }
            baseVec.push_back(base);
            c.next();
        }

        for (const double base : baseVec)
        {
            double r = sign;
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