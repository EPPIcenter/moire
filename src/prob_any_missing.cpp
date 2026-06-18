#include "prob_any_missing.h"

#include <bit>
#include <cmath>
#include <cstdint>
#include <vector>

double probAnyMissingFunctor::operator()(const std::vector<float>& eventProbs,
                                         int numEvents)
{
    const int totalEvents = static_cast<int>(eventProbs.size());
    if (numEvents < totalEvents) {
        return 1.0;
    }

    double prob = 0.0;

    int sign = -1;
    for (int i = 1; i <= totalEvents; ++i) {
        sign = -sign;
        c.reset(totalEvents, i);
        while (!c.completed) {
            double base = 1.0;

            for (const auto j : c.curr) {
                base -= eventProbs[j];
            }
            c.next();

            double r = sign;
            int multCounter = static_cast<signed>(numEvents);
            while (multCounter > 0) {
                if (multCounter & 1) {
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
    const std::vector<float>& eventProbs, unsigned int numEvents)
{
    return vectorized(eventProbs, 1, numEvents);
}

std::vector<double> probAnyMissingFunctor::vectorized(
    const std::vector<float>& eventProbs, unsigned int minNumEvents,
    unsigned int maxNumEvents)
{
    const std::size_t totalEvents = eventProbs.size();

    std::vector<double> probVec(maxNumEvents - minNumEvents + 1, 0.0);

    if (maxNumEvents < totalEvents) {
        std::fill_n(probVec.begin(), maxNumEvents - minNumEvents + 1, 1.0);
        return probVec;
    }

    std::fill_n(probVec.begin(), totalEvents - 1, 1.0);

    if (minNumEvents == 1) {
        const std::size_t n = totalEvents;
        if (n > 0 && n <= 31) {
            const std::uint32_t numMasks =
                (n == 31) ? 0x7FFFFFFFu : (static_cast<std::uint32_t>(1) << n) - 1;
            double running_sum = 0.0;
            std::uint32_t prev_mask = 0;

            for (std::uint32_t i = 1; i <= numMasks; ++i) {
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

    int sign = -1;
    for (std::size_t i = minNumEvents; i <= totalEvents; ++i) {
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
