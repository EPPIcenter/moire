#include "prob_any_missing.h"

#include <cmath>
#include <vector>

double probAnyMissingFunctor::operator()(const std::vector<double> &eventProbs,
                                         int numEvents)
{
    int totalEvents = eventProbs.size();
    int multCounter;
    double r;
    double prob = 0.0;

    if (numEvents < totalEvents)
    {
        return 1.0;
    }

    // Calculate via inclusion-exclusion principle
    int sign = -1;
    for (int i = 1; i <= totalEvents; ++i)
    {
        sign = -sign;
        c.reset(totalEvents, i);
        while (!c.completed)
        {
            base = 1.0;
            multCounter = numEvents;

            for (const auto j : c.curr)
            {
                base -= eventProbs[j];
            }

            r = sign;

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
            c.next();
        }
    }

    return prob;
}
