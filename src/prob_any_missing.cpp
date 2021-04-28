#include "prob_any_missing.h"

#include <cmath>
#include <vector>

double probAnyMissingFunctor::operator()(const std::vector<double> &eventProbs,
                                         int numEvents)
{
    int totalEvents = eventProbs.size();
    int multCounter;
    double r;

    if (numEvents < totalEvents)
    {
        return 1.0;
    }

    prob = 0.0;

    // Calculate via inclusion-exclusion principle
    int sign = -1;
    for (int i = 1; i <= totalEvents; ++i)
    {
        sign = -sign;
        c.reset(totalEvents, i);
        while (!c.completed)
        {
            eventCombo = 0.0;
            multCounter = numEvents;

            for (const auto j : c.curr)
            {
                eventCombo += eventProbs[j];
            }

            r = sign;

            while (multCounter > 0)
            {
                r *= (1 - eventCombo);
                --multCounter;
            }

            prob += r;
            c.next();
        }
    }

    return prob;
}
