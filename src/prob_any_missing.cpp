#include "prob_any_missing.h"

#include <algorithm>
#include <cmath>
#include <vector>

namespace
{
// Hot-path threshold: float inclusion-exclusion for small latent genotypes.
// Larger k uses EGF product (stable, O(k n^2), much faster than 2^k).
constexpr int k_ie_max = 10;

float inclusion_exclusion_float(probAnyMissingFunctor &self,
                                const std::vector<float> &eventProbs,
                                int numEvents)
{
    const int totalEvents = static_cast<int>(eventProbs.size());
    float prob = 0.0f;
    int sign = -1;
    for (int i = 1; i <= totalEvents; ++i)
    {
        sign = -sign;
        self.c.reset(totalEvents, i);
        self.baseVec.resize(0);
        while (!self.c.completed)
        {
            float base = 1.0f;
            for (const auto j : self.c.curr)
            {
                base -= eventProbs[j];
            }
            self.baseVec.push_back(base);
            self.c.next();
        }

        for (float base : self.baseVec)
        {
            float r = static_cast<float>(sign);
            int multCounter = numEvents;
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

// a_m = m! [x^m] prod_i (e^{p_i x} - 1); coverage(n) = a_n.
// Recurrence: new_a[m] = sum_{j=1..m} C(m,j) a[m-j] p^j, new_a[0] = 0.
void fill_egf_coverage(probAnyMissingFunctor &self,
                       const std::vector<float> &eventProbs, int max_n)
{
    const int k = static_cast<int>(eventProbs.size());
    self.egf_a.assign(static_cast<size_t>(max_n) + 1, 0.0);
    self.egf_b.assign(static_cast<size_t>(max_n) + 1, 0.0);
    self.egf_a[0] = 1.0;

    self.egf_p_pow.resize(static_cast<size_t>(max_n) + 1);

    for (int allele = 0; allele < k; ++allele)
    {
        const double p = eventProbs[static_cast<size_t>(allele)];
        std::fill(self.egf_b.begin(), self.egf_b.end(), 0.0);

        self.egf_p_pow[0] = 0.0;
        if (max_n >= 1)
        {
            self.egf_p_pow[1] = p;
        }
        for (int j = 2; j <= max_n; ++j)
        {
            self.egf_p_pow[static_cast<size_t>(j)] =
                self.egf_p_pow[static_cast<size_t>(j - 1)] * p;
        }

        for (int m = 1; m <= max_n; ++m)
        {
            double comb = static_cast<double>(m);  // C(m,1)
            double sum = 0.0;
            for (int j = 1; j <= m; ++j)
            {
                sum += comb * self.egf_a[static_cast<size_t>(m - j)] *
                       self.egf_p_pow[static_cast<size_t>(j)];
                comb *= static_cast<double>(m - j) / static_cast<double>(j + 1);
            }
            self.egf_b[static_cast<size_t>(m)] = sum;
        }
        self.egf_a.swap(self.egf_b);
    }
}

float pam_from_coverage(double coverage)
{
    // Keep result in [0,1]; tiny negatives from roundoff → 0 coverage overshoot.
    const double pam = 1.0 - coverage;
    if (pam <= 0.0)
    {
        return 0.0f;
    }
    if (pam >= 1.0)
    {
        return 1.0f;
    }
    return static_cast<float>(pam);
}

float exact_pam_n_equals_k(const std::vector<float> &eventProbs)
{
    const int k = static_cast<int>(eventProbs.size());
    double cov = 1.0;
    for (float p : eventProbs)
    {
        cov *= p;
    }
    for (int i = 2; i <= k; ++i)
    {
        cov *= i;
    }
    return pam_from_coverage(cov);
}

}  // namespace

float probAnyMissingFunctor::operator()(const std::vector<float> &eventProbs,
                                        int numEvents)
{
    const int k = static_cast<int>(eventProbs.size());
    if (numEvents < k)
    {
        return 1.0f;
    }
    if (k == 0)
    {
        return 0.0f;
    }
    if (numEvents == k)
    {
        return exact_pam_n_equals_k(eventProbs);
    }
    if (k <= k_ie_max)
    {
        return inclusion_exclusion_float(*this, eventProbs, numEvents);
    }

    fill_egf_coverage(*this, eventProbs, numEvents);
    return pam_from_coverage(egf_a[static_cast<size_t>(numEvents)]);
}

std::vector<float> probAnyMissingFunctor::vectorized(
    const std::vector<float> &eventProbs, unsigned int numEvents)
{
    return vectorized(eventProbs, 1, numEvents);
}

std::vector<float> probAnyMissingFunctor::vectorized(
    const std::vector<float> &eventProbs, unsigned int minNumEvents,
    unsigned int maxNumEvents)
{
    const std::size_t totalEvents = eventProbs.size();
    const int k = static_cast<int>(totalEvents);

    std::vector<float> probVec(maxNumEvents - minNumEvents + 1, 0.0f);

    if (maxNumEvents < static_cast<unsigned int>(k))
    {
        std::fill(probVec.begin(), probVec.end(), 1.0f);
        return probVec;
    }

    // Small-k: keep the existing float IE vectorized path (hot relatedness path).
    if (k <= k_ie_max)
    {
        std::fill_n(probVec.begin(), totalEvents - 1, 1.0f);

        int sign = -1;
        for (std::size_t i = minNumEvents; i <= totalEvents; ++i)
        {
            sign = -sign;
            c.reset(totalEvents, static_cast<int>(i));
            baseVec.clear();
            baseVec.reserve(c.numCombinations);

            for (std::size_t kk = 0; kk < c.numCombinations; ++kk)
            {
                float base = 1.0f;
                for (const auto j : c.curr)
                {
                    base -= eventProbs[j];
                }
                baseVec.push_back(base);
                c.next();
            }

            for (const float base : baseVec)
            {
                float r = static_cast<float>(sign);
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

    // Large-k: one EGF build gives coverage for all n <= maxNumEvents.
    fill_egf_coverage(*this, eventProbs, static_cast<int>(maxNumEvents));
    for (unsigned int n = minNumEvents; n <= maxNumEvents; ++n)
    {
        const size_t out_idx = n - minNumEvents;
        if (static_cast<int>(n) < k)
        {
            probVec[out_idx] = 1.0f;
        }
        else
        {
            probVec[out_idx] =
                pam_from_coverage(egf_a[static_cast<size_t>(n)]);
        }
    }
    return probVec;
}
