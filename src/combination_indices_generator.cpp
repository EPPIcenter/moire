#include "combination_indices_generator.h"

#include <numeric>

CombinationIndicesGenerator::CombinationIndicesGenerator(int n, int r)
    : completed(n < 1 or r > n or r == 0), n_(n), r_(r)
{
    curr.resize(r_);
    std::iota(curr.begin(), curr.end(), 0);
    calculateNumCombinations();
}

void CombinationIndicesGenerator::reset(int n, int r)
{
    completed = n < 1 or r > n or r == 0;
    generated = 1;

    n_ = n;
    r_ = r;

    curr.resize(r_);
    std::iota(curr.begin(), curr.end(), 0);
    calculateNumCombinations();
}

void CombinationIndicesGenerator::next() noexcept
{
    completed = true;
    const int lim = n_ - r_;
    for (int i = r_ - 1; i >= 0; --i)
        if (curr[i] < lim + i)
        {
            char j = curr[i] + 1;
            while (i < r_)
            {
                curr[i++] = j++;
            }
            completed = false;
            generated++;
            break;
        }
}

CombinationIndicesGenerator::CombinationIndicesGenerator()
{
    completed = true;
    n_ = 0;
    r_ = 0;
}

void CombinationIndicesGenerator::calculateNumCombinations() noexcept
{
    int tmp = r_;
    if (tmp > n_)
    {
        numCombinations = 0;
        return;
    }

    if (tmp * 2 > n_)
    {
        tmp = n_ - tmp;
    }

    if (tmp == 0)
    {
        numCombinations = 1;
        return;
    }

    numCombinations = n_;
    for (int i = 2; i <= tmp; ++i)
    {
        numCombinations *= (n_ - i + 1);
        numCombinations /= i;
    }
}