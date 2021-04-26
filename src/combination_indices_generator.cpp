#include "combination_indices_generator.h"

#include <numeric>

CombinationIndicesGenerator::CombinationIndicesGenerator(int n, int r)
    : completed(n < 1 or r > n or r == 0), n_(n), r_(r)
{
    curr.resize(r_);
    std::iota(curr.begin(), curr.end(), 0);
}

void CombinationIndicesGenerator::reset(int n, int r)
{
    completed = n < 1 or r > n or r == 0;
    generated = 1;

    n_ = n;
    r_ = r;

    curr.resize(r_);
    std::iota(curr.begin(), curr.end(), 0);
}

void CombinationIndicesGenerator::next() noexcept
{
    completed = true;
    const int lim = n_ - r_;
    for (int i = r_ - 1; i >= 0; --i)
        if (curr[i] < lim + i)
        {
            int j = curr[i] + 1;
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
