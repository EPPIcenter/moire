#ifndef COMBINATIONINDICESGENERATOR_H
#define COMBINATIONINDICESGENERATOR_H

#include <vector>

// Adapted from https://stackoverflow.com/a/9432150/2755374
struct CombinationIndicesGenerator
{
    using combination_t = std::vector<int>;

    bool completed;
    unsigned long generated = 1;
    unsigned long numCombinations = 0;

    /**
     * Generate a sequences of indices representing n choose r element
     * combinations.
     * @param n number of elements
     * @param r number of choices
     */
    CombinationIndicesGenerator(int n, int r);

    CombinationIndicesGenerator();

    void reset(int n, int r);

    void next() noexcept;

    combination_t curr{};

   private:
    int n_;
    int r_;

    void calculateNumCombinations() noexcept;
};

#endif /* COMBINATIONINDICESGENERATOR_H */
