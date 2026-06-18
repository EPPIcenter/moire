#pragma once

#include <vector>

struct LatentGenotype {
    std::vector<int> value;
    float log_prob;
};
