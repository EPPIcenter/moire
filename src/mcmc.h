#pragma once

#ifndef MCMC_H_
#define MCMC_H_

#include "chain.h"
#include "genotyping_data.h"
#include "parameters.h"

#include <Rcpp.h>
#include <progress.hpp>

class MCMC
{
   private:
   public:
    GenotypingData genotyping_data;
    Parameters params;
    std::vector<Chain> chains{};

    std::vector<std::vector<int>> m_store{};
    std::vector<std::vector<std::vector<double>>> p_store{};
    std::vector<std::vector<double>> eps_pos_store{};
    std::vector<std::vector<double>> eps_neg_store{};
    std::vector<std::vector<double>> r_store{};
    std::vector<double> mean_coi_store{};
    std::vector<double> swap_store{};

    std::vector<double> llik_burnin{};
    std::vector<double> llik_sample{};
    std::vector<size_t> swap_indices{};
    size_t num_swaps = 0;

    void burnin(int step);
    void sample(int step);
    void swap_chains();
    double get_llik();

    MCMC(GenotypingData genotyping_data, Parameters params);
};

#endif  // MCMC_H_
