#pragma once

#ifndef MCMC_H_
#define MCMC_H_

#include "chain.h"
#include "genotyping_data.h"
#include "parameters.h"

#include <Rcpp.h>
#include <numeric>
#include <progress.hpp>

class MCMC
{
   private:
   public:
    GenotypingData genotyping_data;
    Parameters params;
    std::vector<Chain> chains{};

    std::vector<std::vector<int>> m_store{};
    // indexed by population, locus, step
    std::vector<std::vector<std::vector<std::vector<float>>>> p_store{};
    std::vector<std::vector<std::vector<std::vector<int>>>> latent_genotypes_store{};
    std::vector<std::vector<float>> data_llik_store{};
    std::vector<std::vector<float>> eps_pos_store{};
    // Per-locus pooled false-positive draws (marginal_ecoi mode); indexed by
    // locus, step.
    std::vector<std::vector<float>> eps_pos_locus_store{};
    std::vector<std::vector<float>> eps_neg_store{};
    std::vector<std::vector<float>> r_store{};
    // Effective COI draws (marginal_ecoi mode); indexed by sample, step.
    std::vector<std::vector<float>> eff_coi_store{};
    // Per-population e-hierarchy parameter draws (marginal_ecoi mode); each entry
    // is the per-population vector at a stored step (indexed by step, population).
    std::vector<std::vector<float>> ecoi_mu_plus_store{};
    std::vector<std::vector<float>> ecoi_k_store{};
    std::vector<float> population_coi_p_store{};
    std::vector<float> population_coi_r_store{};
    std::vector<std::vector<float>> population_responsibility_store{};
    // indexed by sample, step
    std::vector<std::vector<std::vector<float>>> population_assignment_store{};
    std::vector<int> swap_store{};

    std::vector<float> llik_burnin{};
    std::vector<float> llik_sample{};
    std::vector<float> prior_burnin{};
    std::vector<float> prior_sample{};
    std::vector<float> posterior_burnin{};
    std::vector<float> posterior_sample{};

    std::vector<size_t> swap_indices{};
    std::vector<size_t> swap_acceptances{};
    std::vector<float> swap_barriers{};
    std::vector<float> temp_gradient{};
    size_t num_swaps = 0;
    bool even_swap = true;

    void burnin(int step);
    void sample(int step);
    // void add_chain(float temp);
    void swap_chains(int step, bool burnin = false);
    void finalize();
    void adapt_temp();
    float get_llik();
    float get_prior();
    float get_posterior();
    int get_hot_chain();

    MCMC(GenotypingData genotyping_data, Parameters params);
};

#endif  // MCMC_H_
