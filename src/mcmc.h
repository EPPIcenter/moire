#pragma once

#ifndef MCMC_H_
#define MCMC_H_

#include "chain.h"
#include "genotyping_data.h"
#include "lookup.h"
#include "parameters.h"

#include <Rcpp.h>

class MCMC
{
   private:
    Chain chain;

   public:
    GenotypingData genotyping_data;
    Lookup lookup;
    Parameters params;

    std::vector<std::vector<int>> m_store{};
    std::vector<std::vector<std::vector<double>>> p_store{};
    std::vector<std::vector<double>> eps_pos_store{};
    std::vector<std::vector<double>> eps_neg_store{};
    std::vector<double> mean_coi_store{};

    std::vector<double> llik_burnin{};
    std::vector<double> llik_sample{};

    void burnin();
    void sample();

    MCMC(GenotypingData genotyping_data, Lookup lookup, Parameters params);
};

#endif  // MCMC_H_
