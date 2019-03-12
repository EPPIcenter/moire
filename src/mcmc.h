#pragma once

#ifndef MCMC_H_
#define MCMC_H_

#include <Rcpp.h>
#include "genotyping_data.h"
#include "parameters.h"
#include "lookup.h"
#include "chain.h"

class MCMC {
private:


public:
    GenotypingData genotyping_data;
    Lookup lookup;
    Parameters params;

    std::vector<Chain> chains;

    std::vector<std::vector<int> > m_store;
    std::vector<std::vector<std::vector<double > > > p_store;
    std::vector<double> eps_pos_store;
    std::vector<double> eps_neg_store;

    std::vector<std::vector<double > > llik_burnin;
    std::vector<std::vector<double > > llik_sample;
    
    void burnin();
    void sample();

    MCMC(GenotypingData genotyping_data, Lookup lookup, Parameters params);
};


#endif // MCMC_H_