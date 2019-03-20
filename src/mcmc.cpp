
#include <Rcpp.h>

#include "mcmc.h"
#include "chain.h"
#include "mcmc_utils.h"

void MCMC::burnin() {
    UtilFunctions::print("Beginning Burnin");
    UtilFunctions::print("Running Burnin over", chains.size(), "Chains");
    for(size_t i = 0; i < params.num_chains; i++){
        UtilFunctions::print("Running Burnin -- Chain", i);
        std::vector<double> chain_llik(params.burnin);
        llik_burnin.push_back(chain_llik);
        for(int j = 0; j < params.burnin; j++){
            if((j + 1) % 10 == 0) {
                UtilFunctions::print("Burnin Iteration", j + 1);
                Rcpp::checkUserInterrupt();
                UtilFunctions::print("Log Likelihood:", chains[i].get_llik());
            }
            chains[i].update_eps_neg(j + 1);
            chains[i].update_eps_pos(j + 1);
            chains[i].update_p();
            chains[i].update_m(j + 1);
            llik_burnin[i][j] = chains[i].get_llik();
        }
    }
    
}

void MCMC::sample() {
    for(size_t i = 0; i < params.num_chains; i++){
        UtilFunctions::print("Sampling -- Chain", i);
        std::vector<double> chain_llik(params.samples);
        llik_sample.push_back(chain_llik);
        for(int j = 0; j < params.samples; j++){
            if((j + 1) % 10 == 0) {
                UtilFunctions::print("Sampling Iteration", j + 1);
                Rcpp::checkUserInterrupt();
                UtilFunctions::print("Log Likelihood:", chains[i].get_llik());
            }
            chains[i].update_eps_neg(params.burnin + j + 1);
            chains[i].update_eps_pos(params.burnin + j + 1);
            chains[i].update_p();
            chains[i].update_m(params.burnin + j + 1);

            if(params.thin == 0 || j % params.thin == 0) {
                m_store.push_back(chains[i].m);
                p_store.push_back(chains[i].p);
                eps_neg_store.push_back(chains[i].eps_neg);
                eps_pos_store.push_back(chains[i].eps_pos);
                llik_sample[i][j] = chains[i].get_llik();
            }
        }

    }
}


MCMC::MCMC(GenotypingData genotyping_data, Lookup lookup, Parameters params) : 
    genotyping_data(genotyping_data), 
    lookup(lookup), 
    params(params)
    {
        UtilFunctions::print("Starting MCMC...");
        for(size_t i = 0; i < params.num_chains; i++) {
            chains.push_back(Chain(genotyping_data, lookup, params));
        }

        UtilFunctions::print("Created", chains.size(), "Chains");
    };