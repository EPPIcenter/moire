
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
            UtilFunctions::print("Iteration", j + 1);
            chains[i].update_m();
            chains[i].update_p();
            chains[i].update_eps_neg();
            chains[i].update_eps_pos();
            chains[i].calculate_llik();
            llik_burnin[i][j] = chains[i].get_llik();
            // chain_llik[j] = chains[i].get_llik();
            // m_store.push_back(chains[i].m);
            // p_store.push_back(chains[i].p);
            // eps_neg_store.push_back(chains[i].eps_neg);
            // eps_pos_store.push_back(chains[i].eps_pos);
        }
    }
    
}

void MCMC::sample() {
    for(size_t i = 0; i < params.num_chains; i++){
        UtilFunctions::print("Sampling -- Chain", i);
        std::vector<double> chain_llik(params.burnin);
        llik_sample.push_back(chain_llik);
        for(int j = 0; j < params.samples; j++){
            UtilFunctions::print("Iteration:", j + 1);
            chains[i].update_m();
            chains[i].update_p();
            chains[i].update_eps_neg();
            chains[i].update_eps_pos();
            chains[i].calculate_llik();
            // chain_llik[j] = chain.get_llik();
            m_store.push_back(chains[i].m);
            p_store.push_back(chains[i].p);
            eps_neg_store.push_back(chains[i].eps_neg);
            eps_pos_store.push_back(chains[i].eps_pos);
            llik_sample[i][j] = chains[i].get_llik();
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