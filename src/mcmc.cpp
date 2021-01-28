
#include <Rcpp.h>

#include "mcmc.h"
#include "chain.h"
#include "mcmc_utils.h"

void MCMC::burnin() {
    UtilFunctions::print("Beginning Burnin");
    UtilFunctions::print("Running Burnin over", chains.size(), "Chains");
    for(size_t i = 0; i < params.num_chains; i++){
        // UtilFunctions::print("Running Burnin -- Chain", i);
        std::vector<double> chain_llik;
        llik_burnin.push_back(chain_llik);
        for(int j = 0; j < params.burnin; j++){
            Rcpp::checkUserInterrupt();
            if((j + 1) % 10 == 0) {
                UtilFunctions::print("Burnin Iteration", j + 1);
                UtilFunctions::print("Log Likelihood:", chains[i].get_llik());
            }
            chains[i].update_eps_neg(j + 1);
            chains[i].update_eps_pos(j + 1);
            // UtilFunctions::print("Eps Pos: ", chains[i].eps_pos);
            // UtilFunctions::print("Eps Neg: ", chains[i].eps_neg);
            // chains[i].update_eps(j + 1);
            chains[i].update_p(j + 1);
            chains[i].update_m(j + 1);
            chains[i].update_mean_coi(j + 1);
            llik_burnin[i].push_back(chains[i].get_llik());
        }
    }

}

void MCMC::sample() {
    for(size_t i = 0; i < params.num_chains; i++){
        // UtilFunctions::print("Sampling -- Chain", i);
        std::vector<double> chain_llik;
        llik_sample.push_back(chain_llik);
        for(int j = 0; j < params.samples; j++){
            Rcpp::checkUserInterrupt();
            if((j + 1) % 10 == 0) {
                UtilFunctions::print("Sampling Iteration", j + 1);
                UtilFunctions::print("Log Likelihood:", chains[i].get_llik());
            }
            chains[i].update_eps_neg(params.burnin + j + 1);
            chains[i].update_eps_pos(params.burnin + j + 1);
            // UtilFunctions::print("Eps Pos: ", chains[i].eps_pos);
            // UtilFunctions::print("Eps Neg: ", chains[i].eps_neg);
            // chains[i].update_eps(j + 1);
            chains[i].update_p(params.burnin + j + 1);
            chains[i].update_m(params.burnin + j + 1);
            chains[i].update_mean_coi(j + 1);

            if(params.thin == 0 || j % params.thin == 0) {
                m_store.push_back(chains[i].m);
                p_store.push_back(chains[i].p);
                eps_neg_store.push_back(chains[i].eps_neg);
                eps_pos_store.push_back(chains[i].eps_pos);
                mean_coi_store.push_back(chains[i].mean_coi);
                llik_sample[i].push_back(chains[i].get_llik());
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
