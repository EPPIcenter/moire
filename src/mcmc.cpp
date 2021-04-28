#include "mcmc.h"

#include "chain.h"
#include "mcmc_utils.h"

#include <Rcpp.h>
#include <algorithm>

void MCMC::burnin()
{
    if (params.verbose)
    {
        UtilFunctions::print("Beginning burnin");
        UtilFunctions::print("Running burnin over", params.burnin, "samples");
    }

    for (int j = 0; j < params.burnin; j++)
    {
        Rcpp::checkUserInterrupt();
        if (params.verbose)
        {
            UtilFunctions::print("Burnin iteration", j + 1);
            UtilFunctions::print("Log likelihood:", chain.get_llik());
        }

        chain.update_eps_neg(j + 1);
        chain.update_eps_pos(j + 1);
        // chain.update_eps(j + 1);
        chain.update_p(j + 1);
        chain.update_m(j + 1);
        chain.update_individual_parameters(j + 1);
        // chain.update_mean_coi(j + 1);
        llik_burnin.push_back(chain.get_llik());
    }
}

void MCMC::sample()
{
    if (params.verbose)
    {
        UtilFunctions::print("Beginning sampling");
        UtilFunctions::print("Sampling for ", params.samples, "samples");
    }
    for (int j = 0; j < params.samples; j++)
    {
        Rcpp::checkUserInterrupt();
        if (params.verbose)
        {
            UtilFunctions::print("Sampling Iteration", j + 1);
            UtilFunctions::print("Log Likelihood:", chain.get_llik());
        }
        chain.update_eps_neg(params.burnin + j + 1);
        chain.update_eps_pos(params.burnin + j + 1);
        // chain.update_eps(j + 1);
        chain.update_p(params.burnin + j + 1);
        chain.update_m(params.burnin + j + 1);
        chain.update_individual_parameters(params.burnin + j + 1);
        // chain.update_mean_coi(j + 1);

        if (params.thin == 0 || j % params.thin == 0)
        {
            m_store.push_back(chain.m);
            p_store.push_back(chain.p);
            eps_neg_store.push_back(chain.eps_neg);
            eps_pos_store.push_back(chain.eps_pos);
            mean_coi_store.push_back(chain.mean_coi);
            llik_sample.push_back(chain.get_llik());
        }
    }
}

MCMC::MCMC(GenotypingData genotyping_data, Lookup lookup, Parameters params)
    : genotyping_data(genotyping_data), lookup(lookup), params(params)
{
    chain = Chain(genotyping_data, lookup, params);
};
