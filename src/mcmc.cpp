#include "mcmc.h"

#include "chain.h"
#include "mcmc_utils.h"

#include <Rcpp.h>

void MCMC::burnin(int step)
{
    chain.update_eps_neg(step);
    chain.update_eps_pos(step);
    chain.update_p(step);
    chain.update_m(step);
    chain.update_individual_parameters(step);
    llik_burnin.push_back(chain.get_llik());
}

void MCMC::sample(int step)
{
    chain.update_eps_neg(params.burnin + step);
    chain.update_eps_pos(params.burnin + step);
    chain.update_p(params.burnin + step);
    chain.update_m(params.burnin + step);
    chain.update_individual_parameters(params.burnin + step);

    if (params.thin == 0 || step % params.thin == 0)
    {
        m_store.push_back(chain.m);
        p_store.push_back(chain.p);
        eps_neg_store.push_back(chain.eps_neg);
        eps_pos_store.push_back(chain.eps_pos);
        mean_coi_store.push_back(chain.mean_coi);
        llik_sample.push_back(chain.get_llik());
    }
}

MCMC::MCMC(GenotypingData genotyping_data, Lookup lookup, Parameters params)
    : genotyping_data(genotyping_data), lookup(lookup), params(params)
{
    chain = Chain(genotyping_data, lookup, params);
};
