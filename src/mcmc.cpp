#include "mcmc.h"

#include "chain.h"
#include "mcmc_utils.h"

#include <Rcpp.h>

MCMC::MCMC(GenotypingData genotyping_data, Parameters params)
    : genotyping_data(genotyping_data),
      params(params),
      chain(genotyping_data, params)
{
    p_store.resize(genotyping_data.num_loci);
    m_store.resize(genotyping_data.num_samples);
    eps_neg_store.resize(genotyping_data.num_samples);
    eps_pos_store.resize(genotyping_data.num_samples);
};

void MCMC::burnin(int step)
{
    chain.update_eps_neg(step);
    chain.update_eps_pos(step);
    chain.update_p(step);
    chain.update_m(step);
    llik_burnin.push_back(chain.get_llik());
}

void MCMC::sample(int step)
{
    chain.update_eps_neg(params.burnin + step);
    chain.update_eps_pos(params.burnin + step);
    chain.update_p(params.burnin + step);
    chain.update_m(params.burnin + step);

    if (params.thin == 0 || step % params.thin == 0)
    {
        for (size_t ii = 0; ii < genotyping_data.num_loci; ++ii)
        {
            p_store[ii].push_back(chain.p[ii]);
        }

        for (size_t jj = 0; jj < genotyping_data.num_samples; ++jj)
        {
            m_store[jj].push_back(chain.m[jj]);
            eps_neg_store[jj].push_back(chain.eps_neg[jj]);
            eps_pos_store[jj].push_back(chain.eps_pos[jj]);
        }
        llik_sample.push_back(chain.get_llik());
    }
}

double MCMC::get_llik() { return chain.get_llik(); }
