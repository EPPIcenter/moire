#include "mcmc.h"

#include "chain.h"
#include "mcmc_utils.h"

#include <Rcpp.h>

#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_num_threads() 1
#define omp_get_thread_num() 0
#define omp_get_max_threads() 1
#define omp_get_thread_limit() 1
#define omp_get_num_procs() 1
#define omp_set_nested(a)
#define omp_set_num_threads(a)
#define omp_get_wtime() 0
#endif

MCMC::MCMC(GenotypingData genotyping_data, Parameters params)
    : genotyping_data(genotyping_data), params(params)
{
    p_store.resize(genotyping_data.num_loci);
    m_store.resize(genotyping_data.num_samples);
    r_store.resize(genotyping_data.num_samples);
    eps_neg_store.resize(genotyping_data.num_samples);
    eps_pos_store.resize(genotyping_data.num_samples);
    omp_set_num_threads(params.pt_num_threads);

    for (size_t i = 0; i < params.pt_chains; ++i)
    {
        chains.emplace_back(genotyping_data, params,
                            1.0 / std::pow(params.pt_grad, i));

        bool ill_conditioned = std::isnan(chains.back().get_llik());
        int max_tries = 1000;

        while (ill_conditioned and max_tries != 0)
        {
            chains.pop_back();
            chains.emplace_back(genotyping_data, params);
            max_tries--;
            ill_conditioned = std::isnan(chains.back().get_llik());
        }

        if (ill_conditioned)
        {
            Rcpp::stop("Error: Initial Llik is NaN");
        }

        swap_indices.push_back(i);
    }
};

void MCMC::burnin(int step)
{
#pragma omp parallel for
    for (auto &chain : chains)
    {
        chain.update_eps_neg(step);
        chain.update_eps_pos(step);
        chain.update_p(step);
        chain.update_m(step);
        if (params.allow_relatedness)
        {
            chain.update_r(step);
            chain.update_m_r(step);
        }
        chain.update_samples(step);
    }
    llik_burnin.push_back(get_llik());
    swap_chains();
}

void MCMC::swap_chains()
{
    for (size_t i = 0 + num_swaps % 2; i < chains.size() - 1; i += 2)
    {
        auto &chain_a = chains[swap_indices[i]];
        auto &chain_b = chains[swap_indices[i + 1]];
        double curr_llik_a = chain_a.get_llik();
        double temp_a = chain_a.get_temp();
        double curr_llik_b = chain_b.get_llik();
        double temp_b = chain_b.get_temp();

        double prop_llik_a = curr_llik_a / temp_a * temp_b;
        double prop_llik_b = curr_llik_b / temp_b * temp_a;

        long double acceptanceRatio =
            prop_llik_a + prop_llik_b - curr_llik_a - curr_llik_b;
        if (acceptanceRatio > 0 || log(R::runif(0, 1)) < acceptanceRatio)
        {
            std::swap(swap_indices[i], swap_indices[i + 1]);
            chain_a.set_temp(temp_b);
            chain_a.set_llik(prop_llik_a);
            chain_b.set_temp(temp_a);
            chain_b.set_llik(prop_llik_b);
        }
    }
    swap_store.push_back(swap_indices[0]);
    num_swaps++;
}

void MCMC::sample(int step)
{
#pragma omp parallel for
    for (auto &chain : chains)
    {
        chain.update_eps_neg(params.burnin + step);
        chain.update_eps_pos(params.burnin + step);
        chain.update_p(params.burnin + step);
        chain.update_m(params.burnin + step);

        if (params.allow_relatedness)
        {
            chain.update_r(params.burnin + step);
            chain.update_m_r(params.burnin + step);
        }
        chain.update_samples(params.burnin + step);

        if ((params.thin == 0 or step % params.thin == 0) and
            chain.get_temp() == 1.0)
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
                r_store[jj].push_back(chain.r[jj]);
            }
        }
    }
    llik_sample.push_back(get_llik());
    swap_chains();
}

double MCMC::get_llik() { return chains[swap_indices[0]].get_llik(); }
