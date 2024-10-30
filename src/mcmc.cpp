#include "mcmc.h"

#include "include/spline/src/spline.h"
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
    latent_genotypes_store.resize(genotyping_data.num_samples);
    for (size_t i = 0; i < genotyping_data.num_samples; ++i)
    {
        latent_genotypes_store[i].resize(genotyping_data.num_loci);
    }
    m_store.resize(genotyping_data.num_samples);
    r_store.resize(genotyping_data.num_samples);
    eps_neg_store.resize(genotyping_data.num_samples);
    eps_pos_store.resize(genotyping_data.num_samples);
    data_llik_store.resize(genotyping_data.num_samples);
    swap_acceptances.resize(params.pt_chains.size() - 1, 0);
    swap_barriers.resize(params.pt_chains.size() - 1, 0.0);
    swap_indices.resize(params.pt_chains.size(), 0);

    omp_set_num_threads(params.pt_num_threads);

    
    chains.resize(params.pt_chains.size());
    std::vector<bool> any_ill_conditioned(params.pt_chains.size(), false);
    temp_gradient = params.pt_chains;
    #pragma omp parallel for
    for (size_t i = 0; i < params.pt_chains.size(); i++)
    {
        if (std::any_of(any_ill_conditioned.begin(), any_ill_conditioned.end(), [](bool v) { return v; }))
        {
            continue;
        }
        float temp = params.pt_chains[i];
        chains[i] = Chain(genotyping_data, params, temp);
        
        bool ill_conditioned = std::isnan(chains[i].get_llik());
        int max_tries = params.max_initialization_tries;

        while (ill_conditioned and max_tries > 0)
        {
            chains[i].initialize_parameters();
            max_tries--;
            ill_conditioned = std::isnan(chains[i].get_llik());
        }

        any_ill_conditioned[i] = ill_conditioned;
    }

    if (std::any_of(any_ill_conditioned.begin(), any_ill_conditioned.end(), [](bool v) { return v; }))
    {
        Rcpp::stop("Initialization failed due to numerical instability, please check your input data or try increasing the maximum number of initialization tries. If you have a large number of alleles in some loci, consider removing extremely rare alleles.");
    }

    std::iota(swap_indices.begin(), swap_indices.end(), 0);
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
            // chain.update_m_r(step);
            chain.update_eff_coi(step);
        }
        chain.update_samples(step);
        chain.update_mean_coi(step);
    }
    llik_burnin.push_back(get_llik());
    prior_burnin.push_back(get_prior());
    posterior_burnin.push_back(get_posterior());

    if (chains.size() > 1)
    {
        swap_chains(step, true);
        if (num_swaps % params.temp_adapt_steps == 0 and
            step > params.pre_adapt_steps and params.adapt_temp)
        {
            adapt_temp();
        }
    }
}

void MCMC::swap_chains(int step, bool burnin)
{
    for (size_t i = even_swap; i < chains.size() - 1; i += 2)
    {
        auto &chain_a = chains[swap_indices[i]];
        auto &chain_b = chains[swap_indices[i + 1]];

        float V_a = -chain_a.get_llik();
        float temp_a = chain_a.get_temp();

        float V_b = -chain_b.get_llik();
        float temp_b = chain_b.get_temp();

        float acceptanceRatio = (temp_b - temp_a) * (V_b - V_a);

        float acceptanceRate = std::min(1.0f, (float)std::exp(acceptanceRatio));

        if (burnin and step > params.pre_adapt_steps)
        {
            swap_barriers[i] += 1.0 - acceptanceRate;
        }

        float u = log(R::runif(0, 1));

        if ((acceptanceRatio > 0 || u < acceptanceRatio) and
            !std::isnan(acceptanceRatio))
        {
            std::swap(swap_indices[i], swap_indices[i + 1]);
            chain_a.set_temp(temp_b);
            chain_b.set_temp(temp_a);

            if (i == 0)
            {
                chain_b.set_hot(true);
                chain_a.set_hot(false);
            }

            if (!burnin)
            {
                swap_acceptances[i]++;
            }
        }
    }
    swap_store.push_back(swap_indices[0]);

    if (burnin and step > params.pre_adapt_steps)
    {
        num_swaps++;
    }
    even_swap = !even_swap;
}

int MCMC::get_hot_chain() { return swap_indices[0]; }

void MCMC::adapt_temp()
{
    // swap rate starts at t = 1 so we need to reverse the swap rate
    const std::vector<float> reversed_swap_barriers{swap_barriers.rbegin(),
                                                    swap_barriers.rend()};

    // cumulative swap rate starts from t = 0
    std::vector<float> cumulative_swap_rate(params.pt_chains.size(), 0.0);
    for (size_t i = 1; i < cumulative_swap_rate.size(); i++)
    {
        cumulative_swap_rate[i] = cumulative_swap_rate[i - 1] + 
            reversed_swap_barriers[i - 1] / (static_cast<float>(num_swaps) / 2.0f);
    }

    // gradient starts at t = 1 so we need to reverse the gradient
    const std::vector<float> reversed_gradient(temp_gradient.rbegin(), temp_gradient.rend());

    const tk::spline s(reversed_gradient, cumulative_swap_rate, tk::spline::cspline, true);

    // target swap rates
    std::vector<float> cumulative_swap_grid = std::vector<float>(cumulative_swap_rate.size(), 0.0f);
    
    cumulative_swap_grid.back() = cumulative_swap_rate.back();
    float step = cumulative_swap_rate.back() / (cumulative_swap_grid.size() - 1);

    for (size_t i = 1; i < cumulative_swap_grid.size() - 1; i++)
    {
        cumulative_swap_grid[i] = cumulative_swap_grid[0] + i * step;
    }

    std::vector<float> new_temp_gradient(temp_gradient.size());
    new_temp_gradient[0] = temp_gradient.back();
    new_temp_gradient[temp_gradient.size() - 1] = 1.0;
    for (size_t i = 1; i < temp_gradient.size() - 1; i++)
    {
        new_temp_gradient[i] = s.solve_one(cumulative_swap_grid[i]);
    }

    std::reverse(new_temp_gradient.begin(), new_temp_gradient.end());

    // check if any new temperatures are NAN, bail on the update if so and keep updating the swap barriers
    for (size_t i = 0; i < new_temp_gradient.size(); i++)
    {
        if (std::isnan(new_temp_gradient[i]))
        {
            return;
        }
    }

    // update the temperatures using the new gradient
    chains[swap_indices[0]].set_hot(true);
    for (size_t i = 1; i < chains.size() - 1; i++)
    {
        chains[swap_indices[i]].set_temp(new_temp_gradient[i]);
        chains[swap_indices[i]].set_hot(false);
    }

    temp_gradient = new_temp_gradient;
    num_swaps = 0;
    std::fill(swap_barriers.begin(), swap_barriers.end(), 0);
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
            chain.update_eff_coi(params.burnin + step);
        }
        chain.update_samples(params.burnin + step);
        chain.update_mean_coi(params.burnin + step);

        if ((params.thin == 0 or step % params.thin == 0) and
            (chain.get_hot() or chains.size() == 1))
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
                data_llik_store[jj].push_back(chain.get_llik(jj));

                if (params.record_latent_genotypes) {
                    for (size_t kk = 0; kk < genotyping_data.num_loci; ++kk)
                    {
                        latent_genotypes_store[jj][kk].push_back(chain.latent_genotypes_new[kk][jj]);
                    }
                }
            }
            mean_coi_store.push_back(chain.mean_coi);
        }
    }
    llik_sample.push_back(get_llik());
    prior_sample.push_back(get_prior());
    posterior_sample.push_back(get_posterior());

    swap_chains(step, false);
}

void MCMC::finalize()
{
    for (size_t i = 0; i < swap_barriers.size(); ++i)
    {
        swap_barriers[i] /= (num_swaps / 2);
    }
}

float MCMC::get_llik() { return chains[swap_indices[0]].get_llik(); }
float MCMC::get_prior() { return chains[swap_indices[0]].get_prior(); }
float MCMC::get_posterior() { return chains[swap_indices[0]].get_posterior(); }

// [[Rcpp::export]]
SEXP openmp_enabled()
{
    #ifdef _OPENMP
        return Rcpp::wrap(true);
    #else
        return Rcpp::wrap(false);
    #endif
}
