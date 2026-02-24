#include "mcmc.h"

#include "include/spline/src/spline.h"
#include "mcmc_utils.h"
#include "multivector.h"

#include <Rcpp.h>
#include "profiler.h"

MCMC::MCMC(GenotypingData genotyping_data, Parameters params)
    : genotyping_data(genotyping_data), params(params)
{
    p_store.resize(params.num_populations);
    for (size_t i = 0; i < params.num_populations; ++i) {
        p_store[i].resize(genotyping_data.num_loci);
    }
    population_assignment_store.resize(genotyping_data.num_samples);
    
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
    
    // Reserve space for chains first, but don't construct default objects
    chains.clear();
    chains.reserve(params.pt_chains.size());
    std::vector<bool> any_ill_conditioned(params.pt_chains.size(), false);
    temp_gradient = params.pt_chains;
    
    // Construct chains sequentially first to avoid thread-safety issues during construction
    // (Construction involves shared resources like the static mutex)
    for (size_t i = 0; i < params.pt_chains.size(); ++i) {
        float temp = params.pt_chains[i];
        chains.emplace_back(genotyping_data, params, temp);
    }
    
    // Then initialize/retry in parallel if needed
    moire_parallel::parallel_for_always(0, chains.size(), [&](size_t i)
    {
        if (std::any_of(any_ill_conditioned.begin(), any_ill_conditioned.end(), [](bool v) { return v; }))
        {
            return;
        }

        bool ill_conditioned = !std::isfinite(chains[i].get_llik());
        int max_tries = params.max_initialization_tries;

        while (ill_conditioned and max_tries > 0)
        {
            Rcpp::checkUserInterrupt();
            chains[i].initialize_parameters();
            max_tries--;
            ill_conditioned = !std::isfinite(chains[i].get_llik());
        }

        any_ill_conditioned[i] = ill_conditioned;
    });

    if (std::any_of(any_ill_conditioned.begin(), any_ill_conditioned.end(), [](bool v) { return v; }))
    {
        Rcpp::stop("Initialization failed due to numerical instability, please check your input data or try increasing the maximum number of initialization tries. If you have a large number of alleles in some loci, consider removing extremely rare alleles.");
    }

    std::iota(swap_indices.begin(), swap_indices.end(), 0);
};

void MCMC::burnin(int step)
{
    ProfileScope scope("MCMC::burnin");
    // Strategic parallelism to avoid C stack overflow:
    // - If multiple chains: parallelize across chains, inner operations run sequentially (no nested parallelism)
    // - If single chain: sequential at top level, inner operations parallelize
    // NOTE: Nested parallelism causes C stack overflow in R, so we disable inner parallelism
    //       when multiple chains are active
    if (chains.size() > 1) {
        // Parallel tempering: parallelize across chains
        // Disable inner parallelism to avoid nested parallelism and C stack overflow
        // Set thread-local flag inside lambda so each worker thread has its own copy
        moire_parallel::parallel_for_always(0, chains.size(), [&](size_t i)
        {
            // Disable nested parallelism for this thread to prevent C stack overflow
            moire_parallel::disable_nested_parallelism = true;
            auto &chain = chains[i];
            ProfileScope s_update("Chain::updates (burnin)");
            chain.update_eps_neg(step);
            chain.update_eps_pos(step);
            chain.update_p(step);
            chain.update_m(step);
            if (params.allow_relatedness)
            {
                chain.update_r(step);
                chain.update_eff_coi(step);
            }
            chain.update_samples(step);
            chain.update_population_coi_p(step);
            chain.update_population_coi_r(step);
            chain.update_population_responsibility_vector(step);
            // Reset flag for this thread (though not strictly necessary since lambda ends)
            moire_parallel::disable_nested_parallelism = false;
        });
    } else {
        // Single chain: sequential at top level, inner operations parallelize
        auto &chain = chains[0];
        ProfileScope s_update("Chain::updates (burnin)");
        chain.update_eps_neg(step);
        chain.update_eps_pos(step);
        chain.update_p(step);
        chain.update_m(step);
        if (params.allow_relatedness)
        {
            chain.update_r(step);
            chain.update_eff_coi(step);
        }
        chain.update_samples(step);
        chain.update_population_coi_p(step);
        chain.update_population_coi_r(step);
        chain.update_population_responsibility_vector(step);
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
            std::isfinite(acceptanceRatio))
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
        if (!std::isfinite(new_temp_gradient[i]))
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
    ProfileScope scope("MCMC::sample");
    // Strategic parallelism to avoid C stack overflow:
    // - If multiple chains: parallelize across chains, inner operations run sequentially (no nested parallelism)
    // - If single chain: sequential at top level, inner operations parallelize
    // NOTE: Nested parallelism causes C stack overflow in R, so we disable inner parallelism
    //       when multiple chains are active
    if (chains.size() > 1) {
        // Parallel tempering: parallelize across chains
        // Disable inner parallelism to avoid nested parallelism and C stack overflow
        // Set thread-local flag inside lambda so each worker thread has its own copy
        moire_parallel::parallel_for_always(0, chains.size(), [&](size_t i)
        {
            // Disable nested parallelism for this thread to prevent C stack overflow
            moire_parallel::disable_nested_parallelism = true;
            auto &chain = chains[i];
            ProfileScope s_update("Chain::updates (sample)");
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
            chain.update_population_coi_p(params.burnin + step);
            chain.update_population_coi_r(params.burnin + step);
            chain.update_population_responsibility_vector(params.burnin + step);

        if ((params.thin == 0 or step % params.thin == 0) and
            (chain.is_hot() or chains.size() == 1))
        {
            ProfileScope s_store("MCMC::store_samples");
            for (size_t pop_idx = 0; pop_idx < params.num_populations; ++pop_idx) {
                for (size_t locus_idx = 0; locus_idx < genotyping_data.num_loci; ++locus_idx)
                {
                    const auto [begin, end] = chain.p.inner_iterators({pop_idx, locus_idx});
                    p_store[pop_idx][locus_idx].push_back(std::vector(begin, end));
                }
            }

            for (size_t sample_idx = 0; sample_idx < genotyping_data.num_samples; ++sample_idx)
            {
                m_store[sample_idx].push_back(chain.m.at({sample_idx}));
                eps_neg_store[sample_idx].push_back(chain.eps_neg.at({sample_idx}));
                eps_pos_store[sample_idx].push_back(chain.eps_pos.at({sample_idx}));
                r_store[sample_idx].push_back(chain.r.at({sample_idx}));
                // data_llik_store[sample_idx].push_back(chain.get_llik(sample_idx));

                if (params.record_latent_genotypes) {
                    for (size_t sample_locus_idx = 0; sample_locus_idx < genotyping_data.num_loci; ++sample_locus_idx)
                    {
                        const auto [begin, end] = chain.latent_genotypes_new.inner_iterators({sample_idx, sample_locus_idx});
                        latent_genotypes_store[sample_idx][sample_locus_idx].push_back(std::vector(begin, end));
                    }
                }

            }
            // MultiVector methods now automatically choose parallel vs sequential based on workload size
            const auto population_assignment_vec = (chain.transmission_llik_new.sum() + chain.coi_prior_new)
                    .element_add(chain.population_responsibility_vector.log().as_span())
                    .softmax(true);
            for (size_t sample_idx = 0; sample_idx < genotyping_data.num_samples; ++sample_idx) {
                const auto [begin, end] = population_assignment_vec.inner_iterators({sample_idx});
                population_assignment_store[sample_idx].push_back(std::vector(begin, end));
            }
            population_coi_p_store.push_back(chain.population_coi_p);
            population_coi_r_store.push_back(chain.population_coi_r);
            population_responsibility_store.push_back(chain.population_responsibility_vector.data());
            // Reset flag for this thread (though not strictly necessary since lambda ends)
            moire_parallel::disable_nested_parallelism = false;
        }
        });
    } else {
        // Single chain: sequential at top level, inner operations can parallelize
        auto &chain = chains[0];
        ProfileScope s_update("Chain::updates (sample)");
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
        chain.update_population_coi_p(params.burnin + step);
        chain.update_population_coi_r(params.burnin + step);
        chain.update_population_responsibility_vector(params.burnin + step);

        if ((params.thin == 0 or step % params.thin == 0) and chain.is_hot())
        {
            ProfileScope s_store("MCMC::store_samples");
            for (size_t pop_idx = 0; pop_idx < params.num_populations; ++pop_idx) {
                for (size_t locus_idx = 0; locus_idx < genotyping_data.num_loci; ++locus_idx)
                {
                    const auto [begin, end] = chain.p.inner_iterators({pop_idx, locus_idx});
                    p_store[pop_idx][locus_idx].push_back(std::vector(begin, end));
                }
            }

            for (size_t sample_idx = 0; sample_idx < genotyping_data.num_samples; ++sample_idx)
            {
                m_store[sample_idx].push_back(chain.m.at({sample_idx}));
                eps_neg_store[sample_idx].push_back(chain.eps_neg.at({sample_idx}));
                eps_pos_store[sample_idx].push_back(chain.eps_pos.at({sample_idx}));
                r_store[sample_idx].push_back(chain.r.at({sample_idx}));

                if (params.record_latent_genotypes) {
                    for (size_t sample_locus_idx = 0; sample_locus_idx < genotyping_data.num_loci; ++sample_locus_idx)
                    {
                        const auto [begin, end] = chain.latent_genotypes_new.inner_iterators({sample_idx, sample_locus_idx});
                        latent_genotypes_store[sample_idx][sample_locus_idx].push_back(std::vector(begin, end));
                    }
                }
            }
            const auto population_assignment_vec = (chain.transmission_llik_new.sum() + chain.coi_prior_new)
                    .element_add(chain.population_responsibility_vector.log().as_span())
                    .softmax(true);
            for (size_t sample_idx = 0; sample_idx < genotyping_data.num_samples; ++sample_idx) {
                const auto [begin, end] = population_assignment_vec.inner_iterators({sample_idx});
                population_assignment_store[sample_idx].push_back(std::vector(begin, end));
            }
            population_coi_p_store.push_back(chain.population_coi_p);
            population_coi_r_store.push_back(chain.population_coi_r);
            population_responsibility_store.push_back(chain.population_responsibility_vector.data());
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
SEXP tbb_enabled()
{
#ifdef MOIRE_HAVE_PARALLEL
    return Rcpp::wrap(true);
#else
    return Rcpp::wrap(false);
#endif
}
