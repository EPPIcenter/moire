#include "chain.h"

#include "prob_any_missing.h"
#include "pam_fast_paths.h"
#include "multivector_fused.h"
#include "mcmc_utils.h"
#include "sampler.h"
#include "profiler.h"

#include <cmath>
#include <algorithm>

#include "env_defs.h"

#include <cstdlib>
#include <limits>
#include <map>
#include <numeric>
#include <span>

constexpr float min_sampled = std::numeric_limits<float>::min();

void Chain::update_m(int iteration)
{
    ProfileScope scope("Chain::update_m");
    auto indices = std::vector<std::size_t>(genotyping_data.num_samples);
    std::iota(indices.begin(), indices.end(), 0);
    sampler.shuffle_vec(indices);

    for (const auto sample_idx : indices)
    {
        const std::size_t prop_m = m.at({sample_idx}) + sampler.sample_coi_delta(2);
        if (prop_m > 0 and prop_m <= params.max_coi)
        {
            const int prev_m = m.at({sample_idx});
            const float prev_r = r.at({sample_idx});
            m.at({sample_idx}) = prop_m;
            calculate_coi_likelihood(sample_idx);

            float adj_ratio = 0;

            // First, sample latent genotypes for pop_idx == 0 (sequential, uses sampler)
            for (std::size_t locus_idx = 0; locus_idx < genotyping_data.num_loci; ++locus_idx)
            {
                const auto& lg = observation_model_->sample_latent_genotype(
                    sampler,
                    genotyping_data.get_observed_alleles(sample_idx, locus_idx), 
                    m.at({sample_idx}), 
                    eps_pos.at({sample_idx}), 
                    eps_neg.at({sample_idx})
                );
                assign_latent_genotype_new(sample_idx, locus_idx, lg.value);
                lg_adj_new.at({sample_idx, locus_idx}) = lg.log_prob;
                adj_ratio = adj_ratio + lg_adj_new.at({sample_idx, locus_idx}) - lg_adj_old.at({sample_idx, locus_idx});
            }
            moire_parallel::recalc_parallel_for(0, genotyping_data.num_loci, [&](std::size_t locus_idx) {
                calculate_observation_likelihood(sample_idx, locus_idx);
            });
            sync_obs_sum_for_sample(sample_idx);

            float new_llik;
            if (pam_fast_paths::tx_opts_enabled()) {
                recalculate_transmission_for_sample_incremental(
                    sample_idx, prev_m, prev_r);
                new_llik = obs_llik_sum_new_ + tx_llik_sum_new;
            } else {
                invalidate_transmission_llik_cache();
                moire_parallel::recalc_parallel_for_2d(0, params.num_populations, 0, genotyping_data.num_loci,
                    [&](std::size_t pop_idx, std::size_t locus_idx) {
                        calculate_transmission_likelihood(pop_idx, sample_idx, locus_idx);
                    });
                new_llik = calc_new_likelihood();
            }
            const float new_prior = calc_new_prior();
            const float new_post = new_llik * temp + new_prior;
            const float alpha = sampler.sample_log_mh_acceptance();
            const float mh_ratio = new_post - get_posterior() + adj_ratio;
            if (std::isfinite(new_post) && alpha <= mh_ratio)
            // Accept
            {
                llik = new_llik;
                prior = new_prior;

                save_coi_likelihood(sample_idx);
                // Save observation likelihoods for pop_idx == 0 (parallel over loci)
                moire_parallel::parallel_for(0, genotyping_data.num_loci, [&](std::size_t locus_idx) {
                    save_observation_likelihood(sample_idx, locus_idx);
                    lg_adj_old.at({sample_idx, locus_idx}) = lg_adj_new.at({sample_idx, locus_idx});
                    auto [begin, end] = latent_genotypes_new.inner_iterators({sample_idx, locus_idx});
                    std::copy(begin, end, latent_genotypes_old.inner_begin({sample_idx, locus_idx}));
                });
                // Save transmission likelihoods across all populations/loci
                for (std::size_t pop_idx = 0; pop_idx < params.num_populations; ++pop_idx) {
                    for (std::size_t locus_idx = 0; locus_idx < genotyping_data.num_loci; ++locus_idx) {
                        save_transmission_likelihood(pop_idx, sample_idx, locus_idx);
                    }
                }
                m_accept.at({sample_idx}) = m_accept.at({sample_idx}) + 1;
            }
            else
            // Reject
            {
                m.at({sample_idx}) = prev_m;
                restore_coi_likelihood(sample_idx);
                // Restore observation likelihoods for pop_idx == 0 (parallel over loci)
                moire_parallel::parallel_for(0, genotyping_data.num_loci, [&](std::size_t locus_idx) {
                    restore_observation_likelihood(sample_idx, locus_idx);
                    lg_adj_new.at({sample_idx, locus_idx}) = lg_adj_old.at({sample_idx, locus_idx});
                    restore_latent_genotype_new(sample_idx, locus_idx);
                });
                sync_obs_sum_for_sample(sample_idx);
                if (pam_fast_paths::tx_opts_enabled()) {
                    restore_transmission_for_sample_incremental(sample_idx);
                } else {
                    for (std::size_t pop_idx = 0; pop_idx < params.num_populations; ++pop_idx) {
                        for (std::size_t locus_idx = 0; locus_idx < genotyping_data.num_loci; ++locus_idx) {
                            restore_transmission_likelihood(sample_idx, pop_idx, locus_idx);
                        }
                    }
                    invalidate_transmission_llik_cache();
                }
            }
        }
    }
}

void Chain::update_eff_coi(int iteration)
{
    ProfileScope scope("Chain::update_eff_coi");
    auto sample_indices = std::vector<std::size_t>(genotyping_data.num_samples);
    std::iota(sample_indices.begin(), sample_indices.end(), 0);
    sampler.shuffle_vec(sample_indices);

    for (const auto sample_idx : sample_indices)
    {
        const float curr_eff_coi = (m.at({sample_idx}) - 1) * (1.0f - r.at({sample_idx})) + 1.0f;
        const auto prop_adj = sampler.sample_constrained(
            curr_eff_coi, m_r_var.at({sample_idx}), 1, params.max_coi);

        const float prop_eff_coi = std::get<0>(prop_adj);
        const float adj = std::get<1>(prop_adj);

        const std::size_t prop_m = m.at({sample_idx}) + sampler.sample_coi_delta(2);
        const float prop_r = 1.0f - (prop_eff_coi - 1.0f) / (prop_m - 1.0f);

        if (prop_m <= 0 || prop_m > params.max_coi || prop_r > .99999 || prop_r < .00001 || !std::isfinite(prop_r)) {
            continue;
        }

        const int prev_m = m.at({sample_idx});
        const float prev_r = r.at({sample_idx});
        m.at({sample_idx}) = prop_m;
        r.at({sample_idx}) = prop_r;


        float adj_ratio = adj;
        calculate_relatedness_likelihood(sample_idx);
        calculate_coi_likelihood(sample_idx);

        for (std::size_t locus_idx = 0; locus_idx < genotyping_data.num_loci; ++locus_idx)
        {
            const auto &lg = observation_model_->sample_latent_genotype(
                sampler,
                genotyping_data.get_observed_alleles(sample_idx, locus_idx),
                m.at({sample_idx}),
                eps_pos.at({sample_idx}),
                eps_neg.at({sample_idx})
            );
            assign_latent_genotype_new(sample_idx, locus_idx, lg.value);
            lg_adj_new.at({sample_idx, locus_idx}) = lg.log_prob;
            adj_ratio = adj_ratio + lg_adj_new.at({sample_idx, locus_idx}) - lg_adj_old.at({sample_idx, locus_idx});
        }
        moire_parallel::parallel_for(0, genotyping_data.num_loci, [&](std::size_t locus_idx) {
            calculate_observation_likelihood(sample_idx, locus_idx);
        });
        sync_obs_sum_for_sample(sample_idx);

        float new_llik;
        if (pam_fast_paths::tx_opts_enabled()) {
            recalculate_transmission_for_sample_incremental(
                sample_idx, prev_m, prev_r);
            new_llik = obs_llik_sum_new_ + tx_llik_sum_new;
        } else {
            invalidate_transmission_llik_cache();
            moire_parallel::recalc_parallel_for_2d(0, params.num_populations, 0, genotyping_data.num_loci,
                [&](std::size_t pop_idx, std::size_t locus_idx) {
                    calculate_transmission_likelihood(pop_idx, sample_idx, locus_idx);
                });
            new_llik = calc_new_likelihood();
        }

        const float new_prior = calc_new_prior();
        const float new_post = new_llik * temp + new_prior;
        const double mh_ratio = new_post - get_posterior() + adj_ratio;
        const double alpha = sampler.sample_log_mh_acceptance();

        // Reject
        if (!std::isfinite(new_post) || !std::isfinite(mh_ratio) || alpha > mh_ratio)
        {
            m.at({sample_idx}) = prev_m;
            r.at({sample_idx}) = prev_r;
            restore_relatedness_likelihood(sample_idx);
            restore_coi_likelihood(sample_idx);
            for (std::size_t locus_idx = 0; locus_idx < genotyping_data.num_loci; ++locus_idx)
            {
                lg_adj_new.at({sample_idx, locus_idx}) = lg_adj_old.at({sample_idx, locus_idx});
                restore_latent_genotype_new(sample_idx, locus_idx);
                restore_observation_likelihood(sample_idx, locus_idx);
            }
            sync_obs_sum_for_sample(sample_idx);
            if (pam_fast_paths::tx_opts_enabled()) {
                restore_transmission_for_sample_incremental(sample_idx);
            } else {
                for (std::size_t pop_idx = 0; pop_idx < params.num_populations; ++pop_idx) {
                    for (std::size_t locus_idx = 0; locus_idx < genotyping_data.num_loci; ++locus_idx) {
                        restore_transmission_likelihood(sample_idx, pop_idx, locus_idx);
                    }
                }
                invalidate_transmission_llik_cache();
            }
        }
        // Accept
        else
        {
            llik = new_llik;
            prior = new_prior;
            save_relatedness_likelihood(sample_idx);
            save_coi_likelihood(sample_idx);
            // Save observation likelihoods for pop_idx == 0 (sequential)
            for (std::size_t locus_idx = 0; locus_idx < genotyping_data.num_loci; ++locus_idx)
            {
                lg_adj_old.at({sample_idx, locus_idx}) = lg_adj_new.at({sample_idx, locus_idx});
                const auto [begin, end] = latent_genotypes_new.inner_iterators({sample_idx, locus_idx});
                std::copy(begin, end, latent_genotypes_old.inner_begin({sample_idx, locus_idx}));
                save_observation_likelihood(sample_idx, locus_idx);
            }
            // Save transmission likelihoods across all populations/loci
            for (std::size_t pop_idx = 0; pop_idx < params.num_populations; ++pop_idx) {
                for (std::size_t locus_idx = 0; locus_idx < genotyping_data.num_loci; ++locus_idx) {
                    save_transmission_likelihood(pop_idx, sample_idx, locus_idx);
                }
            }
            ++m_r_accept.at({sample_idx});
        }

        if (iteration < params.burnin and
            iteration > 15)  // don't start adapting until there are
                             // at least a few samples
        {
            const float acceptanceRate = m_r_accept.at({sample_idx}) / float(iteration);
            const float update =
                (acceptanceRate - .23) / std::pow(iteration + 1, .5);
            m_r_var.at({sample_idx}) = std::max(m_r_var.at({sample_idx}) + update, .01f);
        }
    }
}

void Chain::update_r(int iteration)
{
    ProfileScope scope("Chain::update_r");
    // Observation likelihood is unchanged across proposals in this update; its
    // running sum is maintained incrementally elsewhere.
    const float obs_llik = obs_llik_sum_new_;
    if (!tx_llik_cache_valid_) {
        rebuild_transmission_llik_cache();
    }

    auto sample_indices = std::vector<size_t>(genotyping_data.num_samples);
    std::iota(sample_indices.begin(), sample_indices.end(), 0);
    sampler.shuffle_vec(sample_indices);

    for (const auto sample_idx : sample_indices)
    {
        const auto prop_adj =
            sampler.sample_constrained(r.at({sample_idx}), r_var.at({sample_idx}), .00001, .99999);
        const float prop_r = std::get<0>(prop_adj);
        const float adj = std::get<1>(prop_adj);

        const float prev_r = r.at({sample_idx});
        r.at({sample_idx}) = prop_r;
        calculate_relatedness_likelihood(sample_idx);

        if (!pam_fast_paths::tx_opts_enabled()) {
            invalidate_transmission_llik_cache();
            for (std::size_t locus_idx = 0; locus_idx < genotyping_data.num_loci; ++locus_idx) {
                for (std::size_t pop_idx = 0; pop_idx < params.num_populations; ++pop_idx) {
                    calculate_transmission_likelihood(pop_idx, sample_idx, locus_idx);
                }
            }
            const float new_llik = calc_new_likelihood();
            const float new_prior = calc_new_prior();
            const float new_post = new_llik * temp + new_prior;
            if (!std::isfinite(new_post) or sampler.sample_log_mh_acceptance() >
                                            (new_post - get_posterior() + adj)) {
                r.at({sample_idx}) = prev_r;
                restore_relatedness_likelihood(sample_idx);
                for (std::size_t locus_idx = 0; locus_idx < genotyping_data.num_loci; ++locus_idx) {
                    for (std::size_t pop_idx = 0; pop_idx < params.num_populations; ++pop_idx) {
                        restore_transmission_likelihood(sample_idx, pop_idx, locus_idx);
                    }
                }
            } else {
                llik = new_llik;
                prior = new_prior;
                save_relatedness_likelihood(sample_idx);
                for (std::size_t locus_idx = 0; locus_idx < genotyping_data.num_loci; ++locus_idx) {
                    for (std::size_t pop_idx = 0; pop_idx < params.num_populations; ++pop_idx) {
                        save_transmission_likelihood(pop_idx, sample_idx, locus_idx);
                    }
                }
                ++r_accept.at({sample_idx});
            }
        } else {
            ProfileScope scope("Chain::update_r::recalc_transmission");
            moire_parallel::parallel_for_2d(
                0, params.num_populations, 0, genotyping_data.num_loci,
                [&](std::size_t pop_idx, std::size_t locus_idx) {
                    calculate_transmission_likelihood_after_r_change(
                        pop_idx, sample_idx, locus_idx);
                });
            for (std::size_t pop_idx = 0; pop_idx < params.num_populations; ++pop_idx) {
                for (std::size_t locus_idx = 0; locus_idx < genotyping_data.num_loci;
                     ++locus_idx) {
                    const float old_cell =
                        transmission_llik_old.unchecked_at({sample_idx, pop_idx, locus_idx});
                    const float new_cell =
                        transmission_llik_new.unchecked_at({sample_idx, pop_idx, locus_idx});
                    apply_transmission_cell_change(sample_idx, pop_idx, old_cell, new_cell);
                }
            }

            float new_llik;
            float new_prior;
            float new_post;
            {
                ProfileScope scope("Chain::update_r::calc_post");
                new_llik = obs_llik + tx_llik_sum_new;
                new_prior = calc_new_prior();
                new_post = new_llik * temp + new_prior;
            }

            if (!std::isfinite(new_post) or sampler.sample_log_mh_acceptance() >
                                            (new_post - get_posterior() + adj))
            {
                ProfileScope scope("Chain::update_r::reject_restore");
                r.at({sample_idx}) = prev_r;
                restore_relatedness_likelihood(sample_idx);
                for (std::size_t pop_idx = 0; pop_idx < params.num_populations; ++pop_idx) {
                    for (std::size_t locus_idx = 0; locus_idx < genotyping_data.num_loci;
                         ++locus_idx) {
                        const float proposed =
                            transmission_llik_new.unchecked_at({sample_idx, pop_idx, locus_idx});
                        const float old_cell =
                            transmission_llik_old.unchecked_at({sample_idx, pop_idx, locus_idx});
                        apply_transmission_cell_change(
                            sample_idx, pop_idx, proposed, old_cell);
                        transmission_llik_new.unchecked_at({sample_idx, pop_idx, locus_idx}) =
                            old_cell;
                    }
                }
            }
            else
            {
                ProfileScope scope("Chain::update_r::accept_save");
                llik = new_llik;
                prior = new_prior;
                save_relatedness_likelihood(sample_idx);
                moire_parallel::parallel_for_2d(
                    0, params.num_populations, 0, genotyping_data.num_loci,
                    [&](std::size_t pop_idx, std::size_t locus_idx) {
                        save_transmission_likelihood(pop_idx, sample_idx, locus_idx);
                    });
                ++r_accept.at({sample_idx});
            }
        }

        if (iteration < params.burnin and
            iteration > 15)  // don't start adapting until there are
                             // at least a few samples
        {
            const float acceptanceRate = r_accept.at({sample_idx}) / float(iteration);
            const float update =
                (acceptanceRate - .23) / std::pow(iteration + 1, .5);
            r_var.at({sample_idx}) = std::max(r_var.at({sample_idx}) + update, .0001f);
        }
    }
}

void Chain::update_m_r(int iteration)
{
    ProfileScope scope("Chain::update_m_r");
    auto sample_indices = std::vector<std::size_t>(genotyping_data.num_samples);
    std::iota(sample_indices.begin(), sample_indices.end(), 0);
    sampler.shuffle_vec(sample_indices);

    for (const auto sample_idx : sample_indices)
    {
        const auto prop_adj =
            sampler.sample_constrained(r.at({sample_idx}), m_r_var.at({sample_idx}), .00001, .99999);
        const float prop_r = std::get<0>(prop_adj);
        const float adj = std::get<1>(prop_adj);

        const std::size_t prop_m = m.at({sample_idx}) + sampler.sample_coi_delta(1);

        if (prop_m <= 0 or prop_m > params.max_coi) {
            continue;
        }

        const float prev_m = m.at({sample_idx});
        const float prev_r = r.at({sample_idx});
        r.at({sample_idx}) = prop_r;
        calculate_relatedness_likelihood(sample_idx);
        m.at({sample_idx}) = prop_m;

        float adj_ratio = adj;

        calculate_coi_likelihood(sample_idx);
        // First, sample latent genotypes for pop_idx == 0 (sequential, uses sampler)
        for (std::size_t locus_idx = 0; locus_idx < genotyping_data.num_loci; ++locus_idx)
        {
            const auto &lg = observation_model_->sample_latent_genotype(
                sampler,
                genotyping_data.get_observed_alleles(sample_idx, locus_idx), 
                m.at({sample_idx}), 
                eps_pos.at({sample_idx}),
                eps_neg.at({sample_idx})
            );
            assign_latent_genotype_new(sample_idx, locus_idx, lg.value);
            lg_adj_new.at({sample_idx, locus_idx}) = lg.log_prob;
            adj_ratio = adj_ratio + lg_adj_new.at({sample_idx, locus_idx}) - lg_adj_old.at({sample_idx, locus_idx});
        }
        moire_parallel::recalc_parallel_for(0, genotyping_data.num_loci, [&](std::size_t locus_idx) {
            calculate_observation_likelihood(sample_idx, locus_idx);
        });
        sync_obs_sum_for_sample(sample_idx);

        float new_llik;
        if (pam_fast_paths::tx_opts_enabled()) {
            recalculate_transmission_for_sample_incremental(
                sample_idx, static_cast<int>(prev_m), prev_r);
            new_llik = obs_llik_sum_new_ + tx_llik_sum_new;
        } else {
            invalidate_transmission_llik_cache();
            moire_parallel::recalc_parallel_for_2d(0, params.num_populations, 0, genotyping_data.num_loci,
                [&](std::size_t pop_idx, std::size_t locus_idx) {
                    calculate_transmission_likelihood(pop_idx, sample_idx, locus_idx);
                });
            new_llik = calc_new_likelihood();
        }

        const float new_prior = calc_new_prior();
        const float new_post = new_llik * temp + new_prior;
        const float alpha = sampler.sample_log_mh_acceptance();
        const float mh_ratio = new_post - get_posterior() + adj_ratio;

        // Reject
        if (!std::isfinite(new_post) or alpha > mh_ratio)
        {
            r.at({sample_idx}) = prev_r;
            restore_relatedness_likelihood(sample_idx);
            m.at({sample_idx}) = prev_m;
            restore_coi_likelihood(sample_idx);
            moire_parallel::parallel_for(0, genotyping_data.num_loci, [&](std::size_t locus_idx) {
                restore_observation_likelihood(sample_idx, locus_idx);
                lg_adj_new.at({sample_idx, locus_idx}) = lg_adj_old.at({sample_idx, locus_idx});
                restore_latent_genotype_new(sample_idx, locus_idx);
            });
            sync_obs_sum_for_sample(sample_idx);
            if (pam_fast_paths::tx_opts_enabled()) {
                restore_transmission_for_sample_incremental(sample_idx);
            } else {
                for (std::size_t pop_idx = 0; pop_idx < params.num_populations; ++pop_idx) {
                    for (std::size_t locus_idx = 0; locus_idx < genotyping_data.num_loci; ++locus_idx) {
                        restore_transmission_likelihood(sample_idx, pop_idx, locus_idx);
                    }
                }
                invalidate_transmission_llik_cache();
            }
        }
        else
        {
            llik = new_llik;
            prior = new_prior;
            save_relatedness_likelihood(sample_idx);
            save_coi_likelihood(sample_idx);
            // Save observation likelihoods for pop_idx == 0 (parallel over loci)
            moire_parallel::parallel_for(0, genotyping_data.num_loci, [&](std::size_t locus_idx) {
                save_observation_likelihood(sample_idx, locus_idx);
                lg_adj_old.at({sample_idx, locus_idx}) = lg_adj_new.at({sample_idx, locus_idx});
                const auto [begin, end] = latent_genotypes_new.inner_iterators({sample_idx, locus_idx});
                std::copy(begin, end, latent_genotypes_old.inner_begin({sample_idx, locus_idx}));
            });
            // Save transmission likelihoods across all populations/loci
            for (std::size_t pop_idx = 0; pop_idx < params.num_populations; ++pop_idx) {
                for (std::size_t locus_idx = 0; locus_idx < genotyping_data.num_loci; ++locus_idx) {
                    save_transmission_likelihood(pop_idx, sample_idx, locus_idx);
                }
            }
            ++m_r_accept.at({sample_idx});
        }

        if (iteration < params.burnin and
            iteration > 15)  // don't start adapting until there are
                             // at least a few samples
        {
            const float acceptanceRate = m_r_accept.at({sample_idx}) / float(iteration);
            const float update =
                (acceptanceRate - .23) / std::pow(iteration + 1, .5);
            m_r_var.at({sample_idx}) = std::max(m_r_var.at({sample_idx}) + update, .0001f);
        }
    }
}

/*
 * SALT Sampler approach.
 * https://doi.org/10.1080/00949655.2017.1376063
 */
void Chain::update_p(int iteration)
{
    ProfileScope scope("Chain::update_p");
    // Observation and prior are unchanged across proposals within this update.
    const float obs_llik = obs_llik_sum_new_;
    if (!tx_llik_cache_valid_) {
        rebuild_transmission_llik_cache();
    }

    auto locus_indices = std::vector<std::size_t>(genotyping_data.num_loci);
    auto pop_indices = std::vector<std::size_t>(params.num_populations);
    std::iota(locus_indices.begin(), locus_indices.end(), 0);
    std::iota(pop_indices.begin(), pop_indices.end(), 0);
    sampler.shuffle_vec(locus_indices);
    sampler.shuffle_vec(pop_indices);

    for (const std::size_t pop_idx : pop_indices) {
        for (const std::size_t locus_idx : locus_indices)
        {
            const std::size_t num_alleles = p.ragged_dimensions(locus_idx);
            // Alternative: rep = num_alleles (fixed rep=3 used for tuning) 
            std::size_t rep = 3;
            const auto [begin, end] = p.inner_iterators({pop_idx, locus_idx});
            while (rep-- > 0)
            {
                const size_t allele_idx = sampler.sample_random_int(0, num_alleles - 1);

                ++p_attempt.at({pop_idx, locus_idx, allele_idx});
                {
                    ProfileScope scope("Chain::update_p::build_logit");
                    // build logit of current simplex
                }
                auto logitPropP = UtilFunctions::logitVec(begin, end);

                const float logitCurr = logitPropP[allele_idx];
                const float logitProp =
                    sampler.sample_epsilon(logitCurr, p_prop_var.at({pop_idx, locus_idx, allele_idx}));

                const auto currLogPQ = UtilFunctions::log_pq(logitCurr);
                const auto propLogPQ = UtilFunctions::log_pq(logitProp);

                logitPropP.erase(logitPropP.begin() + allele_idx);

                const float ls =
                    propLogPQ.second - UtilFunctions::logitSum(logitPropP);
                logitPropP = UtilFunctions::logitScale(logitPropP, ls);
                logitPropP.insert(logitPropP.begin() + allele_idx, logitProp);

                const float logAdj =
                    (currLogPQ.first - propLogPQ.first) +
                    (num_alleles - 1) * (currLogPQ.second - propLogPQ.second);

                {
                    ProfileScope scope("Chain::update_p::expit");
                }
                const auto prop_p = UtilFunctions::expitVec(logitPropP);
                // check to make sure the proposed simplex is within a bounded range
                bool sub_threshold_flag = false;
                for (const auto el : prop_p)
                {
                    // below a very small threshold that can cause numerical
                    // instability
                    if (el < 1e-4 || !std::isfinite(el))
                    {
                        sub_threshold_flag = true;
                        break;
                    }
                }

                if (sub_threshold_flag)
                {
                    break;
                }

                update_p_prev_p_ws_.assign(begin, end);
                p.inner_fill({pop_idx, locus_idx}, prop_p);

                {
                    const std::size_t n_samples = genotyping_data.num_samples;
                    if (pam_fast_paths::tx_opts_enabled()) {
                        recalculate_transmission_at_locus_after_p_change(
                            pop_idx, locus_idx,
                            std::span<const float>(update_p_prev_p_ws_.data(),
                                                   update_p_prev_p_ws_.size()),
                            allele_idx);
                    } else {
                        ProfileScope scope("Chain::update_p::recalc_transmission");
                        moire_parallel::recalc_parallel_for(0, n_samples, [&](std::size_t sample_idx) {
                            calculate_transmission_likelihood(pop_idx, sample_idx, locus_idx);
                        });
                    }
                    apply_transmission_column_change(pop_idx, locus_idx);
                }

                float new_llik;
                float new_post;
                {
                    ProfileScope scope("Chain::update_p::calc_post");
                    new_llik = obs_llik + tx_llik_sum_new;
                    new_post = new_llik * temp + prior;
                }

                const float acceptanceRatio = new_post - get_posterior() + logAdj;

                if (!std::isfinite(new_post) or
                    sampler.sample_log_mh_acceptance() > acceptanceRatio)
                {
                    p.inner_fill({pop_idx, locus_idx}, update_p_prev_p_ws_);
                    restore_transmission_column_change(pop_idx, locus_idx);
                }
                else
                {
                    ProfileScope scope("Chain::update_p::accept_save");
                    llik = new_llik;
                    moire_parallel::parallel_for(0, genotyping_data.num_samples, [&](std::size_t sample_idx) {
                        save_transmission_likelihood(pop_idx, sample_idx, locus_idx);
                    });
                    ++p_accept.at({pop_idx, locus_idx, allele_idx});
                }

                // don't start adapting until there are at least a few samples
                if (iteration < params.burnin and p_attempt.at({pop_idx, locus_idx, allele_idx}) > 15)
                {
                    const float acceptanceRate =
                        (p_accept.at({pop_idx, locus_idx, allele_idx}) + 1) / (float(p_attempt.at({pop_idx, locus_idx, allele_idx})) + 1);
                    p_prop_var.at({pop_idx, locus_idx, allele_idx}) += (acceptanceRate - .23) /
                                        std::pow(p_attempt.at({pop_idx, locus_idx, allele_idx}) + 1, .5);
                    p_prop_var.at({pop_idx, locus_idx, allele_idx}) = std::max(p_prop_var.at({pop_idx, locus_idx, allele_idx}), .01f);
                }
            }
        }
    }
}

void Chain::update_eps_pos(int iteration)
{
    ProfileScope scope("Chain::update_eps_pos");
    auto sample_indices = std::vector<std::size_t>(genotyping_data.num_samples);
    std::iota(sample_indices.begin(), sample_indices.end(), 0);
    sampler.shuffle_vec(sample_indices);

    for (const auto sample_idx : sample_indices)
    {
        const auto prop_adj = sampler.sample_constrained(
            eps_pos.at({sample_idx}), eps_pos_var.at({sample_idx}), min_sampled, 1);
        const float prop_eps_pos = std::get<0>(prop_adj);
        const float adj = std::get<1>(prop_adj);

        if (prop_eps_pos < 1 && prop_eps_pos > 1e-32)
        {
            const float prev_eps_pos = eps_pos.at({sample_idx});
            eps_pos.at({sample_idx}) = prop_eps_pos;
            calculate_eps_pos_likelihood(sample_idx);

            moire_parallel::parallel_for(0, genotyping_data.num_loci, [&](std::size_t locus_idx) {
                calculate_observation_likelihood(sample_idx, locus_idx);
            });
            sync_obs_sum_for_sample(sample_idx);

            const float new_llik = calc_new_likelihood();
            const float new_prior = calc_new_prior();
            const float new_post = new_llik * temp + new_prior;

            // Reject
            if (!std::isfinite(new_post) or sampler.sample_log_mh_acceptance() >
                                            (new_post - get_posterior() + adj))
            {
                eps_pos.at({sample_idx}) = prev_eps_pos;
                restore_eps_pos_likelihood(sample_idx);
                for (std::size_t locus_idx = 0; locus_idx < genotyping_data.num_loci; ++locus_idx)
                {
                    restore_observation_likelihood(sample_idx, locus_idx);
                }
                sync_obs_sum_for_sample(sample_idx);
            }
            else
            {
                llik = new_llik;
                prior = new_prior;
                save_eps_pos_likelihood(sample_idx);
                for (std::size_t locus_idx = 0; locus_idx < genotyping_data.num_loci; ++locus_idx)
                {
                    save_observation_likelihood(sample_idx, locus_idx);
                }
                ++eps_pos_accept.at({sample_idx});
            }

            if (iteration < params.burnin and
                iteration > 15)  // don't start adapting until there are
                                 // at least a few samples
            {
                const float acceptanceRate =
                    eps_pos_accept.at({sample_idx}) / float(iteration);
                const float update =
                    (acceptanceRate - .23) / std::pow(iteration + 1, .5);
                eps_pos_var.at({sample_idx}) = std::max(eps_pos_var.at({sample_idx}) + update, .0001f);
            }
        }
    }
}

void Chain::update_eps_neg(int iteration)
{
    ProfileScope scope("Chain::update_eps_neg");
    auto sample_indices = std::vector<std::size_t>(genotyping_data.num_samples);
    std::iota(sample_indices.begin(), sample_indices.end(), 0);
    sampler.shuffle_vec(sample_indices);

    for (const auto sample_idx : sample_indices)
    {
        const auto prop_adj = sampler.sample_constrained(
            eps_neg.at({sample_idx}), eps_neg_var.at({sample_idx}), min_sampled, 1);
        const float prop_eps_neg = std::get<0>(prop_adj);
        const float adj = std::get<1>(prop_adj);

        if (prop_eps_neg < 1 && prop_eps_neg > 1e-32)
        {
            const float prev_eps_neg = eps_neg.at({sample_idx});
            eps_neg.at({sample_idx}) = prop_eps_neg;
            calculate_eps_neg_likelihood(sample_idx);

            moire_parallel::parallel_for(0, genotyping_data.num_loci, [&](std::size_t locus_idx) {
                calculate_observation_likelihood(sample_idx, locus_idx);
            });
            sync_obs_sum_for_sample(sample_idx);

            const float new_llik = calc_new_likelihood();
            const float new_prior = calc_new_prior();
            const float new_post = new_llik * temp + new_prior;

            // Reject
            if (!std::isfinite(new_post) or sampler.sample_log_mh_acceptance() >
                                            (new_post - get_posterior() + adj))
            {
                eps_neg.at({sample_idx}) = prev_eps_neg;
                restore_eps_neg_likelihood(sample_idx);
                for (std::size_t locus_idx = 0; locus_idx < genotyping_data.num_loci; ++locus_idx)
                {
                    restore_observation_likelihood(sample_idx, locus_idx);
                }
                sync_obs_sum_for_sample(sample_idx);
            }
            else
            {
                llik = new_llik;
                prior = new_prior;
                save_eps_neg_likelihood(sample_idx);
                for (std::size_t locus_idx = 0; locus_idx < genotyping_data.num_loci; ++locus_idx)
                {
                    save_observation_likelihood(sample_idx, locus_idx);
                }
                ++eps_neg_accept.at({sample_idx});
            }

            if (iteration < params.burnin and
                iteration > 15)  // don't start adapting until there are
                                 // at least a few samples
            {
                const float acceptanceRate =
                    eps_neg_accept.at({sample_idx}) / float(iteration);
                const float update =
                    (acceptanceRate - .23) / std::pow(iteration + 1, .5);
                eps_neg_var.at({sample_idx}) = std::max(eps_neg_var.at({sample_idx}) + update, .0001f);
            }
        }
    }
}

void Chain::update_samples(int iteration)
{
    ProfileScope scope("Chain::update_samples");
    auto sample_indices = std::vector<std::size_t>(genotyping_data.num_samples);
    std::iota(sample_indices.begin(), sample_indices.end(), 0);
    sampler.shuffle_vec(sample_indices);

    for (const auto sample_idx : sample_indices)
    {
        const auto eps_neg_prop_adj = sampler.sample_constrained(
            eps_neg.at({sample_idx}), eps_neg_var.at({sample_idx}), min_sampled, 1);
        const float prop_eps_neg = std::get<0>(eps_neg_prop_adj);
        const float eps_neg_adj = std::get<1>(eps_neg_prop_adj);
        const bool valid_prop_eps_neg =
            prop_eps_neg < 1 && prop_eps_neg > 1e-32;

        const auto eps_pos_prop_adj = sampler.sample_constrained(
            eps_pos.at({sample_idx}), eps_pos_var.at({sample_idx}), min_sampled, 1);
        const float prop_eps_pos = std::get<0>(eps_pos_prop_adj);
        const float eps_pos_adj = std::get<1>(eps_pos_prop_adj);
        const bool valid_prop_eps_pos =
            prop_eps_pos < 1 && prop_eps_pos > 1e-32;

        float prop_r = 0;
        float r_adj = 0;
        bool valid_prop_r = true;
        if (params.allow_relatedness)
        {
            auto r_prop_adj =
                sampler.sample_constrained(r.at({sample_idx}), r_var.at({sample_idx}), min_sampled, .99);
            prop_r = std::get<0>(r_prop_adj);
            r_adj = std::get<1>(r_prop_adj);
            valid_prop_r = prop_r < 1 && prop_r > 1e-32;
        }

        const std::size_t prop_m = m.at({sample_idx}) + sampler.sample_coi_delta(2);
        const bool valid_prop_m = prop_m > 0 && prop_m <= params.max_coi;

        if (valid_prop_eps_neg && valid_prop_eps_pos && valid_prop_r &&
            valid_prop_m)
        {
            const float prev_eps_pos = eps_pos.at({sample_idx});
            eps_pos.at({sample_idx}) = prop_eps_pos;
            calculate_eps_pos_likelihood(sample_idx);

            const float prev_eps_neg = eps_neg.at({sample_idx});
            eps_neg.at({sample_idx}) = prop_eps_neg;
            calculate_eps_neg_likelihood(sample_idx);

            const float prev_r = r.at({sample_idx});
            r.at({sample_idx}) = prop_r;
            calculate_relatedness_likelihood(sample_idx);

            const float prev_m = m.at({sample_idx});
            m.at({sample_idx}) = prop_m;

            float adj_ratio = eps_neg_adj + eps_pos_adj + r_adj;

            calculate_coi_likelihood(sample_idx);
            // First, sample latent genotypes for pop_idx == 0 (sequential, uses sampler)
            for (std::size_t locus_idx = 0; locus_idx < genotyping_data.num_loci; ++locus_idx)
            {
                auto lg = observation_model_->sample_latent_genotype(
                    sampler,
                    genotyping_data.get_observed_alleles(sample_idx, locus_idx), m.at({sample_idx}),
                    eps_pos.at({sample_idx}), eps_neg.at({sample_idx}));
                assign_latent_genotype_new(sample_idx, locus_idx, lg.value);
                lg_adj_new.at({sample_idx, locus_idx}) = lg.log_prob;
                adj_ratio = adj_ratio + lg_adj_new.at({sample_idx, locus_idx}) - lg_adj_old.at({sample_idx, locus_idx});
            }
            moire_parallel::recalc_parallel_for(0, genotyping_data.num_loci, [&](std::size_t locus_idx) {
                calculate_observation_likelihood(sample_idx, locus_idx);
            });
            sync_obs_sum_for_sample(sample_idx);

            float new_llik;
            if (pam_fast_paths::tx_opts_enabled()) {
                recalculate_transmission_for_sample_incremental(
                    sample_idx, static_cast<int>(prev_m), prev_r);
                new_llik = obs_llik_sum_new_ + tx_llik_sum_new;
            } else {
                invalidate_transmission_llik_cache();
                moire_parallel::recalc_parallel_for_2d(0, params.num_populations, 0, genotyping_data.num_loci,
                    [&](std::size_t pop_idx, std::size_t locus_idx) {
                        calculate_transmission_likelihood(pop_idx, sample_idx, locus_idx);
                    });
                new_llik = calc_new_likelihood();
            }

            const float new_prior = calc_new_prior();
            const float new_post = new_llik * temp + new_prior;
            const float alpha = sampler.sample_log_mh_acceptance();

            if (std::isfinite(new_post) and
                alpha <= (new_post - get_posterior() + adj_ratio))
            {
                llik = new_llik;
                prior = new_prior;
                save_eps_neg_likelihood(sample_idx);
                save_eps_pos_likelihood(sample_idx);
                save_relatedness_likelihood(sample_idx);
                save_coi_likelihood(sample_idx);
                // Save observation likelihoods for pop_idx == 0 (parallel over loci)
                moire_parallel::parallel_for(0, genotyping_data.num_loci, [&](std::size_t locus_idx) {
                    save_observation_likelihood(sample_idx, locus_idx);
                    lg_adj_old.at({sample_idx, locus_idx}) = lg_adj_new.at({sample_idx, locus_idx});
                    auto [begin, end] = latent_genotypes_new.inner_iterators({sample_idx, locus_idx});
                    std::copy(begin, end, latent_genotypes_old.inner_begin({sample_idx, locus_idx}));
                });
                // Save transmission likelihoods across all populations/loci
                for (std::size_t pop_idx = 0; pop_idx < params.num_populations; ++pop_idx) {
                    for (std::size_t locus_idx = 0; locus_idx < genotyping_data.num_loci; ++locus_idx) {
                        save_transmission_likelihood(pop_idx, sample_idx, locus_idx);
                    }
                }
                ++sample_accept.at({sample_idx});
            }
            else
            {
                m.at({sample_idx}) = prev_m;
                eps_pos.at({sample_idx}) = prev_eps_pos;
                eps_neg.at({sample_idx}) = prev_eps_neg;
                r.at({sample_idx}) = prev_r;
                restore_eps_neg_likelihood(sample_idx);
                restore_eps_pos_likelihood(sample_idx);
                restore_relatedness_likelihood(sample_idx);
                restore_coi_likelihood(sample_idx);
                // Restore observation likelihoods for pop_idx == 0 (parallel over loci)
                moire_parallel::parallel_for(0, genotyping_data.num_loci, [&](std::size_t locus_idx) {
                    restore_observation_likelihood(sample_idx, locus_idx);
                    lg_adj_new.at({sample_idx, locus_idx}) = lg_adj_old.at({sample_idx, locus_idx});
                    restore_latent_genotype_new(sample_idx, locus_idx);
                });
                sync_obs_sum_for_sample(sample_idx);
                if (pam_fast_paths::tx_opts_enabled()) {
                    restore_transmission_for_sample_incremental(sample_idx);
                } else {
                    for (std::size_t pop_idx = 0; pop_idx < params.num_populations; ++pop_idx) {
                        for (std::size_t locus_idx = 0; locus_idx < genotyping_data.num_loci; ++locus_idx) {
                            restore_transmission_likelihood(sample_idx, pop_idx, locus_idx);
                        }
                    }
                    invalidate_transmission_llik_cache();
                }
            }
        }
    }
}

void Chain::update_population_coi_p(int iteration)
{
    const auto prop_adj = sampler.sample_constrained(population_coi_p, population_coi_p_sampling_variance,
                                                    0.01, 0.99);
    const float prop_p = std::get<0>(prop_adj);
    const float adj = std::get<1>(prop_adj);

    const float prev_p = population_coi_p;
    population_coi_p = prop_p;
    calculate_population_coi_p_likelihood();
    for (std::size_t sample_idx = 0; sample_idx < genotyping_data.num_samples; ++sample_idx)
    {
        calculate_coi_likelihood(sample_idx);
    }
    refresh_all_samples_tx_logsumexp();
    const float new_llik = obs_llik_sum_new_ + tx_llik_sum_new;
    const float new_prior = calc_new_prior();
    const float new_post = new_llik * temp + new_prior;

    const float alpha = sampler.sample_log_mh_acceptance();
    if (std::isfinite(new_post) and alpha <= (new_post - get_posterior() + adj))
    {
        llik = new_llik;
        prior = new_prior;
        save_population_coi_p_likelihood();
        for (std::size_t sample_idx = 0; sample_idx < genotyping_data.num_samples; ++sample_idx)
        {
            save_coi_likelihood(sample_idx);
        }
        ++population_coi_p_accept;
    }
    else
    {
        population_coi_p = prev_p;
        restore_population_coi_p_likelihood();
        for (std::size_t sample_idx = 0; sample_idx < genotyping_data.num_samples; ++sample_idx)
        {
            restore_coi_likelihood(sample_idx);
        }
        refresh_all_samples_tx_logsumexp();
    }

    if (iteration < params.burnin and iteration > 15)
    {
        const float acceptance_rate = population_coi_p_accept / (iteration + 1);
        const float update =
            (acceptance_rate - .23) / std::pow(iteration + 1, .5);
        population_coi_p_sampling_variance = std::max(population_coi_p_sampling_variance + update, .0001f);
    }
}

void Chain::update_population_coi_r(int iteration)
{
    const auto prop_adj = sampler.sample_constrained(population_coi_r, population_coi_r_sampling_variance,
                                                    .01, 100);
    const float prop_r = std::get<0>(prop_adj);
    const float adj = std::get<1>(prop_adj);

    const float prev_r = population_coi_r;
    population_coi_r = prop_r;
    calculate_population_coi_r_likelihood();
    for (std::size_t sample_idx = 0; sample_idx < genotyping_data.num_samples; ++sample_idx)
    {
        calculate_coi_likelihood(sample_idx);
    }
    refresh_all_samples_tx_logsumexp();
    const float new_llik = obs_llik_sum_new_ + tx_llik_sum_new;
    const float new_prior = calc_new_prior();
    const float new_post = new_llik * temp + new_prior;

    const float alpha = sampler.sample_log_mh_acceptance();

    if (std::isfinite(new_post) and alpha <= (new_post - get_posterior() + adj))
    {
        llik = new_llik;
        prior = new_prior;
        save_population_coi_r_likelihood();
        for (std::size_t sample_idx = 0; sample_idx < genotyping_data.num_samples; ++sample_idx)
        {
            save_coi_likelihood(sample_idx);
        }
        ++population_coi_r_accept;
    }
    else
    {
        population_coi_r = prev_r;
        restore_population_coi_r_likelihood();
        for (std::size_t sample_idx = 0; sample_idx < genotyping_data.num_samples; ++sample_idx)
        {
            restore_coi_likelihood(sample_idx);
        }
        refresh_all_samples_tx_logsumexp();
    }

    if (iteration < params.burnin and iteration > 15)
    {
        const float acceptance_rate = population_coi_r_accept / (iteration + 1);
        const float update = (acceptance_rate - .23) / std::pow(iteration + 1, .5);
        population_coi_r_sampling_variance = std::max(population_coi_r_sampling_variance + update, .0001f);
    }
}

void Chain::update_population_responsibility_vector(int iteration) 
{
    const std::size_t num_populations = params.num_populations;
    if (num_populations == 1) {
        return;
    }
    std::size_t rep = 2;
    const auto [begin, end] = population_responsibility_vector.inner_iterators();
    while (rep-- > 0) 
    {
        const size_t pop_idx = sampler.sample_random_int(0, num_populations - 1);
        ++population_responsibility_vector_attempt.at({pop_idx});

        auto logitPropP = UtilFunctions::logitVec(begin, end);

        const float logitCurr = logitPropP[pop_idx];
        const float logitProp = sampler.sample_epsilon(logitCurr, population_responsibility_vector_prop_var.at({pop_idx}));

        const auto currLogPQ = UtilFunctions::log_pq(logitCurr);
        const auto propLogPQ = UtilFunctions::log_pq(logitProp);

        logitPropP.erase(logitPropP.begin() + pop_idx);

        const float ls = propLogPQ.second - UtilFunctions::logitSum(logitPropP);
        logitPropP = UtilFunctions::logitScale(logitPropP, ls);
        logitPropP.insert(logitPropP.begin() + pop_idx, logitProp);

        const float logAdj = (currLogPQ.first - propLogPQ.first) + (num_populations - 1) * (currLogPQ.second - propLogPQ.second);
        auto prop_p = UtilFunctions::expitVec(logitPropP);

        // bool sub_threshold_flag = false;
        float sum = 0.0f;
        for (auto &el : prop_p)
        {
            el = std::max(el, 1e-6f);
            sum += el;
        }
        for (auto &el : prop_p) {
            el = el / sum;
        }

        // sort the proposed population responsibility vector 
        auto prop_p_begin = prop_p.begin();
        auto prop_p_end = prop_p.end();
        auto [var_begin, var_end] = population_responsibility_vector_prop_var.inner_iterators();
        // // Create a vector of indices
        std::vector<size_t> indices(params.num_populations);
        std::iota(indices.begin(), indices.end(), 0);
        
        // // Sort indices based on population responsibility values
        std::sort(indices.begin(), indices.end(), 
            [prop_p_begin](size_t i1, size_t i2) { return prop_p_begin[i1] > prop_p_begin[i2]; });
        
        std::vector<float> sorted_resp(params.num_populations);
        // sorted_var (alternative formulation, kept for reference)

        for (size_t i = 0; i < params.num_populations; ++i) {
            sorted_resp[i] = prop_p_begin[indices[i]];
        }

        const auto prev_p = std::vector<float>(begin, end);
        // population_responsibility_vector.inner_fill(prop_p);
        population_responsibility_vector.inner_fill(sorted_resp);
        // Invalidate cached log when vector changes
        population_responsibility_vector_log_valid_ = false;
        calculate_population_responsibility_vector_likelihood();
        refresh_all_samples_tx_logsumexp();
        const float new_llik = obs_llik_sum_new_ + tx_llik_sum_new;
        const float new_prior = calc_new_prior();
        const float new_post = new_llik * temp + new_prior;

        const float acceptance_ratio = new_post - get_posterior() + logAdj;

        bool reject = sampler.sample_log_mh_acceptance() > acceptance_ratio;

        if (!std::isfinite(new_post) or reject)
        {
            population_responsibility_vector.inner_fill(prev_p);
            // Invalidate cached log when reverting
            population_responsibility_vector_log_valid_ = false;
            restore_population_responsibility_vector_likelihood();
            refresh_all_samples_tx_logsumexp();
        }
        else
        {
            llik = new_llik;
            prior = new_prior;
            save_population_responsibility_vector_likelihood();
            ++population_responsibility_vector_accept.at({pop_idx});
            // for (size_t i = 0; i < params.num_populations; ++i) {
            //     sorted_var[i] = var_begin[indices[i]];
            // }
            // population_responsibility_vector_prop_var.inner_fill(sorted_var);
        }

        if (iteration < params.burnin and iteration > 15)
        {
            const float acceptance_rate = population_responsibility_vector_accept.at({pop_idx}) / (float(population_responsibility_vector_attempt.at({pop_idx})) + 1);
            const float update = (acceptance_rate - .23) / std::pow(float(population_responsibility_vector_attempt.at({pop_idx})) + 1, .5);
            population_responsibility_vector_prop_var.at({pop_idx}) = std::max(population_responsibility_vector_prop_var.at({pop_idx}) + update, .0001f);
        }
    }
        
}
