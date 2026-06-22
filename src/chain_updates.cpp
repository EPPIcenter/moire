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

namespace {
// Per-worker RNG stream for the parallel per-sample moves (collapsed eCOI moves
// and the standard COI / eff-COI / sample moves). Samples are independent given
// the shared parameters held fixed during each move, so each worker thread draws
// all of its proposals (jumps, latent-genotype reproposals) from its own stream.
// This lets the proposal draws -- previously a serial pre-draw phase -- run
// inside the parallel-over-samples region with no shared-RNG contention: each
// stream is touched by one thread at a time, so it stays race-free, and every
// proposal is still an independent draw from the proposal kernel (valid
// Metropolis-Hastings; correctness does not depend on which stream produces a
// given sample's draws). At a single thread this collapses to one hot sampler
// used sequentially -- the same cost profile as the old shared-sampler pre-draw.
// Like the rest of the C++ sampler, it is NOT seeded by R's set.seed, so runs
// are nondeterministic regardless of thread count.
Sampler &worker_rng()
{
    static thread_local Sampler rng;
    return rng;
}
}  // namespace

void Chain::update_m(int iteration)
{
    ProfileScope scope("Chain::update_m");
    const std::size_t N = genotyping_data.num_samples;
    const std::size_t n_loci = genotyping_data.num_loci;

    // Samples are independent given the shared parameters held fixed here: a
    // COI move on sample s changes only s's observation row, s's transmission
    // term, and s's COI prior (which is folded into the transmission logsumexp,
    // not calc_new_prior), so the MH ratio reduces to
    //   temp * (delta_obs_s + delta_tx_s) + latent-proposal adjustment,
    // with the global llik / prior cancelling. Each sample evaluates this from
    // its own per-sample slots using a private RNG stream, so the whole move --
    // proposal draws, latent reproposals, likelihood recompute, accept/reject --
    // runs inside one parallel-over-samples region. The deferred global sums and
    // update_p caches are reconciled serially afterward.
    prewarm_caches_for_parallel_sample_move();
    in_parallel_sample_region_ = true;
    moire_parallel::recalc_parallel_for(0, N, [&](std::size_t sample_idx) {
        Sampler &rng = worker_rng();
        const std::size_t prop_m =
            m.at({sample_idx}) + rng.sample_coi_delta(2);
        if (prop_m == 0 || prop_m > params.max_coi) {
            return;
        }

        const int prev_m = m.at({sample_idx});
        const float prev_r = r.at({sample_idx});
        const float old_contrib =
            obs_row_sum_new_[sample_idx] + tx_sample_logsumexp_new.at({sample_idx});

        m.at({sample_idx}) = prop_m;
        calculate_coi_likelihood(sample_idx);

        float adj_ratio = 0.f;
        for (std::size_t locus_idx = 0; locus_idx < n_loci; ++locus_idx) {
            const auto &lg = observation_model_->sample_latent_genotype(
                rng,
                genotyping_data.get_observed_alleles(sample_idx, locus_idx),
                m.at({sample_idx}), eps_pos.at({sample_idx}),
                eps_neg.at({sample_idx}));
            assign_latent_genotype_new(sample_idx, locus_idx, lg.value);
            lg_adj_new.at({sample_idx, locus_idx}) = lg.log_prob;
            adj_ratio += lg_adj_new.at({sample_idx, locus_idx}) -
                         lg_adj_old.at({sample_idx, locus_idx});
            calculate_observation_likelihood(sample_idx, locus_idx);
        }
        sync_obs_sum_for_sample(sample_idx);
        recalculate_transmission_for_sample_incremental(sample_idx, prev_m,
                                                        prev_r);

        const float new_contrib =
            obs_row_sum_new_[sample_idx] + tx_sample_logsumexp_new.at({sample_idx});
        const float mh_ratio = temp * (new_contrib - old_contrib) + adj_ratio;
        const float alpha = rng.sample_log_mh_acceptance();

        if (std::isfinite(new_contrib) && alpha <= mh_ratio) {
            // Accept: promote the per-sample new state into the old buffers.
            save_coi_likelihood(sample_idx);
            for (std::size_t locus_idx = 0; locus_idx < n_loci; ++locus_idx) {
                save_observation_likelihood(sample_idx, locus_idx);
                lg_adj_old.at({sample_idx, locus_idx}) =
                    lg_adj_new.at({sample_idx, locus_idx});
                auto [begin, end] =
                    latent_genotypes_new.inner_iterators({sample_idx, locus_idx});
                std::copy(begin, end,
                          latent_genotypes_old.inner_begin({sample_idx, locus_idx}));
            }
            for (std::size_t pop_idx = 0; pop_idx < params.num_populations; ++pop_idx) {
                for (std::size_t locus_idx = 0; locus_idx < n_loci; ++locus_idx) {
                    save_transmission_likelihood(pop_idx, sample_idx, locus_idx);
                }
            }
            m_accept.at({sample_idx}) += 1;
        } else {
            // Reject: restore the per-sample old state.
            m.at({sample_idx}) = prev_m;
            restore_coi_likelihood(sample_idx);
            for (std::size_t locus_idx = 0; locus_idx < n_loci; ++locus_idx) {
                restore_observation_likelihood(sample_idx, locus_idx);
                lg_adj_new.at({sample_idx, locus_idx}) =
                    lg_adj_old.at({sample_idx, locus_idx});
                restore_latent_genotype_new(sample_idx, locus_idx);
            }
            sync_obs_sum_for_sample(sample_idx);
            restore_transmission_for_sample_incremental(sample_idx);
        }
    });
    in_parallel_sample_region_ = false;
    reduce_likelihood_sums_after_parallel_sample_move();
}

void Chain::update_eff_coi(int iteration)
{
    ProfileScope scope("Chain::update_eff_coi");
    const std::size_t N = genotyping_data.num_samples;
    const std::size_t n_loci = genotyping_data.num_loci;

    // Joint (m, r) move via a continuous effective-COI proposal. As with
    // update_m the samples are independent: the move changes only sample s's
    // observation row, transmission term (which folds in s's COI prior) and s's
    // relatedness prior, so the MH ratio is
    //   temp * (delta_obs_s + delta_tx_s) + delta_relatedness_prior_s + adj,
    // entirely from s's own slots. Run in parallel over samples with a private
    // RNG; reconcile the deferred global sums serially afterward.
    prewarm_caches_for_parallel_sample_move();
    in_parallel_sample_region_ = true;
    moire_parallel::recalc_parallel_for(0, N, [&](std::size_t sample_idx) {
        Sampler &rng = worker_rng();
        const float curr_eff_coi =
            (m.at({sample_idx}) - 1) * (1.0f - r.at({sample_idx})) + 1.0f;
        const auto prop_adj = rng.sample_constrained(
            curr_eff_coi, m_r_var.at({sample_idx}), 1, params.max_coi);

        const float prop_eff_coi = std::get<0>(prop_adj);
        const float adj = std::get<1>(prop_adj);

        const std::size_t prop_m = m.at({sample_idx}) + rng.sample_coi_delta(2);
        const float prop_r = 1.0f - (prop_eff_coi - 1.0f) / (prop_m - 1.0f);

        if (prop_m <= 0 || prop_m > params.max_coi || prop_r > .99999 ||
            prop_r < .00001 || !std::isfinite(prop_r)) {
            return;
        }

        const int prev_m = m.at({sample_idx});
        const float prev_r = r.at({sample_idx});
        const float old_lik =
            obs_row_sum_new_[sample_idx] + tx_sample_logsumexp_new.at({sample_idx});
        const float old_rel_prior = relatedness_prior_new.at({sample_idx});

        m.at({sample_idx}) = prop_m;
        r.at({sample_idx}) = prop_r;

        float adj_ratio = adj;
        calculate_relatedness_likelihood(sample_idx);
        calculate_coi_likelihood(sample_idx);

        for (std::size_t locus_idx = 0; locus_idx < n_loci; ++locus_idx) {
            const auto &lg = observation_model_->sample_latent_genotype(
                rng,
                genotyping_data.get_observed_alleles(sample_idx, locus_idx),
                m.at({sample_idx}), eps_pos.at({sample_idx}),
                eps_neg.at({sample_idx}));
            assign_latent_genotype_new(sample_idx, locus_idx, lg.value);
            lg_adj_new.at({sample_idx, locus_idx}) = lg.log_prob;
            adj_ratio += lg_adj_new.at({sample_idx, locus_idx}) -
                         lg_adj_old.at({sample_idx, locus_idx});
            calculate_observation_likelihood(sample_idx, locus_idx);
        }
        sync_obs_sum_for_sample(sample_idx);
        recalculate_transmission_for_sample_incremental(sample_idx, prev_m,
                                                        prev_r);

        const float new_lik =
            obs_row_sum_new_[sample_idx] + tx_sample_logsumexp_new.at({sample_idx});
        const float new_rel_prior = relatedness_prior_new.at({sample_idx});
        const float mh_ratio = temp * (new_lik - old_lik) +
                               (new_rel_prior - old_rel_prior) + adj_ratio;
        const float alpha = rng.sample_log_mh_acceptance();

        if (!std::isfinite(new_lik) || !std::isfinite(mh_ratio) ||
            alpha > mh_ratio) {
            // Reject: restore per-sample state.
            m.at({sample_idx}) = prev_m;
            r.at({sample_idx}) = prev_r;
            restore_relatedness_likelihood(sample_idx);
            restore_coi_likelihood(sample_idx);
            for (std::size_t locus_idx = 0; locus_idx < n_loci; ++locus_idx) {
                lg_adj_new.at({sample_idx, locus_idx}) =
                    lg_adj_old.at({sample_idx, locus_idx});
                restore_latent_genotype_new(sample_idx, locus_idx);
                restore_observation_likelihood(sample_idx, locus_idx);
            }
            sync_obs_sum_for_sample(sample_idx);
            restore_transmission_for_sample_incremental(sample_idx);
        } else {
            // Accept: promote the per-sample new state into the old buffers.
            save_relatedness_likelihood(sample_idx);
            save_coi_likelihood(sample_idx);
            for (std::size_t locus_idx = 0; locus_idx < n_loci; ++locus_idx) {
                lg_adj_old.at({sample_idx, locus_idx}) =
                    lg_adj_new.at({sample_idx, locus_idx});
                const auto [begin, end] =
                    latent_genotypes_new.inner_iterators({sample_idx, locus_idx});
                std::copy(begin, end,
                          latent_genotypes_old.inner_begin({sample_idx, locus_idx}));
                save_observation_likelihood(sample_idx, locus_idx);
            }
            for (std::size_t pop_idx = 0; pop_idx < params.num_populations; ++pop_idx) {
                for (std::size_t locus_idx = 0; locus_idx < n_loci; ++locus_idx) {
                    save_transmission_likelihood(pop_idx, sample_idx, locus_idx);
                }
            }
            ++m_r_accept.at({sample_idx});
        }

        if (iteration < params.burnin && iteration > 15) {
            const float acceptanceRate =
                m_r_accept.at({sample_idx}) / float(iteration);
            const float update =
                (acceptanceRate - .23) / std::pow(iteration + 1, .5);
            m_r_var.at({sample_idx}) =
                std::max(m_r_var.at({sample_idx}) + update, .01f);
        }
    });
    in_parallel_sample_region_ = false;
    reduce_likelihood_sums_after_parallel_sample_move();
}

// Fully-collapsed effective-COI move. The discrete COI is integrated out, so the
// proposal acts on the continuous e in (1, max_coi] with latent genotypes held
// fixed; acceptance compares the marginalized transmission term M(G_s, e) at the
// proposed vs current e. The continuous e-prior f_p(e) is folded inside M, so no
// separate prior term enters.
//
// Each sample gets two sub-moves: (a) a local random-walk step for fine mixing,
// and (b) an independence step proposing e' from the population-0 continuous
// prior 1 + Gamma(k_0, k_0 / mu_plus_0). The independence step is essential: the
// per-sample e posterior under uniform relatedness has a secondary mode near
// e = 1 (low observed diversity explained as a few highly-related strains), and a
// pure random walk gets trapped there. Proposing from the prior covers the whole
// (1, max_coi] range so trapped samples can escape; for a single population the
// prior-densities cancel and acceptance is the pure marginal-likelihood ratio.
void Chain::update_ecoi(int iteration)
{
    ProfileScope scope("Chain::update_ecoi");
    const std::size_t N = genotyping_data.num_samples;
    const float k0 = ecoi_k.empty() ? 2.0f : ecoi_k[0];
    const float mu0 = ecoi_mu_plus.empty() ? 3.0f : ecoi_mu_plus[0];

    // Samples are independent in the collapsed marginal: each sample's MH ratio
    // reduces to temp*(new_marg - old_marg) + adjustment, with the global llik
    // and prior cancelling out. Each sample owns a private RNG stream, so its
    // proposal draws and the (expensive) marginal evaluations both run on the
    // worker thread -- no serial pre-draw phase, no shared-RNG contention.
    moire_parallel::recalc_parallel_for(0, N, [&](std::size_t s) {
        Sampler &rng = worker_rng();
        float curr_e = eff_coi.at({s});
        double cur_marg = static_cast<double>(ecoi_sample_marg_[s]);
        int accepted = 0;

        // (a) local random-walk step.
        {
            const auto [prop_e, adj_a] = rng.sample_constrained(
                curr_e, eff_coi_var.at({s}), 1, params.max_coi);
            const float u_a = rng.sample_log_mh_acceptance();
            if ((prop_e > 1.0f) && prop_e <= params.max_coi &&
                std::isfinite(prop_e)) {
                const double new_marg = ecoi_sample_rel_marginal_llik(
                    s, static_cast<double>(prop_e));
                if (std::isfinite(new_marg)) {
                    const float mh_ratio =
                        static_cast<float>(temp * (new_marg - cur_marg)) + adj_a;
                    if (u_a <= mh_ratio) {
                        curr_e = prop_e;
                        cur_marg = new_marg;
                        ++accepted;
                    }
                }
            }
        }

        // (b) independence step from the population-0 continuous prior.
        {
            // rgamma(shape, scale) draws Gamma(k0, mu0/k0) == the thread-unsafe
            // R::rgamma(k0, scale = mu0/k0) the serial version used, but from this
            // sample's own engine so it is safe inside the parallel region.
            const float prop_e = 1.0f + rng.rgamma(k0, mu0 / k0);
            const float u_b = rng.sample_log_mh_acceptance();
            if ((prop_e > 1.0f) && prop_e <= params.max_coi &&
                std::isfinite(prop_e)) {
                const double new_marg = ecoi_sample_rel_marginal_llik(
                    s, static_cast<double>(prop_e));
                // Hastings correction for the asymmetric prior proposal q = f_0.
                const double logq_new = ecoi_log_f(0, static_cast<double>(prop_e));
                const double logq_old = ecoi_log_f(0, static_cast<double>(curr_e));
                if (std::isfinite(new_marg) && std::isfinite(logq_new)) {
                    const float mh_ratio =
                        static_cast<float>(temp * (new_marg - cur_marg)) +
                        static_cast<float>(logq_old - logq_new);
                    if (u_b <= mh_ratio) {
                        curr_e = prop_e;
                        cur_marg = new_marg;
                        ++accepted;
                    }
                }
            }
        }

        eff_coi.at({s}) = curr_e;
        ecoi_sample_marg_[s] = static_cast<float>(cur_marg);
        if (accepted > 0) {
            eff_coi_accept.at({s}) += accepted;
        }

        if (iteration < params.burnin && iteration > 15) {
            const float acceptanceRate = eff_coi_accept.at({s}) / float(iteration);
            const float update = (acceptanceRate - .23) / std::pow(iteration + 1, .5);
            eff_coi_var.at({s}) =
                std::max(eff_coi_var.at({s}) + update, .01f);
        }
    });

    // Reduce the per-sample marginals into the running totals (serial; the
    // parallel region only touched per-sample state).
    double marg_sum = 0.0;
    for (std::size_t s = 0; s < N; ++s) {
        marg_sum += static_cast<double>(ecoi_sample_marg_[s]);
    }
    ecoi_marg_sum_ = marg_sum;
    llik = obs_llik_sum_new_ + static_cast<float>(marg_sum);

    // Optional audit: recompute each sample's marginal at its committed e and
    // confirm the parallel move left ecoi_sample_marg_ / ecoi_marg_sum_ exactly
    // consistent (catches any cross-sample corruption from the parallel region).
    if (std::getenv("MOIRE_ECOI_CHECK_MOVES") != nullptr) {
        double recomputed = 0.0;
        double max_abs = 0.0;
        for (std::size_t s = 0; s < N; ++s) {
            const double tr = ecoi_sample_rel_marginal_llik(
                s, static_cast<double>(eff_coi.at({s})));
            max_abs = std::max(
                max_abs, std::abs(tr - static_cast<double>(ecoi_sample_marg_[s])));
            recomputed += tr;
        }
        Rcpp::Rcout << "[ecoi-move-audit] update_ecoi max|stored-recompute|="
                    << max_abs << " sum_diff="
                    << std::abs(recomputed - ecoi_marg_sum_) << "\n";
    }
}

// Fully-collapsed latent-genotype move. Genotypes are reproposed from the
// COI-independent proposal (whose density is evaluable for the reverse move) and
// accepted against the observation likelihood plus the marginalized transmission
// term. Effective COI and error rates are held fixed.
void Chain::ecoi_check_latent_move_consistency(const char *label)
{
    const std::size_t N = genotyping_data.num_samples;
    double rec_marg = 0.0, rec_obs = 0.0, max_marg = 0.0, max_obs = 0.0;
    for (std::size_t s = 0; s < N; ++s) {
        const double tr = ecoi_sample_marginal_llik(s);
        max_marg = std::max(
            max_marg, std::abs(tr - static_cast<double>(ecoi_sample_marg_[s])));
        rec_marg += tr;
        const auto [begin, end] = observation_llik_new.inner_iterators({s});
        double row = 0.0;
        for (auto it = begin; it != end; ++it) {
            row += static_cast<double>(*it);
        }
        max_obs = std::max(
            max_obs, std::abs(row - static_cast<double>(obs_row_sum_new_[s])));
        rec_obs += row;
    }
    Rcpp::Rcout << "[ecoi-move-audit] " << label
                << " max|marg|=" << max_marg << " max|obsrow|=" << max_obs
                << " marg_sum_diff=" << std::abs(rec_marg - ecoi_marg_sum_)
                << " obs_sum_diff="
                << std::abs(rec_obs - static_cast<double>(obs_llik_sum_new_))
                << "\n";
}

void Chain::update_latent_marginal(int iteration)
{
    ProfileScope scope("Chain::update_latent_marginal");
    (void)iteration;
    const std::size_t N = genotyping_data.num_samples;
    const std::size_t L = genotyping_data.num_loci;

    // Samples are independent in the collapsed likelihood, so each sample's MH
    // ratio reduces to its own (obs + marginal) delta (global llik and prior
    // cancel). Each sample owns a private RNG stream, so the genotype-proposal
    // draws (previously a serial pre-draw phase) and the obs/marginal recompute
    // both run on the worker thread. Every write is to per-sample storage
    // (genotypes, obs-llik rows, obs_row_sum_new_, ecoi_sample_marg_,
    // sample_accept) or the sample's own support-cache column.
    moire_parallel::recalc_parallel_for(0, N, [&](std::size_t s) {
        Sampler &rng = worker_rng();
        const float eps_neg_s = eps_neg.at({s});

        // Propose new latent genotypes for this sample from its own RNG stream.
        float adj_ratio = 0.0f;
        for (std::size_t l = 0; l < L; ++l) {
            if (genotyping_data.is_missing(s, l)) {
                continue;
            }
            const float eps_pos_l = eps_pos_at(s, l);
            const auto observed = genotyping_data.get_observed_alleles(s, l);
            const auto old_support = latent_allele_support(s, l);
            const float rev_logq = observation_model_->latent_genotype_log_prob_marginal(
                old_support, observed, eps_pos_l, eps_neg_s);

            const auto lg = observation_model_->propose_latent_genotype_marginal(
                rng, observed, eps_pos_l, eps_neg_s);
            assign_latent_genotype_new(s, l, lg.value);
            lg_adj_new.at({s, l}) = lg.log_prob;
            adj_ratio += rev_logq - lg.log_prob;
        }
        const float accept_u = rng.sample_log_mh_acceptance();

        // Recompute obs likelihood + marginal for the proposal and accept/reject.
        ecoi_support_cache_invalidate_sample(s);
        for (std::size_t l = 0; l < L; ++l) {
            calculate_observation_likelihood(s, l);
        }
        double new_row = 0.0;
        {
            const auto [begin, end] = observation_llik_new.inner_iterators({s});
            for (auto it = begin; it != end; ++it) {
                new_row += static_cast<double>(*it);
            }
        }
        const double old_row = static_cast<double>(obs_row_sum_new_[s]);
        // e and the monoclonal indicator are fixed here, so the per-sample
        // marginal is recomputed at the current state and the e-prior is fixed.
        const double new_marg = ecoi_sample_marginal_llik(s);
        const double old_marg = static_cast<double>(ecoi_sample_marg_[s]);
        const float mh_ratio = static_cast<float>(
                                   temp * ((new_row - old_row) +
                                           (new_marg - old_marg))) +
                               adj_ratio;

        if (std::isfinite(new_marg) && std::isfinite(new_row) &&
            accept_u <= mh_ratio) {
            for (std::size_t l = 0; l < L; ++l) {
                if (genotyping_data.is_missing(s, l)) {
                    continue;
                }
                lg_adj_old.at({s, l}) = lg_adj_new.at({s, l});
                const auto [begin, end] =
                    latent_genotypes_new.inner_iterators({s, l});
                std::copy(begin, end, latent_genotypes_old.inner_begin({s, l}));
                save_observation_likelihood(s, l);
            }
            obs_row_sum_new_[s] = static_cast<float>(new_row);
            ecoi_sample_marg_[s] = static_cast<float>(new_marg);
            ++sample_accept.at({s});
        } else {
            for (std::size_t l = 0; l < L; ++l) {
                if (genotyping_data.is_missing(s, l)) {
                    continue;
                }
                lg_adj_new.at({s, l}) = lg_adj_old.at({s, l});
                restore_latent_genotype_new(s, l);
                restore_observation_likelihood(s, l);
            }
            // Supports restored; rebuild the cache from the restored supports.
            ecoi_support_cache_invalidate_sample(s);
        }
    });

    // Reduce running totals (serial; the parallel region only touched per-sample
    // state, so obs_llik_sum_new_ is rebuilt from the committed row sums).
    double obs_sum = 0.0;
    double marg_sum = 0.0;
    for (std::size_t s = 0; s < N; ++s) {
        obs_sum += static_cast<double>(obs_row_sum_new_[s]);
        marg_sum += static_cast<double>(ecoi_sample_marg_[s]);
    }
    obs_llik_sum_new_ = static_cast<float>(obs_sum);
    ecoi_marg_sum_ = marg_sum;
    llik = obs_llik_sum_new_ + static_cast<float>(marg_sum);

    if (std::getenv("MOIRE_ECOI_CHECK_MOVES") != nullptr) {
        ecoi_check_latent_move_consistency("latent_marginal");
    }
}

// Joint (effective COI, latent-genotype) move. Per sample we propose e' from the
// population-0 continuous prior (independence) and simultaneously repropose all of
// the sample's latent genotypes from the COI-independent proposal, accepting the
// pair in one step. Because the COI-independent proposal draws supports sized to
// the observed diversity (independent of e), a large e jump arrives with supports
// that can actually support it -- so a sample trapped in the e ~ 1 / small-support
// joint mode can escape, which a fixed-G e move cannot do. The e-prior f_0 cancels
// against the proposal density for a single population, leaving the marginal +
// observation likelihood ratio and the genotype-proposal Hastings term.
void Chain::update_ecoi_latent_joint(int iteration)
{
    ProfileScope scope("Chain::update_ecoi_latent_joint");
    (void)iteration;
    const std::size_t N = genotyping_data.num_samples;
    const std::size_t L = genotyping_data.num_loci;
    const float k0 = ecoi_k.empty() ? 2.0f : ecoi_k[0];
    const float mu0 = ecoi_mu_plus.empty() ? 3.0f : ecoi_mu_plus[0];

    // Per sample (parallel, each on its own RNG stream): draw the joint
    // (e', latent-genotype) proposal and accept/reject against this sample's own
    // (obs + marginal) delta plus the genotype/e Hastings terms. Samples with an
    // invalid e' proposal are skipped entirely (no genotype proposal, no accept
    // draw), matching the sequential version. Every write is to per-sample
    // storage or the sample's own support-cache column.
    moire_parallel::recalc_parallel_for(0, N, [&](std::size_t s) {
        Sampler &rng = worker_rng();
        const float curr_e = eff_coi.at({s});
        // rgamma(shape, scale) draws Gamma(k0, mu0/k0) == the thread-unsafe
        // R::rgamma(k0, scale = mu0/k0) the serial version used.
        const float pe = 1.0f + rng.rgamma(k0, mu0 / k0);
        if (!(pe > 1.0f) || pe > params.max_coi || !std::isfinite(pe)) {
            return;
        }
        const double logq_old_e = ecoi_log_f(0, static_cast<double>(curr_e));
        const double logq_new_e = ecoi_log_f(0, static_cast<double>(pe));
        if (!std::isfinite(logq_new_e)) {
            return;
        }
        const float e_hastings = static_cast<float>(logq_old_e - logq_new_e);

        const float eps_neg_s = eps_neg.at({s});
        float adj_ratio = 0.0f;
        for (std::size_t l = 0; l < L; ++l) {
            if (genotyping_data.is_missing(s, l)) {
                continue;
            }
            const float eps_pos_l = eps_pos_at(s, l);
            const auto observed = genotyping_data.get_observed_alleles(s, l);
            const auto old_support = latent_allele_support(s, l);
            const float rev_logq = observation_model_->latent_genotype_log_prob_marginal(
                old_support, observed, eps_pos_l, eps_neg_s);

            const auto lg = observation_model_->propose_latent_genotype_marginal(
                rng, observed, eps_pos_l, eps_neg_s);
            assign_latent_genotype_new(s, l, lg.value);
            lg_adj_new.at({s, l}) = lg.log_prob;
            adj_ratio += rev_logq - lg.log_prob;
        }
        const float accept_u = rng.sample_log_mh_acceptance();

        ecoi_support_cache_invalidate_sample(s);
        for (std::size_t l = 0; l < L; ++l) {
            calculate_observation_likelihood(s, l);
        }
        double new_row = 0.0;
        {
            const auto [begin, end] = observation_llik_new.inner_iterators({s});
            for (auto it = begin; it != end; ++it) {
                new_row += static_cast<double>(*it);
            }
        }
        const double old_row = static_cast<double>(obs_row_sum_new_[s]);
        const double new_marg =
            ecoi_sample_rel_marginal_llik(s, static_cast<double>(pe));
        const double old_marg = static_cast<double>(ecoi_sample_marg_[s]);
        const float mh_ratio = static_cast<float>(
                                   temp * ((new_row - old_row) +
                                           (new_marg - old_marg))) +
                               adj_ratio + e_hastings;

        if (std::isfinite(new_marg) && std::isfinite(new_row) &&
            accept_u <= mh_ratio) {
            for (std::size_t l = 0; l < L; ++l) {
                if (genotyping_data.is_missing(s, l)) {
                    continue;
                }
                lg_adj_old.at({s, l}) = lg_adj_new.at({s, l});
                const auto [begin, end] =
                    latent_genotypes_new.inner_iterators({s, l});
                std::copy(begin, end, latent_genotypes_old.inner_begin({s, l}));
                save_observation_likelihood(s, l);
            }
            eff_coi.at({s}) = pe;
            obs_row_sum_new_[s] = static_cast<float>(new_row);
            ecoi_sample_marg_[s] = static_cast<float>(new_marg);
            ++sample_accept.at({s});
        } else {
            for (std::size_t l = 0; l < L; ++l) {
                if (genotyping_data.is_missing(s, l)) {
                    continue;
                }
                lg_adj_new.at({s, l}) = lg_adj_old.at({s, l});
                restore_latent_genotype_new(s, l);
                restore_observation_likelihood(s, l);
            }
            ecoi_support_cache_invalidate_sample(s);
        }
    });

    // Reduce running totals (serial).
    double obs_sum = 0.0;
    double marg_sum = 0.0;
    for (std::size_t s = 0; s < N; ++s) {
        obs_sum += static_cast<double>(obs_row_sum_new_[s]);
        marg_sum += static_cast<double>(ecoi_sample_marg_[s]);
    }
    obs_llik_sum_new_ = static_cast<float>(obs_sum);
    ecoi_marg_sum_ = marg_sum;
    llik = obs_llik_sum_new_ + static_cast<float>(marg_sum);

    if (std::getenv("MOIRE_ECOI_CHECK_MOVES") != nullptr) {
        ecoi_check_latent_move_consistency("ecoi_latent_joint");
    }
}

// Population-e hierarchy move. With the per-population continuous e-prior folded
// inside each sample's population mixture, the per-population (mu_plus_p, k_p)
// enter the marginal likelihood (through f_p), so changing them shifts every
// sample's marginal. For each population we take a log-scale random-walk
// Metropolis step and accept against the full collapsed marginal sum plus the
// hyperpriors (a responsibility-consistent move: a sample's pull on population p
// is exactly its posterior weight in the LSE, no z augmentation needed).
void Chain::update_population_e(int iteration)
{
    ProfileScope scope("Chain::update_population_e");
    (void)iteration;
    const std::size_t n = genotyping_data.num_samples;
    if (n == 0) {
        return;
    }
    const float obs_llik = obs_llik_sum_new_;

    auto pop_indices = std::vector<std::size_t>(params.num_populations);
    std::iota(pop_indices.begin(), pop_indices.end(), 0);
    sampler.shuffle_vec(pop_indices);

    for (const std::size_t pop : pop_indices)
    {
        const float cur_mu = ecoi_mu_plus[pop];
        const float cur_k = ecoi_k[pop];

        const float mu_prop = cur_mu * std::exp(sampler.sample_epsilon(0.0f, ecoi_mu_log_sd));
        const float k_prop = cur_k * std::exp(sampler.sample_epsilon(0.0f, ecoi_k_log_sd));
        if (!(mu_prop > 0.0f) || !(k_prop > 0.0f) || !std::isfinite(mu_prop) ||
            !std::isfinite(k_prop)) {
            continue;
        }

        const double saved_hyper = ecoi_hyperprior_;
        ecoi_mu_plus[pop] = mu_prop;
        ecoi_k[pop] = k_prop;

        const double new_marg_sum =
            recompute_collapsed_marginals(ecoi_sample_marg_scratch_);
        ecoi_hyperprior_ = ecoi_population_hyperprior();
        const float new_llik = obs_llik + static_cast<float>(new_marg_sum);
        const float new_prior = calc_new_prior();
        const float new_post = new_llik * temp + new_prior;

        // Jacobian of the log-scale random walk.
        const double log_jac = std::log(static_cast<double>(mu_prop)) -
                               std::log(static_cast<double>(cur_mu)) +
                               std::log(static_cast<double>(k_prop)) -
                               std::log(static_cast<double>(cur_k));
        const double mh_ratio = new_post - get_posterior() + log_jac;
        ++ecoi_pop_attempt;

        if (std::isfinite(new_post) &&
            sampler.sample_log_mh_acceptance() <= static_cast<float>(mh_ratio))
        {
            std::swap(ecoi_sample_marg_, ecoi_sample_marg_scratch_);
            ecoi_marg_sum_ = new_marg_sum;
            llik = new_llik;
            prior = new_prior;
            ++ecoi_pop_accept;
        }
        else
        {
            ecoi_mu_plus[pop] = cur_mu;
            ecoi_k[pop] = cur_k;
            ecoi_hyperprior_ = saved_hyper;
        }
    }
}

void Chain::update_r(int iteration)
{
    ProfileScope scope("Chain::update_r");
    const std::size_t N = genotyping_data.num_samples;

    // Relatedness enters only the transmission process (which folds in s's COI
    // prior); the observation row is unchanged and cancels. Per-sample MH ratio
    // is temp*(delta_tx_s) + delta_relatedness_prior_s + adj, all from s's own
    // slots, so samples update independently in parallel. The incremental
    // transmission recompute uses s's current COI (unchanged here) and the new
    // relatedness value.
    prewarm_caches_for_parallel_sample_move();
    in_parallel_sample_region_ = true;
    moire_parallel::recalc_parallel_for(0, N, [&](std::size_t sample_idx) {
        Sampler &rng = worker_rng();
        const auto prop_adj = rng.sample_constrained(
            r.at({sample_idx}), r_var.at({sample_idx}), .00001, .99999);
        const float prop_r = std::get<0>(prop_adj);
        const float adj = std::get<1>(prop_adj);

        const float old_tx = tx_sample_logsumexp_new.at({sample_idx});
        const float old_rel_prior = relatedness_prior_new.at({sample_idx});

        const float prev_r = r.at({sample_idx});
        r.at({sample_idx}) = prop_r;
        calculate_relatedness_likelihood(sample_idx);
        recalculate_transmission_for_sample_incremental(
            sample_idx, static_cast<int>(m.at({sample_idx})), prev_r);

        const float new_tx = tx_sample_logsumexp_new.at({sample_idx});
        const float mh_ratio =
            temp * (new_tx - old_tx) +
            (relatedness_prior_new.at({sample_idx}) - old_rel_prior) + adj;

        if (!std::isfinite(new_tx) || rng.sample_log_mh_acceptance() > mh_ratio) {
            r.at({sample_idx}) = prev_r;
            restore_relatedness_likelihood(sample_idx);
            restore_transmission_for_sample_incremental(sample_idx);
        } else {
            save_relatedness_likelihood(sample_idx);
            for (std::size_t pop_idx = 0; pop_idx < params.num_populations; ++pop_idx) {
                for (std::size_t locus_idx = 0; locus_idx < genotyping_data.num_loci;
                     ++locus_idx) {
                    save_transmission_likelihood(pop_idx, sample_idx, locus_idx);
                }
            }
            ++r_accept.at({sample_idx});
        }

        if (iteration < params.burnin && iteration > 15) {
            const float acceptanceRate = r_accept.at({sample_idx}) / float(iteration);
            const float update =
                (acceptanceRate - .23) / std::pow(iteration + 1, .5);
            r_var.at({sample_idx}) =
                std::max(r_var.at({sample_idx}) + update, .0001f);
        }
    });
    in_parallel_sample_region_ = false;
    reduce_likelihood_sums_after_parallel_sample_move();
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
void Chain::update_p_standard_parallel(int iteration)
{
    ProfileScope scope("Chain::update_p");
    const std::size_t pop_idx = 0;
    const std::size_t n_loci = genotyping_data.num_loci;
    const std::size_t n_samples = genotyping_data.num_samples;

    // Prewarm shared caches on the dispatch thread: the per-locus group cache is
    // built lazily otherwise, and the transmission cache must be valid before the
    // region so no worker triggers a rebuild mid-flight.
    if (!tx_llik_cache_valid_) {
        rebuild_transmission_llik_cache();
    }
    for (std::size_t locus_idx = 0; locus_idx < n_loci; ++locus_idx) {
        ensure_update_p_locus_group_cache(locus_idx);
    }

    // Each locus owns its transmission column [*, 0, locus] and its (pop, locus)
    // p column, group cache and PAM slots; the 3 reps for a locus run serially on
    // one worker exactly like the serial path. A proposal changes only locus L's
    // column, and with one population sample s's logsumexp == coi_prior_s +
    // loci_sum_s, so the change to the total transmission llik is exactly the
    // column-sum delta. Score each proposal against that delta (obs and prior are
    // invariant under p moves), defer the global tx-sum bookkeeping, and rebuild
    // it once after the region.
    in_parallel_sample_region_ = true;
    moire_parallel::recalc_parallel_for(0, n_loci, [&](std::size_t locus_idx) {
        // Already parallel over loci; keep the inner per-sample transmission
        // recompute serial to avoid nested spawning.
        const bool prev_disable = moire_parallel::disable_nested_parallelism;
        moire_parallel::disable_nested_parallelism = true;
        Sampler &rng = worker_rng();

        const std::size_t num_alleles = p.ragged_dimensions(locus_idx);
        const auto [begin, end] = p.inner_iterators({pop_idx, locus_idx});

        thread_local std::vector<float> prev_p_ws;
        thread_local std::vector<float> old_col;
        old_col.assign(n_samples, 0.f);
        float col_sum = 0.f;
        for (std::size_t s = 0; s < n_samples; ++s) {
            old_col[s] = transmission_llik_new.unchecked_at({s, pop_idx, locus_idx});
            col_sum += old_col[s];
        }

        std::size_t rep = 3;
        while (rep-- > 0) {
            const std::size_t allele_idx =
                rng.sample_random_int(0, num_alleles - 1);
            ++p_attempt.at({pop_idx, locus_idx, allele_idx});

            auto logitPropP = UtilFunctions::logitVec(begin, end);
            const float logitCurr = logitPropP[allele_idx];
            const float logitProp = rng.sample_epsilon(
                logitCurr, p_prop_var.at({pop_idx, locus_idx, allele_idx}));

            const auto currLogPQ = UtilFunctions::log_pq(logitCurr);
            const auto propLogPQ = UtilFunctions::log_pq(logitProp);

            logitPropP.erase(logitPropP.begin() + allele_idx);
            const float ls = propLogPQ.second - UtilFunctions::logitSum(logitPropP);
            logitPropP = UtilFunctions::logitScale(logitPropP, ls);
            logitPropP.insert(logitPropP.begin() + allele_idx, logitProp);

            const float logAdj =
                (currLogPQ.first - propLogPQ.first) +
                (num_alleles - 1) * (currLogPQ.second - propLogPQ.second);

            const auto prop_p = UtilFunctions::expitVec(logitPropP);
            bool sub_threshold_flag = false;
            for (const auto el : prop_p) {
                if (el < 1e-4 || !std::isfinite(el)) {
                    sub_threshold_flag = true;
                    break;
                }
            }
            if (sub_threshold_flag) {
                break;
            }

            prev_p_ws.assign(begin, end);
            p.inner_fill({pop_idx, locus_idx}, prop_p);
            recalculate_transmission_at_locus_after_p_change(
                pop_idx, locus_idx,
                std::span<const float>(prev_p_ws.data(), prev_p_ws.size()),
                allele_idx);

            float new_col_sum = 0.f;
            for (std::size_t s = 0; s < n_samples; ++s) {
                new_col_sum +=
                    transmission_llik_new.unchecked_at({s, pop_idx, locus_idx});
            }
            const float mh_ratio = temp * (new_col_sum - col_sum) + logAdj;

            if (!std::isfinite(new_col_sum) ||
                rng.sample_log_mh_acceptance() > mh_ratio) {
                // Reject: restore p and the transmission column.
                p.inner_fill({pop_idx, locus_idx}, prev_p_ws);
                for (std::size_t s = 0; s < n_samples; ++s) {
                    transmission_llik_new.unchecked_at({s, pop_idx, locus_idx}) =
                        old_col[s];
                }
            } else {
                // Accept: column already holds the new values.
                col_sum = new_col_sum;
                for (std::size_t s = 0; s < n_samples; ++s) {
                    old_col[s] =
                        transmission_llik_new.unchecked_at({s, pop_idx, locus_idx});
                }
                ++p_accept.at({pop_idx, locus_idx, allele_idx});
            }

            if (iteration < params.burnin &&
                p_attempt.at({pop_idx, locus_idx, allele_idx}) > 15) {
                const float acceptanceRate =
                    (p_accept.at({pop_idx, locus_idx, allele_idx}) + 1) /
                    (float(p_attempt.at({pop_idx, locus_idx, allele_idx})) + 1);
                p_prop_var.at({pop_idx, locus_idx, allele_idx}) +=
                    (acceptanceRate - .23) /
                    std::pow(p_attempt.at({pop_idx, locus_idx, allele_idx}) + 1, .5);
                p_prop_var.at({pop_idx, locus_idx, allele_idx}) =
                    std::max(p_prop_var.at({pop_idx, locus_idx, allele_idx}), .01f);
            }
        }

        // Sync the committed column into the old buffer so subsequent moves'
        // per-sample restore sees a consistent baseline.
        for (std::size_t s = 0; s < n_samples; ++s) {
            transmission_llik_old.unchecked_at({s, pop_idx, locus_idx}) =
                transmission_llik_new.unchecked_at({s, pop_idx, locus_idx});
        }
        moire_parallel::disable_nested_parallelism = prev_disable;
    });
    in_parallel_sample_region_ = false;

    // Rebuild the global transmission sums from the per-locus columns and refresh
    // the chain likelihood. Prior is invariant under allele-frequency moves.
    rebuild_transmission_llik_cache();
    llik = calc_new_likelihood();
}

void Chain::update_p(int iteration)
{
    ProfileScope scope("Chain::update_p");
    // Single-population, non-marginal allele frequencies parallelize cleanly over
    // loci (the transmission llik is additively separable across loci when there
    // is one population), so use the loci-parallel path there.
    if (!params.marginal_ecoi && params.num_populations == 1 &&
        pam_fast_paths::tx_opts_enabled()) {
        update_p_standard_parallel(iteration);
        return;
    }
    // Observation and prior are unchanged across proposals within this update.
    const float obs_llik = obs_llik_sum_new_;
    if (!tx_llik_cache_valid_) {
        rebuild_transmission_llik_cache();
    }

    // Single-population eCOI fast path: build each sample's per-m transmission
    // decomposition once, then update only the changed locus' contribution per
    // proposal (O(max_coi) per sample) instead of rebuilding every locus.
    const bool ecoi_fast = params.marginal_ecoi && params.num_populations == 1 &&
                           ecoi_fast_update_p_;
    if (ecoi_fast) {
        ecoi_up_build_decomp();
        // One-time consistency guard: the rebuilt decomposition must reproduce
        // the persistent marginal sum. If not, fall back to the safe full
        // recompute for the rest of the run.
        double check_sum = 0.0;
        for (std::size_t s = 0; s < genotyping_data.num_samples; ++s) {
            check_sum += up_L_[s];
        }
        const double tol = 1e-2 * (1.0 + std::abs(ecoi_marg_sum_));
        if (!std::isfinite(check_sum) || std::abs(check_sum - ecoi_marg_sum_) > tol) {
            ecoi_fast_update_p_ = false;
        }
    }
    const bool ecoi_fast_active = ecoi_fast && ecoi_fast_update_p_;
    const std::size_t up_stride = up_m_stride_;

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
            // Fast path: the OLD per-(sample, m) transmission terms for this locus
            // depend only on the locus' current p (not on which allele a proposal
            // perturbs) and on the fixed e_s / latent support, so they are
            // identical across this locus' reps. Derive them ONCE here rather than
            // recomputing them on every proposal -- the dominant per-locus kernel.
            // After an accepted proposal the freshly computed new terms (up_tnew_)
            // are the current ones, so the accept branch copies them back into
            // up_told_ to keep this invariant across reps.
            if (ecoi_fast_active) {
                ecoi_up_locus_terms(locus_idx, up_told_);
            }
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

                if (params.marginal_ecoi) {
                    // Only this (pop, locus) column of allele frequencies changed;
                    // invalidate just those cached LocusSupport entries.
                    ecoi_support_cache_invalidate_locus(pop_idx, locus_idx);
                }

                if (!params.marginal_ecoi) {
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
                double collapsed_sum = 0.0;
                {
                    ProfileScope scope("Chain::update_p::calc_post");
                    if (ecoi_fast_active) {
                        // Incremental: only locus_idx's transmission term changed.
                        // Fused new-locus-term + marginal reassembly: each sample
                        // rebuilds its new LocusSupport for the changed locus,
                        // computes the new per-m transmission terms (up_tnew_),
                        // and folds the (new - old) delta straight into its
                        // marginal LSE -- all in ONE parallel region instead of
                        // two (locus_terms map + reassembly map). At small N the
                        // per-proposal spawn/join overhead dominates this inner
                        // loop (~3 spawns per proposal, ~hundreds of proposals per
                        // sweep), so halving the spawns here is the win. Each
                        // sample writes only its own up_tnew_ / up_newL_ slices
                        // (thread_local scratch); the result is identical to the
                        // unfused two-pass version.
                        const int max_coi = static_cast<int>(params.max_coi);
                        const std::size_t n_samples =
                            genotyping_data.num_samples;
                        const std::size_t locus = locus_idx;
                        moire_parallel::recalc_parallel_for(
                            0, n_samples, [&](std::size_t s) {
                                static thread_local std::vector<double> terms;
                                terms.clear();
                                const ecoi_marginal::LocusSupport *ls =
                                    ecoi_cached_support(s, 0, locus);
                                const double e =
                                    static_cast<double>(eff_coi.at({s}));
                                for (int m = up_mlo_[s]; m <= max_coi; ++m) {
                                    const std::size_t idx =
                                        s * up_stride + static_cast<std::size_t>(m);
                                    // New per-locus transmission term (matches
                                    // ecoi_up_locus_terms exactly).
                                    double tnew = 0.0;
                                    if (ls != nullptr) {
                                        const double r = ecoi_marginal::r_of(m, e);
                                        if (r > 0.0 && r < 1.0) {
                                            const double t =
                                                ecoi_marginal::tx_loglik_locus(*ls, m, r);
                                            tnew = std::isfinite(t) ? t : 0.0;
                                        }
                                    }
                                    up_tnew_[idx] = tnew;
                                    if (up_ninf_[idx] > 0 ||
                                        !std::isfinite(up_a_[idx])) {
                                        terms.push_back(ecoi_marginal::kNegInf);
                                        continue;
                                    }
                                    terms.push_back(
                                        up_a_[idx] + up_slog_[idx] +
                                        (tnew - up_told_[idx]));
                                }
                                up_newL_[s] =
                                    up_foff_[s] +
                                    ecoi_marginal::log_sum_exp(
                                        std::span<const double>(terms));
                            });
                        collapsed_sum = 0.0;
                        for (std::size_t s = 0; s < n_samples; ++s) {
                            collapsed_sum += up_newL_[s];
                        }
                        new_llik = obs_llik + static_cast<float>(collapsed_sum);

                        // Optional exact-correctness audit: compare the
                        // incremental sum to a full recompute (same float inputs,
                        // different summation order -> expect ~1e-5).
                        static const bool ecoi_audit =
                            std::getenv("MOIRE_ECOI_CHECK_FAST") != nullptr;
                        if (ecoi_audit) {
                            const double full =
                                recompute_collapsed_marginals(ecoi_sample_marg_scratch_);
                            const double d = std::abs(full - collapsed_sum);
                            if (d > ecoi_fast_audit_max_) {
                                ecoi_fast_audit_max_ = d;
                            }
                        }
                    } else if (params.marginal_ecoi) {
                        // p changes shift every sample's marginal; obs is fixed.
                        collapsed_sum = recompute_collapsed_marginals(ecoi_sample_marg_scratch_);
                        new_llik = obs_llik + static_cast<float>(collapsed_sum);
                    } else {
                        new_llik = obs_llik + tx_llik_sum_new;
                    }
                    new_post = new_llik * temp + prior;
                }

                const float acceptanceRatio = new_post - get_posterior() + logAdj;

                if (!std::isfinite(new_post) or
                    sampler.sample_log_mh_acceptance() > acceptanceRatio)
                {
                    p.inner_fill({pop_idx, locus_idx}, update_p_prev_p_ws_);
                    if (params.marginal_ecoi) {
                        // Restored p; rebuild this column from the restored values.
                        ecoi_support_cache_invalidate_locus(pop_idx, locus_idx);
                    }
                    if (!params.marginal_ecoi) {
                        restore_transmission_column_change(pop_idx, locus_idx);
                    }
                }
                else
                {
                    ProfileScope scope("Chain::update_p::accept_save");
                    llik = new_llik;
                    if (ecoi_fast_active) {
                        // Commit the locus delta into the persistent decomposition
                        // and the per-sample marginal cache.
                        const int max_coi = static_cast<int>(params.max_coi);
                        for (std::size_t s = 0; s < genotyping_data.num_samples; ++s) {
                            for (int m = up_mlo_[s]; m <= max_coi; ++m) {
                                const std::size_t idx =
                                    s * up_stride + static_cast<std::size_t>(m);
                                up_slog_[idx] += (up_tnew_[idx] - up_told_[idx]);
                                // Accepted: the new locus terms are now current,
                                // so they become up_told_ for this locus' next rep.
                                up_told_[idx] = up_tnew_[idx];
                            }
                            up_L_[s] = up_newL_[s];
                            ecoi_sample_marg_[s] = static_cast<float>(up_newL_[s]);
                        }
                        ecoi_marg_sum_ = collapsed_sum;
                    } else if (params.marginal_ecoi) {
                        std::swap(ecoi_sample_marg_, ecoi_sample_marg_scratch_);
                        ecoi_marg_sum_ = collapsed_sum;
                    } else {
                        moire_parallel::parallel_for(0, genotyping_data.num_samples, [&](std::size_t sample_idx) {
                            save_transmission_likelihood(pop_idx, sample_idx, locus_idx);
                        });
                    }
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

    if (ecoi_fast_active && std::getenv("MOIRE_ECOI_CHECK_FAST") != nullptr) {
        Rcpp::Rcout << "[ecoi-fast-audit] running max|incremental - full| = "
                    << ecoi_fast_audit_max_ << "\n";
    }
}

void Chain::update_eps_pos(int iteration)
{
    ProfileScope scope("Chain::update_eps_pos");
    if (params.marginal_ecoi)
    {
        update_eps_pos_locus(iteration);
        return;
    }
    const std::size_t N = genotyping_data.num_samples;
    const std::size_t n_loci = genotyping_data.num_loci;

    // eps_pos enters only the observation model, so the move changes sample s's
    // observation row and its eps_pos prior; the transmission term is untouched
    // and cancels. Per-sample MH ratio is temp*(delta_obs_s) + delta_prior_s +
    // adj, all from s's own slots, so samples update independently in parallel.
    prewarm_caches_for_parallel_sample_move();
    in_parallel_sample_region_ = true;
    moire_parallel::recalc_parallel_for(0, N, [&](std::size_t sample_idx) {
        Sampler &rng = worker_rng();
        const auto prop_adj = rng.sample_constrained(
            eps_pos.at({sample_idx}), eps_pos_var.at({sample_idx}), min_sampled, 1);
        const float prop_eps_pos = std::get<0>(prop_adj);
        const float adj = std::get<1>(prop_adj);

        if (prop_eps_pos < 1 && prop_eps_pos > 1e-32) {
            const float old_obs = obs_row_sum_new_[sample_idx];
            const float old_prior = eps_pos_prior_new.at({sample_idx});

            const float prev_eps_pos = eps_pos.at({sample_idx});
            eps_pos.at({sample_idx}) = prop_eps_pos;
            calculate_eps_pos_likelihood(sample_idx);
            for (std::size_t locus_idx = 0; locus_idx < n_loci; ++locus_idx) {
                calculate_observation_likelihood(sample_idx, locus_idx);
            }
            sync_obs_sum_for_sample(sample_idx);

            const float new_obs = obs_row_sum_new_[sample_idx];
            const float mh_ratio =
                temp * (new_obs - old_obs) +
                (eps_pos_prior_new.at({sample_idx}) - old_prior) + adj;

            if (!std::isfinite(new_obs) ||
                rng.sample_log_mh_acceptance() > mh_ratio) {
                eps_pos.at({sample_idx}) = prev_eps_pos;
                restore_eps_pos_likelihood(sample_idx);
                for (std::size_t locus_idx = 0; locus_idx < n_loci; ++locus_idx) {
                    restore_observation_likelihood(sample_idx, locus_idx);
                }
                sync_obs_sum_for_sample(sample_idx);
            } else {
                save_eps_pos_likelihood(sample_idx);
                for (std::size_t locus_idx = 0; locus_idx < n_loci; ++locus_idx) {
                    save_observation_likelihood(sample_idx, locus_idx);
                }
                ++eps_pos_accept.at({sample_idx});
            }

            if (iteration < params.burnin && iteration > 15) {
                const float acceptanceRate =
                    eps_pos_accept.at({sample_idx}) / float(iteration);
                const float update =
                    (acceptanceRate - .23) / std::pow(iteration + 1, .5);
                eps_pos_var.at({sample_idx}) =
                    std::max(eps_pos_var.at({sample_idx}) + update, .0001f);
            }
        }
    });
    in_parallel_sample_region_ = false;
    reduce_likelihood_sums_after_parallel_sample_move();
}

void Chain::update_eps_pos_locus(int iteration)
{
    ProfileScope scope("Chain::update_eps_pos_locus");
    const std::size_t N = genotyping_data.num_samples;
    auto locus_indices = std::vector<std::size_t>(genotyping_data.num_loci);
    std::iota(locus_indices.begin(), locus_indices.end(), 0);
    sampler.shuffle_vec(locus_indices);

    for (const auto locus_idx : locus_indices)
    {
        const auto prop_adj = sampler.sample_constrained(
            eps_pos_locus.at({locus_idx}), eps_pos_locus_var.at({locus_idx}), min_sampled, 1);
        const float prop_eps_pos = std::get<0>(prop_adj);
        const float adj = std::get<1>(prop_adj);

        if (!(prop_eps_pos < 1) || !(prop_eps_pos > 1e-32))
        {
            continue;
        }

        const float prev_eps_pos = eps_pos_locus.at({locus_idx});
        eps_pos_locus.at({locus_idx}) = prop_eps_pos;
        calculate_eps_pos_locus_likelihood(locus_idx);

        // The false-positive rate enters only the observation model, so changing
        // it touches this locus' column across every sample; the marginalized
        // transmission term (ecoi_marg_sum_) is unaffected.
        moire_parallel::parallel_for(0, N, [&](std::size_t sample_idx) {
            calculate_observation_likelihood(sample_idx, locus_idx);
        });
        for (std::size_t sample_idx = 0; sample_idx < N; ++sample_idx)
        {
            sync_obs_sum_for_sample(sample_idx);
        }

        const float new_llik = calc_new_likelihood();
        const float new_prior = calc_new_prior();
        const float new_post = new_llik * temp + new_prior;

        if (!std::isfinite(new_post) or sampler.sample_log_mh_acceptance() >
                                        (new_post - get_posterior() + adj))
        {
            eps_pos_locus.at({locus_idx}) = prev_eps_pos;
            restore_eps_pos_locus_likelihood(locus_idx);
            for (std::size_t sample_idx = 0; sample_idx < N; ++sample_idx)
            {
                restore_observation_likelihood(sample_idx, locus_idx);
            }
            for (std::size_t sample_idx = 0; sample_idx < N; ++sample_idx)
            {
                sync_obs_sum_for_sample(sample_idx);
            }
        }
        else
        {
            llik = new_llik;
            prior = new_prior;
            save_eps_pos_locus_likelihood(locus_idx);
            for (std::size_t sample_idx = 0; sample_idx < N; ++sample_idx)
            {
                save_observation_likelihood(sample_idx, locus_idx);
            }
            ++eps_pos_locus_accept.at({locus_idx});
        }

        if (iteration < params.burnin and iteration > 15)
        {
            const float acceptanceRate =
                eps_pos_locus_accept.at({locus_idx}) / float(iteration);
            const float update =
                (acceptanceRate - .23) / std::pow(iteration + 1, .5);
            eps_pos_locus_var.at({locus_idx}) =
                std::max(eps_pos_locus_var.at({locus_idx}) + update, .0001f);
        }
    }
}

void Chain::update_eps_neg(int iteration)
{
    ProfileScope scope("Chain::update_eps_neg");
    const std::size_t N = genotyping_data.num_samples;
    const std::size_t n_loci = genotyping_data.num_loci;

    // eps_neg enters only the observation model, so the move changes sample s's
    // observation row and its eps_neg prior; the transmission term is untouched
    // and cancels. Per-sample MH ratio is temp*(delta_obs_s) + delta_prior_s +
    // adj, all from s's own slots, so samples update independently in parallel.
    prewarm_caches_for_parallel_sample_move();
    in_parallel_sample_region_ = true;
    moire_parallel::recalc_parallel_for(0, N, [&](std::size_t sample_idx) {
        Sampler &rng = worker_rng();
        const auto prop_adj = rng.sample_constrained(
            eps_neg.at({sample_idx}), eps_neg_var.at({sample_idx}), min_sampled, 1);
        const float prop_eps_neg = std::get<0>(prop_adj);
        const float adj = std::get<1>(prop_adj);

        if (prop_eps_neg < 1 && prop_eps_neg > 1e-32) {
            const float old_obs = obs_row_sum_new_[sample_idx];
            const float old_prior = eps_neg_prior_new.at({sample_idx});

            const float prev_eps_neg = eps_neg.at({sample_idx});
            eps_neg.at({sample_idx}) = prop_eps_neg;
            calculate_eps_neg_likelihood(sample_idx);
            for (std::size_t locus_idx = 0; locus_idx < n_loci; ++locus_idx) {
                calculate_observation_likelihood(sample_idx, locus_idx);
            }
            sync_obs_sum_for_sample(sample_idx);

            const float new_obs = obs_row_sum_new_[sample_idx];
            const float mh_ratio =
                temp * (new_obs - old_obs) +
                (eps_neg_prior_new.at({sample_idx}) - old_prior) + adj;

            if (!std::isfinite(new_obs) ||
                rng.sample_log_mh_acceptance() > mh_ratio) {
                eps_neg.at({sample_idx}) = prev_eps_neg;
                restore_eps_neg_likelihood(sample_idx);
                for (std::size_t locus_idx = 0; locus_idx < n_loci; ++locus_idx) {
                    restore_observation_likelihood(sample_idx, locus_idx);
                }
                sync_obs_sum_for_sample(sample_idx);
            } else {
                save_eps_neg_likelihood(sample_idx);
                for (std::size_t locus_idx = 0; locus_idx < n_loci; ++locus_idx) {
                    save_observation_likelihood(sample_idx, locus_idx);
                }
                ++eps_neg_accept.at({sample_idx});
            }

            if (iteration < params.burnin && iteration > 15) {
                const float acceptanceRate =
                    eps_neg_accept.at({sample_idx}) / float(iteration);
                const float update =
                    (acceptanceRate - .23) / std::pow(iteration + 1, .5);
                eps_neg_var.at({sample_idx}) =
                    std::max(eps_neg_var.at({sample_idx}) + update, .0001f);
            }
        }
    });
    in_parallel_sample_region_ = false;
    reduce_likelihood_sums_after_parallel_sample_move();
}

void Chain::update_samples(int iteration)
{
    ProfileScope scope("Chain::update_samples");
    const std::size_t N = genotyping_data.num_samples;
    const std::size_t n_loci = genotyping_data.num_loci;

    // Joint per-sample move over (eps_neg, eps_pos, r, m). Independent across
    // samples: the move changes only sample s's observation row, transmission
    // term (folds in s's COI prior) and s's eps_neg / eps_pos / relatedness
    // priors, so the MH ratio is
    //   temp * (delta_obs_s + delta_tx_s)
    //        + delta_eps_neg_prior_s + delta_eps_pos_prior_s
    //        + delta_relatedness_prior_s + adj,
    // from s's own slots. Parallel over samples with a private RNG; reconcile
    // the deferred global sums serially afterward.
    prewarm_caches_for_parallel_sample_move();
    in_parallel_sample_region_ = true;
    moire_parallel::recalc_parallel_for(0, N, [&](std::size_t sample_idx) {
        Sampler &rng = worker_rng();
        const auto eps_neg_prop_adj = rng.sample_constrained(
            eps_neg.at({sample_idx}), eps_neg_var.at({sample_idx}), min_sampled, 1);
        const float prop_eps_neg = std::get<0>(eps_neg_prop_adj);
        const float eps_neg_adj = std::get<1>(eps_neg_prop_adj);
        const bool valid_prop_eps_neg = prop_eps_neg < 1 && prop_eps_neg > 1e-32;

        const auto eps_pos_prop_adj = rng.sample_constrained(
            eps_pos.at({sample_idx}), eps_pos_var.at({sample_idx}), min_sampled, 1);
        const float prop_eps_pos = std::get<0>(eps_pos_prop_adj);
        const float eps_pos_adj = std::get<1>(eps_pos_prop_adj);
        const bool valid_prop_eps_pos = prop_eps_pos < 1 && prop_eps_pos > 1e-32;

        float prop_r = 0;
        float r_adj = 0;
        bool valid_prop_r = true;
        if (params.allow_relatedness) {
            auto r_prop_adj = rng.sample_constrained(
                r.at({sample_idx}), r_var.at({sample_idx}), min_sampled, .99);
            prop_r = std::get<0>(r_prop_adj);
            r_adj = std::get<1>(r_prop_adj);
            valid_prop_r = prop_r < 1 && prop_r > 1e-32;
        }

        const std::size_t prop_m = m.at({sample_idx}) + rng.sample_coi_delta(2);
        const bool valid_prop_m = prop_m > 0 && prop_m <= params.max_coi;

        if (!(valid_prop_eps_neg && valid_prop_eps_pos && valid_prop_r &&
              valid_prop_m)) {
            return;
        }

        const float old_lik =
            obs_row_sum_new_[sample_idx] + tx_sample_logsumexp_new.at({sample_idx});
        const float old_eps_neg_prior = eps_neg_prior_new.at({sample_idx});
        const float old_eps_pos_prior = eps_pos_prior_new.at({sample_idx});
        const float old_rel_prior = relatedness_prior_new.at({sample_idx});

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

        for (std::size_t locus_idx = 0; locus_idx < n_loci; ++locus_idx) {
            auto lg = observation_model_->sample_latent_genotype(
                rng, genotyping_data.get_observed_alleles(sample_idx, locus_idx),
                m.at({sample_idx}), eps_pos.at({sample_idx}),
                eps_neg.at({sample_idx}));
            assign_latent_genotype_new(sample_idx, locus_idx, lg.value);
            lg_adj_new.at({sample_idx, locus_idx}) = lg.log_prob;
            adj_ratio += lg_adj_new.at({sample_idx, locus_idx}) -
                         lg_adj_old.at({sample_idx, locus_idx});
            calculate_observation_likelihood(sample_idx, locus_idx);
        }
        sync_obs_sum_for_sample(sample_idx);
        recalculate_transmission_for_sample_incremental(
            sample_idx, static_cast<int>(prev_m), prev_r);

        const float new_lik =
            obs_row_sum_new_[sample_idx] + tx_sample_logsumexp_new.at({sample_idx});
        const float prior_delta =
            (eps_neg_prior_new.at({sample_idx}) - old_eps_neg_prior) +
            (eps_pos_prior_new.at({sample_idx}) - old_eps_pos_prior) +
            (relatedness_prior_new.at({sample_idx}) - old_rel_prior);
        const float mh_ratio =
            temp * (new_lik - old_lik) + prior_delta + adj_ratio;
        const float alpha = rng.sample_log_mh_acceptance();

        if (std::isfinite(new_lik) && alpha <= mh_ratio) {
            // Accept.
            save_eps_neg_likelihood(sample_idx);
            save_eps_pos_likelihood(sample_idx);
            save_relatedness_likelihood(sample_idx);
            save_coi_likelihood(sample_idx);
            for (std::size_t locus_idx = 0; locus_idx < n_loci; ++locus_idx) {
                save_observation_likelihood(sample_idx, locus_idx);
                lg_adj_old.at({sample_idx, locus_idx}) =
                    lg_adj_new.at({sample_idx, locus_idx});
                auto [begin, end] =
                    latent_genotypes_new.inner_iterators({sample_idx, locus_idx});
                std::copy(begin, end,
                          latent_genotypes_old.inner_begin({sample_idx, locus_idx}));
            }
            for (std::size_t pop_idx = 0; pop_idx < params.num_populations; ++pop_idx) {
                for (std::size_t locus_idx = 0; locus_idx < n_loci; ++locus_idx) {
                    save_transmission_likelihood(pop_idx, sample_idx, locus_idx);
                }
            }
            ++sample_accept.at({sample_idx});
        } else {
            // Reject.
            m.at({sample_idx}) = prev_m;
            eps_pos.at({sample_idx}) = prev_eps_pos;
            eps_neg.at({sample_idx}) = prev_eps_neg;
            r.at({sample_idx}) = prev_r;
            restore_eps_neg_likelihood(sample_idx);
            restore_eps_pos_likelihood(sample_idx);
            restore_relatedness_likelihood(sample_idx);
            restore_coi_likelihood(sample_idx);
            for (std::size_t locus_idx = 0; locus_idx < n_loci; ++locus_idx) {
                restore_observation_likelihood(sample_idx, locus_idx);
                lg_adj_new.at({sample_idx, locus_idx}) =
                    lg_adj_old.at({sample_idx, locus_idx});
                restore_latent_genotype_new(sample_idx, locus_idx);
            }
            sync_obs_sum_for_sample(sample_idx);
            restore_transmission_for_sample_incremental(sample_idx);
        }
    });
    in_parallel_sample_region_ = false;
    reduce_likelihood_sums_after_parallel_sample_move();
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
    float new_llik;
    double collapsed_sum = 0.0;
    if (params.marginal_ecoi) {
        // population_coi_p enters w(m|e); the whole marginal shifts.
        collapsed_sum = recompute_collapsed_marginals(ecoi_sample_marg_scratch_);
        new_llik = obs_llik_sum_new_ + static_cast<float>(collapsed_sum);
    } else {
        refresh_all_samples_tx_logsumexp();
        new_llik = obs_llik_sum_new_ + tx_llik_sum_new;
    }
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
        if (params.marginal_ecoi) {
            std::swap(ecoi_sample_marg_, ecoi_sample_marg_scratch_);
            ecoi_marg_sum_ = collapsed_sum;
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
        if (!params.marginal_ecoi) {
            refresh_all_samples_tx_logsumexp();
        }
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
    float new_llik;
    double collapsed_sum = 0.0;
    if (params.marginal_ecoi) {
        collapsed_sum = recompute_collapsed_marginals(ecoi_sample_marg_scratch_);
        new_llik = obs_llik_sum_new_ + static_cast<float>(collapsed_sum);
    } else {
        refresh_all_samples_tx_logsumexp();
        new_llik = obs_llik_sum_new_ + tx_llik_sum_new;
    }
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
        if (params.marginal_ecoi) {
            std::swap(ecoi_sample_marg_, ecoi_sample_marg_scratch_);
            ecoi_marg_sum_ = collapsed_sum;
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
        if (!params.marginal_ecoi) {
            refresh_all_samples_tx_logsumexp();
        }
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
        float new_llik;
        double collapsed_sum = 0.0;
        if (params.marginal_ecoi) {
            // population weights enter the per-sample population mixture in M.
            collapsed_sum = recompute_collapsed_marginals(ecoi_sample_marg_scratch_);
            new_llik = obs_llik_sum_new_ + static_cast<float>(collapsed_sum);
        } else {
            refresh_all_samples_tx_logsumexp();
            new_llik = obs_llik_sum_new_ + tx_llik_sum_new;
        }
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
            if (!params.marginal_ecoi) {
                refresh_all_samples_tx_logsumexp();
            }
        }
        else
        {
            llik = new_llik;
            prior = new_prior;
            save_population_responsibility_vector_likelihood();
            if (params.marginal_ecoi) {
                std::swap(ecoi_sample_marg_, ecoi_sample_marg_scratch_);
                ecoi_marg_sum_ = collapsed_sum;
            }
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
