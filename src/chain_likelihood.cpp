#include "chain.h"

#include "prob_any_missing.h"
#include "prob_any_missing_cache.h"
#include "transmission_incremental.h"
#include "pam_fast_paths.h"
#include "multivector_fused.h"
#include "parallel_backend.h"
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

const PamCachedVectors& compute_and_cache_pam_vector(
    probAnyMissingFunctor& functor,
    std::span<const float> q,
    unsigned min_events,
    unsigned max_events)
{
    thread_local std::vector<double> fast_scratch;
    if (pam_fast_paths::try_fill_pam_vector(q, min_events, max_events, fast_scratch)) {
        ProfileScope scope("Chain::pam_vec_fast_low_k");
        return pam_vector_cache().store_pam(
            q, min_events, max_events, std::span<const double>(fast_scratch));
    }

    ProfileScope scope("Chain::pam_vec_gray_code");
    thread_local std::vector<float> q_copy;
    q_copy.assign(q.begin(), q.end());
    thread_local std::vector<double> gray_scratch;
    gray_scratch = functor.vectorized(q_copy, min_events, max_events);
    return pam_vector_cache().store_pam(
        q, min_events, max_events, std::span<const double>(gray_scratch));
}

const PamCachedVectors& cached_pam_vector(
    probAnyMissingFunctor& functor,
    const std::vector<float>& event_probs,
    unsigned min_events,
    unsigned max_events)
{
    const auto q = std::span<const float>(event_probs);
    auto& cache = pam_vector_cache();
    if (const PamCachedVectors* hit = cache.lookup(q, min_events, max_events)) {
        ProfileScope scope("Chain::pam_vec_cache_hit");
        if (pam_cache::Config::instance().verify) {
            thread_local std::vector<double> exact;
            if (pam_fast_paths::try_fill_pam_vector(q, min_events, max_events, exact)) {
                pam_cache_verify_hit(q, min_events, max_events, *hit, exact);
            } else {
                exact = functor.vectorized(event_probs, min_events, max_events);
                pam_cache_verify_hit(q, min_events, max_events, *hit, exact);
            }
        }
        return *hit;
    }
    ProfileScope scope("Chain::pam_vec_cache_miss");
    ++pam_cache::stats().misses;
    return compute_and_cache_pam_vector(functor, q, min_events, max_events);
}

float finish_transmission_log(const PamCachedVectors& pam,
                              Sampler& sampler,
                              bool allow_relatedness,
                              int coi,
                              std::size_t total_alleles,
                              float relatedness,
                              float log_constrained_set_total_prob)
{
    if (!allow_relatedness) {
        return transmission_incremental::transmission_log_no_relatedness(
            std::span<const float>(pam.log_one_minus_pam), coi,
            log_constrained_set_total_prob);
    }
    ProfileScope scope_loop("Chain::calc_transmission_process::loop");
    return transmission_incremental::transmission_log_with_relatedness(
        std::span<const double>(pam.pam), std::span<const float>(pam.log_one_minus_pam),
        sampler, coi, total_alleles, relatedness, log_constrained_set_total_prob);
}

} // namespace

float Chain::calc_transmission_process(
    std::span<int const> full_allele_index_vec,
    std::span<float const> allele_frequencies, int coi, float relatedness)
{
    ProfileScope scope_all("Chain::calc_transmission_process::all");
    // transmission process - prob that after "coi" number of draws, all
    // alleles are drawn at least once conditional on all draws come
    // from the constrained set, where the constrained set is the set of
    // positive alleles in the latent genotype

    // get span up to first instance of -1
    auto first_instance = std::find(full_allele_index_vec.begin(), full_allele_index_vec.end(), -1);
    std::span<int const> allele_index_vec = full_allele_index_vec.subspan(0, first_instance - full_allele_index_vec.begin());

    if (allele_index_vec.size() > coi) 
    {
        // Invalid state: more alleles than COI is impossible
        static int log_count = 0;
        if (log_count < 10) {
            UtilFunctions::print("WARNING: calc_transmission_process: allele_index_vec.size() > coi:",
                                 "size =", allele_index_vec.size(), "coi =", coi);
            log_count++;
        }
        return -std::numeric_limits<float>::infinity();
    }

    float constrained_set_total_prob = 0;
    thread_local std::vector<float> prVec_;
    prVec_.clear();
    const size_t total_alleles = allele_index_vec.size();
    if (total_alleles == 0) {
        return -std::numeric_limits<float>::infinity();
    }

    prVec_.reserve(total_alleles);
    for (const auto &e : allele_index_vec)
    {
        // Bounds check: ensure e is valid index into allele_frequencies
        if (static_cast<size_t>(e) >= allele_frequencies.size() || e < 0) {
            // Invalid allele index - return -infinity instead of throwing
            // (throwing exceptions from parallel code causes issues with R's error handling)
            static int log_count = 0;
            if (log_count < 10) {
                UtilFunctions::print("WARNING: calc_transmission_process: invalid allele index:",
                                     "e =", e, "allele_frequencies.size() =", allele_frequencies.size());
                log_count++;
            }
            return -std::numeric_limits<float>::infinity();
        }
        prVec_.push_back(allele_frequencies[e]);
        constrained_set_total_prob += prVec_.back();
    }

    // Check for zero total probability (shouldn't happen but safety check)
    if (constrained_set_total_prob <= 0.0f || !std::isfinite(constrained_set_total_prob)) {
        static int log_count = 0;
        if (log_count < 10) {
            UtilFunctions::print("WARNING: calc_transmission_process: invalid constrained_set_total_prob:",
                                 "value =", constrained_set_total_prob, "total_alleles =", total_alleles);
            log_count++;
        }
        return -std::numeric_limits<float>::infinity();
    }

    const float log_constrained_set_total_prob = std::log(constrained_set_total_prob);

    // normalize the vector
    for (float &k : prVec_)
    {
        k = k / constrained_set_total_prob;
    }

    if (total_alleles == 1 && pam_fast_paths::tx_opts_enabled()) {
        if (!params.allow_relatedness) {
            return log_constrained_set_total_prob * static_cast<float>(coi);
        }
        return transmission_incremental::transmission_log_with_relatedness_k1(
            sampler, coi, relatedness, log_constrained_set_total_prob);
    }

    thread_local probAnyMissingFunctor local_probAnyMissing;
    const PamCachedVectors& pam =
        cached_pam_vector(local_probAnyMissing, prVec_, 1u, static_cast<unsigned>(coi));

    if (params.allow_relatedness) {
        return finish_transmission_log(
            pam, sampler, true, coi, total_alleles, relatedness,
            log_constrained_set_total_prob);
    }

    return transmission_incremental::transmission_log_no_relatedness(
        std::span<const float>(pam.log_one_minus_pam), coi, log_constrained_set_total_prob);
}

float Chain::calc_transmission_process_after_r_change(
    std::span<int const> full_allele_index_vec,
    std::span<float const> allele_frequencies,
    int coi,
    float relatedness)
{
    ProfileScope scope("Chain::calc_transmission_process::r_update");

    const auto first_instance =
        std::find(full_allele_index_vec.begin(), full_allele_index_vec.end(), -1);
    const std::span<int const> allele_index_vec =
        full_allele_index_vec.subspan(0, first_instance - full_allele_index_vec.begin());

    if (allele_index_vec.empty()) {
        return -std::numeric_limits<float>::infinity();
    }
    if (allele_index_vec.size() > static_cast<std::size_t>(coi)) {
        return -std::numeric_limits<float>::infinity();
    }

    if (!pam_fast_paths::tx_opts_enabled()) {
        return calc_transmission_process(
            full_allele_index_vec, allele_frequencies, coi, relatedness);
    }

    thread_local std::vector<float> q;
    float sum = 0.f;
    if (!transmission_incremental::build_constrained_q(
            allele_index_vec, allele_frequencies, q, sum)) {
        return -std::numeric_limits<float>::infinity();
    }

    const float log_sum = std::log(sum);
    if (allele_index_vec.size() == 1) {
        if (!params.allow_relatedness) {
            return log_sum * static_cast<float>(coi);
        }
        return transmission_incremental::transmission_log_with_relatedness_k1(
            sampler, coi, relatedness, log_sum);
    }

    ProfileScope scope_inc("Chain::calc_transmission_process::r_inc");
    thread_local probAnyMissingFunctor local_probAnyMissing;
    const PamCachedVectors& pam =
        cached_pam_vector(local_probAnyMissing, q, 1u, static_cast<unsigned>(coi));

    if (!params.allow_relatedness) {
        return transmission_incremental::transmission_log_no_relatedness(
            std::span<const float>(pam.log_one_minus_pam), coi, log_sum);
    }
    return transmission_incremental::transmission_log_with_relatedness(
        std::span<const double>(pam.pam), std::span<const float>(pam.log_one_minus_pam),
        sampler, coi, allele_index_vec.size(), relatedness, log_sum);
}

void Chain::invalidate_transmission_llik_cache() {
    tx_llik_cache_valid_ = false;
}

void Chain::rebuild_transmission_llik_cache() {
    if (!population_responsibility_vector_log_valid_) {
        population_responsibility_vector_log_ = population_responsibility_vector.log();
        population_responsibility_vector_log_valid_ = true;
    }

    const std::size_t num_samples = genotyping_data.num_samples;
    const std::size_t num_populations = params.num_populations;
    const std::size_t num_loci = genotyping_data.num_loci;

    tx_loci_sum_new.resize({num_samples, num_populations});
    tx_sample_logsumexp_new.resize({num_samples});

    std::vector<float> logits(num_populations);
    const auto pop_log = population_responsibility_vector_log_.as_span();
    const std::size_t tx_stride0 = transmission_llik_new.strides()[0];
    const std::size_t tx_stride1 = transmission_llik_new.strides()[1];
    const std::size_t coi_stride0 = coi_prior_new.strides()[0];
    const float* tx_data = transmission_llik_new.data().data();
    const float* coi_data = coi_prior_new.data().data();

    float total = 0.f;
    for (std::size_t s = 0; s < num_samples; ++s) {
        for (std::size_t p = 0; p < num_populations; ++p) {
            const std::size_t start_idx = s * tx_stride0 + p * tx_stride1;
            const float loci_sum = moire_fused::detail::sum_loci_slice<float>(
                tx_data, start_idx, num_loci);
            tx_loci_sum_new.at({s, p}) = loci_sum;
            logits[p] = loci_sum + coi_data[s * coi_stride0 + p] + pop_log[p];
        }
        const float sample_lse = moire_fused::detail::population_logits_logsumexp(
            logits.data(), num_populations);
        tx_sample_logsumexp_new.at({s}) = sample_lse;
        total += sample_lse;
    }
    tx_llik_sum_new = total;
    tx_llik_cache_valid_ = true;
}

float Chain::apply_transmission_cell_change(
    std::size_t sample_idx, std::size_t pop_idx, float old_val, float new_val)
{
    if (old_val == new_val) {
        return 0.f;
    }
    if (!tx_llik_cache_valid_) {
        rebuild_transmission_llik_cache();
    }
    if (!population_responsibility_vector_log_valid_) {
        population_responsibility_vector_log_ = population_responsibility_vector.log();
        population_responsibility_vector_log_valid_ = true;
    }

    const float pop_log = population_responsibility_vector_log_.as_span()[pop_idx];
    const float coi_term = coi_prior_new.at({sample_idx, pop_idx});
    const float old_loci_sum = tx_loci_sum_new.at({sample_idx, pop_idx});
    const float old_logit_k = old_loci_sum + coi_term + pop_log;
    const float old_L = tx_sample_logsumexp_new.at({sample_idx});

    const float new_loci_sum = old_loci_sum - old_val + new_val;
    tx_loci_sum_new.at({sample_idx, pop_idx}) = new_loci_sum;
    const float new_logit_k = new_loci_sum + coi_term + pop_log;

    float new_L;
    if (params.num_populations == 1) {
        new_L = new_logit_k;
    } else {
        new_L = moire_fused::logsumexp_replace_one(old_L, old_logit_k, new_logit_k);
    }

    tx_sample_logsumexp_new.at({sample_idx}) = new_L;
    const float delta = new_L - old_L;
    tx_llik_sum_new += delta;
    return delta;
}

float Chain::calc_transmission_llik_sum() {
    if (!tx_llik_cache_valid_) {
        rebuild_transmission_llik_cache();
    }
    return tx_llik_sum_new;
}

float Chain::calc_new_likelihood() {
    ProfileScope scope("Chain::calc_new_likelihood");
    // observation_llik_new sum is maintained incrementally (see
    // sync_obs_sum_for_sample); avoids an O(N*L) reduction per proposal.
    const float observation_llik = obs_llik_sum_new_;

    float transmission_llik;
    {
        ProfileScope s_tx("Chain::calc_new_likelihood::transmission_sum");
        transmission_llik = calc_transmission_llik_sum();
    }

    return observation_llik + transmission_llik;
}

float Chain::calc_new_prior() {
    ProfileScope scope("Chain::calc_new_prior");
    // Per-sample prior sums are maintained incrementally in calculate_*/restore_*
    // to avoid an O(N) reduction per proposal.
    return eps_neg_prior_sum_new_ +
                eps_pos_prior_sum_new_ +
                relatedness_prior_sum_new_ +
                population_coi_p_hyper_prior_new + 
                population_coi_r_hyper_prior_new + 
                population_responsibility_vector_prior_new;
}

float Chain::get_llik() { return llik; }
float Chain::get_prior() { return prior; }
float Chain::get_posterior() { return llik * temp + prior; }

void Chain::calculate_observation_likelihood(std::size_t sample_idx, std::size_t locus_idx) {
    ProfileScope scope("Chain::calculate_observation_likelihood");
    if (genotyping_data.is_missing(sample_idx, locus_idx)) {
        observation_llik_new.at({sample_idx, locus_idx}) = 0;
    } else {
        auto [latent_genotypes_begin, latent_genotypes_end] = latent_genotypes_new.inner_iterators({sample_idx, locus_idx});
        const float obs_prob = observation_model_->log_likelihood(
            std::span(latent_genotypes_begin, latent_genotypes_end),
            genotyping_data.get_observed_alleles(sample_idx, locus_idx),
            eps_neg.at({sample_idx}), eps_pos.at({sample_idx})
        );
        observation_llik_new.at({sample_idx, locus_idx}) = obs_prob;
    }
}

void Chain::sync_obs_sum_for_sample(std::size_t sample_idx)
{
    // Recompute the row sum of observation_llik_new for this sample and fold the
    // delta into the running total. Must be called sequentially (outside the
    // parallel-over-loci recompute/restore loops) to avoid races on the scalar.
    const auto [begin, end] = observation_llik_new.inner_iterators({sample_idx});
    float row = 0.f;
    for (auto it = begin; it != end; ++it) {
        row += *it;
    }
    obs_llik_sum_new_ += row - obs_row_sum_new_[sample_idx];
    obs_row_sum_new_[sample_idx] = row;
}

void Chain::calculate_transmission_likelihood(std::size_t population_idx, std::size_t sample_idx, std::size_t locus_idx) {
    ProfileScope scope("Chain::calculate_transmission_likelihood");
    if (genotyping_data.is_missing(sample_idx, locus_idx)) {
        transmission_llik_new.unchecked_at({sample_idx, population_idx, locus_idx}) = 0;
    } else {
        const auto [latent_genotypes_begin, latent_genotypes_end] = latent_genotypes_new.inner_iterators({sample_idx, locus_idx});
        const auto [p_begin, p_end] = p.inner_iterators({population_idx, locus_idx});
        const auto latent_genotypes = std::span(latent_genotypes_begin, latent_genotypes_end);
        const auto p_span = std::span(p_begin, p_end);
        const float transmission_prob = calc_transmission_process(
            latent_genotypes,
            p_span, m.at({sample_idx}), r.at({sample_idx}));
        transmission_llik_new.unchecked_at({sample_idx, population_idx, locus_idx}) = transmission_prob;
    }
}

void Chain::calculate_transmission_likelihood_after_r_change(
    std::size_t population_idx,
    std::size_t sample_idx,
    std::size_t locus_idx)
{
    ProfileScope scope("Chain::calculate_transmission_likelihood_after_r_change");
    if (genotyping_data.is_missing(sample_idx, locus_idx)) {
        transmission_llik_new.unchecked_at({sample_idx, population_idx, locus_idx}) = 0;
    } else {
        const auto [latent_genotypes_begin, latent_genotypes_end] =
            latent_genotypes_new.inner_iterators({sample_idx, locus_idx});
        const auto [p_begin, p_end] = p.inner_iterators({population_idx, locus_idx});
        const auto latent_genotypes =
            std::span(latent_genotypes_begin, latent_genotypes_end);
        const auto p_span = std::span<const float>(p_begin, p_end);
        const float transmission_prob = calc_transmission_process_after_r_change(
            latent_genotypes, p_span, m.at({sample_idx}), r.at({sample_idx}));
        transmission_llik_new.unchecked_at({sample_idx, population_idx, locus_idx}) =
            transmission_prob;
    }
}

void Chain::calculate_eps_neg_likelihood(std::size_t sample_idx)
{
    ProfileScope scope("Chain::calculate_eps_neg_likelihood");
    const float val = sampler.get_beta_log_prior(
        eps_neg.at({sample_idx}), params.eps_neg_alpha, params.eps_neg_beta);
    eps_neg_prior_sum_new_ += val - eps_neg_prior_new.at({sample_idx});
    eps_neg_prior_new.at({sample_idx}) = val;
}

void Chain::calculate_eps_pos_likelihood(std::size_t sample_idx)
{
    ProfileScope scope("Chain::calculate_eps_pos_likelihood");
    const float val = sampler.get_beta_log_prior(
        eps_pos.at({sample_idx}), params.eps_pos_alpha, params.eps_pos_beta);
    eps_pos_prior_sum_new_ += val - eps_pos_prior_new.at({sample_idx});
    eps_pos_prior_new.at({sample_idx}) = val;
}

void Chain::calculate_relatedness_likelihood(std::size_t sample_idx)
{
    ProfileScope scope("Chain::calculate_relatedness_likelihood");
    const float val = sampler.get_relatedness_log_prior(
        r.at({sample_idx}), params.r_alpha, params.r_beta);
    relatedness_prior_sum_new_ += val - relatedness_prior_new.at({sample_idx});
    relatedness_prior_new.at({sample_idx}) = val;
}

void Chain::calculate_coi_likelihood(std::size_t sample_idx)
{
    ProfileScope scope("Chain::calculate_coi_likelihood");
    coi_prior_new.at({sample_idx}) =
        sampler.get_coi_log_prior(m.at({sample_idx}), population_coi_p, population_coi_r);
}

void Chain::calculate_population_coi_p_likelihood()
{
    ProfileScope scope("Chain::calculate_population_coi_p_likelihood");
    population_coi_p_hyper_prior_new = sampler.get_beta_log_prior(
        population_coi_p, params.population_coi_p_alpha, params.population_coi_p_beta);
}

void Chain::calculate_population_coi_r_likelihood()
{
    ProfileScope scope("Chain::calculate_population_coi_r_likelihood");
    population_coi_r_hyper_prior_new = sampler.get_gamma_log_prior(
        population_coi_r, params.population_coi_r_shape, params.population_coi_r_rate);
}

void Chain::calculate_population_responsibility_vector_likelihood()
{
    ProfileScope scope("Chain::calculate_population_responsibility_vector_likelihood");
    population_responsibility_vector_prior_new = sampler.unnormalized_dirichlet_log_prior(
        population_responsibility_vector, params.population_responsibility_vector_alpha
    );
}

void Chain::save_population_responsibility_vector_likelihood()
{
    population_responsibility_vector_prior_old = population_responsibility_vector_prior_new;
}

void Chain::initialize_likelihood()
{
    ProfileScope scope("Chain::initialize_likelihood");
    const std::size_t num_samples = genotyping_data.num_samples;
    const std::size_t num_loci = genotyping_data.num_loci;
    const std::size_t num_populations = params.num_populations;

    observation_llik_new.resize({num_samples, num_loci});
    observation_llik_old.resize({num_samples, num_loci});

    transmission_llik_new.resize({num_samples, num_populations, num_loci});
    transmission_llik_old.resize({num_samples, num_populations, num_loci});

    eps_neg_prior_new.resize({num_samples});
    eps_neg_prior_old.resize({num_samples});

    eps_pos_prior_new.resize({num_samples});
    eps_pos_prior_old.resize({num_samples});

    relatedness_prior_new.resize({num_samples});
    relatedness_prior_old.resize({num_samples});
    
    coi_prior_new.resize({num_samples, num_populations});
    coi_prior_old.resize({num_samples, num_populations});

    // population_coi_mean_hyper_prior resize (optional formulation, kept for reference)
    // population_coi_mean_hyper_prior_new.resize({num_populations});
    // population_coi_mean_hyper_prior_old.resize({num_populations});

    // Calculate likelihoods of sample level parameters
    for (std::size_t sample_idx = 0; sample_idx < genotyping_data.num_samples; ++sample_idx) {
        calculate_eps_neg_likelihood(sample_idx);
        save_eps_neg_likelihood(sample_idx);
        
        calculate_eps_pos_likelihood(sample_idx);
        save_eps_pos_likelihood(sample_idx);
        
        calculate_relatedness_likelihood(sample_idx);
        save_relatedness_likelihood(sample_idx);
    }
    // Observation likelihoods: independent across (sample, locus), parallelize
    moire_parallel::parallel_for_2d(0, genotyping_data.num_samples, 0, genotyping_data.num_loci,
        [&](std::size_t sample_idx, std::size_t locus_idx) {
            calculate_observation_likelihood(sample_idx, locus_idx);
            save_observation_likelihood(sample_idx, locus_idx);
        });

    // Calculate likelihoods of sample specific population varying parameters
    for (std::size_t sample_idx = 0; sample_idx < genotyping_data.num_samples; ++sample_idx)
    {
        calculate_coi_likelihood(sample_idx);
        save_coi_likelihood(sample_idx);
        for (std::size_t pop_idx = 0; pop_idx < params.num_populations; ++pop_idx) {
            for (std::size_t locus_idx = 0; locus_idx < genotyping_data.num_loci; ++locus_idx)
            {
                calculate_transmission_likelihood(pop_idx, sample_idx, locus_idx);
                save_transmission_likelihood(pop_idx, sample_idx, locus_idx);
            }
        }
    }

    // Calculate likelihoods of population mean coi
    calculate_population_coi_p_likelihood();
    save_population_coi_p_likelihood();
    calculate_population_coi_r_likelihood();
    save_population_coi_r_likelihood();

    // Calculate likelihood of population responsibility vector
    calculate_population_responsibility_vector_likelihood();
    save_population_responsibility_vector_likelihood();

    // Rebuild incremental running sums from scratch. Done here (rather than
    // trusting the incremental updates above) so repeated initialize_parameters()
    // retries always start from a consistent state.
    eps_neg_prior_sum_new_ = eps_neg_prior_new.full_sum();
    eps_pos_prior_sum_new_ = eps_pos_prior_new.full_sum();
    relatedness_prior_sum_new_ = relatedness_prior_new.full_sum();

    obs_row_sum_new_.assign(num_samples, 0.f);
    float obs_total = 0.f;
    for (std::size_t sample_idx = 0; sample_idx < num_samples; ++sample_idx) {
        const auto [begin, end] = observation_llik_new.inner_iterators({sample_idx});
        float row = 0.f;
        for (auto it = begin; it != end; ++it) {
            row += *it;
        }
        obs_row_sum_new_[sample_idx] = row;
        obs_total += row;
    }
    obs_llik_sum_new_ = obs_total;

    llik = calc_new_likelihood();
    prior = calc_new_prior();
}

void Chain::save_observation_likelihood(std::size_t sample_idx, std::size_t locus_idx)
{
    observation_llik_old.at({sample_idx, locus_idx}) = observation_llik_new.at({sample_idx, locus_idx});
}

void Chain::save_transmission_likelihood(std::size_t population_idx, std::size_t sample_idx, std::size_t locus_idx)
{
    transmission_llik_old.unchecked_at({sample_idx, population_idx, locus_idx}) =
        transmission_llik_new.unchecked_at({sample_idx, population_idx, locus_idx});
}

void Chain::save_eps_neg_likelihood(std::size_t sample_idx)
{
    eps_neg_prior_old.at({sample_idx}) = eps_neg_prior_new.at({sample_idx});
}

void Chain::save_eps_pos_likelihood(std::size_t sample_idx)
{
    eps_pos_prior_old.at({sample_idx}) = eps_pos_prior_new.at({sample_idx});
}

void Chain::save_relatedness_likelihood(std::size_t sample_idx)
{
    relatedness_prior_old.at({sample_idx}) = relatedness_prior_new.at({sample_idx});
}

void Chain::save_coi_likelihood(std::size_t sample_idx)
{
    coi_prior_old.at({sample_idx}) = coi_prior_new.at({sample_idx});
}

void Chain::save_population_coi_p_likelihood()
{
    population_coi_p_hyper_prior_old = population_coi_p_hyper_prior_new;
}

void Chain::save_population_coi_r_likelihood()
{
    population_coi_r_hyper_prior_old = population_coi_r_hyper_prior_new;
}

void Chain::restore_observation_likelihood(std::size_t sample_idx, std::size_t locus_idx)
{
    observation_llik_new.at({sample_idx, locus_idx}) = observation_llik_old.at({sample_idx, locus_idx});
}

void Chain::restore_transmission_likelihood(std::size_t sample_idx, std::size_t population_idx, std::size_t locus_idx)
{
    transmission_llik_new.unchecked_at({sample_idx, population_idx, locus_idx}) =
        transmission_llik_old.unchecked_at({sample_idx, population_idx, locus_idx});
}

void Chain::restore_eps_neg_likelihood(std::size_t sample_idx)
{
    const float old_val = eps_neg_prior_old.at({sample_idx});
    eps_neg_prior_sum_new_ += old_val - eps_neg_prior_new.at({sample_idx});
    eps_neg_prior_new.at({sample_idx}) = old_val;
}

void Chain::restore_eps_pos_likelihood(std::size_t sample_idx)
{
    const float old_val = eps_pos_prior_old.at({sample_idx});
    eps_pos_prior_sum_new_ += old_val - eps_pos_prior_new.at({sample_idx});
    eps_pos_prior_new.at({sample_idx}) = old_val;
}

void Chain::restore_relatedness_likelihood(std::size_t sample_idx)
{
    const float old_val = relatedness_prior_old.at({sample_idx});
    relatedness_prior_sum_new_ += old_val - relatedness_prior_new.at({sample_idx});
    relatedness_prior_new.at({sample_idx}) = old_val;
}

void Chain::restore_coi_likelihood(std::size_t sample_idx)
{
    coi_prior_new.at({sample_idx}) = coi_prior_old.at({sample_idx});
}

void Chain::restore_population_coi_p_likelihood()
{
    population_coi_p_hyper_prior_new = population_coi_p_hyper_prior_old;
}

void Chain::restore_population_coi_r_likelihood()
{
    population_coi_r_hyper_prior_new = population_coi_r_hyper_prior_old;
}

void Chain::restore_population_responsibility_vector_likelihood()
{
    population_responsibility_vector_prior_new = population_responsibility_vector_prior_old;
}
