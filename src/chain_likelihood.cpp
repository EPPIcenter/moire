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

const PamCachedVectors& compute_pam_inline_low_k(std::span<const float> q,
                                               unsigned min_events,
                                               unsigned max_events)
{
    ProfileScope scope("Chain::pam_vec_fast_low_k_inline");
    thread_local PamCachedVectors scratch;
    thread_local std::vector<double> pam_buf;
    pam_fast_paths::fill_pam_vector_low_k(q, min_events, max_events, pam_buf);
    scratch.pam.assign(pam_buf.begin(), pam_buf.end());
    fill_log_one_minus_pam(scratch.pam, scratch.log_one_minus_pam);
    return scratch;
}

const PamCachedVectors& compute_and_cache_pam_vector(
    probAnyMissingFunctor& functor,
    std::span<const float> q,
    unsigned min_events,
    unsigned max_events)
{
    ProfileScope scope("Chain::pam_vec_gray_code");
    thread_local std::vector<float> q_copy;
    q_copy.assign(q.begin(), q.end());
    thread_local std::vector<double> gray_scratch;
    gray_scratch = functor.vectorized(q_copy, min_events, max_events);
    return pam_vector_cache().store_pam(
        q, min_events, max_events, std::span<const double>(gray_scratch));
}

const PamCachedVectors& compute_pam_from_vector(std::span<const double> pam_vec)
{
    thread_local PamCachedVectors scratch;
    scratch.pam.assign(pam_vec.begin(), pam_vec.end());
    fill_log_one_minus_pam(scratch.pam, scratch.log_one_minus_pam);
    return scratch;
}

bool q_spans_equal(std::span<const float> a, std::span<const float> b) noexcept
{
    return a.size() == b.size() &&
           std::equal(a.begin(), a.end(), b.begin());
}

const PamCachedVectors& cached_pam_vector(
    probAnyMissingFunctor& functor,
    std::span<const float> q,
    unsigned min_events,
    unsigned max_events,
    unsigned prev_max_events = 0,
    std::span<const float> q_prev = {})
{
    auto& cache = pam_vector_cache();
    if (prev_max_events > 0 && prev_max_events != max_events &&
        pam_fast_paths::tx_opts_enabled() &&
        q.size() <= pam_fast_paths::kLowKMaxSupport) {
        const int delta =
            static_cast<int>(max_events) - static_cast<int>(prev_max_events);
        if (delta > 0 && delta <= 2) {
            thread_local std::vector<double> prev_pam;
            pam_fast_paths::fill_pam_vector_low_k(
                q, min_events, prev_max_events, prev_pam);
            thread_local std::vector<double> extended;
            if (pam_fast_paths::extend_pam_vector_low_k(
                    q, min_events, prev_max_events, max_events, prev_pam, extended)) {
                ProfileScope scope("Chain::pam_vec_coi_extend");
                ++pam_cache::stats().coi_extend;
                return compute_pam_from_vector(
                    std::span<const double>(extended.data(), extended.size()));
            }
        }
    }

    if (!q_prev.empty() && q_prev.size() == q.size() &&
        pam_fast_paths::tx_opts_enabled() &&
        q.size() <= pam_fast_paths::kLowKMaxSupport) {
        if (q_spans_equal(q, q_prev)) {
            if (const PamCachedVectors* hit = cache.lookup(q, min_events, max_events)) {
                ProfileScope scope("Chain::pam_vec_p_q_unchanged");
                ++pam_cache::stats().p_q_unchanged;
                return *hit;
            }
        } else {
            thread_local std::vector<double> pam_buf;
            if (pam_fast_paths::try_fill_pam_vector_from_q_change(
                    q_prev, q, min_events, max_events, pam_buf)) {
                ProfileScope scope("Chain::pam_vec_p_q_one_step");
                ++pam_cache::stats().p_q_one_step;
                return compute_pam_from_vector(
                    std::span<const double>(pam_buf.data(), pam_buf.size()));
            }
        }
    }

    if (const PamCachedVectors* hit = cache.lookup(q, min_events, max_events)) {
        ProfileScope scope("Chain::pam_vec_cache_hit");
        if (pam_cache::Config::instance().verify) {
            thread_local std::vector<double> exact;
            if (pam_fast_paths::try_fill_pam_vector(q, min_events, max_events, exact)) {
                pam_cache_verify_hit(q, min_events, max_events, *hit, exact);
            } else {
                thread_local std::vector<float> q_copy;
                q_copy.assign(q.begin(), q.end());
                exact = functor.vectorized(q_copy, min_events, max_events);
                pam_cache_verify_hit(q, min_events, max_events, *hit, exact);
            }
        }
        return *hit;
    }

    // k <= kLowKMaxSupport: inclusion-exclusion is cheaper than cache store on miss.
    if (pam_fast_paths::tx_opts_enabled() &&
        q.size() <= pam_fast_paths::kLowKMaxSupport) {
        ++pam_cache::stats().low_k_inline;
        return compute_pam_inline_low_k(q, min_events, max_events);
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

float finish_transmission_from_group(const PamCachedVectors* pam,
                                     Sampler& sampler,
                                     bool allow_relatedness,
                                     int coi,
                                     std::size_t total_alleles,
                                     float relatedness,
                                     float log_sum,
                                     bool pam_valid)
{
    constexpr float kNegInf = -std::numeric_limits<float>::infinity();
    if (total_alleles == 0 || !std::isfinite(log_sum)) {
        return kNegInf;
    }
    if (total_alleles > static_cast<std::size_t>(coi)) {
        return kNegInf;
    }
    if (total_alleles == 1 && pam_fast_paths::tx_opts_enabled()) {
        if (!allow_relatedness) {
            return log_sum * static_cast<float>(coi);
        }
        return transmission_incremental::transmission_log_with_relatedness_k1(
            sampler, coi, relatedness, log_sum);
    }
    if (!pam_valid || pam == nullptr) {
        return kNegInf;
    }
    if (!allow_relatedness) {
        return transmission_incremental::transmission_log_no_relatedness(
            std::span<const float>(pam->log_one_minus_pam), coi, log_sum);
    }
    return transmission_incremental::transmission_log_with_relatedness(
        std::span<const double>(pam->pam), std::span<const float>(pam->log_one_minus_pam),
        sampler, coi, total_alleles, relatedness, log_sum);
}

void copy_pam_cached(const PamCachedVectors& src, PamCachedVectors& dst)
{
    dst.pam = src.pam;
    dst.log_one_minus_pam = src.log_one_minus_pam;
}

bool support_spans_equal(std::span<const int> a, std::span<const int> b) noexcept
{
    return a.size() == b.size() &&
           std::equal(a.begin(), a.end(), b.begin());
}

struct PChangePamGroup {
    std::vector<int> support;
    int coi{0};
    PamCachedVectors pam;
    float log_sum{0.f};
    std::size_t total_alleles{0};
    bool pam_valid{false};
};

} // namespace

void Chain::refresh_latent_support_k(std::size_t sample_idx, std::size_t locus_idx)
{
    const auto [begin, end] = latent_genotypes_new.inner_iterators({sample_idx, locus_idx});
    std::size_t k = 0;
    const std::size_t cap = static_cast<std::size_t>(end - begin);
    while (k < cap && begin[k] != -1) {
        ++k;
    }
    latent_support_k_.at({sample_idx, locus_idx}) = static_cast<std::uint16_t>(k);
}

void Chain::assign_latent_genotype_new(std::size_t sample_idx, std::size_t locus_idx,
                                       std::span<const int> value)
{
    latent_genotypes_new.inner_fill({sample_idx, locus_idx}, -1);
    latent_genotypes_new.inner_fill({sample_idx, locus_idx}, value);
    refresh_latent_support_k(sample_idx, locus_idx);
}

void Chain::restore_latent_genotype_new(std::size_t sample_idx, std::size_t locus_idx)
{
    const auto [begin, end] =
        latent_genotypes_old.inner_iterators({sample_idx, locus_idx});
    std::copy(begin, end, latent_genotypes_new.inner_begin({sample_idx, locus_idx}));
    refresh_latent_support_k(sample_idx, locus_idx);
}

std::span<const int> Chain::latent_allele_support(std::size_t sample_idx,
                                                  std::size_t locus_idx) const
{
    const auto [begin, end] = latent_genotypes_new.inner_iterators({sample_idx, locus_idx});
    const std::size_t k = latent_support_k_.at({sample_idx, locus_idx});
    (void)end;
    return std::span<const int>(begin, k);
}

std::span<const int> Chain::latent_allele_support_old(std::size_t sample_idx,
                                                      std::size_t locus_idx) const
{
    const auto [begin, end] = latent_genotypes_old.inner_iterators({sample_idx, locus_idx});
    std::size_t k = 0;
    const std::size_t cap = static_cast<std::size_t>(end - begin);
    while (k < cap && begin[k] != -1) {
        ++k;
    }
    return std::span<const int>(begin, k);
}

bool Chain::locus_tx_inputs_unchanged(std::size_t sample_idx,
                                      std::size_t locus_idx,
                                      int prev_coi,
                                      float prev_r) const
{
    if (m.at({sample_idx}) != prev_coi) {
        return false;
    }
    if (r.at({sample_idx}) != prev_r) {
        return false;
    }
    return support_spans_equal(latent_allele_support(sample_idx, locus_idx),
                               latent_allele_support_old(sample_idx, locus_idx));
}

float Chain::calc_transmission_process(
    std::span<int const> allele_index_vec,
    std::span<float const> allele_frequencies, int coi, float relatedness)
{
    ProfileScope scope_all("Chain::calc_transmission_process::all");
    // transmission process - prob that after "coi" number of draws, all
    // alleles are drawn at least once conditional on all draws come
    // from the constrained set, where the constrained set is the set of
    // positive alleles in the latent genotype

    if (allele_index_vec.size() > static_cast<std::size_t>(coi))
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
    std::span<int const> allele_index_vec,
    std::span<float const> allele_frequencies,
    int coi,
    float relatedness)
{
    ProfileScope scope("Chain::calc_transmission_process::r_update");

    if (allele_index_vec.empty()) {
        return -std::numeric_limits<float>::infinity();
    }
    if (allele_index_vec.size() > static_cast<std::size_t>(coi)) {
        return -std::numeric_limits<float>::infinity();
    }

    if (!pam_fast_paths::tx_opts_enabled()) {
        return calc_transmission_process(
            allele_index_vec, allele_frequencies, coi, relatedness);
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

float Chain::calc_transmission_process_after_p_change(
    std::span<int const> allele_index_vec,
    std::span<float const> allele_frequencies,
    int coi,
    float relatedness)
{
    ProfileScope scope("Chain::calc_transmission_process::p_update");

    if (allele_index_vec.empty()) {
        return -std::numeric_limits<float>::infinity();
    }
    if (allele_index_vec.size() > static_cast<std::size_t>(coi)) {
        return -std::numeric_limits<float>::infinity();
    }

    if (!pam_fast_paths::tx_opts_enabled()) {
        return calc_transmission_process(
            allele_index_vec, allele_frequencies, coi, relatedness);
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

    ProfileScope scope_inc("Chain::calc_transmission_process::p_inc");
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

void Chain::refresh_sample_tx_after_coi_change(std::size_t sample_idx)
{
    if (!tx_llik_cache_valid_) {
        rebuild_transmission_llik_cache();
        return;
    }
    if (!population_responsibility_vector_log_valid_) {
        population_responsibility_vector_log_ = population_responsibility_vector.log();
        population_responsibility_vector_log_valid_ = true;
    }

    const float old_L = tx_sample_logsumexp_new.at({sample_idx});
    const std::size_t num_populations = params.num_populations;
    thread_local std::vector<float> logits;
    if (logits.size() < num_populations) {
        logits.resize(num_populations);
    }
    const auto pop_log = population_responsibility_vector_log_.as_span();
    for (std::size_t p = 0; p < num_populations; ++p) {
        logits[p] = tx_loci_sum_new.at({sample_idx, p}) +
                    coi_prior_new.at({sample_idx, p}) + pop_log[p];
    }
    const float new_L = moire_fused::detail::population_logits_logsumexp(
        logits.data(), num_populations);
    tx_sample_logsumexp_new.at({sample_idx}) = new_L;
    tx_llik_sum_new += new_L - old_L;
}

void Chain::refresh_all_samples_tx_logsumexp()
{
    ProfileScope scope("Chain::refresh_all_samples_tx_logsumexp");
    if (!tx_llik_cache_valid_) {
        rebuild_transmission_llik_cache();
        return;
    }
    const std::size_t num_samples = genotyping_data.num_samples;
    for (std::size_t sample_idx = 0; sample_idx < num_samples; ++sample_idx) {
        refresh_sample_tx_after_coi_change(sample_idx);
    }
}

void Chain::recalculate_transmission_for_sample_incremental(std::size_t sample_idx,
                                                          int prev_coi,
                                                          float prev_r)
{
    ProfileScope scope("Chain::recalculate_transmission_for_sample");
    if (!tx_llik_cache_valid_) {
        rebuild_transmission_llik_cache();
    }
    refresh_sample_tx_after_coi_change(sample_idx);

    const int coi = m.at({sample_idx});
    const float relatedness = r.at({sample_idx});
    const std::size_t n_pops = params.num_populations;
    const std::size_t n_loci = genotyping_data.num_loci;
    constexpr float kNegInf = -std::numeric_limits<float>::infinity();

    // Group loci by (support, coi) so support indices are resolved once per group.
    std::vector<PChangePamGroup> support_groups;
    std::vector<int> locus_support_group(n_loci, -1);
    for (std::size_t locus_idx = 0; locus_idx < n_loci; ++locus_idx) {
        if (genotyping_data.is_missing(sample_idx, locus_idx)) {
            continue;
        }
        const auto support = latent_allele_support(sample_idx, locus_idx);
        int group_idx = -1;
        for (std::size_t g = 0; g < support_groups.size(); ++g) {
            if (support_groups[g].coi == coi &&
                support_spans_equal(support, support_groups[g].support)) {
                group_idx = static_cast<int>(g);
                break;
            }
        }
        if (group_idx < 0) {
            PChangePamGroup group;
            group.support.assign(support.begin(), support.end());
            group.coi = coi;
            group.total_alleles = group.support.size();
            support_groups.push_back(std::move(group));
            group_idx = static_cast<int>(support_groups.size() - 1);
        }
        locus_support_group[locus_idx] = group_idx;
    }

    struct LocusPopPamGroup {
        std::vector<float> q;
        float log_sum{0.f};
        PamCachedVectors pam;
        bool pam_valid{false};
        std::vector<std::size_t> pops;
    };

    moire_parallel::recalc_parallel_for(0, n_loci, [&](std::size_t locus_idx) {
        if (genotyping_data.is_missing(sample_idx, locus_idx)) {
            for (std::size_t pop_idx = 0; pop_idx < n_pops; ++pop_idx) {
                transmission_llik_new.unchecked_at({sample_idx, pop_idx, locus_idx}) = 0.f;
            }
            return;
        }

        if (locus_tx_inputs_unchanged(sample_idx, locus_idx, prev_coi, prev_r)) {
            return;
        }

        const int support_group_idx = locus_support_group[locus_idx];
        const PChangePamGroup& support_group =
            support_groups[static_cast<std::size_t>(support_group_idx)];
        const auto support = std::span<const int>(support_group.support);
        const std::size_t total_alleles = support_group.total_alleles;
        if (total_alleles == 0 || total_alleles > static_cast<std::size_t>(coi)) {
            for (std::size_t pop_idx = 0; pop_idx < n_pops; ++pop_idx) {
                transmission_llik_new.unchecked_at({sample_idx, pop_idx, locus_idx}) = kNegInf;
            }
            return;
        }

        if (!pam_fast_paths::tx_opts_enabled()) {
            for (std::size_t pop_idx = 0; pop_idx < n_pops; ++pop_idx) {
                calculate_transmission_likelihood(pop_idx, sample_idx, locus_idx);
            }
            return;
        }

        thread_local probAnyMissingFunctor functor;
        thread_local std::vector<float> q_scratch;
        thread_local std::vector<LocusPopPamGroup> pop_groups;
        pop_groups.clear();

        for (std::size_t pop_idx = 0; pop_idx < n_pops; ++pop_idx) {
            const auto [p_begin, p_end] = p.inner_iterators({pop_idx, locus_idx});
            const std::span<const float> p_span(p_begin, p_end);

            q_scratch.clear();
            float sum = 0.f;
            if (!transmission_incremental::build_constrained_q(
                    support, p_span, q_scratch, sum)) {
                transmission_llik_new.unchecked_at({sample_idx, pop_idx, locus_idx}) =
                    kNegInf;
                continue;
            }

            int group_idx = -1;
            for (std::size_t g = 0; g < pop_groups.size(); ++g) {
                if (q_spans_equal(q_scratch, pop_groups[g].q)) {
                    group_idx = static_cast<int>(g);
                    break;
                }
            }
            if (group_idx < 0) {
                LocusPopPamGroup group;
                group.q = q_scratch;
                group.log_sum = std::log(sum);
                pop_groups.push_back(std::move(group));
                group_idx = static_cast<int>(pop_groups.size() - 1);
            }
            pop_groups[static_cast<std::size_t>(group_idx)].pops.push_back(pop_idx);
        }

        {
            ProfileScope scope_groups("Chain::sample_tx::pam_groups");
            for (auto& group : pop_groups) {
                if (total_alleles == 1) {
                    group.pam_valid = false;
                    continue;
                }
                const unsigned prev_coi_u =
                    (coi != prev_coi) ? static_cast<unsigned>(prev_coi) : 0u;
                const PamCachedVectors& pam_ref = cached_pam_vector(
                    functor, group.q, 1u, static_cast<unsigned>(coi), prev_coi_u);
                copy_pam_cached(pam_ref, group.pam);
                group.pam_valid = true;
            }
        }

        for (const auto& group : pop_groups) {
            const float tx = finish_transmission_from_group(
                group.pam_valid ? &group.pam : nullptr,
                sampler,
                params.allow_relatedness,
                coi,
                total_alleles,
                relatedness,
                group.log_sum,
                group.pam_valid);
            for (const std::size_t pop_idx : group.pops) {
                transmission_llik_new.unchecked_at({sample_idx, pop_idx, locus_idx}) = tx;
            }
        }
    });

    // Fold cell deltas sequentially — tx_sample_logsumexp is shared per sample.
    for (std::size_t pop_idx = 0; pop_idx < n_pops; ++pop_idx) {
        for (std::size_t locus_idx = 0; locus_idx < n_loci; ++locus_idx) {
            const float old_cell =
                transmission_llik_old.unchecked_at({sample_idx, pop_idx, locus_idx});
            const float new_cell =
                transmission_llik_new.unchecked_at({sample_idx, pop_idx, locus_idx});
            apply_transmission_cell_change(sample_idx, pop_idx, old_cell, new_cell);
        }
    }
}

void Chain::restore_transmission_for_sample_incremental(std::size_t sample_idx)
{
    moire_parallel::parallel_for_2d(
        0, params.num_populations, 0, genotyping_data.num_loci,
        [&](std::size_t pop_idx, std::size_t locus_idx) {
            const float proposed =
                transmission_llik_new.unchecked_at({sample_idx, pop_idx, locus_idx});
            const float old_cell =
                transmission_llik_old.unchecked_at({sample_idx, pop_idx, locus_idx});
            apply_transmission_cell_change(sample_idx, pop_idx, proposed, old_cell);
            transmission_llik_new.unchecked_at({sample_idx, pop_idx, locus_idx}) = old_cell;
        });
    refresh_sample_tx_after_coi_change(sample_idx);
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
        const auto [p_begin, p_end] = p.inner_iterators({population_idx, locus_idx});
        const auto allele_support = latent_allele_support(sample_idx, locus_idx);
        const auto p_span = std::span(p_begin, p_end);
        const float transmission_prob = calc_transmission_process(
            allele_support,
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
        const auto [p_begin, p_end] = p.inner_iterators({population_idx, locus_idx});
        const auto allele_support = latent_allele_support(sample_idx, locus_idx);
        const auto p_span = std::span<const float>(p_begin, p_end);
        const float transmission_prob = calc_transmission_process_after_r_change(
            allele_support, p_span, m.at({sample_idx}), r.at({sample_idx}));
        transmission_llik_new.unchecked_at({sample_idx, population_idx, locus_idx}) =
            transmission_prob;
    }
}

void Chain::recalculate_transmission_at_locus_after_p_change(
    std::size_t pop_idx,
    std::size_t locus_idx,
    std::span<const float> p_old_span)
{
    ProfileScope scope("Chain::update_p::recalc_transmission");
    const std::size_t n_samples = genotyping_data.num_samples;
    const auto [p_begin, p_end] = p.inner_iterators({pop_idx, locus_idx});
    const std::span<const float> p_new_span(p_begin, p_end);

    std::vector<PChangePamGroup> groups;
    std::vector<int> sample_group;
    groups.clear();
    sample_group.assign(n_samples, -1);

    for (std::size_t s = 0; s < n_samples; ++s) {
        if (genotyping_data.is_missing(s, locus_idx)) {
            transmission_llik_new.unchecked_at({s, pop_idx, locus_idx}) = 0.f;
            continue;
        }

        const auto support = latent_allele_support(s, locus_idx);
        const int coi = m.at({s});
        if (support.empty() || support.size() > static_cast<std::size_t>(coi)) {
            transmission_llik_new.unchecked_at({s, pop_idx, locus_idx}) =
                -std::numeric_limits<float>::infinity();
            continue;
        }

        int group_idx = -1;
        for (std::size_t g = 0; g < groups.size(); ++g) {
            if (groups[g].coi == coi && support_spans_equal(support, groups[g].support)) {
                group_idx = static_cast<int>(g);
                break;
            }
        }
        if (group_idx < 0) {
            PChangePamGroup group;
            group.support.assign(support.begin(), support.end());
            group.coi = coi;
            group.total_alleles = group.support.size();
            groups.push_back(std::move(group));
            group_idx = static_cast<int>(groups.size() - 1);
        }
        sample_group[s] = group_idx;
    }

    {
        ProfileScope scope_groups("Chain::update_p::pam_groups");
        thread_local probAnyMissingFunctor functor;
        thread_local std::vector<float> q;
        thread_local std::vector<float> q_prev;
        for (auto& group : groups) {
            q.clear();
            q_prev.clear();
            float sum = 0.f;
            float sum_prev = 0.f;
            if (!transmission_incremental::build_constrained_q(
                    std::span<const int>(group.support), p_new_span, q, sum) ||
                !transmission_incremental::build_constrained_q(
                    std::span<const int>(group.support), p_old_span, q_prev, sum_prev)) {
                group.log_sum = -std::numeric_limits<float>::infinity();
                group.pam_valid = false;
                continue;
            }
            group.log_sum = std::log(sum);
            group.total_alleles = group.support.size();

            if (group.total_alleles == 1 && pam_fast_paths::tx_opts_enabled()) {
                group.pam_valid = false;
                continue;
            }

            const PamCachedVectors& pam_ref = cached_pam_vector(
                functor, q, 1u, static_cast<unsigned>(group.coi), 0u,
                std::span<const float>(q_prev));
            copy_pam_cached(pam_ref, group.pam);
            group.pam_valid = true;
        }
    }

    moire_parallel::recalc_parallel_for(0, n_samples, [&](std::size_t s) {
        const int group_idx = sample_group[s];
        if (group_idx < 0) {
            return;
        }
        const PChangePamGroup& group = groups[static_cast<std::size_t>(group_idx)];
        const float tx = finish_transmission_from_group(
            group.pam_valid ? &group.pam : nullptr,
            sampler,
            params.allow_relatedness,
            m.at({s}),
            group.total_alleles,
            r.at({s}),
            group.log_sum,
            group.pam_valid);
        transmission_llik_new.unchecked_at({s, pop_idx, locus_idx}) = tx;
    });
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
