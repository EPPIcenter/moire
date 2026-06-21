#include "chain.h"

#include "prob_any_missing.h"
#include "prob_any_missing_cache.h"
#include "transmission_incremental.h"
#include "ecoi_marginal.h"
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
                                     bool pam_valid,
                                     std::span<const float> log_binom_w = {})
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
        if (!log_binom_w.empty()) {
            return transmission_incremental::transmission_log_with_relatedness_k1(
                log_binom_w, coi, log_sum);
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
    if (!log_binom_w.empty()) {
        return transmission_incremental::transmission_log_with_relatedness(
            std::span<const double>(pam->pam),
            std::span<const float>(pam->log_one_minus_pam),
            log_binom_w, coi, total_alleles, log_sum);
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

bool support_contains_allele(std::span<const int> support,
                             std::size_t allele_idx) noexcept
{
    for (int allele : support) {
        if (static_cast<std::size_t>(allele) == allele_idx) {
            return true;
        }
    }
    return false;
}

struct PChangePamGroup {
    std::vector<int> support;
    int coi{0};
    PamCachedVectors pam;
    float log_sum{0.f};
    std::size_t total_alleles{0};
    bool pam_valid{false};
};

struct LocusPopPamGroup {
    std::vector<float> q;
    float log_sum{0.f};
    std::vector<std::size_t> pops;
};

/// Dedup (q, coi) PAM vectors across loci in one sample recalculation pass.
struct SamplePamCache {
    struct Entry {
        std::vector<float> q;
        unsigned coi{0};
        unsigned prev_coi_u{0};
        PamCachedVectors pam;
        bool pam_valid{false};
    };

    std::vector<Entry> entries;

    const Entry* find(std::span<const float> q,
                      unsigned coi,
                      unsigned prev_coi_u) const noexcept
    {
        for (const Entry& entry : entries) {
            if (entry.coi == coi && entry.prev_coi_u == prev_coi_u &&
                q_spans_equal(q, entry.q)) {
                return &entry;
            }
        }
        return nullptr;
    }

    int find_index(std::span<const float> q,
                   unsigned coi,
                   unsigned prev_coi_u) const noexcept
    {
        for (std::size_t i = 0; i < entries.size(); ++i) {
            const Entry& entry = entries[i];
            if (entry.coi == coi && entry.prev_coi_u == prev_coi_u &&
                q_spans_equal(q, entry.q)) {
                return static_cast<int>(i);
            }
        }
        return -1;
    }

    const Entry& find_or_add(std::span<const float> q,
                             unsigned coi,
                             unsigned prev_coi_u,
                             probAnyMissingFunctor& functor)
    {
        if (const Entry* hit = find(q, coi, prev_coi_u)) {
            return *hit;
        }
        Entry entry;
        entry.q.assign(q.begin(), q.end());
        entry.coi = coi;
        entry.prev_coi_u = prev_coi_u;
        if (q.size() <= 1) {
            entry.pam_valid = false;
        } else {
            const PamCachedVectors& pam_ref =
                cached_pam_vector(functor, q, 1u, coi, prev_coi_u);
            copy_pam_cached(pam_ref, entry.pam);
            entry.pam_valid = true;
        }
        entries.push_back(std::move(entry));
        return entries.back();
    }
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
    invalidate_update_p_locus_group_cache(locus_idx);
}

void Chain::restore_latent_genotype_new(std::size_t sample_idx, std::size_t locus_idx)
{
    const auto [begin, end] =
        latent_genotypes_old.inner_iterators({sample_idx, locus_idx});
    std::copy(begin, end, latent_genotypes_new.inner_begin({sample_idx, locus_idx}));
    refresh_latent_support_k(sample_idx, locus_idx);
    invalidate_update_p_locus_group_cache(locus_idx);
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

void Chain::fold_sample_tx_after_cell_updates(std::size_t sample_idx)
{
    ProfileScope scope("Chain::fold_sample_tx");
    if (!tx_llik_cache_valid_) {
        rebuild_transmission_llik_cache();
        return;
    }
    const std::size_t n_pops = params.num_populations;
    const std::size_t n_loci = genotyping_data.num_loci;
    const std::size_t tx_stride0 = transmission_llik_new.strides()[0];
    const std::size_t tx_stride1 = transmission_llik_new.strides()[1];
    const float* tx_data = transmission_llik_new.data().data();
    for (std::size_t pop_idx = 0; pop_idx < n_pops; ++pop_idx) {
        const std::size_t start_idx = sample_idx * tx_stride0 + pop_idx * tx_stride1;
        tx_loci_sum_new.at({sample_idx, pop_idx}) =
            moire_fused::detail::sum_loci_slice<float>(tx_data, start_idx, n_loci);
    }
    refresh_sample_tx_after_coi_change(sample_idx);
}

void Chain::apply_transmission_column_change(std::size_t pop_idx, std::size_t locus_idx)
{
    ProfileScope scope("Chain::apply_transmission_column_change");
    const std::size_t n_samples = genotyping_data.num_samples;
    moire_parallel::parallel_for(0, n_samples, [&](std::size_t sample_idx) {
        const float old_cell =
            transmission_llik_old.unchecked_at({sample_idx, pop_idx, locus_idx});
        const float new_cell =
            transmission_llik_new.unchecked_at({sample_idx, pop_idx, locus_idx});
        apply_transmission_cell_change(sample_idx, pop_idx, old_cell, new_cell);
    });
}

void Chain::restore_transmission_column_change(std::size_t pop_idx, std::size_t locus_idx)
{
    ProfileScope scope("Chain::update_p::reject_restore");
    const std::size_t n_samples = genotyping_data.num_samples;
    moire_parallel::parallel_for(0, n_samples, [&](std::size_t sample_idx) {
        const float proposed =
            transmission_llik_new.unchecked_at({sample_idx, pop_idx, locus_idx});
        const float old_cell =
            transmission_llik_old.unchecked_at({sample_idx, pop_idx, locus_idx});
        apply_transmission_cell_change(sample_idx, pop_idx, proposed, old_cell);
        transmission_llik_new.unchecked_at({sample_idx, pop_idx, locus_idx}) = old_cell;
    });
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

    const unsigned prev_coi_u =
        (coi != prev_coi) ? static_cast<unsigned>(prev_coi) : 0u;
    const unsigned coi_u = static_cast<unsigned>(coi);

    std::vector<std::size_t> dirty_loci;
    dirty_loci.clear();
    dirty_loci.reserve(n_loci);
    for (std::size_t locus_idx = 0; locus_idx < n_loci; ++locus_idx) {
        if (genotyping_data.is_missing(sample_idx, locus_idx)) {
            for (std::size_t pop_idx = 0; pop_idx < n_pops; ++pop_idx) {
                transmission_llik_new.unchecked_at({sample_idx, pop_idx, locus_idx}) = 0.f;
            }
            continue;
        }
        if (!locus_tx_inputs_unchanged(sample_idx, locus_idx, prev_coi, prev_r)) {
            dirty_loci.push_back(locus_idx);
        }
    }

    if (dirty_loci.empty()) {
        return;
    }

    SamplePamCache sample_pam;
    if (!dirty_loci.empty() && pam_fast_paths::tx_opts_enabled()) {
        ProfileScope scope_prefill("Chain::sample_tx::pam_prefill");
        thread_local probAnyMissingFunctor prefill_functor;
        thread_local std::vector<float> q_prefill;
        for (const std::size_t locus_idx : dirty_loci) {
            const int support_group_idx = locus_support_group[locus_idx];
            if (support_group_idx < 0) {
                continue;
            }
            const PChangePamGroup& support_group =
                support_groups[static_cast<std::size_t>(support_group_idx)];
            const auto support = std::span<const int>(support_group.support);
            const std::size_t total_alleles = support_group.total_alleles;
            if (total_alleles <= 1 || total_alleles > static_cast<std::size_t>(coi)) {
                continue;
            }
            for (std::size_t pop_idx = 0; pop_idx < n_pops; ++pop_idx) {
                const auto [p_begin, p_end] = p.inner_iterators({pop_idx, locus_idx});
                q_prefill.clear();
                float sum = 0.f;
                if (!transmission_incremental::build_constrained_q(
                        support, std::span<const float>(p_begin, p_end), q_prefill, sum)) {
                    continue;
                }
                sample_pam.find_or_add(q_prefill, coi_u, prev_coi_u, prefill_functor);
            }
        }
    }

    std::vector<float> log_binom_w;
    if (params.allow_relatedness && coi >= 1) {
        log_binom_w.resize(static_cast<std::size_t>(coi));
        for (std::size_t i = 0; i < log_binom_w.size(); ++i) {
            log_binom_w[i] =
                sampler.dbinom(static_cast<int>(i), coi - 1, relatedness);
        }
    }
    const std::span<const float> log_binom_w_span(log_binom_w);

    // dirty_loci must be an ordinary local vector — not thread_local. TBB workers
    // each get their own thread_local storage; a prior parallel attempt read empty
    // lists and segfaulted. Per-locus work here is too small for recalc_parallel_for
    // to beat serial (TBB overhead dominates typical dirty-locus counts).
    const SamplePamCache* const pam_table = &sample_pam;
    const auto process_dirty_locus = [&](std::size_t dirty_idx) {
        const std::size_t locus_idx = dirty_loci[dirty_idx];

        const int support_group_idx = locus_support_group[locus_idx];
        if (support_group_idx < 0) {
            return;
        }
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

        for (const auto& group : pop_groups) {
            thread_local PamCachedVectors pam_local;
            const PamCachedVectors* pam_ptr = nullptr;
            bool pam_valid = false;
            if (total_alleles > 1) {
                const int pam_idx =
                    pam_table->find_index(group.q, coi_u, prev_coi_u);
                if (pam_idx >= 0) {
                    const SamplePamCache::Entry& entry =
                        pam_table->entries[static_cast<std::size_t>(pam_idx)];
                    if (entry.pam_valid) {
                        copy_pam_cached(entry.pam, pam_local);
                        pam_ptr = &pam_local;
                        pam_valid = true;
                    }
                }
            }
            const float tx = finish_transmission_from_group(
                pam_ptr,
                sampler,
                params.allow_relatedness,
                coi,
                total_alleles,
                relatedness,
                group.log_sum,
                pam_valid,
                log_binom_w_span);
            for (const std::size_t pop_idx : group.pops) {
                transmission_llik_new.unchecked_at({sample_idx, pop_idx, locus_idx}) = tx;
            }
        }
    };

    for (std::size_t dirty_idx = 0; dirty_idx < dirty_loci.size(); ++dirty_idx) {
        process_dirty_locus(dirty_idx);
    }

    fold_sample_tx_after_cell_updates(sample_idx);
}

void Chain::restore_transmission_for_sample_incremental(std::size_t sample_idx)
{
    moire_parallel::parallel_for_2d(
        0, params.num_populations, 0, genotyping_data.num_loci,
        [&](std::size_t pop_idx, std::size_t locus_idx) {
            transmission_llik_new.unchecked_at({sample_idx, pop_idx, locus_idx}) =
                transmission_llik_old.unchecked_at({sample_idx, pop_idx, locus_idx});
        });
    fold_sample_tx_after_cell_updates(sample_idx);
}

float Chain::calc_transmission_llik_sum() {
    if (!tx_llik_cache_valid_) {
        rebuild_transmission_llik_cache();
    }
    return tx_llik_sum_new;
}

// eCOI mode: full marginalized transmission log-likelihood, summed over samples,
// recomputed from scratch from the current latent supports, p, and eff_coi. The
// discrete COI (and relatedness, via r(m,e)) is integrated analytically. This is
// used *inside* the collapsed eCOI move (update_ecoi) to score effective-COI
// proposals; it is NOT the persistent chain likelihood (which stays the legacy
// single-(m,r) form, because the latent-genotype proposal keeps COI as an
// auxiliary variable whose proposal density is only available at draw time).
float Chain::calc_marginal_transmission_llik_sum() {
    ProfileScope scope("Chain::calc_marginal_transmission_llik_sum");
    double total = 0.0;
    for (std::size_t s = 0; s < genotyping_data.num_samples; ++s) {
        total += ecoi_sample_marginal_transmission_llik(
            s, static_cast<double>(eff_coi.at({s})));
    }
    return static_cast<float>(total);
}

float Chain::calc_new_likelihood() {
    ProfileScope scope("Chain::calc_new_likelihood");
    // observation_llik_new sum is maintained incrementally (see
    // sync_obs_sum_for_sample); avoids an O(N*L) reduction per proposal.
    const float observation_llik = obs_llik_sum_new_;

    // Fully-collapsed (eCOI) mode: transmission is the marginalized term tracked
    // in ecoi_marg_sum_. Moves that change the marginal (p, eCOI, latent
    // genotypes, population params) update ecoi_marg_sum_ themselves; moves that
    // only touch the observation model (eps) read it unchanged here.
    if (params.marginal_ecoi) {
        return observation_llik + static_cast<float>(ecoi_marg_sum_);
    }

    float transmission_llik;
    {
        ProfileScope s_tx("Chain::calc_new_likelihood::transmission_sum");
        transmission_llik = calc_transmission_llik_sum();
    }

    return observation_llik + transmission_llik;
}

double Chain::recompute_collapsed_marginals(std::vector<float> &out) {
    out.resize(genotyping_data.num_samples);
    double sum = 0.0;
    for (std::size_t s = 0; s < genotyping_data.num_samples; ++s) {
        const double m = ecoi_sample_marginal_llik(s);
        out[s] = static_cast<float>(m);
        sum += m;
    }
    return sum;
}

void Chain::initialize_collapsed_marginal_cache() {
    ecoi_marg_sum_ = recompute_collapsed_marginals(ecoi_sample_marg_);
    ecoi_sample_marg_scratch_.assign(genotyping_data.num_samples, 0.0f);
    llik = obs_llik_sum_new_ + static_cast<float>(ecoi_marg_sum_);
}

float Chain::calc_new_prior() {
    ProfileScope scope("Chain::calc_new_prior");
    // Per-sample prior sums are maintained incrementally in calculate_*/restore_*
    // to avoid an O(N) reduction per proposal.
    // In the population-e (eCOI) hierarchy the relatedness weight and the
    // per-population continuous e-prior log f_p(e_s) are both folded into the
    // marginal likelihood (ecoi_marg_sum_, inside the population mixture), so the
    // only eCOI prior term left here is the hyperprior on the per-population
    // (mu_plus_p, k_p), tracked in ecoi_hyperprior_.
    if (params.marginal_ecoi) {
        // False positives are pooled per-locus in eCOI mode, so the FP prior is
        // the per-locus sum rather than the (unused) per-sample sum.
        return eps_neg_prior_sum_new_ +
               eps_pos_locus_prior_sum_new_ +
               static_cast<float>(ecoi_hyperprior_) +
               population_responsibility_vector_prior_new;
    }
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
            eps_neg.at({sample_idx}), eps_pos_at(sample_idx, locus_idx)
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

void Chain::invalidate_update_p_locus_group_cache(std::size_t locus_idx)
{
    if (locus_idx < update_p_locus_group_cache_.size()) {
        update_p_locus_group_cache_[locus_idx].valid = false;
        for (std::size_t pop_idx = 0; pop_idx < update_p_pam_slots_.size(); ++pop_idx) {
            update_p_pam_slots_[pop_idx][locus_idx].clear();
        }
    }
}

void Chain::ensure_update_p_locus_group_cache(std::size_t locus_idx)
{
    UpdatePLocusGroupCache& cache = update_p_locus_group_cache_[locus_idx];
    if (cache.valid) {
        return;
    }

    ProfileScope scope("Chain::update_p::build_group_cache");
    const std::size_t n_samples = genotyping_data.num_samples;
    cache.groups.clear();
    cache.sample_group.assign(n_samples, UpdatePLocusGroupCache::kSampleMissing);

    for (std::size_t s = 0; s < n_samples; ++s) {
        if (genotyping_data.is_missing(s, locus_idx)) {
            continue;
        }

        const auto support = latent_allele_support(s, locus_idx);
        const int coi = m.at({s});
        if (support.empty() || support.size() > static_cast<std::size_t>(coi)) {
            cache.sample_group[s] = UpdatePLocusGroupCache::kSampleInvalid;
            continue;
        }

        int group_idx = -1;
        for (std::size_t g = 0; g < cache.groups.size(); ++g) {
            if (cache.groups[g].coi == coi &&
                support_spans_equal(support, cache.groups[g].support)) {
                group_idx = static_cast<int>(g);
                break;
            }
        }
        if (group_idx < 0) {
            UpdatePLocusGroupCache::Group group;
            group.support.assign(support.begin(), support.end());
            group.coi = coi;
            group.total_alleles = group.support.size();
            cache.groups.push_back(std::move(group));
            group_idx = static_cast<int>(cache.groups.size() - 1);
        }
        cache.sample_group[s] = group_idx;
    }
    cache.valid = true;
}

void Chain::recalculate_transmission_at_locus_after_p_change(
    std::size_t pop_idx,
    std::size_t locus_idx,
    std::span<const float> p_old_span,
    std::size_t changed_allele_idx)
{
    ProfileScope scope("Chain::update_p::recalc_transmission");
    const std::size_t n_samples = genotyping_data.num_samples;
    const auto [p_begin, p_end] = p.inner_iterators({pop_idx, locus_idx});
    const std::span<const float> p_new_span(p_begin, p_end);

    ensure_update_p_locus_group_cache(locus_idx);
    const UpdatePLocusGroupCache& locus_cache = update_p_locus_group_cache_[locus_idx];
    const std::vector<int>& sample_group = locus_cache.sample_group;

    constexpr float kNegInf = -std::numeric_limits<float>::infinity();
    for (std::size_t s = 0; s < n_samples; ++s) {
        const int group_idx = sample_group[s];
        if (group_idx == UpdatePLocusGroupCache::kSampleMissing) {
            transmission_llik_new.unchecked_at({s, pop_idx, locus_idx}) = 0.f;
        } else if (group_idx == UpdatePLocusGroupCache::kSampleInvalid) {
            transmission_llik_new.unchecked_at({s, pop_idx, locus_idx}) = kNegInf;
        }
    }

    std::vector<PChangePamGroup> groups;
    groups.reserve(locus_cache.groups.size());
    for (const UpdatePLocusGroupCache::Group& cached : locus_cache.groups) {
        PChangePamGroup group;
        group.support = cached.support;
        group.coi = cached.coi;
        group.total_alleles = cached.total_alleles;
        groups.push_back(std::move(group));
    }

    std::vector<UpdatePPamSlot>& pam_slots = update_p_pam_slots_[pop_idx][locus_idx];
    if (pam_slots.size() < groups.size()) {
        pam_slots.resize(groups.size());
    }

    {
        ProfileScope scope_groups("Chain::update_p::pam_groups");
        thread_local probAnyMissingFunctor functor;
        thread_local std::vector<float> q;
        thread_local std::vector<float> q_prev;
        for (std::size_t gi = 0; gi < groups.size(); ++gi) {
            PChangePamGroup& group = groups[gi];
            UpdatePPamSlot& slot = pam_slots[gi];

            if (changed_allele_idx < p_new_span.size() &&
                !support_contains_allele(
                    std::span<const int>(group.support), changed_allele_idx) &&
                slot.pam_valid) {
                group.log_sum = slot.log_sum;
                group.pam_valid = true;
                copy_pam_cached(slot.pam, group.pam);
                continue;
            }

            q.clear();
            q_prev.clear();
            float sum = 0.f;
            float sum_prev = 0.f;
            if (!transmission_incremental::build_constrained_q(
                    std::span<const int>(group.support), p_new_span, q, sum) ||
                !transmission_incremental::build_constrained_q(
                    std::span<const int>(group.support), p_old_span, q_prev, sum_prev)) {
                group.log_sum = kNegInf;
                group.pam_valid = false;
                slot.pam_valid = false;
                continue;
            }
            group.log_sum = std::log(sum);
            group.total_alleles = group.support.size();

            if (group.total_alleles == 1 && pam_fast_paths::tx_opts_enabled()) {
                group.pam_valid = false;
                slot.pam_valid = false;
                continue;
            }

            if (q_spans_equal(q, q_prev) && slot.pam_valid && q_spans_equal(q, slot.q)) {
                group.pam_valid = true;
                copy_pam_cached(slot.pam, group.pam);
                continue;
            }

            if (q_spans_equal(q, slot.q) && slot.pam_valid) {
                group.pam_valid = true;
                copy_pam_cached(slot.pam, group.pam);
                slot.log_sum = group.log_sum;
                continue;
            }

            const PamCachedVectors& pam_ref = cached_pam_vector(
                functor, q, 1u, static_cast<unsigned>(group.coi), 0u,
                std::span<const float>(q_prev));
            copy_pam_cached(pam_ref, group.pam);
            group.pam_valid = true;
            slot.q.assign(q.begin(), q.end());
            copy_pam_cached(group.pam, slot.pam);
            slot.log_sum = group.log_sum;
            slot.pam_valid = true;
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

void Chain::calculate_eps_pos_locus_likelihood(std::size_t locus_idx)
{
    const float val = sampler.get_beta_log_prior(
        eps_pos_locus.at({locus_idx}), params.eps_pos_locus_alpha,
        params.eps_pos_locus_beta);
    eps_pos_locus_prior_sum_new_ += val - eps_pos_locus_prior_new.at({locus_idx});
    eps_pos_locus_prior_new.at({locus_idx}) = val;
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

    if (params.marginal_ecoi)
    {
        eps_pos_locus_prior_new.resize({num_loci});
        eps_pos_locus_prior_old.resize({num_loci});
    }

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
    if (params.marginal_ecoi)
    {
        for (std::size_t locus_idx = 0; locus_idx < num_loci; ++locus_idx)
        {
            calculate_eps_pos_locus_likelihood(locus_idx);
            save_eps_pos_locus_likelihood(locus_idx);
        }
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
    if (params.marginal_ecoi)
    {
        eps_pos_locus_prior_sum_new_ = eps_pos_locus_prior_new.full_sum();
    }
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

    update_p_locus_group_cache_.assign(num_loci, UpdatePLocusGroupCache{});
    update_p_pam_slots_.assign(
        num_populations,
        std::vector<std::vector<UpdatePPamSlot>>(num_loci));

    if (params.marginal_ecoi) {
        const double diff = ecoi_transmission_selfcheck();
        if (params.verbose) {
            UtilFunctions::print(
                "eCOI core vs production transmission self-check: max|diff| =", diff);
        }
        if (params.allow_relatedness) {
            const double adiff = ecoi_assembly_selfcheck();
            if (params.verbose) {
                UtilFunctions::print(
                    "eCOI marginal assembly vs production brute-force: max|diff| =",
                    adiff);
            }
        }
        // Switch the persistent likelihood onto the collapsed marginal.
        initialize_collapsed_marginal_cache();
    }
}

// Diagnostic (Stage 3a validation gate): confirm the eCOI marginal core's
// per-locus transmission term (inference/ecoi_marginal.h) reproduces the
// production calc_transmission_process on the chain's *real* latent supports and
// allele frequencies, across several (COI, relatedness) values. Returns the max
// absolute log-likelihood discrepancy (inf if one path is finite and the other
// is not). This is the non-circular check that the eCOI implementation matches
// the existing transmission machinery before it replaces it in the hot path.
double Chain::ecoi_transmission_selfcheck()
{
    double max_abs = 0.0;
    std::vector<float> support_p;
    const int max_coi = static_cast<int>(params.max_coi);

    for (std::size_t s = 0; s < genotyping_data.num_samples; ++s) {
        for (std::size_t pop = 0; pop < params.num_populations; ++pop) {
            for (std::size_t l = 0; l < genotyping_data.num_loci; ++l) {
                if (genotyping_data.is_missing(s, l)) {
                    continue;
                }
                const auto support = latent_allele_support(s, l);
                if (support.empty()) {
                    continue;
                }
                const auto [p_begin, p_end] = p.inner_iterators({pop, l});
                const std::span<const float> p_span(p_begin, p_end);

                support_p.clear();
                bool ok = true;
                for (int a : support) {
                    if (a < 0 || static_cast<std::size_t>(a) >= p_span.size()) {
                        ok = false;
                        break;
                    }
                    support_p.push_back(p_span[a]);
                }
                if (!ok) {
                    continue;
                }

                const int k = static_cast<int>(support.size());
                const auto pre = ecoi_marginal::precompute_locus(
                    std::span<const float>(support_p.data(), support_p.size()), max_coi);

                const int m_candidates[] = {k, k + 1, std::min(k + 3, max_coi), max_coi};
                for (int m : m_candidates) {
                    if (m < k || m > max_coi) {
                        continue;
                    }
                    // The production Sampler::dbinom is only well-defined for
                    // relatedness strictly inside (0, 1) (the sampler clamps r to
                    // (1e-5, 1-1e-5)); r == 0 yields 0*log(0)=NaN there. So when
                    // relatedness is enabled, probe interior r; otherwise exercise
                    // the no-relatedness path (production ignores r, core uses r=0).
                    std::vector<float> r_values;
                    if (params.allow_relatedness) {
                        r_values = {0.001f, 0.3f, 0.7f};
                    } else {
                        r_values = {0.0f};
                    }
                    for (float rr : r_values) {
                        const float prod =
                            calc_transmission_process(support, p_span, m, rr);
                        const double core = ecoi_marginal::tx_loglik_locus(
                            pre, m, static_cast<double>(rr));
                        const bool fp = std::isfinite(prod);
                        const bool fc = std::isfinite(core);
                        if (fp && fc) {
                            max_abs = std::max(
                                max_abs,
                                std::fabs(static_cast<double>(prod) - core));
                        } else if (fp != fc) {
                            max_abs = std::numeric_limits<double>::infinity();
                        }
                    }
                }
            }
        }
    }
    return max_abs;
}

// Per-sample marginalized transmission log-likelihood at effective COI e:
//   LSE_p[ log pi_p + LSE_m( log w(m|e) + sum_l t_l(s,p,m,r(m,e)) ) ]
// The COI prior folds into w(m|e); missing/empty-support loci contribute 0
// (they neither constrain m nor add to the product). Continuous part (e > 1).
double Chain::ecoi_sample_marginal_transmission_llik(std::size_t sample_idx, double e)
{
    const ecoi_marginal::Hyperparams hp{
        static_cast<double>(population_coi_p), static_cast<double>(population_coi_r),
        static_cast<double>(params.r_alpha), static_cast<double>(params.r_beta)};
    const int max_coi = static_cast<int>(params.max_coi);
    const std::size_t num_pop = params.num_populations;

    std::vector<double> per_pop(num_pop, ecoi_marginal::kNegInf);
    std::vector<ecoi_marginal::LocusSupport> loci;
    std::vector<float> support_p;

    for (std::size_t pop = 0; pop < num_pop; ++pop) {
        loci.clear();
        for (std::size_t l = 0; l < genotyping_data.num_loci; ++l) {
            if (genotyping_data.is_missing(sample_idx, l)) {
                continue;
            }
            const auto support = latent_allele_support(sample_idx, l);
            if (support.empty()) {
                continue;
            }
            const auto [pb, pe] = p.inner_iterators({pop, l});
            const std::span<const float> p_span(pb, pe);
            support_p.clear();
            bool ok = true;
            for (int a : support) {
                if (a < 0 || static_cast<std::size_t>(a) >= p_span.size()) {
                    ok = false;
                    break;
                }
                support_p.push_back(p_span[a]);
            }
            if (!ok) {
                continue;
            }
            loci.push_back(ecoi_marginal::precompute_locus(
                std::span<const float>(support_p.data(), support_p.size()), max_coi));
        }
        const double pop_log =
            std::log(static_cast<double>(population_responsibility_vector.at({pop})));
        const double inner = ecoi_marginal::log_marginal_e(
            e, std::span<const ecoi_marginal::LocusSupport>(loci.data(), loci.size()),
            hp, max_coi);
        per_pop[pop] = pop_log + inner;
    }
    return ecoi_marginal::log_sum_exp(std::span<const double>(per_pop));
}

// ---- Population-e hierarchy: COI-prior-free per-sample marginals -----------

void Chain::build_sample_pop_loci(
    std::size_t sample_idx, int max_coi,
    std::vector<std::vector<ecoi_marginal::LocusSupport>> &out) const
{
    const std::size_t num_pop = params.num_populations;
    out.assign(num_pop, {});
    std::vector<float> support_p;
    for (std::size_t pop = 0; pop < num_pop; ++pop) {
        std::vector<ecoi_marginal::LocusSupport> &loci = out[pop];
        for (std::size_t l = 0; l < genotyping_data.num_loci; ++l) {
            if (genotyping_data.is_missing(sample_idx, l)) {
                continue;
            }
            const auto support = latent_allele_support(sample_idx, l);
            if (support.empty()) {
                continue;
            }
            const auto [pb, pe] = p.inner_iterators({pop, l});
            const std::span<const float> p_span(pb, pe);
            support_p.clear();
            bool ok = true;
            for (int a : support) {
                if (a < 0 || static_cast<std::size_t>(a) >= p_span.size()) {
                    ok = false;
                    break;
                }
                support_p.push_back(p_span[a]);
            }
            if (!ok) {
                continue;
            }
            loci.push_back(ecoi_marginal::precompute_locus(
                std::span<const float>(support_p.data(), support_p.size()), max_coi));
        }
    }
}

// ---- per-(sample, pop, locus) LocusSupport cache --------------------------

void Chain::ecoi_support_cache_init()
{
    const std::size_t N = genotyping_data.num_samples;
    const std::size_t P = params.num_populations;
    const std::size_t L = genotyping_data.num_loci;
    ecoi_cache_pop_stride_ = L;
    ecoi_cache_sample_stride_ = P * L;
    const std::size_t total = N * P * L;
    ecoi_support_cache_.assign(total, ecoi_marginal::LocusSupport{});
    ecoi_support_valid_.assign(total, 0);
    ecoi_support_contributes_.assign(total, 0);
}

void Chain::ecoi_support_cache_invalidate_locus(std::size_t pop_idx,
                                                std::size_t locus_idx)
{
    if (ecoi_support_valid_.empty()) {
        return;
    }
    const std::size_t N = genotyping_data.num_samples;
    for (std::size_t s = 0; s < N; ++s) {
        const std::size_t idx =
            s * ecoi_cache_sample_stride_ + pop_idx * ecoi_cache_pop_stride_ + locus_idx;
        ecoi_support_valid_[idx] = 0;
    }
}

void Chain::ecoi_support_cache_invalidate_sample(std::size_t sample_idx)
{
    if (ecoi_support_valid_.empty()) {
        return;
    }
    const std::size_t base = sample_idx * ecoi_cache_sample_stride_;
    for (std::size_t i = 0; i < ecoi_cache_sample_stride_; ++i) {
        ecoi_support_valid_[base + i] = 0;
    }
}

const ecoi_marginal::LocusSupport *Chain::ecoi_cached_support(
    std::size_t sample_idx, std::size_t pop_idx, std::size_t locus_idx)
{
    const std::size_t idx = sample_idx * ecoi_cache_sample_stride_ +
                            pop_idx * ecoi_cache_pop_stride_ + locus_idx;
    if (!ecoi_support_valid_[idx]) {
        bool contributes = false;
        if (!genotyping_data.is_missing(sample_idx, locus_idx)) {
            const auto support = latent_allele_support(sample_idx, locus_idx);
            if (!support.empty()) {
                const auto [pb, pe] = p.inner_iterators({pop_idx, locus_idx});
                const std::span<const float> p_span(pb, pe);
                ecoi_support_p_scratch_.clear();
                bool ok = true;
                for (int a : support) {
                    if (a < 0 || static_cast<std::size_t>(a) >= p_span.size()) {
                        ok = false;
                        break;
                    }
                    ecoi_support_p_scratch_.push_back(p_span[a]);
                }
                if (ok) {
                    ecoi_support_cache_[idx] = ecoi_marginal::precompute_locus(
                        std::span<const float>(ecoi_support_p_scratch_.data(),
                                               ecoi_support_p_scratch_.size()),
                        static_cast<int>(params.max_coi));
                    contributes = true;
                }
            }
        }
        ecoi_support_contributes_[idx] = contributes ? 1 : 0;
        ecoi_support_valid_[idx] = 1;
    }
    return ecoi_support_contributes_[idx] ? &ecoi_support_cache_[idx] : nullptr;
}

// log of the continuous per-population e-prior density at e:
//   log Gamma(e - 1; shape = k_p, rate = k_p / mu_plus_p).
double Chain::ecoi_log_f(std::size_t pop_idx, double e) const
{
    const double y = e - 1.0;
    if (y <= 0.0) {
        return ecoi_marginal::kNegInf;
    }
    const double k = static_cast<double>(ecoi_k[pop_idx]);
    const double mu = static_cast<double>(ecoi_mu_plus[pop_idx]);
    const double rate = k / mu;
    return k * std::log(rate) - std::lgamma(k) + (k - 1.0) * std::log(y) -
           rate * y;
}

// Per-sample marginal log-likelihood at effective COI e > 1. The discrete COI is
// integrated out under uniform relatedness (the only weight is the Jacobian
// 1/(m-1)), and the per-population continuous e-prior log f_p(e) is folded inside
// the population mixture so that the population that explains the data also
// supplies the e-prior:
//   LSE_p[ log pi_p + log f_p(e) + LSE_m( -log(m-1) + sum_l t_l(s,p,m,r(m,e)) ) ].
// Reads the cached per-(sample,pop,locus) LocusSupport (built lazily).
double Chain::ecoi_sample_rel_marginal_llik(std::size_t sample_idx, double e)
{
    const int max_coi = static_cast<int>(params.max_coi);
    const std::size_t L = genotyping_data.num_loci;
    const std::size_t P = params.num_populations;

    std::vector<double> per_pop(P, ecoi_marginal::kNegInf);
    const int m_lo = std::max(2, static_cast<int>(std::ceil(e)));

    for (std::size_t pop = 0; pop < P; ++pop) {
        const double pop_log = std::log(
            static_cast<double>(population_responsibility_vector.at({pop})));
        if (e <= 1.0 || m_lo > max_coi) {
            per_pop[pop] = pop_log + ecoi_marginal::kNegInf;
            continue;
        }
        ecoi_contrib_scratch_.clear();
        for (std::size_t l = 0; l < L; ++l) {
            const ecoi_marginal::LocusSupport *ls =
                ecoi_cached_support(sample_idx, pop, l);
            if (ls != nullptr) {
                ecoi_contrib_scratch_.push_back(ls);
            }
        }
        ecoi_terms_scratch_.clear();
        for (int m = m_lo; m <= max_coi; ++m) {
            const double r = ecoi_marginal::r_of(m, e);
            if (r <= 0.0 || r >= 1.0) {
                ecoi_terms_scratch_.push_back(ecoi_marginal::kNegInf);
                continue;
            }
            double tx = 0.0;
            bool ok = true;
            for (const ecoi_marginal::LocusSupport *ls : ecoi_contrib_scratch_) {
                const double t = ecoi_marginal::tx_loglik_locus(*ls, m, r);
                if (!std::isfinite(t)) {
                    ok = false;
                    break;
                }
                tx += t;
            }
            ecoi_terms_scratch_.push_back(
                ok ? -std::log(static_cast<double>(m - 1)) + tx
                   : ecoi_marginal::kNegInf);
        }
        const double inner =
            ecoi_marginal::log_sum_exp(std::span<const double>(ecoi_terms_scratch_));
        per_pop[pop] = pop_log + ecoi_log_f(pop, e) + inner;
    }
    return ecoi_marginal::log_sum_exp(std::span<const double>(per_pop));
}

double Chain::ecoi_sample_rel_transmission_at_e(std::size_t sample_idx, double e)
{
    const int max_coi = static_cast<int>(params.max_coi);
    const std::size_t L = genotyping_data.num_loci;
    const std::size_t P = params.num_populations;

    std::vector<double> per_pop(P, ecoi_marginal::kNegInf);
    const int m_lo = std::max(2, static_cast<int>(std::ceil(e)));

    for (std::size_t pop = 0; pop < P; ++pop) {
        const double pop_log = std::log(
            static_cast<double>(population_responsibility_vector.at({pop})));
        if (e <= 1.0 || m_lo > max_coi) {
            per_pop[pop] = pop_log + ecoi_marginal::kNegInf;
            continue;
        }
        ecoi_contrib_scratch_.clear();
        for (std::size_t l = 0; l < L; ++l) {
            const ecoi_marginal::LocusSupport *ls =
                ecoi_cached_support(sample_idx, pop, l);
            if (ls != nullptr) {
                ecoi_contrib_scratch_.push_back(ls);
            }
        }
        ecoi_terms_scratch_.clear();
        for (int m = m_lo; m <= max_coi; ++m) {
            const double r = ecoi_marginal::r_of(m, e);
            if (r <= 0.0 || r >= 1.0) {
                ecoi_terms_scratch_.push_back(ecoi_marginal::kNegInf);
                continue;
            }
            double tx = 0.0;
            bool ok = true;
            for (const ecoi_marginal::LocusSupport *ls : ecoi_contrib_scratch_) {
                const double t = ecoi_marginal::tx_loglik_locus(*ls, m, r);
                if (!std::isfinite(t)) {
                    ok = false;
                    break;
                }
                tx += t;
            }
            ecoi_terms_scratch_.push_back(
                ok ? -std::log(static_cast<double>(m - 1)) + tx
                   : ecoi_marginal::kNegInf);
        }
        const double inner =
            ecoi_marginal::log_sum_exp(std::span<const double>(ecoi_terms_scratch_));
        per_pop[pop] = pop_log + inner;
    }
    return ecoi_marginal::log_sum_exp(std::span<const double>(per_pop));
}

double Chain::ecoi_sample_marginal_llik(std::size_t sample_idx)
{
    return ecoi_sample_rel_marginal_llik(
        sample_idx, static_cast<double>(eff_coi.at({sample_idx})));
}

double Chain::ecoi_population_hyperprior() const
{
    // Shared hyperpriors applied to each population's (mu_plus_p, k_p):
    //   mu_plus_p ~ Exponential(ecoi_mu_rate),
    //   k_p       ~ Gamma(shape = ecoi_k_shape, rate = ecoi_k_rate).
    const double mu_rate = static_cast<double>(params.ecoi_mu_rate);
    const double ks = static_cast<double>(params.ecoi_k_shape);
    const double kr = static_cast<double>(params.ecoi_k_rate);
    const double ks_const = ks * std::log(kr) - std::lgamma(ks);
    double lp = 0.0;
    for (std::size_t pop = 0; pop < params.num_populations; ++pop) {
        lp += std::log(mu_rate) - mu_rate * static_cast<double>(ecoi_mu_plus[pop]);
        lp += ks_const + (ks - 1.0) * std::log(static_cast<double>(ecoi_k[pop])) -
              kr * static_cast<double>(ecoi_k[pop]);
    }
    return lp;
}

void Chain::recompute_ecoi_prior()
{
    ecoi_hyperprior_ = ecoi_population_hyperprior();
}

// ---- update_p incremental transmission decomposition (single population) ---

void Chain::ecoi_up_build_decomp()
{
    const int max_coi = static_cast<int>(params.max_coi);
    up_m_stride_ = static_cast<std::size_t>(max_coi) + 1;
    const std::size_t stride = up_m_stride_;
    const std::size_t N = genotyping_data.num_samples;
    const std::size_t L = genotyping_data.num_loci;
    const double pop_log = std::log(
        static_cast<double>(population_responsibility_vector.at({0})));

    up_mlo_.assign(N, 0);
    up_foff_.assign(N, 0.0);
    up_a_.assign(N * stride, ecoi_marginal::kNegInf);
    up_slog_.assign(N * stride, 0.0);
    up_ninf_.assign(N * stride, 0);
    up_L_.assign(N, ecoi_marginal::kNegInf);
    up_newL_.assign(N, ecoi_marginal::kNegInf);
    up_told_.assign(N * stride, 0.0);
    up_tnew_.assign(N * stride, 0.0);

    std::vector<const ecoi_marginal::LocusSupport *> contrib;
    std::vector<double> terms;

    for (std::size_t s = 0; s < N; ++s) {
        const double e = static_cast<double>(eff_coi.at({s}));
        const int m_lo = std::max(2, static_cast<int>(std::ceil(e)));
        up_mlo_[s] = m_lo;
        // Population-0 mixing weight + folded continuous e-prior; constant in p
        // for the duration of this update_p call (e_s is fixed).
        up_foff_[s] = pop_log + ecoi_log_f(0, e);

        contrib.clear();
        for (std::size_t l = 0; l < L; ++l) {
            const ecoi_marginal::LocusSupport *ls = ecoi_cached_support(s, 0, l);
            if (ls != nullptr) {
                contrib.push_back(ls);
            }
        }

        terms.clear();
        for (int m = m_lo; m <= max_coi; ++m) {
            const std::size_t idx = s * stride + static_cast<std::size_t>(m);
            const double r = ecoi_marginal::r_of(m, e);
            if (r <= 0.0 || r >= 1.0) {
                up_a_[idx] = ecoi_marginal::kNegInf;
                terms.push_back(ecoi_marginal::kNegInf);
                continue;
            }
            up_a_[idx] = -std::log(static_cast<double>(m - 1));  // uniform relatedness
            double sfin = 0.0;
            int ninf = 0;
            for (const ecoi_marginal::LocusSupport *ls : contrib) {
                const double t = ecoi_marginal::tx_loglik_locus(*ls, m, r);
                if (std::isfinite(t)) {
                    sfin += t;
                } else {
                    ++ninf;
                }
            }
            up_slog_[idx] = sfin;
            up_ninf_[idx] = ninf;
            terms.push_back(ninf == 0 ? up_a_[idx] + sfin : ecoi_marginal::kNegInf);
        }
        up_L_[s] = up_foff_[s] +
                   ecoi_marginal::log_sum_exp(std::span<const double>(terms));
    }
}

// Per-locus transmission term t_l(m, r(m, e_s)) for population 0, written into
// `out` over each sample's m range. Non-finite (incompatible) and
// non-contributing entries are written as 0 so that the proposal delta
// tnew - told is exactly 0 for them (the incompatibility of a locus does not
// change with p, since the latent support is fixed).
void Chain::ecoi_up_locus_terms(std::size_t locus_idx, std::vector<double> &out)
{
    const int max_coi = static_cast<int>(params.max_coi);
    const std::size_t stride = up_m_stride_;
    const std::size_t N = genotyping_data.num_samples;

    for (std::size_t s = 0; s < N; ++s) {
        const ecoi_marginal::LocusSupport *ls = ecoi_cached_support(s, 0, locus_idx);
        const double e = static_cast<double>(eff_coi.at({s}));
        const int m_lo = up_mlo_[s];
        for (int m = m_lo; m <= max_coi; ++m) {
            const std::size_t idx = s * stride + static_cast<std::size_t>(m);
            if (ls == nullptr) {
                out[idx] = 0.0;
                continue;
            }
            const double r = ecoi_marginal::r_of(m, e);
            if (r <= 0.0 || r >= 1.0) {
                out[idx] = 0.0;
                continue;
            }
            const double t = ecoi_marginal::tx_loglik_locus(*ls, m, r);
            out[idx] = std::isfinite(t) ? t : 0.0;
        }
    }
}

// Stage 3b validation gate: confirm the in-Chain marginal assembly (population
// mixture + induced weight w(m|e) + LSE over m, all via ecoi_marginal.h) matches
// an independent brute-force over m that uses the *production*
// calc_transmission_process for the per-locus term. Only meaningful with
// relatedness enabled (eCOI == COI otherwise). Returns max abs discrepancy.
double Chain::ecoi_assembly_selfcheck()
{
    const ecoi_marginal::Hyperparams hp{
        static_cast<double>(population_coi_p), static_cast<double>(population_coi_r),
        static_cast<double>(params.r_alpha), static_cast<double>(params.r_beta)};
    const int max_coi = static_cast<int>(params.max_coi);
    double max_abs = 0.0;

    for (std::size_t s = 0; s < genotyping_data.num_samples; ++s) {
        const double e = static_cast<double>(eff_coi.at({s}));
        const double assembled = ecoi_sample_marginal_transmission_llik(s, e);

        const int m_lo = std::max(2, static_cast<int>(std::ceil(e)));
        std::vector<double> per_pop(params.num_populations, ecoi_marginal::kNegInf);
        for (std::size_t pop = 0; pop < params.num_populations; ++pop) {
            std::vector<double> terms;
            for (int mm = m_lo; mm <= max_coi; ++mm) {
                const double rr = ecoi_marginal::r_of(mm, e);
                if (rr <= 0.0 || rr >= 1.0) {
                    terms.push_back(ecoi_marginal::kNegInf);
                    continue;
                }
                double tx = 0.0;
                bool ok = true;
                for (std::size_t l = 0; l < genotyping_data.num_loci; ++l) {
                    if (genotyping_data.is_missing(s, l)) {
                        continue;
                    }
                    const auto support = latent_allele_support(s, l);
                    if (support.empty()) {
                        continue;
                    }
                    const auto [pb, pe] = p.inner_iterators({pop, l});
                    const float prod = calc_transmission_process(
                        support, std::span<const float>(pb, pe), mm,
                        static_cast<float>(rr));
                    if (!std::isfinite(prod)) {
                        ok = false;
                        break;
                    }
                    tx += static_cast<double>(prod);
                }
                const double logw = ecoi_marginal::dztnbinom_log(mm, hp.coi_p, hp.coi_r) +
                                    ecoi_marginal::dbeta_log(rr, hp.r_alpha, hp.r_beta) -
                                    std::log(static_cast<double>(mm - 1));
                terms.push_back(ok ? logw + tx : ecoi_marginal::kNegInf);
            }
            const double pop_log = std::log(
                static_cast<double>(population_responsibility_vector.at({pop})));
            per_pop[pop] = pop_log + ecoi_marginal::log_sum_exp(
                                         std::span<const double>(terms));
        }
        const double brute = ecoi_marginal::log_sum_exp(std::span<const double>(per_pop));

        const bool fa = std::isfinite(assembled);
        const bool fb = std::isfinite(brute);
        if (fa && fb) {
            max_abs = std::max(max_abs, std::fabs(assembled - brute));
        } else if (fa != fb) {
            max_abs = std::numeric_limits<double>::infinity();
        }
    }
    return max_abs;
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

void Chain::save_eps_pos_locus_likelihood(std::size_t locus_idx)
{
    eps_pos_locus_prior_old.at({locus_idx}) = eps_pos_locus_prior_new.at({locus_idx});
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

void Chain::restore_eps_pos_locus_likelihood(std::size_t locus_idx)
{
    const float old_val = eps_pos_locus_prior_old.at({locus_idx});
    eps_pos_locus_prior_sum_new_ += old_val - eps_pos_locus_prior_new.at({locus_idx});
    eps_pos_locus_prior_new.at({locus_idx}) = old_val;
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
