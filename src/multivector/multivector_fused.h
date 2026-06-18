#pragma once

// Fused MultiVector kernels for MCMC hot paths (no intermediate allocations).
// Layout: transmission[sample, population, locus], coi_prior[sample, population],
// population_log[population].

#include "multivector.h"
#include "multivector_algorithms.h"
#include "multivector_ops.h"

#include <array>
#include <cmath>
#include <limits>
#include <span>
#include <vector>

namespace moire_fused {

namespace detail {

template<typename T>
inline T sum_loci_slice(
    const T* data, size_t start_idx, size_t dim2) {
    return moire_kernels::reduce_slice(
        data + start_idx, data + start_idx + dim2, T{}, std::plus<T>());
}

template<typename T>
inline T population_logits_logsumexp(
    const T* logits, size_t dim1) {
    return moire_kernels::logsumexp_slice(logits, logits + dim1);
}

template<typename T>
inline void population_logits_log_softmax(
    const T* logits, size_t dim1, T* out) {
    if (dim1 == 0) return;
    const T log_norm = moire_kernels::logsumexp_slice(logits, logits + dim1);
    for (size_t j = 0; j < dim1; ++j) {
        out[j] = logits[j] - log_norm;
    }
}

template<typename T>
inline void fill_population_logits(
    const MultiVector<T, 3>& transmission,
    const MultiVector<T, 2>& coi_prior,
    std::span<const T> population_log,
    size_t sample_idx,
    T* logits) {
    const size_t dim1 = transmission.dimensions()[1];
    const size_t dim2 = transmission.dimensions()[2];
    const size_t stride0 = transmission.strides()[0];
    const size_t stride1 = transmission.strides()[1];
    const size_t coi_stride0 = coi_prior.strides()[0];
    const T* tx = transmission.data().data();
    const T* coi = coi_prior.data().data();

    for (size_t j = 0; j < dim1; ++j) {
        const size_t start_idx = sample_idx * stride0 + j * stride1;
        logits[j] = sum_loci_slice<T>(tx, start_idx, dim2)
            + coi[sample_idx * coi_stride0 + j]
            + population_log[j];
    }
}

} // namespace detail

/// Update logsumexp when a single logit changes from a_old to a_new.
template<typename T>
inline T logsumexp_replace_one(T L, T a_old, T a_new) {
    constexpr T neg_inf = -std::numeric_limits<T>::infinity();
    if (a_old == a_new) {
        return L;
    }
    if (!std::isfinite(L)) {
        return a_new;
    }
    const T m = std::max({L, a_old, a_new});
    const T sum = std::exp(L - m) - std::exp(a_old - m) + std::exp(a_new - m);
    if (sum <= T{}) {
        return neg_inf;
    }
    return m + std::log(sum);
}

/// sum over loci, add coi_prior, broadcast-add population_log, logsumexp over
/// populations, sum over samples — zero intermediate MultiVector allocations.
template<typename T>
inline T transmission_loglikelihood_sum(
    const MultiVector<T, 3>& transmission,
    const MultiVector<T, 2>& coi_prior,
    const MultiVector<T, 1>& population_log) {
    const size_t dim0 = transmission.dimensions()[0];
    const size_t dim1 = transmission.dimensions()[1];
    if (population_log.dimensions()[0] != dim1) {
        throw std::invalid_argument("population_log size must match population dimension");
    }

    std::vector<T> logits(dim1);
    const auto pop_log = population_log.as_span();
    T total = T{};
    for (size_t i = 0; i < dim0; ++i) {
        detail::fill_population_logits(transmission, coi_prior, pop_log, i, logits.data());
        total += detail::population_logits_logsumexp(logits.data(), dim1);
    }
    return total;
}

/// Fused sum over loci + coi + broadcast population_log + log-softmax over populations.
/// One output allocation instead of three intermediate MultiVectors.
template<typename T>
inline MultiVector<T, 2> population_assignment_log_softmax(
    const MultiVector<T, 3>& transmission,
    const MultiVector<T, 2>& coi_prior,
    const MultiVector<T, 1>& population_log) {
    const size_t dim0 = transmission.dimensions()[0];
    const size_t dim1 = transmission.dimensions()[1];
    if (population_log.dimensions()[0] != dim1) {
        throw std::invalid_argument("population_log size must match population dimension");
    }

    MultiVector<T, 2> result({dim0, dim1});
    std::vector<T> logits(dim1);
    const auto pop_log = population_log.as_span();
    for (size_t i = 0; i < dim0; ++i) {
        detail::fill_population_logits(transmission, coi_prior, pop_log, i, logits.data());
        const size_t out_start = i * result.strides()[0];
        detail::population_logits_log_softmax(
            logits.data(), dim1, result.data_.data() + out_start);
    }
    return result;
}

} // namespace moire_fused
