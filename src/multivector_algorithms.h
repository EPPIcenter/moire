#pragma once

// Reusable slice-level kernels for MultiVector algorithms.
// Used by multivector.h; requires <execution>, <numeric>, <cmath>.

#include <execution>
#include <numeric>
#include <cmath>

namespace moire_kernels {

/// Reduce [first, last) with init and binary_op.
template<typename It, typename T, typename BinaryOp>
inline T reduce_slice(It first, It last, const T& init, BinaryOp op, std::execution::sequenced_policy) {
    return std::reduce(std::execution::seq, first, last, init, op);
}

template<typename It, typename T, typename BinaryOp, typename ExecutionPolicy>
inline T reduce_slice(It first, It last, const T& init, BinaryOp op, ExecutionPolicy policy) {
    return std::reduce(policy, first, last, init, op);
}

/// Log-sum-exp over [first, last) for numerical stability.
template<typename It>
inline auto logsumexp_slice(It first, It last) -> typename std::iterator_traits<It>::value_type {
    using T = typename std::iterator_traits<It>::value_type;
    if (first == last) return T{};
    T max_val = *std::max_element(std::execution::seq, first, last);
    T sum = T{};
    for (It it = first; it != last; ++it)
        sum += std::exp(*it - max_val);
    return std::log(sum) + max_val;
}

/// Softmax (or log-softmax) over [first, last), write to result_first.
/// If is_log_values, write log(softmax(x)) = x - logsumexp(x); else write softmax(x).
template<typename ItIn, typename ItOut>
inline void softmax_slice(ItIn first, ItIn last, ItOut result_first, bool is_log_values) {
    using T = typename std::iterator_traits<ItIn>::value_type;
    if (first == last) return;
    T max_val = *std::max_element(std::execution::seq, first, last);
    T sum = T{};
    for (ItIn it = first; it != last; ++it)
        sum += std::exp(*it - max_val);
    T log_sum = std::log(sum) + max_val;
    if (is_log_values) {
        for (; first != last; ++first, ++result_first)
            *result_first = *first - log_sum;
    } else {
        for (; first != last; ++first, ++result_first)
            *result_first = std::exp(*first - max_val) / sum;
    }
}

/// Unary transform: result[i] = op(input[i]) for [first, last) -> result_first.
template<typename ItIn, typename ItOut, typename UnaryOp>
inline void transform_slice_unary(ItIn first, ItIn last, ItOut result_first, UnaryOp op, std::execution::sequenced_policy) {
    std::transform(std::execution::seq, first, last, result_first, op);
}

template<typename ItIn, typename ItOut, typename UnaryOp, typename ExecutionPolicy>
inline void transform_slice_unary(ItIn first, ItIn last, ItOut result_first, UnaryOp op, ExecutionPolicy policy) {
    std::transform(policy, first, last, result_first, op);
}

/// Binary transform with scalar: result[i] = op(input[i], value).
template<typename ItIn, typename ItOut, typename T, typename BinaryOp>
inline void transform_slice_binary(ItIn first, ItIn last, ItOut result_first, const T& value, BinaryOp op) {
    std::transform(std::execution::seq, first, last, result_first, [&](const T& x) { return op(x, value); });
}

} // namespace moire_kernels
