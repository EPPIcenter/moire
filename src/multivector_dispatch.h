#pragma once

// Single source of truth for parallel/sequential backend dispatch.
// Include multivector_ops.h before using; this header only adds dispatch helpers
// used by multivector algorithms (for_each_slice_1d, for_each_slice_2d, reduce_range).

#include <execution>
#include <numeric>

// Requires: multivector_ops.h already included (defines should_parallelize, moire_parallel::*; backend in parallel_backend.h).

namespace moire_dispatch {

/// Run a loop over 1D slice indices [0, n). Uses TBB or std::execution::par when enabled and threshold met; else sequential.
template<typename Func>
inline void for_each_slice_1d(size_t n, size_t threshold, Func&& func) {
    if (!should_parallelize(n, threshold)) {
        for (size_t i = 0; i < n; ++i) func(i);
        return;
    }
    moire_parallel::parallel_for(size_t(0), n, std::forward<Func>(func));
}

/// Run a loop over 2D slice indices [0,dim0) x [0,dim1). Uses TBB 2D or std::execution when enabled and threshold met; else sequential.
template<typename Func>
inline void for_each_slice_2d(size_t dim0, size_t dim1, size_t threshold, Func&& func) {
    const size_t total = dim0 * dim1;
    if (!should_parallelize(total, threshold)) {
        for (size_t i = 0; i < dim0; ++i)
            for (size_t j = 0; j < dim1; ++j) func(i, j);
        return;
    }
    moire_parallel::parallel_for_2d(size_t(0), dim0, size_t(0), dim1, std::forward<Func>(func));
}

/// Reduce range [first, last) with init and binary_op. Uses parallel backend when MOIRE_ENABLE_PARALLEL and threshold met.
template<typename It, typename T, typename BinaryOp>
inline T reduce_range(It first, It last, T init, BinaryOp op, size_t size_hint, size_t threshold) {
    if (!should_parallelize(size_hint, threshold)) {
        return std::reduce(std::execution::seq, first, last, init, op);
    }
    return moire_parallel::reduce(first, last, init, op);
}

/// Reduce range [first, last) with init and binary_op using given execution policy (no threshold).
template<typename It, typename T, typename BinaryOp, typename ExecutionPolicy>
inline T reduce_range(It first, It last, T init, BinaryOp op, ExecutionPolicy policy) {
    return std::reduce(policy, first, last, init, op);
}

} // namespace moire_dispatch
