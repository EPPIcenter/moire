#pragma once

#include "parallel_backend.h"

// Parallelization configuration and dispatch (included by multivector.h after its includes).
// To enable automatic parallelization in parallel_* methods, define MOIRE_ENABLE_PARALLEL.
// Fallback: TBB -> std::execution::par -> sequential.

#ifndef MOIRE_PARALLEL_ELEMENT_THRESHOLD
#define MOIRE_PARALLEL_ELEMENT_THRESHOLD 100000  // For elementwise operations (e.g., transform, add)
#endif
#ifndef MOIRE_PARALLEL_SLICE_THRESHOLD
#define MOIRE_PARALLEL_SLICE_THRESHOLD 10000     // For reduction operations over slices
#endif
#ifndef MOIRE_PARALLEL_FOR_THRESHOLD
#define MOIRE_PARALLEL_FOR_THRESHOLD 1000        // For parallel_for loops
#endif
#ifndef MOIRE_PARALLEL_FOR_2D_THRESHOLD
#define MOIRE_PARALLEL_FOR_2D_THRESHOLD 500      // For parallel_for_2d loops (populations × loci)
#endif

constexpr bool should_parallelize(size_t workload_size, size_t threshold) {
#ifdef MOIRE_ENABLE_PARALLEL
    return workload_size >= threshold;
#else
    return false;
#endif
}

namespace moire_parallel {
    thread_local inline bool disable_nested_parallelism = false;

    template<typename Func>
    inline void parallel_for(size_t begin, size_t end, Func&& func) {
        const size_t workload = end - begin;
        if (disable_nested_parallelism) {
            for (size_t i = begin; i < end; ++i) func(i);
            return;
        }
        if (!should_parallelize(workload, MOIRE_PARALLEL_FOR_THRESHOLD)) {
            for (size_t i = begin; i < end; ++i) func(i);
            return;
        }
#ifdef MOIRE_ENABLE_PARALLEL
#ifdef TBB_AVAILABLE
        tbb::parallel_for(tbb::blocked_range<size_t>(begin, end),
            [&](const tbb::blocked_range<size_t>& range) {
                for (size_t i = range.begin(); i != range.end(); ++i) func(i);
            });
        return;
#elif defined(HAS_EXECUTION)
        std::vector<size_t> indices(workload);
        std::iota(indices.begin(), indices.end(), begin);
        std::for_each(std::execution::par, indices.begin(), indices.end(), func);
        return;
#endif
#endif
        for (size_t i = begin; i < end; ++i) func(i);
    }

    template<typename Func>
    inline void parallel_for_always(size_t begin, size_t end, Func&& func) {
        const size_t workload = end - begin;
        if (workload <= 1) {
            if (begin < end) func(begin);
            return;
        }
#ifdef MOIRE_ENABLE_PARALLEL
#ifdef TBB_AVAILABLE
        tbb::parallel_for(tbb::blocked_range<size_t>(begin, end),
            [&](const tbb::blocked_range<size_t>& range) {
                for (size_t i = range.begin(); i != range.end(); ++i) func(i);
            });
        return;
#elif defined(HAS_EXECUTION)
        std::vector<size_t> indices(workload);
        std::iota(indices.begin(), indices.end(), begin);
        std::for_each(std::execution::par, indices.begin(), indices.end(), func);
        return;
#endif
#endif
        for (size_t i = begin; i < end; ++i) func(i);
    }

    template<typename Func>
    inline void parallel_for_2d(size_t dim0_begin, size_t dim0_end,
                                size_t dim1_begin, size_t dim1_end, Func&& func) {
        const size_t dim0_size = dim0_end - dim0_begin;
        const size_t dim1_size = dim1_end - dim1_begin;
        const size_t total_work = dim0_size * dim1_size;
        if (disable_nested_parallelism) {
            for (size_t i = dim0_begin; i < dim0_end; ++i)
                for (size_t j = dim1_begin; j < dim1_end; ++j) func(i, j);
            return;
        }
        if (!should_parallelize(total_work, MOIRE_PARALLEL_FOR_2D_THRESHOLD)) {
            for (size_t i = dim0_begin; i < dim0_end; ++i)
                for (size_t j = dim1_begin; j < dim1_end; ++j) func(i, j);
            return;
        }
#ifdef MOIRE_ENABLE_PARALLEL
#ifdef TBB_AVAILABLE
        tbb::parallel_for(
            tbb::blocked_range2d<size_t>(dim0_begin, dim0_end, dim1_begin, dim1_end),
            [&](const tbb::blocked_range2d<size_t>& range) {
                for (size_t i = range.rows().begin(); i != range.rows().end(); ++i)
                    for (size_t j = range.cols().begin(); j != range.cols().end(); ++j)
                        func(i, j);
            });
        return;
#elif defined(HAS_EXECUTION)
        std::vector<std::pair<size_t, size_t>> indices;
        indices.reserve(total_work);
        for (size_t i = dim0_begin; i < dim0_end; ++i)
            for (size_t j = dim1_begin; j < dim1_end; ++j)
                indices.emplace_back(i, j);
        std::for_each(std::execution::par, indices.begin(), indices.end(),
            [&](const std::pair<size_t, size_t>& idx) { func(idx.first, idx.second); });
        return;
#endif
#endif
        for (size_t i = dim0_begin; i < dim0_end; ++i)
            for (size_t j = dim1_begin; j < dim1_end; ++j) func(i, j);
    }

    template<typename InputIt, typename OutputIt, typename UnaryOp>
    inline void transform(InputIt first, InputIt last, OutputIt result, UnaryOp op) {
#ifdef MOIRE_ENABLE_PARALLEL
#ifdef TBB_AVAILABLE
        const size_t size = std::distance(first, last);
        tbb::parallel_for(tbb::blocked_range<size_t>(0, size),
            [&](const tbb::blocked_range<size_t>& range) {
                auto it_in = first + range.begin();
                auto it_out = result + range.begin();
                for (size_t i = range.begin(); i < range.end(); ++i, ++it_in, ++it_out)
                    *it_out = op(*it_in);
            });
        return;
#elif defined(HAS_EXECUTION)
        std::transform(std::execution::par, first, last, result, op);
        return;
#endif
#endif
        std::transform(first, last, result, op);
    }

    template<typename InputIt1, typename InputIt2, typename OutputIt, typename BinaryOp>
    inline void transform(InputIt1 first1, InputIt1 last1, InputIt2 first2, OutputIt result, BinaryOp op) {
#ifdef MOIRE_ENABLE_PARALLEL
#ifdef TBB_AVAILABLE
        const size_t size = std::distance(first1, last1);
        tbb::parallel_for(tbb::blocked_range<size_t>(0, size),
            [&](const tbb::blocked_range<size_t>& range) {
                auto it1 = first1 + range.begin();
                auto it2 = first2 + range.begin();
                auto it_out = result + range.begin();
                for (size_t i = range.begin(); i < range.end(); ++i, ++it1, ++it2, ++it_out)
                    *it_out = op(*it1, *it2);
            });
        return;
#elif defined(HAS_EXECUTION)
        std::transform(std::execution::par, first1, last1, first2, result, op);
        return;
#endif
#endif
        std::transform(first1, last1, first2, result, op);
    }

    template<typename InputIt, typename T, typename BinaryOp>
    inline T reduce(InputIt first, InputIt last, T init, BinaryOp op) {
#ifdef MOIRE_ENABLE_PARALLEL
#ifdef TBB_AVAILABLE
        return tbb::parallel_reduce(
            tbb::blocked_range<InputIt>(first, last),
            init,
            [&](const tbb::blocked_range<InputIt>& range, T local_init) {
                return std::reduce(std::execution::seq, range.begin(), range.end(), local_init, op);
            },
            op);
#elif defined(HAS_EXECUTION)
        return std::reduce(std::execution::par, first, last, init, op);
#else
        return std::reduce(std::execution::seq, first, last, init, op);
#endif
#else
        return std::reduce(std::execution::seq, first, last, init, op);
#endif
    }

    template<typename InputIt, typename T, typename BinaryOp, typename UnaryOp>
    inline T transform_reduce(InputIt first, InputIt last, T init, BinaryOp reduce_op, UnaryOp transform_op) {
#ifdef MOIRE_ENABLE_PARALLEL
#ifdef TBB_AVAILABLE
        return tbb::parallel_reduce(
            tbb::blocked_range<InputIt>(first, last),
            init,
            [&](const tbb::blocked_range<InputIt>& range, T local_init) {
                for (auto it = range.begin(); it != range.end(); ++it) {
                    local_init = reduce_op(std::move(local_init), transform_op(*it));
                }
                return local_init;
            },
            reduce_op);
#elif defined(HAS_EXECUTION)
        return std::transform_reduce(std::execution::par, first, last, init, reduce_op, transform_op);
#else
        return std::transform_reduce(std::execution::seq, first, last, init, reduce_op, transform_op);
#endif
#else
        return std::transform_reduce(std::execution::seq, first, last, init, reduce_op, transform_op);
#endif
    }

    /// Reduce with unseq if available (vectorization), else seq.
    template<typename InputIt, typename T, typename BinaryOp>
    inline T reduce_unseq_or_seq(InputIt first, InputIt last, T init, BinaryOp op) {
#if defined(HAS_EXECUTION)
        return std::reduce(std::execution::unseq, first, last, init, op);
#else
        return std::reduce(std::execution::seq, first, last, init, op);
#endif
    }

    /// Transform-reduce with unseq if available, else seq. Use for e.g. logSumExp.
    template<typename InputIt, typename T, typename BinaryOp, typename UnaryOp>
    inline T transform_reduce_unseq_or_seq(InputIt first, InputIt last, T init, BinaryOp reduce_op, UnaryOp transform_op) {
#if defined(HAS_EXECUTION)
        return std::transform_reduce(std::execution::unseq, first, last, init, reduce_op, transform_op);
#else
        return std::transform_reduce(std::execution::seq, first, last, init, reduce_op, transform_op);
#endif
    }
}
