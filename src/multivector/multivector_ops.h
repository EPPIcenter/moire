#pragma once

#include "parallel_backend.h"

#include <algorithm>
#include <iterator>
#include <numeric>
#include <utility>
#include <vector>

#ifndef MOIRE_PARALLEL_ELEMENT_THRESHOLD
#define MOIRE_PARALLEL_ELEMENT_THRESHOLD 100000
#endif
#ifndef MOIRE_PARALLEL_SLICE_THRESHOLD
#define MOIRE_PARALLEL_SLICE_THRESHOLD 10000
#endif
#ifndef MOIRE_PARALLEL_FOR_THRESHOLD
#define MOIRE_PARALLEL_FOR_THRESHOLD 1000
#endif
#ifndef MOIRE_PARALLEL_FOR_2D_THRESHOLD
#define MOIRE_PARALLEL_FOR_2D_THRESHOLD 500
#endif
#ifndef MOIRE_PARALLEL_SLICE_UNARY_INNER_THRESHOLD
#define MOIRE_PARALLEL_SLICE_UNARY_INNER_THRESHOLD 256
#endif

constexpr bool should_parallelize(size_t workload_size, size_t threshold) {
    return workload_size >= threshold;
}

namespace moire_parallel {
    thread_local inline bool disable_nested_parallelism = false;

    template<typename Func>
    inline void parallel_for(size_t begin, size_t end, Func&& func) {
        const size_t workload = end - begin;
        if (disable_nested_parallelism || !should_parallelize(workload, MOIRE_PARALLEL_FOR_THRESHOLD)) {
            for (size_t i = begin; i < end; ++i) func(i);
            return;
        }
        tbb::parallel_for(tbb::blocked_range<size_t>(begin, end),
            [&](const tbb::blocked_range<size_t>& range) {
                for (size_t i = range.begin(); i != range.end(); ++i) func(i);
            });
    }

    template<typename Func>
    inline void parallel_for_always(size_t begin, size_t end, Func&& func) {
        const size_t workload = end - begin;
        if (workload <= 1) {
            if (begin < end) func(begin);
            return;
        }
        if (disable_nested_parallelism) {
            for (size_t i = begin; i < end; ++i) func(i);
            return;
        }
        tbb::parallel_for(tbb::blocked_range<size_t>(begin, end),
            [&](const tbb::blocked_range<size_t>& range) {
                for (size_t i = range.begin(); i != range.end(); ++i) func(i);
            });
    }

    template<typename Func>
    inline void parallel_for_2d(size_t dim0_begin, size_t dim0_end,
                                size_t dim1_begin, size_t dim1_end, Func&& func) {
        const size_t total_work = (dim0_end - dim0_begin) * (dim1_end - dim1_begin);
        if (disable_nested_parallelism || !should_parallelize(total_work, MOIRE_PARALLEL_FOR_2D_THRESHOLD)) {
            for (size_t i = dim0_begin; i < dim0_end; ++i)
                for (size_t j = dim1_begin; j < dim1_end; ++j) func(i, j);
            return;
        }
        tbb::parallel_for(
            tbb::blocked_range2d<size_t>(dim0_begin, dim0_end, dim1_begin, dim1_end),
            [&](const tbb::blocked_range2d<size_t>& range) {
                for (size_t i = range.rows().begin(); i != range.rows().end(); ++i)
                    for (size_t j = range.cols().begin(); j != range.cols().end(); ++j)
                        func(i, j);
            });
    }

    template<typename Func>
    inline void parallel_for_2d_always(size_t dim0_begin, size_t dim0_end,
                                       size_t dim1_begin, size_t dim1_end, Func&& func) {
        const size_t total_work = (dim0_end - dim0_begin) * (dim1_end - dim1_begin);
        if (disable_nested_parallelism || total_work <= 1) {
            for (size_t i = dim0_begin; i < dim0_end; ++i)
                for (size_t j = dim1_begin; j < dim1_end; ++j) func(i, j);
            return;
        }
        tbb::parallel_for(
            tbb::blocked_range2d<size_t>(dim0_begin, dim0_end, dim1_begin, dim1_end),
            [&](const tbb::blocked_range2d<size_t>& range) {
                for (size_t i = range.rows().begin(); i != range.rows().end(); ++i)
                    for (size_t j = range.cols().begin(); j != range.cols().end(); ++j)
                        func(i, j);
            });
    }

    // Per-sample likelihood recomputes (transmission / observation) are small
    // (O(loci) or O(pop*loci)) but run thousands of times per MCMC iteration.
    // When only one chain is active the worker cores are otherwise idle, so it
    // pays to parallelize these below the generic threshold; with >1 chain these
    // fall back to serial because nested parallelism is disabled. These thin
    // wrappers name that intent at the call sites.
    template<typename Func>
    inline void recalc_parallel_for(size_t begin, size_t end, Func&& func) {
        parallel_for_always(begin, end, std::forward<Func>(func));
    }

    template<typename Func>
    inline void recalc_parallel_for_2d(size_t dim0_begin, size_t dim0_end,
                                       size_t dim1_begin, size_t dim1_end, Func&& func) {
        parallel_for_2d_always(dim0_begin, dim0_end, dim1_begin, dim1_end,
                               std::forward<Func>(func));
    }

    template<typename InputIt, typename OutputIt, typename UnaryOp>
    inline void transform(InputIt first, InputIt last, OutputIt result, UnaryOp op) {
        const size_t size = static_cast<size_t>(std::distance(first, last));
        if (!should_parallelize(size, MOIRE_PARALLEL_ELEMENT_THRESHOLD)) {
            std::transform(first, last, result, op);
            return;
        }
        tbb::parallel_for(tbb::blocked_range<size_t>(0, size),
            [&](const tbb::blocked_range<size_t>& range) {
                auto it_in = first + range.begin();
                auto it_out = result + range.begin();
                for (size_t i = range.begin(); i != range.end(); ++i, ++it_in, ++it_out)
                    *it_out = op(*it_in);
            });
    }

    template<typename InputIt1, typename InputIt2, typename OutputIt, typename BinaryOp>
    inline void transform(InputIt1 first1, InputIt1 last1, InputIt2 first2, OutputIt result, BinaryOp op) {
        const size_t size = static_cast<size_t>(std::distance(first1, last1));
        if (!should_parallelize(size, MOIRE_PARALLEL_ELEMENT_THRESHOLD)) {
            std::transform(first1, last1, first2, result, op);
            return;
        }
        tbb::parallel_for(tbb::blocked_range<size_t>(0, size),
            [&](const tbb::blocked_range<size_t>& range) {
                auto it1 = first1 + range.begin();
                auto it2 = first2 + range.begin();
                auto it_out = result + range.begin();
                for (size_t i = range.begin(); i != range.end(); ++i, ++it1, ++it2, ++it_out)
                    *it_out = op(*it1, *it2);
            });
    }

    template<typename InputIt, typename T, typename BinaryOp>
    inline T reduce(InputIt first, InputIt last, T init, BinaryOp op) {
        const size_t size = static_cast<size_t>(std::distance(first, last));
        if (!should_parallelize(size, MOIRE_PARALLEL_ELEMENT_THRESHOLD)) {
            return std::reduce(first, last, init, op);
        }
        return tbb::parallel_reduce(
            tbb::blocked_range<InputIt>(first, last),
            init,
            [&](const tbb::blocked_range<InputIt>& range, T local_init) {
                return std::reduce(range.begin(), range.end(), local_init, op);
            },
            op);
    }

    template<typename InputIt, typename T, typename BinaryOp, typename UnaryOp>
    inline T transform_reduce(InputIt first, InputIt last, T init, BinaryOp reduce_op, UnaryOp transform_op) {
        const size_t size = static_cast<size_t>(std::distance(first, last));
        if (!should_parallelize(size, MOIRE_PARALLEL_ELEMENT_THRESHOLD)) {
            return std::transform_reduce(first, last, init, reduce_op, transform_op);
        }
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
    }

    template<typename Func>
    inline void for_each_slice_1d(size_t n, size_t threshold, Func&& func) {
        if (!should_parallelize(n, threshold)) {
            for (size_t i = 0; i < n; ++i) func(i);
            return;
        }
        parallel_for(size_t(0), n, std::forward<Func>(func));
    }

    template<typename Func>
    inline void for_each_slice_2d(size_t dim0, size_t dim1, size_t threshold, Func&& func) {
        const size_t total = dim0 * dim1;
        if (!should_parallelize(total, threshold)) {
            for (size_t i = 0; i < dim0; ++i)
                for (size_t j = 0; j < dim1; ++j) func(i, j);
            return;
        }
        parallel_for_2d(size_t(0), dim0, size_t(0), dim1, std::forward<Func>(func));
    }

    template<typename It, typename T, typename BinaryOp>
    inline T reduce_range(It first, It last, T init, BinaryOp op, size_t size_hint, size_t threshold) {
        if (!should_parallelize(size_hint, threshold)) {
            return std::reduce(first, last, init, op);
        }
        return moire_parallel::reduce<It, T, BinaryOp>(first, last, init, op);
    }
}
