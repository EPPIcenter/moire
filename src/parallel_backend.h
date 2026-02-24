#pragma once

// Parallel backend: TBB (via RcppParallel or direct headers) preferred;
// else C++17 std::execution; else sequential.
// All parallel code uses moire_parallel:: in multivector_ops.h only.

// 1) When building as R package, RcppParallel may be linked; we do not set TBB_AVAILABLE here
//    because RcppParallel's TBB may not provide tbb::blocked_range2d (required by parallel_for_2d).
//    TBB_AVAILABLE is only set in (2) when direct TBB headers with blocked_range2d are present.
#if defined(__has_include) && __has_include(<RcppParallel.h>)
#include <RcppParallel.h>
#endif

// 2) Else try direct TBB headers.
#ifndef TBB_AVAILABLE
#if defined(__has_include) && __has_include(<tbb/parallel_for.h>) && __has_include(<tbb/blocked_range2d.h>)
#include <tbb/parallel_for.h>
#include <tbb/blocked_range2d.h>
#if defined(__has_include) && __has_include(<tbb/parallel_reduce.h>)
#include <tbb/parallel_reduce.h>
#endif
#define TBB_AVAILABLE 1
#endif
#endif

// 3) C++17 parallel algorithms. Include <execution> whenever supported (needed for std::execution::seq in TBB and fallbacks).
#if (defined(__cplusplus) && __cplusplus >= 201703L) || (defined(__cpp_lib_execution) && (__cpp_lib_execution >= 201603))
#include <execution>
#if !defined(TBB_AVAILABLE)
#define HAS_EXECUTION 1
#endif
#endif

#if defined(TBB_AVAILABLE) || defined(HAS_EXECUTION)
#define MOIRE_HAVE_PARALLEL 1
#endif
