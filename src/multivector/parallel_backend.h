#pragma once

// TBB via RcppParallel (R package) or direct TBB headers.

#if defined(__has_include) && __has_include(<RcppParallel.h>)
#include <RcppParallel.h>
#endif

#if defined(__has_include) && __has_include(<tbb/parallel_for.h>) && __has_include(<tbb/blocked_range2d.h>)
#include <tbb/parallel_for.h>
#include <tbb/blocked_range2d.h>
#include <tbb/parallel_reduce.h>
#define TBB_AVAILABLE 1
#endif

#ifndef TBB_AVAILABLE
#error "moire requires Intel TBB (install RcppParallel or libtbb)"
#endif
