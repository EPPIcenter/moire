#pragma once

#include "prob_any_missing_cache.h"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <span>
#include <string>
#include <vector>

namespace pam_fast_paths {

inline constexpr std::size_t kLowKMaxSupport = 5;

namespace detail {

inline bool read_tx_opts_from_env() noexcept
{
    const char* env = std::getenv("MOIRE_DISABLE_TX_OPT");
    return !(env && (std::string(env) == "1" || std::string(env) == "true" ||
                       std::string(env) == "yes"));
}

inline bool& tx_opts_flag() noexcept
{
    // Resolved once on first use; refreshed explicitly via refresh_tx_opts_from_env().
    static bool flag = read_tx_opts_from_env();
    return flag;
}

}  // namespace detail

/// Re-read MOIRE_DISABLE_TX_OPT from the environment. Call once per run from the
/// main thread (before any parallel work) so the hot path can read a cached bool
/// instead of hitting getenv + string allocation on every transmission eval.
inline void refresh_tx_opts_from_env() noexcept
{
    detail::tx_opts_flag() = detail::read_tx_opts_from_env();
}

inline bool tx_opts_enabled() noexcept
{
    return detail::tx_opts_flag();
}

/// P(all categories in support seen at least once in exactly n draws), q normalized on support.
inline double prob_all_seen_after_n(std::span<const float> q, unsigned n)
{
    const std::size_t k = q.size();
    if (k == 0 || n == 0) {
        return 0.0;
    }
    if (k == 1) {
        return 1.0;
    }

    double prob_all_seen = 0.0;
    const std::size_t num_masks = static_cast<std::size_t>(1) << k;
    for (std::size_t excluded = 0; excluded < num_masks; ++excluded) {
        double subset_sum = 0.0;
        for (std::size_t i = 0; i < k; ++i) {
            if ((excluded & (static_cast<std::size_t>(1) << i)) == 0) {
                subset_sum += static_cast<double>(q[i]);
            }
        }
        const double term = std::pow(subset_sum, static_cast<double>(n));
        const int popcount = std::popcount(static_cast<unsigned>(excluded));
        if (popcount & 1) {
            prob_all_seen -= term;
        } else {
            prob_all_seen += term;
        }
    }
    return prob_all_seen;
}

inline void fill_pam_vector_low_k(std::span<const float> q,
                                  unsigned min_events,
                                  unsigned max_events,
                                  std::vector<double>& out)
{
    const std::size_t k = q.size();
    const std::size_t len = static_cast<std::size_t>(max_events - min_events + 1);
    out.assign(len, 0.0);

    if (max_events < k) {
        std::fill(out.begin(), out.end(), 1.0);
        return;
    }

    for (std::size_t j = 0; j < k - 1 && j < len; ++j) {
        out[j] = 1.0;
    }

    for (unsigned n = static_cast<unsigned>(k); n <= max_events; ++n) {
        const std::size_t idx = static_cast<std::size_t>(n - min_events);
        if (idx >= len) {
            break;
        }
        const double all_seen = prob_all_seen_after_n(q, n);
        out[idx] = 1.0 - all_seen;
    }
}

inline bool try_fill_pam_vector(std::span<const float> q,
                                unsigned min_events,
                                unsigned max_events,
                                std::vector<double>& out)
{
    if (!tx_opts_enabled() || q.size() > kLowKMaxSupport) {
        return false;
    }
    fill_pam_vector_low_k(q, min_events, max_events, out);
    return true;
}

} // namespace pam_fast_paths
