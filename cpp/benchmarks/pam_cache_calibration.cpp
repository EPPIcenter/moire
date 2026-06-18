/**
 * Calibrate P(any missing) cache quantization: sweep delta and report
 * max |pamVec(q) - pamVec(quantize(q))| and log-likelihood proxy error.
 *
 * Usage: pam_cache_calibration [trials_per_cell]
 * Env: none required
 */
#include "prob_any_missing.h"
#include "prob_any_missing_cache.h"

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <random>
#include <vector>

namespace {

constexpr unsigned SEED = 42;

std::vector<float> random_simplex(std::size_t k, std::mt19937& gen)
{
    std::uniform_real_distribution<float> dis(0.01f, 1.f);
    std::vector<float> q(k);
    float sum = 0.f;
    for (std::size_t i = 0; i < k; ++i) {
        q[i] = dis(gen);
        sum += q[i];
    }
    for (float& v : q) {
        v /= sum;
    }
    return q;
}

double max_pam_diff(const std::vector<double>& a, const std::vector<double>& b)
{
    const std::size_t n = std::min(a.size(), b.size());
    double m = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        m = std::max(m, std::fabs(a[i] - b[i]));
    }
    return m;
}

double max_log_proxy_diff(const std::vector<double>& a,
                          const std::vector<double>& b,
                          unsigned coi)
{
    if (coi == 0 || a.size() < coi || b.size() < coi) {
        return 0.0;
    }
    const std::size_t idx = static_cast<std::size_t>(coi - 1);
    const double log_a = std::log1p(-a[idx]);
    const double log_b = std::log1p(-b[idx]);
    return std::fabs(log_a - log_b);
}

void calibrate_delta(float delta,
                     int trials_per_cell,
                     std::mt19937& gen,
                     double& out_max_pam,
                     double& out_max_log,
                     std::uint64_t& out_pairs)
{
    probAnyMissingFunctor functor;
    out_max_pam = 0.0;
    out_max_log = 0.0;
    out_pairs = 0;

    const std::vector<std::size_t> k_values = {1, 2, 3, 4, 5};
    const std::vector<unsigned> coi_values = {1, 2, 3, 5, 8, 10, 15};
    std::uniform_real_distribution<float> jitter(-0.5f, 0.5f);

    for (std::size_t k : k_values) {
        for (unsigned coi : coi_values) {
            if (coi < k) {
                continue;
            }
            for (int t = 0; t < trials_per_cell; ++t) {
                const std::vector<float> q_stored = random_simplex(k, gen);
                const std::vector<float> key =
                    pam_cache::quantize_simplex(q_stored, delta);

                std::vector<float> q_query = q_stored;
                for (int attempt = 0; attempt < 16; ++attempt) {
                    q_query = q_stored;
                    float sum = 0.f;
                    for (std::size_t i = 0; i < k; ++i) {
                        q_query[i] = std::max(
                            0.f, q_stored[i] + jitter(gen) * delta);
                        sum += q_query[i];
                    }
                    if (sum <= 0.f) {
                        continue;
                    }
                    for (float& v : q_query) {
                        v /= sum;
                    }
                    const std::vector<float> q_key =
                        pam_cache::quantize_simplex(q_query, delta);
                    if (q_key.size() == key.size() &&
                        std::memcmp(q_key.data(), key.data(),
                                    key.size() * sizeof(float)) == 0 &&
                        std::memcmp(q_query.data(), q_stored.data(),
                                    k * sizeof(float)) != 0) {
                        break;
                    }
                }

                if (std::memcmp(q_query.data(), q_stored.data(),
                                k * sizeof(float)) == 0) {
                    continue;
                }
                const std::vector<float> q_key =
                    pam_cache::quantize_simplex(q_query, delta);
                if (q_key.size() != key.size() ||
                    std::memcmp(q_key.data(), key.data(),
                                key.size() * sizeof(float)) != 0) {
                    continue;
                }

                ++out_pairs;
                const auto pam_stored = functor.vectorized(q_stored, 1u, coi);
                const auto pam_query = functor.vectorized(q_query, 1u, coi);
                out_max_pam = std::max(out_max_pam, max_pam_diff(pam_stored, pam_query));
                out_max_log = std::max(
                    out_max_log, max_log_proxy_diff(pam_stored, pam_query, coi));
            }
        }
    }
}

} // namespace

int main(int argc, char** argv)
{
    int trials = 2000;
    if (argc >= 2) {
        trials = std::max(1, std::atoi(argv[1]));
    }

    std::mt19937 gen(SEED);
    const std::vector<float> deltas = {1e-3f, 5e-4f, 1e-4f, 5e-5f, 2e-5f, 1e-5f};

    std::cout << "P(any missing) quantization calibration (trials_per_cell="
              << trials << ", seed=" << SEED << ")\n";
    std::cout << "Measures cache bucket error: |pamVec(q_query) - pamVec(q_stored)|\n";
    std::cout << "Targets: pam_diff < 1e-5, log_proxy_diff < 1e-4 for safe MH\n\n";
    std::cout << std::setw(12) << "delta"
              << std::setw(14) << "max_pam_diff"
              << std::setw(16) << "max_log_diff"
              << std::setw(12) << "pairs"
              << std::setw(8) << "ok"
              << "\n";

    float recommended = 1e-4f;
    bool found = false;
    for (float delta : deltas) {
        double max_pam = 0.0;
        double max_log = 0.0;
        std::uint64_t pairs = 0;
        calibrate_delta(delta, trials, gen, max_pam, max_log, pairs);
        const bool ok = max_pam < 1e-5 && max_log < 1e-4;
        if (ok && !found) {
            recommended = delta;
            found = true;
        }
        std::cout << std::setw(12) << std::scientific << delta
                  << std::setw(14) << max_pam
                  << std::setw(16) << max_log
                  << std::setw(12) << pairs
                  << std::setw(8) << (ok ? "yes" : "no")
                  << "\n";
    }

    std::cout << "\nRecommended delta (largest passing): " << recommended
              << " (set MOIRE_PAM_CACHE_DELTA=" << recommended << ")\n";
    return 0;
}
