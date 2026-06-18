/**
 * Regression-style benchmarks for prob_any_missing (inclusion-exclusion).
 * Sweep over (total_alleles, coi) with fixed seed; optional CSV when
 * MOIRE_BENCHMARK_CSV=1 for before/after and CI.
 */
#include "prob_any_missing.h"
#include "benchmark_framework.hpp"
#include <chrono>
#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <vector>

namespace {

constexpr unsigned SEED = 42;
constexpr int WARMUP = 2;
constexpr int ITERS = 5;
// Calls per timed block; use fewer for large total_alleles so sweep finishes in reasonable time
static int reps_per_iter_for(std::size_t total_alleles) {
    if (total_alleles <= 5) return 100000;
    if (total_alleles <= 10) return 20000;
    if (total_alleles <= 15) return 2000;
    return 200;
}

/// Build normalized event probs of size total_alleles using fixed seed.
std::vector<float> make_event_probs(std::size_t total_alleles, unsigned seed) {
    if (total_alleles == 0) return {};
    std::mt19937 gen(seed);
    std::uniform_real_distribution<float> dis(0.01f, 1.0f);
    std::vector<float> p(total_alleles);
    float sum = 0.0f;
    for (std::size_t i = 0; i < total_alleles; ++i) {
        p[i] = dis(gen);
        sum += p[i];
    }
    if (sum <= 0.0f) { sum = 1.0f; }
    for (std::size_t i = 0; i < total_alleles; ++i)
        p[i] /= sum;
    return p;
}

/// Run warmup + iters timed blocks; each block does REPS_PER_ITER calls. Returns mean time per call (ms).
double run_timed_per_call(int warmup, int iters, int reps_per_iter,
                          std::function<void()> one_block) {
    for (int i = 0; i < warmup; ++i) one_block();
    std::vector<double> block_times_ms;
    block_times_ms.reserve(static_cast<std::size_t>(iters));
    for (int i = 0; i < iters; ++i) {
        auto start = std::chrono::high_resolution_clock::now();
        one_block();
        auto end = std::chrono::high_resolution_clock::now();
        block_times_ms.push_back(
            std::chrono::duration<double, std::milli>(end - start).count());
    }
    double sum = 0.0;
    for (double t : block_times_ms) sum += t;
    double mean_block_ms = sum / static_cast<double>(block_times_ms.size());
    return mean_block_ms / static_cast<double>(reps_per_iter);
}

void emit_prob_any_missing_csv_row(std::ostream& out,
                                    const std::string& scenario,
                                    std::size_t total_alleles, std::size_t coi,
                                    const std::string& operation, double time_ms) {
    out << scenario << "," << total_alleles << "," << coi
        << "," << operation << ","
        << std::fixed << std::setprecision(6) << time_ms << "\n";
}

void run_scenario(const std::string& scenario,
                  std::size_t total_alleles, std::size_t coi,
                  const std::string& operation,
                  std::function<double()> timed_fn) {
    double time_per_call_ms = timed_fn();
    if (moire_bench::output_csv()) {
        emit_prob_any_missing_csv_row(std::cout, scenario,
                                      total_alleles, coi, operation, time_per_call_ms);
    } else {
        std::cout << "prob_any_missing " << total_alleles << " alleles, coi=" << coi
                  << " " << operation << ": "
                  << std::fixed << std::setprecision(6) << time_per_call_ms << " ms/call\n";
    }
}

} // namespace

int main(int argc, char* argv[]) {
    const bool csv = moire_bench::output_csv();

    // Custom run: usage <eventProb1> <eventProb2> ... <numEvents>
    if (argc >= 3) {
        std::vector<float> eventProbs;
        for (int i = 1; i < argc - 1; ++i) {
            float prob;
            std::istringstream(argv[i]) >> prob;
            eventProbs.push_back(prob);
        }
        unsigned int numEvents;
        std::istringstream(argv[argc - 1]) >> numEvents;

        float sum = 0.0f;
        for (float p : eventProbs) sum += p;
        if (sum > 0.0f) {
            for (float& p : eventProbs) p /= sum;
        }

        probAnyMissingFunctor functor;
        const std::size_t total_alleles = eventProbs.size();

        const int reps = 100000;
        if (!csv) {
            std::cout << "prob_any_missing custom (total_alleles=" << total_alleles
                      << ", coi=" << numEvents << ", warmup=" << WARMUP
                      << ", iters=" << ITERS << ", reps=" << reps << ")\n";
        }

        auto time_vectorized = run_timed_per_call(
            WARMUP, ITERS, reps,
            [&]() {
                for (int r = 0; r < reps; ++r)
                    (void)functor.vectorized(eventProbs, numEvents);
            });
        if (csv) {
            emit_prob_any_missing_csv_row(std::cout, "custom", total_alleles,
                                          static_cast<std::size_t>(numEvents),
                                          "vectorized", time_vectorized);
        } else {
            std::cout << "  vectorized: " << std::fixed << std::setprecision(6)
                      << time_vectorized << " ms/call\n";
        }

        auto time_combination = run_timed_per_call(
            WARMUP, ITERS, reps,
            [&]() {
                for (int r = 0; r < reps; ++r)
                    (void)functor.vectorized_combination(eventProbs, numEvents);
            });
        if (csv) {
            emit_prob_any_missing_csv_row(std::cout, "custom", total_alleles,
                                          static_cast<std::size_t>(numEvents),
                                          "vectorized_combination", time_combination);
        } else {
            std::cout << "  vectorized_combination (legacy): " << std::fixed << std::setprecision(6)
                      << time_combination << " ms/call\n";
        }

        auto time_mobius = run_timed_per_call(
            WARMUP, ITERS, reps,
            [&]() {
                for (int r = 0; r < reps; ++r)
                    (void)functor.vectorized_mobius(eventProbs, numEvents);
            });
        if (csv) {
            emit_prob_any_missing_csv_row(std::cout, "custom", total_alleles,
                                          static_cast<std::size_t>(numEvents),
                                          "vectorized_mobius", time_mobius);
        } else {
            std::cout << "  vectorized_mobius: " << std::fixed << std::setprecision(6)
                      << time_mobius << " ms/call\n";
        }

        auto time_operator = run_timed_per_call(
            WARMUP, ITERS, reps,
            [&]() {
                for (int r = 0; r < reps; ++r)
                    (void)functor(eventProbs, static_cast<int>(numEvents));
            });
        if (csv) {
            emit_prob_any_missing_csv_row(std::cout, "custom", total_alleles,
                                          static_cast<std::size_t>(numEvents),
                                          "operator", time_operator);
        } else {
            std::cout << "  operator:   " << std::fixed << std::setprecision(6)
                      << time_operator << " ms/call\n";
        }
        return 0;
    }

    // Default: sweep over (total_alleles, coi)
    if (csv) {
        std::cout << "scenario,total_alleles,coi,operation,time_ms\n";
    } else {
        std::cout << "prob_any_missing regression (seed=" << SEED
                  << ", warmup=" << WARMUP << ", iters=" << ITERS
                  << ", reps_per_iter=adaptive)\n";
    }

    const std::vector<std::size_t> total_alleles_vec = {2, 5, 10};
    const std::vector<std::size_t> coi_vec = {1, 5, 10};

    // Validate Gray-code vs combination (same formula, different iteration order)
    const double tol = 1e-5;
    probAnyMissingFunctor functor_validate;
    for (std::size_t total_alleles : total_alleles_vec) {
        std::vector<float> eventProbs = make_event_probs(total_alleles, SEED);
        for (std::size_t coi : coi_vec) {
            if (coi < total_alleles) continue;
            const unsigned int coi_u = static_cast<unsigned int>(coi);
            auto gray = functor_validate.vectorized(eventProbs, 1u, coi_u);
            auto comb = functor_validate.vectorized_combination(eventProbs, coi_u);
            if (gray.size() != comb.size()) {
                std::cerr << "FAIL: Gray vs combination size mismatch total_alleles=" << total_alleles
                          << " coi=" << coi << " gray=" << gray.size() << " comb=" << comb.size() << "\n";
                return 1;
            }
            for (std::size_t i = 0; i < gray.size(); ++i) {
                const double d = std::fabs(gray[i] - comb[i]);
                if (d > tol) {
                    std::cerr << "FAIL: Gray vs combination diff total_alleles=" << total_alleles
                              << " coi=" << coi << " i=" << i << " gray=" << gray[i]
                              << " comb=" << comb[i] << " |diff|=" << d << " (tol=" << tol << ")\n";
                    return 1;
                }
            }
        }
    }
    if (!csv) std::cout << "Gray vs combination: OK (tol=" << tol << ")\n";

    for (std::size_t total_alleles : total_alleles_vec) {
        std::vector<float> eventProbs = make_event_probs(total_alleles, SEED);
        probAnyMissingFunctor functor;

        for (std::size_t coi : coi_vec) {
            if (coi < total_alleles) continue;
            const unsigned int coi_u = static_cast<unsigned int>(coi);
            const int reps = reps_per_iter_for(total_alleles);

            run_scenario("sweep", total_alleles, coi, "vectorized", [&]() {
                return run_timed_per_call(
                    WARMUP, ITERS, reps,
                    [&]() {
                        for (int r = 0; r < reps; ++r)
                            (void)functor.vectorized(eventProbs, coi_u);
                    });
            });

            run_scenario("sweep", total_alleles, coi, "vectorized_combination", [&]() {
                return run_timed_per_call(
                    WARMUP, ITERS, reps,
                    [&]() {
                        for (int r = 0; r < reps; ++r)
                            (void)functor.vectorized_combination(eventProbs, coi_u);
                    });
            });

            run_scenario("sweep", total_alleles, coi, "vectorized_mobius", [&]() {
                return run_timed_per_call(
                    WARMUP, ITERS, reps,
                    [&]() {
                        for (int r = 0; r < reps; ++r)
                            (void)functor.vectorized_mobius(eventProbs, coi_u);
                    });
            });

            run_scenario("sweep", total_alleles, coi, "operator", [&]() {
                return run_timed_per_call(
                    WARMUP, ITERS, reps,
                    [&]() {
                        for (int r = 0; r < reps; ++r)
                            (void)functor(eventProbs, static_cast<int>(coi));
                    });
            });
        }
    }

    return 0;
}
