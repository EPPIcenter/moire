/**
 * Repeatable regression-oriented MultiVector benchmarks.
 * Uses fixed seed (42) for deterministic data. Emit CSV when MOIRE_BENCHMARK_CSV=1
 * for CI trend tracking and before/after comparisons.
 */
#include "multivector.h"
#include "benchmark_framework.hpp"
#include <array>
#include <chrono>
#include <cmath>
#include <execution>
#include <functional>
#include <iomanip>
#include <iostream>
#include <random>
#include <vector>

namespace {

constexpr unsigned SEED = 42;
constexpr int WARMUP = 2;
constexpr int ITERS = 5;

template<typename T>
void fill_random(MultiVector<T, 3>& mv, unsigned seed) {
    std::mt19937 gen(seed);
    std::uniform_real_distribution<T> dis(0.01, 2.0);
    const auto d = mv.dimensions();
    for (size_t i = 0; i < d[0]; ++i)
        for (size_t j = 0; j < d[1]; ++j)
            for (size_t k = 0; k < d[2]; ++k)
                mv.at({i, j, k}) = dis(gen);
}

double run_timed(int warmup, int iters, std::function<void()> fn) {
    for (int i = 0; i < warmup; ++i) fn();
    std::vector<double> times;
    times.reserve(static_cast<size_t>(iters));
    for (int i = 0; i < iters; ++i) {
        auto start = std::chrono::high_resolution_clock::now();
        fn();
        auto end = std::chrono::high_resolution_clock::now();
        times.push_back(std::chrono::duration<double, std::milli>(end - start).count());
    }
    double sum = 0;
    for (double t : times) sum += t;
    return sum / times.size();
}

void run_scenario(const std::string& scenario, size_t d0, size_t d1, size_t d2,
                  const std::string& operation, const std::string& policy,
                  std::function<void()> fn) {
    double time_ms = run_timed(WARMUP, ITERS, std::move(fn));
    if (moire_bench::output_csv()) {
        moire_bench::RegressionResult r;
        r.scenario = scenario;
        r.dim0 = d0; r.dim1 = d1; r.dim2 = d2;
        r.operation = operation;
        r.policy = policy;
        r.time_ms = time_ms;
        moire_bench::emit_csv_row(std::cout, r);
    } else {
        std::cout << scenario << " " << d0 << "x" << d1 << "x" << d2
                  << " " << operation << " " << policy << ": "
                  << std::fixed << std::setprecision(3) << time_ms << " ms\n";
    }
}

} // namespace

int main(int argc, char** argv) {
    const bool csv = moire_bench::output_csv();
    if (csv) {
        moire_bench::emit_csv_header(std::cout);
    } else {
        std::cout << "MultiVector regression benchmarks (seed=" << SEED
                  << ", warmup=" << WARMUP << ", iters=" << ITERS << ")\n";
    }

    const std::vector<std::array<size_t, 3>> sizes = {
        {10, 10, 10},
        {25, 25, 25},
        {50, 50, 50},
        {100, 100, 100},
    };

    for (const auto& dims : sizes) {
        const size_t d0 = dims[0], d1 = dims[1], d2 = dims[2];
        std::array<size_t, 3> dims_arr = {d0, d1, d2};
        MultiVector<double, 3> mv(dims_arr);
        fill_random(mv, SEED);

        run_scenario("size_sweep", d0, d1, d2, "reduce_sum", "seq",
                     [&]() { (void)mv.reduce(std::plus<double>(), 0.0, std::execution::seq); });
        run_scenario("size_sweep", d0, d1, d2, "reduce_sum", "parallel",
                     [&]() { (void)mv.parallel_reduce(std::plus<double>(), 0.0); });
        run_scenario("size_sweep", d0, d1, d2, "sum", "seq",
                     [&]() { (void)mv.sum(std::execution::seq); });
        run_scenario("size_sweep", d0, d1, d2, "sum", "parallel",
                     [&]() { (void)mv.parallel_sum(); });
        run_scenario("size_sweep", d0, d1, d2, "logsumexp", "seq",
                     [&]() { (void)mv.logsumexp(std::execution::seq); });
        run_scenario("size_sweep", d0, d1, d2, "logsumexp", "parallel",
                     [&]() { (void)mv.parallel_logsumexp(); });
        run_scenario("size_sweep", d0, d1, d2, "transform_sqrt", "seq",
                     [&]() { (void)mv.transform([](double x) { return std::sqrt(x); }, std::execution::seq); });
        run_scenario("size_sweep", d0, d1, d2, "transform_sqrt", "parallel",
                     [&]() { (void)mv.parallel_sqrt(); });
    }

    return 0;
}
