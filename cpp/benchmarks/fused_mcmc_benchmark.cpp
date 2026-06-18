#include "multivector_fused.h"
#include "benchmark_framework.hpp"
#include <chrono>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <numeric>
#include <random>
#include <vector>

namespace {

constexpr unsigned SEED = 42;

template<typename F>
double bench_ms(int warmup, int iters, F fn) {
    for (int i = 0; i < warmup; ++i) fn();
    std::vector<double> times;
    for (int i = 0; i < iters; ++i) {
        auto t0 = std::chrono::high_resolution_clock::now();
        fn();
        auto t1 = std::chrono::high_resolution_clock::now();
        times.push_back(std::chrono::duration<double, std::milli>(t1 - t0).count());
    }
    return std::accumulate(times.begin(), times.end(), 0.0) / times.size();
}

void fill_random(MultiVector<float, 3>& tx, MultiVector<float, 2>& coi, MultiVector<float, 1>& pop_log, unsigned seed) {
    std::mt19937 gen(seed);
    std::uniform_real_distribution<float> dis(-5.f, 0.f);
    const auto d = tx.dimensions();
    for (size_t i = 0; i < d[0]; ++i)
        for (size_t j = 0; j < d[1]; ++j)
            for (size_t k = 0; k < d[2]; ++k)
                tx.unchecked_at({i, j, k}) = dis(gen);
    for (size_t i = 0; i < d[0]; ++i)
        for (size_t j = 0; j < d[1]; ++j)
            coi.unchecked_at({i, j}) = dis(gen);
    for (size_t j = 0; j < d[1]; ++j)
        pop_log.unchecked_at({j}) = dis(gen);
}

void run_scenario(size_t s, size_t p, size_t l) {
    std::array<size_t, 3> dims = {s, p, l};
    MultiVector<float, 3> tx(dims);
    MultiVector<float, 2> coi({s, p});
    MultiVector<float, 1> pop_log({p});
    fill_random(tx, coi, pop_log, SEED);
    MultiVector<float, 1> pop_log_cached = pop_log.log();

    volatile float sink = 0;
    const double chained = bench_ms(3, 20, [&] {
        sink = (tx.sum() + coi)
            .element_add(pop_log_cached.as_span())
            .logsumexp()
            .full_sum();
    });
    const double fused = bench_ms(3, 20, [&] {
        sink = moire_fused::transmission_loglikelihood_sum(
            tx, coi, pop_log_cached);
    });

    std::cout << s << "x" << p << "x" << l
              << " chained=" << std::fixed << std::setprecision(4) << chained << " ms"
              << " fused=" << fused << " ms"
              << " speedup=" << (chained / fused) << "x\n";
}

} // namespace

int main() {
    moire_bench::BenchmarkSuite suite("Fused MCMC transmission benchmarks");
    suite.add("small", [] { run_scenario(50, 3, 20); });
    suite.add("medium", [] { run_scenario(200, 5, 100); });
    suite.add("large", [] { run_scenario(500, 10, 200); });
    suite.run_all();
    return 0;
}
