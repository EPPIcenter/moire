#pragma once

#include <chrono>
#include <cstdlib>
#include <functional>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace moire_bench {

struct BenchmarkResult {
    std::string name;
    double duration_ms;
};

/// Regression-oriented result: scenario name, dimensions, operation, policy, time_ms (for CSV/CI).
struct RegressionResult {
    std::string scenario;
    size_t dim0{0}, dim1{0}, dim2{0};
    std::string operation;
    std::string policy;
    double time_ms{0.0};
};

/// Row for committed baseline files (section=summary|profiler).
struct BaselineRow {
    std::string section;
    std::string preset;
    std::string metric;
    double value{0.0};
    std::string unit;
    long long calls{0};
};

inline void emit_baseline_header(std::ostream& out) {
    out << "section,preset,metric,value,unit,calls\n";
}

inline void emit_baseline_row(std::ostream& out, const BaselineRow& r) {
    out << r.section << "," << r.preset << "," << r.metric << ","
        << std::fixed << std::setprecision(6) << r.value << ","
        << r.unit << "," << r.calls << "\n";
}

template<typename F>
double run_timed_ms(int warmup, int iters, F&& fn) {
    for (int i = 0; i < warmup; ++i) fn();
    double sum = 0.0;
    for (int i = 0; i < iters; ++i) {
        const auto start = std::chrono::high_resolution_clock::now();
        fn();
        const auto end = std::chrono::high_resolution_clock::now();
        sum += std::chrono::duration<double, std::milli>(end - start).count();
    }
    return sum / static_cast<double>(iters);
}

inline bool env_flag(const char* name) {
    const char* env = std::getenv(name);
    return env && (std::string(env) == "1" || std::string(env) == "true" || std::string(env) == "yes");
}

inline bool output_csv() {
    return env_flag("MOIRE_BENCHMARK_CSV");
}

/// When set, benchmarks emit baseline rows (same CSV schema as MCMC baseline).
inline bool output_baseline() {
    return env_flag("MOIRE_BENCHMARK_BASELINE");
}

inline double regression_threshold_pct() {
    const char* env = std::getenv("MOIRE_BENCH_REGRESSION_THRESHOLD");
    if (!env) return 5.0;
    try {
        return std::stod(env);
    } catch (...) {
        return 5.0;
    }
}

inline void emit_csv_header(std::ostream& out) {
    out << "scenario,dim0,dim1,dim2,operation,policy,time_ms\n";
}

inline void emit_csv_row(std::ostream& out, const RegressionResult& r) {
    out << r.scenario << "," << r.dim0 << "," << r.dim1 << "," << r.dim2
        << "," << r.operation << "," << r.policy << ","
        << std::fixed << std::setprecision(6) << r.time_ms << "\n";
}

class BenchmarkSuite {
private:
    std::string suite_name;
    std::vector<std::pair<std::string, std::function<void()>>> benchmarks;
    std::vector<BenchmarkResult> results;

public:
    explicit BenchmarkSuite(std::string name) : suite_name(std::move(name)) {}

    void add(const std::string &name, std::function<void()> fn) {
        benchmarks.emplace_back(name, std::move(fn));
    }

    void run_all() {
        if (!output_csv()) {
            std::cout << "\n⚡ Benchmark suite: " << suite_name << "\n";
            std::cout << std::string(suite_name.size() + 20, '=') << "\n";
        }
        for (auto &b : benchmarks) {
            auto start = std::chrono::high_resolution_clock::now();
            b.second();
            auto end = std::chrono::high_resolution_clock::now();
            const double ms = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0;
            results.push_back({b.first, ms});
            if (!output_csv()) {
                std::cout << std::setw(36) << std::left << (b.first + ":")
                          << std::fixed << std::setprecision(3) << ms << " ms\n";
            }
        }
        if (!output_csv()) print_summary();
    }

    void print_summary() const {
        double total_ms = 0.0;
        for (const auto &r : results) total_ms += r.duration_ms;
        std::cout << "\n📊 Benchmark summary for " << suite_name << ":\n";
        std::cout << "   Count:  " << results.size() << "\n";
        std::cout << "   Total:  " << std::fixed << std::setprecision(3) << total_ms << " ms\n";
        if (!results.empty()) {
            std::cout << "   Avg:    " << (total_ms / results.size()) << " ms\n";
        }
    }

    const std::vector<BenchmarkResult>& get_results() const { return results; }
};

} // namespace moire_bench


