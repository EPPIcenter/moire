#include "profiler.h"
#include <Rcpp.h>

#ifdef MOIRE_ENABLE_PROFILER_REGISTRY
ProfilerRegistry& ProfilerRegistry::instance() {
    static ProfilerRegistry inst;
    return inst;
}

void ProfilerRegistry::add_sample(const std::string& key, long long durationNs) {
    // Fast path: look up without lock; if missing, lock and insert
    {
        std::lock_guard<std::mutex> lock(mutex_);
        auto& s = stats_[key];
        s.totalNanoseconds.fetch_add(durationNs, std::memory_order_relaxed);
        s.numCalls.fetch_add(1, std::memory_order_relaxed);
    }
}

void ProfilerRegistry::reset() {
    std::lock_guard<std::mutex> lock(mutex_);
    stats_.clear();
}

std::vector<ProfilerRegistry::Snapshot> ProfilerRegistry::snapshot() {
    std::lock_guard<std::mutex> lock(mutex_);
    std::vector<Snapshot> out;
    out.reserve(stats_.size());
    for (auto& kv : stats_) {
        Snapshot s;
        s.key = kv.first;
        s.totalNanoseconds = kv.second.totalNanoseconds.load();
        s.numCalls = kv.second.numCalls.load();
        out.emplace_back(std::move(s));
    }
    return out;
}

ProfileScope::~ProfileScope() {
    const auto end = std::chrono::steady_clock::now();
    const auto dur = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start_).count();
    ProfilerRegistry::instance().add_sample(key_, dur);
}

#else
ProfilerRegistry& ProfilerRegistry::instance() {
    static ProfilerRegistry inst;
    return inst;
}
#endif  // MOIRE_ENABLE_PROFILER_REGISTRY

// [[Rcpp::export]]
Rcpp::DataFrame moire_profiler_stats() {
#ifdef MOIRE_ENABLE_PROFILER_REGISTRY
    auto snap = ProfilerRegistry::instance().snapshot();
    std::vector<std::string> key;
    std::vector<double> total_ms;
    std::vector<double> avg_ms;
    std::vector<double> calls;
    key.reserve(snap.size());
    total_ms.reserve(snap.size());
    avg_ms.reserve(snap.size());
    calls.reserve(snap.size());
    for (auto& s : snap) {
        const auto total_ns = static_cast<double>(s.totalNanoseconds);
        const auto n = static_cast<double>(s.numCalls);
        key.emplace_back(s.key);
        total_ms.emplace_back(total_ns / 1e6);
        avg_ms.emplace_back(n > 0 ? (total_ns / 1e6) / n : 0.0);
        calls.emplace_back(n);
    }
    return Rcpp::DataFrame::create(Rcpp::Named("key") = key,
                                   Rcpp::Named("calls") = calls,
                                   Rcpp::Named("total_ms") = total_ms,
                                   Rcpp::Named("avg_ms") = avg_ms);
#else
    // Profiler disabled - return empty data frame
    return Rcpp::DataFrame::create(
        Rcpp::Named("key") = std::vector<std::string>(),
        Rcpp::Named("calls") = std::vector<double>(),
        Rcpp::Named("total_ms") = std::vector<double>(),
        Rcpp::Named("avg_ms") = std::vector<double>()
    );
#endif
}

// [[Rcpp::export]]
void moire_profiler_reset() {
#ifdef MOIRE_ENABLE_PROFILER_REGISTRY
    ProfilerRegistry::instance().reset();
#endif
}

#ifdef ENABLE_PROFILER
#include <gperftools/profiler.h>
#endif

// [[Rcpp::export]]
SEXP start_profiler(SEXP str) {
#ifdef ENABLE_PROFILER
    ProfilerStart(Rcpp::as<const char*>(str));
#else
    Rcpp::Rcerr << "Profiler not enabled. Enable by setting ENABLE_PROFILER=1 when compiling the package." << std::endl;
#endif
    return R_NilValue;
}

// [[Rcpp::export]]
SEXP stop_profiler() {
#ifdef ENABLE_PROFILER
    ProfilerStop();
#else
    Rcpp::Rcerr << "Profiler not enabled. Enable by setting ENABLE_PROFILER=1 when compiling the package." << std::endl;
#endif
  return R_NilValue;
}
