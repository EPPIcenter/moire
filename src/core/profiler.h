#ifndef MOIRE_PROFILER_H_
#define MOIRE_PROFILER_H_

// ProfileScope is disabled by default to avoid lock contention in parallel code
// Define MOIRE_ENABLE_PROFILER_REGISTRY to enable it
#ifdef MOIRE_ENABLE_PROFILER_REGISTRY

#include <atomic>
#include <chrono>
#include <mutex>
#include <string>
#include <unordered_map>
#include <vector>

class ProfilerRegistry {
  public:
    struct Stat {
        std::atomic<long long> totalNanoseconds{0};
        std::atomic<long long> numCalls{0};
    };
    struct Snapshot {
        std::string key;
        long long totalNanoseconds;
        long long numCalls;
    };

    static ProfilerRegistry& instance();

    void add_sample(const std::string& key, long long durationNs);
    void reset();

    // immutable snapshot for reporting
    std::vector<Snapshot> snapshot();

  private:
    ProfilerRegistry() = default;
    std::unordered_map<std::string, Stat> stats_;
    std::mutex mutex_;
};

class ProfileScope {
  public:
    explicit ProfileScope(const char* key)
        : key_(key), start_(std::chrono::steady_clock::now()) {}
    explicit ProfileScope(const std::string& key)
        : key_(key), start_(std::chrono::steady_clock::now()) {}
    ~ProfileScope();

  private:
    std::string key_;
    std::chrono::steady_clock::time_point start_;
};

#else
// ProfileScope and ProfilerRegistry are no-ops unless MOIRE_ENABLE_PROFILER_REGISTRY is set.
#include <string>
#include <vector>

class ProfileScope {
  public:
    explicit ProfileScope(const char* key) {}
    explicit ProfileScope(const std::string& key) {}
    ~ProfileScope() = default;
};

class ProfilerRegistry {
  public:
    struct Snapshot {
        std::string key;
        long long totalNanoseconds;
        long long numCalls;
    };
    static ProfilerRegistry& instance();
    std::vector<Snapshot> snapshot() { return {}; }
    void reset() {}
};

#endif  // MOIRE_ENABLE_PROFILER_REGISTRY

#endif  // MOIRE_PROFILER_H_


