#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <span>
#include <string>
#include <vector>

namespace pam_cache {

inline bool env_flag(const char* name) noexcept
{
    const char* env = std::getenv(name);
    return env && (std::string(env) == "1" || std::string(env) == "true" ||
                   std::string(env) == "yes");
}

struct Config {
    /// When true, bucket keys by quantize_simplex(q, delta).
    bool quantize{true};
    /// Grid step on the simplex; default matches update_p's 1e-4 stability floor.
    float delta{1e-4f};
    /// When true, recompute pamVec on cache hits and track approximation error.
    bool verify{false};

    static Config& instance() noexcept
    {
        static Config cfg = load_from_env();
        return cfg;
    }

    static Config load_from_env() noexcept
    {
        Config cfg;
        if (const char* q = std::getenv("MOIRE_PAM_CACHE_QUANTIZE")) {
            const std::string v(q);
            cfg.quantize = !(v == "0" || v == "false" || v == "no");
        }
        if (const char* d = std::getenv("MOIRE_PAM_CACHE_DELTA")) {
            cfg.delta = std::max(std::strtof(d, nullptr), 1e-8f);
        }
        cfg.verify = env_flag("MOIRE_PAM_CACHE_VERIFY");
        return cfg;
    }
};

struct Stats {
    std::uint64_t exact_hits{0};
    std::uint64_t quant_hits{0};
    std::uint64_t misses{0};
    /// k <= kLowKMaxSupport: computed inline without cache store on miss.
    std::uint64_t low_k_inline{0};
    /// low-k PAM extended from cached vector when COI changes by at most 2.
    std::uint64_t coi_extend{0};
    /// update_p: constrained q unchanged (changed allele not in latent support).
    std::uint64_t p_q_unchanged{0};
    /// update_p: PAM refilled via single-component q delta on low-k path.
    std::uint64_t p_q_one_step{0};
    std::uint64_t verify_checks{0};
    std::uint64_t verify_over_tol{0};
    double verify_max_pam_diff{0.0};
    double verify_max_log_diff{0.0};

    void reset() noexcept { *this = Stats{}; }
};

inline constexpr double kVerifyPamTol = 3e-4;
inline constexpr double kVerifyLogTol = 1e-3;

inline Stats& stats() noexcept
{
    thread_local Stats s;
    return s;
}

/// Round each component to delta grid and renormalize to the simplex.
inline void quantize_simplex_into(std::span<const float> q, float delta, std::vector<float>& out)
{
    out.resize(q.size());
    if (q.empty() || delta <= 0.f) {
        return;
    }
    float sum = 0.f;
    for (std::size_t i = 0; i < q.size(); ++i) {
        out[i] = std::round(q[i] / delta) * delta;
        if (out[i] < 0.f) {
            out[i] = 0.f;
        }
        sum += out[i];
    }
    if (sum <= 0.f) {
        const float uniform = 1.f / static_cast<float>(q.size());
        std::fill(out.begin(), out.end(), uniform);
        return;
    }
    for (float& v : out) {
        v /= sum;
    }
}

inline std::vector<float> quantize_simplex(std::span<const float> q, float delta)
{
    std::vector<float> out;
    quantize_simplex_into(q, delta, out);
    return out;
}

inline std::vector<float> cache_key(std::span<const float> q)
{
    const Config& cfg = Config::instance();
    if (!cfg.quantize) {
        return std::vector<float>(q.begin(), q.end());
    }
    return quantize_simplex(q, cfg.delta);
}

inline void cache_key_into(std::span<const float> q, std::vector<float>& scratch)
{
    const Config& cfg = Config::instance();
    if (!cfg.quantize) {
        scratch.assign(q.begin(), q.end());
        return;
    }
    quantize_simplex_into(q, cfg.delta, scratch);
}

} // namespace pam_cache

struct PamCachedVectors {
    std::vector<double> pam;
    std::vector<float> log_one_minus_pam;
};

inline void fill_log_one_minus_pam(const std::vector<double>& pam,
                                   std::vector<float>& log_one_minus_out)
{
    log_one_minus_out.resize(pam.size());
    for (std::size_t i = 0; i < pam.size(); ++i) {
        log_one_minus_out[i] = static_cast<float>(std::log1p(-pam[i]));
    }
}

inline PamCachedVectors make_pam_cached_vectors(std::vector<double> pam)
{
    PamCachedVectors out;
    out.pam = std::move(pam);
    fill_log_one_minus_pam(out.pam, out.log_one_minus_pam);
    return out;
}

/// Thread-local cache for P(any missing) vectors and log(1 - PAM).
///
/// Backed by a fixed-size, direct-mapped table with a short linear-probe
/// window. Slots are allocated once and overwritten in place, so a steady-state
/// store performs no heap allocation and never triggers a full-table clear --
/// the previous std::unordered_map design re-allocated an Entry (two vectors)
/// on every miss and wiped all 512 entries on overflow, which dominated the
/// update_p hot path at the measured ~80% miss rate.
class PamVectorCache {
public:
    static constexpr std::size_t kCapacity = 2048;  // power of two
    static constexpr std::size_t kProbe = 8;

    void clear() noexcept
    {
        for (Slot& s : slots_) {
            s.occupied = false;
        }
    }

    const PamCachedVectors* lookup(std::span<const float> q,
                                   unsigned min_events,
                                   unsigned max_events) const noexcept
    {
        if (slots_.empty()) {
            return nullptr;
        }
        thread_local std::vector<float> key_scratch;
        pam_cache::cache_key_into(q, key_scratch);
        const std::span<const float> key_q(key_scratch);
        const std::size_t h = hash_key(key_q, min_events, max_events);
        const std::size_t base = h & (kCapacity - 1);
        for (std::size_t i = 0; i < kProbe; ++i) {
            const Slot& s = slots_[(base + i) & (kCapacity - 1)];
            if (!s.occupied) {
                break;  // open-addressing: a match would precede the first gap
            }
            if (s.hash == h && s.min_events == min_events &&
                s.max_events == max_events && floats_equal(s.key_q, key_q)) {
                record_hit(s.key_q, q);
                return &s.vectors;
            }
        }
        return nullptr;
    }

    /// Store the freshly computed PAM vector, copying into reused slot storage.
    /// `pam` is borrowed (caller keeps ownership of its scratch buffer).
    const PamCachedVectors& store_pam(std::span<const float> q,
                                      unsigned min_events,
                                      unsigned max_events,
                                      std::span<const double> pam)
    {
        if (slots_.empty()) {
            slots_.resize(kCapacity);
        }
        thread_local std::vector<float> key_scratch;
        pam_cache::cache_key_into(q, key_scratch);
        const std::span<const float> key_q(key_scratch);
        const std::size_t h = hash_key(key_q, min_events, max_events);

        const std::size_t base = h & (kCapacity - 1);
        Slot* target = nullptr;
        for (std::size_t i = 0; i < kProbe; ++i) {
            Slot& s = slots_[(base + i) & (kCapacity - 1)];
            if (!s.occupied ||
                (s.hash == h && s.min_events == min_events &&
                 s.max_events == max_events && floats_equal(s.key_q, key_q))) {
                target = &s;
                break;
            }
        }
        if (target == nullptr) {
            target = &slots_[base];  // probe window full: evict chain head
        }
        return fill_slot(*target, h, min_events, max_events, key_q, pam);
    }

    const PamCachedVectors& store_and_get(std::span<const float> q,
                                          unsigned min_events,
                                          unsigned max_events,
                                          PamCachedVectors vectors)
    {
        return store_pam(q, min_events, max_events,
                         std::span<const double>(vectors.pam));
    }

    const PamCachedVectors& store_and_get(std::span<const float> q,
                                          unsigned min_events,
                                          unsigned max_events,
                                          std::vector<double> pam)
    {
        return store_pam(q, min_events, max_events,
                         std::span<const double>(pam));
    }

    /// Legacy accessor for tests and verify paths.
    const std::vector<double>& store_and_get_pam_only(std::span<const float> q,
                                                      unsigned min_events,
                                                      unsigned max_events,
                                                      std::vector<double> pam)
    {
        return store_pam(q, min_events, max_events,
                         std::span<const double>(pam))
            .pam;
    }

private:
    struct Slot {
        bool occupied{false};
        std::size_t hash{0};
        unsigned min_events{0};
        unsigned max_events{0};
        std::vector<float> key_q;
        PamCachedVectors vectors;
    };

    static void record_hit(const std::vector<float>& key_q,
                           std::span<const float> q) noexcept
    {
        if (pam_cache::Config::instance().quantize && !floats_equal(key_q, q)) {
            ++pam_cache::stats().quant_hits;
        } else {
            ++pam_cache::stats().exact_hits;
        }
    }

    static const PamCachedVectors& fill_slot(Slot& s, std::size_t h,
                                             unsigned min_events,
                                             unsigned max_events,
                                             std::span<const float> key_q,
                                             std::span<const double> pam)
    {
        s.occupied = true;
        s.hash = h;
        s.min_events = min_events;
        s.max_events = max_events;
        s.key_q.assign(key_q.begin(), key_q.end());
        s.vectors.pam.assign(pam.begin(), pam.end());
        fill_log_one_minus_pam(s.vectors.pam, s.vectors.log_one_minus_pam);
        return s.vectors;
    }

    static std::size_t hash_key(std::span<const float> q,
                                unsigned min_events,
                                unsigned max_events) noexcept
    {
        std::size_t h = static_cast<std::size_t>(min_events) ^
                        (static_cast<std::size_t>(max_events) << 16);
        for (float f : q) {
            std::uint32_t bits = 0;
            std::memcpy(&bits, &f, sizeof(bits));
            h ^= static_cast<std::size_t>(bits) + 0x9e3779b9 + (h << 6) + (h >> 2);
        }
        return h;
    }

    static bool floats_equal(const std::vector<float>& a,
                             std::span<const float> b) noexcept
    {
        if (a.size() != b.size()) {
            return false;
        }
        return std::memcmp(a.data(), b.data(), a.size() * sizeof(float)) == 0;
    }

    std::vector<Slot> slots_;
};

inline PamVectorCache& pam_vector_cache()
{
    thread_local PamVectorCache cache;
    return cache;
}

inline void pam_cache_verify_hit(std::span<const float> q,
                                 unsigned min_events,
                                 unsigned max_events,
                                 const PamCachedVectors& cached,
                                 const std::vector<double>& exact_pam)
{
    pam_cache::Stats& s = pam_cache::stats();
    ++s.verify_checks;
    const std::size_t n = std::min(cached.pam.size(), exact_pam.size());
    for (std::size_t i = 0; i < n; ++i) {
        const double pam_diff = std::fabs(cached.pam[i] - exact_pam[i]);
        if (pam_diff > s.verify_max_pam_diff) {
            s.verify_max_pam_diff = pam_diff;
        }
        if (pam_diff > pam_cache::kVerifyPamTol) {
            ++s.verify_over_tol;
        }
        const double log_diff = std::fabs(
            static_cast<double>(cached.log_one_minus_pam[i]) -
            std::log1p(-exact_pam[i]));
        if (log_diff > s.verify_max_log_diff) {
            s.verify_max_log_diff = log_diff;
        }
        if (log_diff > pam_cache::kVerifyLogTol) {
            ++s.verify_over_tol;
        }
    }
    (void)min_events;
    (void)max_events;
    (void)q;
}
