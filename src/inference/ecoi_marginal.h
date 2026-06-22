#pragma once

// ---------------------------------------------------------------------------
// Effective-COI (eCOI) marginal transmission likelihood.
//
// Reparameterize a sample's (COI m, relatedness r) -> (eCOI e, m) with
//     e = (m - 1)(1 - r) + 1
// and treat e as the sampled parameter, marginalizing the discrete COI m
// analytically against its induced prior. This collapses the (m, r) ridge that
// the (m, r) sampler fights into a single smooth, identified eCOI dimension,
// while still permitting exact post-hoc recovery of the COI and relatedness
// marginals from the conditional pi(m | e, data) (Rao-Blackwellization).
//
// The math here is the C++ port of the validated R oracle in
// inst/scripts/prototype_ecoi_marginal.R (functions coi_log_prior,
// rel_log_prior, tx_loglik_locus_vec, optB_log_marginal_e, optB_log_atom1,
// optB_conditional_m). It is intentionally self-contained: per-locus PAM via
// the header-only pam_fast_paths inclusion-exclusion, and priors via exact
// lgamma. Note this uses the *exact* zero-truncated negative binomial log-pmf
// rather than Sampler::dztnbinom, which indexes a precomputed lgamma LUT with a
// truncated (int)(x + r) and is therefore only approximate; the eCOI path is
// the reference-quality implementation.
// ---------------------------------------------------------------------------

#include "pam_fast_paths.h"

#include <bit>
#include <cmath>
#include <cstddef>
#include <limits>
#include <span>
#include <vector>

namespace ecoi_marginal {

constexpr double kNegInf = -std::numeric_limits<double>::infinity();

// Largest latent-support size for which the closed-form inclusion-exclusion
// transmission term is precomputed and used. The closed form costs O(2^k) pow()
// per (m, r) while the binomial mixture costs ~O(max_coi), so beyond this the
// mixture is cheaper; correctness, however, holds for any k (the cancellation
// gate below makes the closed form fall back to the mixture when it is
// ill-conditioned). Tunable.
inline constexpr int kClosedFormMaxK = 6;

// Cancellation guard for the closed-form signed sum. The inclusion-exclusion
// terms are evaluated in linear space, so when relatedness is high (r -> 1) the
// true value is a tiny difference of O(T) terms and loses precision. Fall back
// to the cancellation-free log-space mixture whenever sum|term| exceeds this
// multiple of |sum| (i.e. more than ~8 of a double's ~15.9 digits are lost).
// Validated in inst/scripts/validate_ecoi_closed_form.R: marginal L(e) error
// < 5e-9 with a ~0.3% fallback rate at this threshold.
inline constexpr double kClosedFormCondLimit = 1e8;

// ---- numerics -------------------------------------------------------------

// log-sum-exp over double values, dropping non-finite entries (matches the R
// oracle's log_sum_exp, which filters !is.finite before reducing).
inline double log_sum_exp(std::span<const double> values)
{
    double max_el = kNegInf;
    for (double v : values) {
        if (std::isfinite(v) && v > max_el) {
            max_el = v;
        }
    }
    if (!std::isfinite(max_el)) {
        return kNegInf;
    }
    double sum = 0.0;
    for (double v : values) {
        if (std::isfinite(v)) {
            sum += std::exp(v - max_el);
        }
    }
    return max_el + std::log(sum);
}

// ---- eCOI map -------------------------------------------------------------

inline double ecoi_of(int m, double r)
{
    return static_cast<double>(m - 1) * (1.0 - r) + 1.0;
}

// Invert e = (m-1)(1-r)+1 at fixed m (valid for m >= 2).
inline double r_of(int m, double e)
{
    return 1.0 - (e - 1.0) / static_cast<double>(m - 1);
}

// ---- model priors (exact lgamma; match the R oracle) ----------------------

// Memoized log-factorial: lgamma_factorial(n) = log(n!) = lgamma(n + 1). The
// transmission inner loop evaluates binomial coefficients O(max_coi^2) times per
// proposal, so caching log(n!) (grown lazily, thread_local to stay race-free)
// removes the dominant lgamma cost from the hot path.
inline double lgamma_factorial(int n)
{
    static thread_local std::vector<double> cache;
    if (n < 0) {
        return 0.0;
    }
    if (static_cast<std::size_t>(n) >= cache.size()) {
        const std::size_t old = cache.size();
        cache.resize(static_cast<std::size_t>(n) + 1);
        for (std::size_t i = old; i < cache.size(); ++i) {
            cache[i] = std::lgamma(static_cast<double>(i) + 1.0);
        }
    }
    return cache[static_cast<std::size_t>(n)];
}

// log Binomial(x; size, prob), with the 0 * log(0) terms defined to 0 so the
// degenerate size == 0 case returns log(1) = 0.
inline double dbinom_log(int x, int size, double prob)
{
    if (x < 0 || x > size) {
        return kNegInf;
    }
    const double lcoef =
        lgamma_factorial(size) - lgamma_factorial(x) - lgamma_factorial(size - x);
    const double lp = (x > 0) ? x * std::log(prob) : 0.0;
    const double lq = (size - x > 0) ? (size - x) * std::log1p(-prob) : 0.0;
    return lcoef + lp + lq;
}

// Zero-truncated negative binomial log-pmf, m >= 1.
inline double dztnbinom_log(int m, double p, double r)
{
    if (m < 1) {
        return kNegInf;
    }
    return std::lgamma(m + r) + r * std::log(p) + m * std::log1p(-p) -
           std::lgamma(r) - std::lgamma(m + 1.0) - std::log1p(-std::pow(p, r));
}

// log Beta(x; alpha, beta) density for x in (0, 1).
inline double dbeta_log(double x, double alpha, double beta)
{
    const double lbeta =
        std::lgamma(alpha) + std::lgamma(beta) - std::lgamma(alpha + beta);
    const double ta = (alpha == 1.0) ? 0.0 : (alpha - 1.0) * std::log(x);
    const double tb = (beta == 1.0) ? 0.0 : (beta - 1.0) * std::log1p(-x);
    return ta + tb - lbeta;
}

// ---- hyperparameters ------------------------------------------------------

struct Hyperparams {
    double coi_p;    // ZTNB success prob
    double coi_r;    // ZTNB dispersion
    double r_alpha;  // Beta relatedness prior
    double r_beta;
};

// ---- per-locus precompute (independent of m, r) ---------------------------

// Holds the (m, r)-independent transmission inputs for one locus' latent
// support: logT = log(sum of allele freqs over the support), and
// logpas[n - 1] = log P(all support alleles seen at least once in n draws),
// for n = 1 .. max_coi. Mirrors precompute_locus() in the R oracle.
struct LocusSupport {
    int k = 0;
    double logT = 0.0;
    std::vector<double> logpas;  // length max_coi (binomial-mixture path)
    // Closed-form inclusion-exclusion terms, indexed by excluded-allele mask S:
    //   subset_V[S] = sum of UNNORMALIZED support freqs NOT in S   (subset_V[0] = T).
    // The per-(m, r) locus term then collapses (binomial theorem) to
    //   t_l(m, r) = sum_S (-1)^popcount(S) * subset_V[S] * (r + (1-r) subset_V[S])^{m-1}.
    // Empty when k > kClosedFormMaxK, in which case tx_loglik_locus uses the
    // binomial mixture instead.
    std::vector<double> subset_V;  // length 2^k when populated
};

// support_p: unnormalized allele frequencies restricted to the latent support.
inline LocusSupport precompute_locus(std::span<const float> support_p, int max_coi)
{
    LocusSupport pre;
    pre.k = static_cast<int>(support_p.size());
    double total = 0.0;
    for (float v : support_p) {
        total += static_cast<double>(v);
    }
    pre.logT = std::log(total);

    // Closed-form inclusion-exclusion terms V_S = sum of (unnormalized) support
    // freqs not in S. compute_subset_sums(support_p) returns exactly this,
    // indexed by the excluded mask (subset_V[0] = total). Only populated when the
    // support is small enough for 2^k to beat the mixture (see kClosedFormMaxK).
    if (pre.k >= 1 && pre.k <= kClosedFormMaxK) {
        pam_fast_paths::compute_subset_sums(support_p, pre.subset_V);
    }

    std::vector<float> q;
    q.reserve(support_p.size());
    for (float v : support_p) {
        q.push_back(static_cast<float>(static_cast<double>(v) / total));
    }

    pre.logpas.assign(static_cast<std::size_t>(std::max(max_coi, 0)), kNegInf);
    if (pre.k <= 1) {
        // single category: all-seen probability is 1 for every n >= 1.
        std::fill(pre.logpas.begin(), pre.logpas.end(), 0.0);
        return pre;
    }

    std::vector<double> subset_sums;
    pam_fast_paths::compute_subset_sums(q, subset_sums);
    for (int n = 1; n <= max_coi; ++n) {
        double pas = pam_fast_paths::prob_all_seen_from_subset_sums(
            std::span<const double>(subset_sums), static_cast<unsigned>(n));
        if (pas < 1e-300) {
            pas = 1e-300;
        }
        pre.logpas[static_cast<std::size_t>(n - 1)] = std::log(pas);
    }
    return pre;
}

// ---- transmission likelihood (mirrors tx_loglik_locus_vec) ----------------

// log t_l(m, r) for one locus:
//   t = sum_{i=0}^{m-k} Binom(i; m-1, r) * P(all seen | m-i draws) * T^{m-i}
inline double tx_loglik_locus(const LocusSupport& pre, int m, double r)
{
    if (m < pre.k || m < 1) {
        return kNegInf;
    }

    // Closed-form path (exact binomial-theorem collapse of the IBD mixture):
    //   t = sum_S (-1)^popcount(S) V_S (r + (1-r) V_S)^{m-1}.
    // O(2^k) pow() vs the O(max_coi) mixture below. Guarded against catastrophic
    // cancellation (r -> 1): if the signed sum has lost too much precision we
    // fall through to the cancellation-free log-space mixture.
    if (!pre.subset_V.empty()) {
        const double exponent = static_cast<double>(m - 1);
        const double one_minus_r = 1.0 - r;
        double sum = 0.0;
        double abs_sum = 0.0;
        const std::size_t n_masks = pre.subset_V.size();
        for (std::size_t mask = 0; mask < n_masks; ++mask) {
            const double V = pre.subset_V[mask];
            if (V <= 0.0) {
                continue;  // V = 0 (full exclusion) contributes nothing
            }
            const double term = V * std::pow(r + one_minus_r * V, exponent);
            abs_sum += term;
            if (std::popcount(static_cast<unsigned>(mask)) & 1) {
                sum -= term;
            } else {
                sum += term;
            }
        }
        if (sum > 0.0 && abs_sum <= kClosedFormCondLimit * sum) {
            return std::log(sum);
        }
        // else: ill-conditioned (or numerically <= 0) -> mixture fallback.
    }

    const int n_terms = m - pre.k + 1;
    // r is constant across the mixture, so the binomial coefficient's
    // i*log(r) + (size-i)*log1p(-r) factors are computed from two hoisted
    // logs instead of re-evaluating dbinom_log (two transcendental calls) per
    // term. The scratch is thread_local so this stays allocation-free and safe
    // to call from parallel per-sample recomputes.
    static thread_local std::vector<double> terms;
    terms.resize(static_cast<std::size_t>(n_terms));
    const int size = m - 1;
    const double lg_size = lgamma_factorial(size);
    const double logr = (r > 0.0) ? std::log(r) : kNegInf;
    const double log1mr = std::log1p(-r);
    for (int i = 0; i < n_terms; ++i) {
        const int n_draws = m - i;  // 1 .. m
        const double logpas = pre.logpas[static_cast<std::size_t>(n_draws - 1)];
        const double lcoef =
            lg_size - lgamma_factorial(i) - lgamma_factorial(size - i);
        const double lp = (i > 0) ? i * logr : 0.0;
        const double lq = (size - i > 0) ? (size - i) * log1mr : 0.0;
        terms[static_cast<std::size_t>(i)] =
            lcoef + lp + lq + logpas + n_draws * pre.logT;
    }
    return log_sum_exp(std::span<const double>(terms));
}

inline double tx_loglik_sample(std::span<const LocusSupport> loci, int m, double r)
{
    double acc = 0.0;
    for (const LocusSupport& pre : loci) {
        const double t = tx_loglik_locus(pre, m, r);
        if (!std::isfinite(t)) {
            return kNegInf;
        }
        acc += t;
    }
    return acc;
}

// ---- Option B: parameter e, marginalize m ---------------------------------

// Induced log weight for a single COI m given eCOI e (continuous part, m >= 2):
//   log w(m | e) = log NegBinom(m) + log Beta(r(m,e)) - log(m - 1) [Jacobian]
//                  + transmission(m, r(m,e))
// Returns -inf when r(m, e) is outside (0, 1).
inline double log_weight_m_given_e(int m, double e,
                                   std::span<const LocusSupport> loci,
                                   const Hyperparams& hp)
{
    const double r = r_of(m, e);
    if (r <= 0.0 || r >= 1.0) {
        return kNegInf;
    }
    return dztnbinom_log(m, hp.coi_p, hp.coi_r) +
           dbeta_log(r, hp.r_alpha, hp.r_beta) -
           std::log(static_cast<double>(m - 1)) +
           tx_loglik_sample(loci, m, r);
}

// log unnormalized marginal density of e (continuous part, e > 1).
inline double log_marginal_e(double e, std::span<const LocusSupport> loci,
                             const Hyperparams& hp, int max_coi)
{
    if (e <= 1.0) {
        return kNegInf;
    }
    const int m_lo = std::max(2, static_cast<int>(std::ceil(e)));
    if (m_lo > max_coi) {
        return kNegInf;
    }
    std::vector<double> terms;
    terms.reserve(static_cast<std::size_t>(max_coi - m_lo + 1));
    for (int m = m_lo; m <= max_coi; ++m) {
        terms.push_back(log_weight_m_given_e(m, e, loci, hp));
    }
    return log_sum_exp(std::span<const double>(terms));
}

// Atom at e = 1 (monoclonal, m = 1). The relatedness prior integrates to 1, so
// only the COI prior and the (relatedness-free) transmission term appear.
inline double log_atom1(std::span<const LocusSupport> loci, const Hyperparams& hp)
{
    return dztnbinom_log(1, hp.coi_p, hp.coi_r) +
           tx_loglik_sample(loci, /*m=*/1, /*r=*/0.0);
}

// ---- COI-prior-free (population-e hierarchy) variants ----------------------
//
// In the population-e design the prior on the COI is *not* part of the model:
// the headline parameter is the effective COI e, and the per-population
// continuous prior f_p(e) = 1 + Gamma(k_p, k_p / mu_plus_p) supplies the prior on
// e directly (folded into the per-sample population mixture by the sampler,
// outside this header; there is no e = 1 atom). The discrete COI m is still
// integrated out, but the only weight carried in that integral is the relatedness
// pdf times the change-of-variables Jacobian 1/(m-1). The marginal_ecoi sampler
// fits under uniform relatedness (Beta(1, 1)); these Beta(r_alpha, r_beta) forms
// are retained for the on-chain self-checks and for post-hoc COI/relatedness
// recovery from the fitted eCOI posterior.

// Relatedness-only induced log weight for a single COI m given eCOI e (m >= 2):
//   log w(m | e) = log Beta(r(m,e); r_alpha, r_beta) - log(m - 1) [Jacobian]
//                  + transmission(m, r(m,e))
// Returns -inf when r(m, e) is outside (0, 1).
inline double log_weight_m_given_e_rel(int m, double e,
                                       std::span<const LocusSupport> loci,
                                       double r_alpha, double r_beta)
{
    const double r = r_of(m, e);
    if (r <= 0.0 || r >= 1.0) {
        return kNegInf;
    }
    return dbeta_log(r, r_alpha, r_beta) -
           std::log(static_cast<double>(m - 1)) + tx_loglik_sample(loci, m, r);
}

// COI-prior-free continuous (polyclonal) marginal log-likelihood at eCOI e > 1.
inline double log_marginal_e_rel(double e, std::span<const LocusSupport> loci,
                                 double r_alpha, double r_beta, int max_coi)
{
    if (e <= 1.0) {
        return kNegInf;
    }
    const int m_lo = std::max(2, static_cast<int>(std::ceil(e)));
    if (m_lo > max_coi) {
        return kNegInf;
    }
    std::vector<double> terms;
    terms.reserve(static_cast<std::size_t>(max_coi - m_lo + 1));
    for (int m = m_lo; m <= max_coi; ++m) {
        terms.push_back(log_weight_m_given_e_rel(m, e, loci, r_alpha, r_beta));
    }
    return log_sum_exp(std::span<const double>(terms));
}

// COI-prior-free monoclonal atom log-likelihood (m = 1, r = 0): just the
// relatedness-free transmission term. The prior mass pi_mono is applied by the
// sampler, not here.
inline double log_atom_rel(std::span<const LocusSupport> loci)
{
    return tx_loglik_sample(loci, /*m=*/1, /*r=*/0.0);
}

// ---- Rao-Blackwell recovery: conditional pi(m | e, data) ------------------

struct ConditionalM {
    std::vector<int> m;
    std::vector<double> r;  // r(m, e) for each m
    std::vector<double> w;  // normalized weights (all zero if degenerate)
};

inline ConditionalM conditional_m(double e, std::span<const LocusSupport> loci,
                                  const Hyperparams& hp, int max_coi)
{
    ConditionalM out;
    const int m_lo = std::max(2, static_cast<int>(std::ceil(e)));
    if (m_lo > max_coi) {
        return out;
    }
    std::vector<double> logw;
    for (int m = m_lo; m <= max_coi; ++m) {
        out.m.push_back(m);
        out.r.push_back(r_of(m, e));
        logw.push_back(log_weight_m_given_e(m, e, loci, hp));
    }
    const double lse = log_sum_exp(std::span<const double>(logw));
    out.w.assign(logw.size(), 0.0);
    if (std::isfinite(lse)) {
        for (std::size_t i = 0; i < logw.size(); ++i) {
            out.w[i] = std::exp(logw[i] - lse);
        }
    }
    return out;
}

// Post-hoc conditional pi(m | e, data) using the relatedness-only weight, for
// recovering the COI/relatedness marginals from the sampled e (population-e
// hierarchy: the COI prior is not part of the model).
inline ConditionalM conditional_m_rel(double e, std::span<const LocusSupport> loci,
                                      double r_alpha, double r_beta, int max_coi)
{
    ConditionalM out;
    const int m_lo = std::max(2, static_cast<int>(std::ceil(e)));
    if (m_lo > max_coi) {
        return out;
    }
    std::vector<double> logw;
    for (int m = m_lo; m <= max_coi; ++m) {
        out.m.push_back(m);
        out.r.push_back(r_of(m, e));
        logw.push_back(log_weight_m_given_e_rel(m, e, loci, r_alpha, r_beta));
    }
    const double lse = log_sum_exp(std::span<const double>(logw));
    out.w.assign(logw.size(), 0.0);
    if (std::isfinite(lse)) {
        for (std::size_t i = 0; i < logw.size(); ++i) {
            out.w[i] = std::exp(logw[i] - lse);
        }
    }
    return out;
}

// ---- per-sample population mixture ----------------------------------------
//
// A sample's COI/relatedness (hence the induced weight w(m|e)) are shared
// across populations; only the per-locus transmission term depends on the
// population's allele frequencies. So the marginalization over m sits *inside*
// the population mixture, and the single-population log_marginal_e / log_atom1
// above are exactly the per-population inner terms:
//
//   L_s(e) = LSE_p [ log pi_p + log_marginal_e(e, supports[p], .) ]
//
// `pops[p]` holds the LocusSupport (built from population p's allele
// frequencies on the sample's latent supports) for every locus.

// Continuous (polyclonal) per-sample log-likelihood at eCOI e (e > 1).
inline double sample_log_marginal_e(std::span<const std::vector<LocusSupport>> pops,
                                    std::span<const double> log_pi, double e,
                                    const Hyperparams& hp, int max_coi)
{
    std::vector<double> terms(pops.size());
    for (std::size_t p = 0; p < pops.size(); ++p) {
        terms[p] = log_pi[p] +
                   log_marginal_e(e, std::span<const LocusSupport>(pops[p]), hp, max_coi);
    }
    return log_sum_exp(std::span<const double>(terms));
}

// Monoclonal (m = 1, e = 1) per-sample log-likelihood. Selected, against the
// continuous part above, by the per-sample monoclonal indicator the sampler
// carries (Stage 3b); both pieces are exposed here so the likelihood core stays
// free of sampler state.
inline double sample_log_atom1(std::span<const std::vector<LocusSupport>> pops,
                               std::span<const double> log_pi,
                               const Hyperparams& hp)
{
    std::vector<double> terms(pops.size());
    for (std::size_t p = 0; p < pops.size(); ++p) {
        terms[p] = log_pi[p] + log_atom1(std::span<const LocusSupport>(pops[p]), hp);
    }
    return log_sum_exp(std::span<const double>(terms));
}

}  // namespace ecoi_marginal
