// Rcpp bridge for the effective-COI (eCOI) marginal likelihood core
// (inference/ecoi_marginal.h). These exports exist purely to validate the C++
// implementation against the R oracle (inst/scripts/prototype_ecoi_marginal.R)
// through the real package toolchain; they are not part of the sampler hot
// path. Inputs mirror the oracle: each locus is given as the *unnormalized*
// allele frequencies restricted to that locus' latent support.

#include <Rcpp.h>

#include "ecoi_marginal.h"

#include <span>
#include <vector>

namespace {

std::vector<ecoi_marginal::LocusSupport> build_population(const Rcpp::List& loci_freqs,
                                                          int max_coi)
{
    std::vector<ecoi_marginal::LocusSupport> loci;
    loci.reserve(static_cast<std::size_t>(loci_freqs.size()));
    for (R_xlen_t l = 0; l < loci_freqs.size(); ++l) {
        Rcpp::NumericVector v = loci_freqs[l];
        std::vector<float> support_p;
        support_p.reserve(static_cast<std::size_t>(v.size()));
        for (R_xlen_t j = 0; j < v.size(); ++j) {
            support_p.push_back(static_cast<float>(v[j]));
        }
        loci.push_back(ecoi_marginal::precompute_locus(
            std::span<const float>(support_p.data(), support_p.size()), max_coi));
    }
    return loci;
}

ecoi_marginal::Hyperparams make_hp(double coi_p, double coi_r, double r_alpha,
                                   double r_beta)
{
    return ecoi_marginal::Hyperparams{coi_p, coi_r, r_alpha, r_beta};
}

// Mirror the Chain plumbing: for each locus take the full population allele
// frequencies p_full and select the latent support entries (0-based indices),
// exactly as Chain::latent_allele_support + p.inner_iterators feed
// calc_transmission_process.
std::vector<ecoi_marginal::LocusSupport> build_population_from_indices(
    const Rcpp::List& p_full, const Rcpp::List& support_idx, int max_coi)
{
    std::vector<ecoi_marginal::LocusSupport> loci;
    loci.reserve(static_cast<std::size_t>(p_full.size()));
    for (R_xlen_t l = 0; l < p_full.size(); ++l) {
        Rcpp::NumericVector p = p_full[l];
        Rcpp::IntegerVector idx = support_idx[l];
        std::vector<float> support_p;
        support_p.reserve(static_cast<std::size_t>(idx.size()));
        for (R_xlen_t j = 0; j < idx.size(); ++j) {
            support_p.push_back(static_cast<float>(p[idx[j]]));
        }
        loci.push_back(ecoi_marginal::precompute_locus(
            std::span<const float>(support_p.data(), support_p.size()), max_coi));
    }
    return loci;
}

}  // namespace

// [[Rcpp::export]]
double ecoi_log_marginal_e_cpp(Rcpp::List loci_freqs, double coi_p, double coi_r,
                               double r_alpha, double r_beta, double e, int max_coi)
{
    const auto hp = make_hp(coi_p, coi_r, r_alpha, r_beta);
    const auto loci = build_population(loci_freqs, max_coi);
    return ecoi_marginal::log_marginal_e(
        e, std::span<const ecoi_marginal::LocusSupport>(loci.data(), loci.size()), hp,
        max_coi);
}

// [[Rcpp::export]]
double ecoi_log_atom1_cpp(Rcpp::List loci_freqs, double coi_p, double coi_r,
                          double r_alpha, double r_beta, int max_coi)
{
    const auto hp = make_hp(coi_p, coi_r, r_alpha, r_beta);
    const auto loci = build_population(loci_freqs, max_coi);
    return ecoi_marginal::log_atom1(
        std::span<const ecoi_marginal::LocusSupport>(loci.data(), loci.size()), hp);
}

// [[Rcpp::export]]
Rcpp::List ecoi_conditional_m_cpp(Rcpp::List loci_freqs, double coi_p, double coi_r,
                                  double r_alpha, double r_beta, double e, int max_coi)
{
    const auto hp = make_hp(coi_p, coi_r, r_alpha, r_beta);
    const auto loci = build_population(loci_freqs, max_coi);
    const auto cond = ecoi_marginal::conditional_m(
        e, std::span<const ecoi_marginal::LocusSupport>(loci.data(), loci.size()), hp,
        max_coi);
    return Rcpp::List::create(Rcpp::Named("m") = Rcpp::wrap(cond.m),
                              Rcpp::Named("r") = Rcpp::wrap(cond.r),
                              Rcpp::Named("w") = Rcpp::wrap(cond.w));
}

// [[Rcpp::export]]
double ecoi_log_marginal_e_from_indices_cpp(Rcpp::List p_full, Rcpp::List support_idx,
                                            double coi_p, double coi_r, double r_alpha,
                                            double r_beta, double e, int max_coi)
{
    const auto hp = make_hp(coi_p, coi_r, r_alpha, r_beta);
    const auto loci = build_population_from_indices(p_full, support_idx, max_coi);
    return ecoi_marginal::log_marginal_e(
        e, std::span<const ecoi_marginal::LocusSupport>(loci.data(), loci.size()), hp,
        max_coi);
}

// [[Rcpp::export]]
double ecoi_sample_log_marginal_e_cpp(Rcpp::List pops, Rcpp::NumericVector log_pi,
                                      double coi_p, double coi_r, double r_alpha,
                                      double r_beta, double e, int max_coi)
{
    const auto hp = make_hp(coi_p, coi_r, r_alpha, r_beta);
    std::vector<std::vector<ecoi_marginal::LocusSupport>> by_pop;
    by_pop.reserve(static_cast<std::size_t>(pops.size()));
    for (R_xlen_t p = 0; p < pops.size(); ++p) {
        by_pop.push_back(build_population(Rcpp::as<Rcpp::List>(pops[p]), max_coi));
    }
    std::vector<double> lp(log_pi.begin(), log_pi.end());
    return ecoi_marginal::sample_log_marginal_e(
        std::span<const std::vector<ecoi_marginal::LocusSupport>>(by_pop.data(),
                                                                  by_pop.size()),
        std::span<const double>(lp.data(), lp.size()), e, hp, max_coi);
}

// [[Rcpp::export]]
double ecoi_sample_log_atom1_cpp(Rcpp::List pops, Rcpp::NumericVector log_pi,
                                 double coi_p, double coi_r, double r_alpha,
                                 double r_beta, int max_coi)
{
    const auto hp = make_hp(coi_p, coi_r, r_alpha, r_beta);
    std::vector<std::vector<ecoi_marginal::LocusSupport>> by_pop;
    by_pop.reserve(static_cast<std::size_t>(pops.size()));
    for (R_xlen_t p = 0; p < pops.size(); ++p) {
        by_pop.push_back(build_population(Rcpp::as<Rcpp::List>(pops[p]), max_coi));
    }
    std::vector<double> lp(log_pi.begin(), log_pi.end());
    return ecoi_marginal::sample_log_atom1(
        std::span<const std::vector<ecoi_marginal::LocusSupport>>(by_pop.data(),
                                                                  by_pop.size()),
        std::span<const double>(lp.data(), lp.size()), hp);
}
