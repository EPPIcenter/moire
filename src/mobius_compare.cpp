#include <Rcpp.h>
#include <vector>
#include <numeric>
#include <cmath>
#include <algorithm>
#include "zeta_mobius.h"
#include "prob_any_missing.h"

// [[Rcpp::export]]
Rcpp::List compare_mobius_inclusion_exclusion(Rcpp::NumericVector p_full,
                                              Rcpp::IntegerVector support_vec,
                                              int max_n)
{
    // Inputs:
    // - p_full: probabilities over all categories (sum <= 1, not necessarily 1)
    // - support_vec: 0/1 for which categories define the constrained set S
    // - max_n: compute for n = 0..max_n
    const std::size_t num_categories = static_cast<std::size_t>(p_full.size());
    std::vector<double> p(num_categories);
    for (std::size_t i = 0; i < num_categories; ++i) p[i] = p_full[i];

    // Build mask and list of indices in S
    std::vector<int> support_mask(num_categories, 0);
    std::vector<int> support_indices;
    support_indices.reserve(num_categories);
    for (std::size_t i = 0; i < num_categories; ++i) {
        const int s = support_vec[i];
        support_mask[i] = s ? 1 : 0;
        if (s) support_indices.push_back(static_cast<int>(i));
    }
    const std::size_t K = support_indices.size();

    // Precompute total probability mass of the constrained set
    double total_mass_S = 0.0;
    for (int idx : support_indices) total_mass_S += p[idx];

    // Normalized probabilities within S
    std::vector<double> q;
    q.reserve(K);
    for (int idx : support_indices) {
        q.push_back(total_mass_S > 0.0 ? p[idx] / total_mass_S : 0.0);
    }

    // Mobius-based exact-support probabilities using the full distribution p
    auto mobius = zeta_mobius::create_multinomial_mobius_transform(p, static_cast<std::size_t>(max_n));
    std::vector<double> mobius_probs = mobius.query_probabilities(support_mask);

    // Mobius-based original-path equivalent: compute P(any missing | restricted to S)
    // use existing support_indices built above
    std::vector<double> p_any_missing_mobius = zeta_mobius::any_missing_restricted_via_mobius(p, support_indices, static_cast<std::size_t>(max_n));

    // Original implementation path: use probAnyMissingFunctor on normalized q within S
    std::vector<float> q_float;
    q_float.reserve(K);
    for (int idx : support_indices) {
        q_float.push_back(static_cast<float>(total_mass_S > 0.0 ? p[idx] / total_mass_S : 0.0));
    }
    probAnyMissingFunctor pam;
    // This returns a vector for n = 1..max_n (original code indexes pamVec[n-1])
    std::vector<double> pam_vec_1_to_N = pam.vectorized(q_float, static_cast<unsigned int>(max_n));

    std::vector<double> p_any_missing(max_n + 1, 0.0);
    std::vector<double> p_support_eq_S_IE(max_n + 1, 0.0);
    // n = 0
    p_any_missing[0] = (K == 0) ? 0.0 : 1.0;
    p_support_eq_S_IE[0] = (K == 0) ? 1.0 : 0.0;
    for (int n = 1; n <= max_n; ++n) {
        const double pam_n = pam_vec_1_to_N[static_cast<std::size_t>(n - 1)];
        p_any_missing[static_cast<std::size_t>(n)] = pam_n;
        p_support_eq_S_IE[static_cast<std::size_t>(n)] = std::pow(total_mass_S, n) * (1.0 - pam_n);
    }

    // Prepare result vectors for n = 0..max_n
    Rcpp::NumericVector nseq(max_n + 1);
    for (int n = 0; n <= max_n; ++n) nseq[n] = n;

    Rcpp::NumericVector mobius_out(max_n + 1);
    for (int n = 0; n <= max_n; ++n) mobius_out[n] = mobius_probs[static_cast<std::size_t>(n)];

    Rcpp::NumericVector ie_support_out(max_n + 1);
    for (int n = 0; n <= max_n; ++n) ie_support_out[n] = p_support_eq_S_IE[static_cast<std::size_t>(n)];

    Rcpp::NumericVector p_any_missing_out(max_n + 1);
    for (int n = 0; n <= max_n; ++n) p_any_missing_out[n] = p_any_missing[static_cast<std::size_t>(n)];

    Rcpp::NumericVector p_any_missing_ie_out(max_n + 1);
    for (int n = 0; n <= max_n; ++n) p_any_missing_ie_out[n] = p_any_missing[static_cast<std::size_t>(n)];

    Rcpp::NumericVector p_any_missing_mobius_out(max_n + 1);
    for (int n = 0; n <= max_n; ++n) p_any_missing_mobius_out[n] = p_any_missing_mobius[static_cast<std::size_t>(n)];

    Rcpp::NumericVector abs_diff_support(max_n + 1);
    for (int n = 0; n <= max_n; ++n) abs_diff_support[n] = std::abs(mobius_out[n] - ie_support_out[n]);

    Rcpp::NumericVector abs_diff_pam(max_n + 1);
    for (int n = 0; n <= max_n; ++n) abs_diff_pam[n] = std::abs(p_any_missing_ie_out[n] - p_any_missing_mobius_out[n]);

    return Rcpp::List::create(
        Rcpp::Named("n") = nseq,
        Rcpp::Named("mobius_support_eq_S") = mobius_out,
        Rcpp::Named("ie_support_eq_S_from_prob_any_missing") = ie_support_out,
        Rcpp::Named("prob_any_missing_restricted_S__ie") = p_any_missing_ie_out,
        Rcpp::Named("prob_any_missing_restricted_S__mobius") = p_any_missing_mobius_out,
        Rcpp::Named("abs_diff_support") = abs_diff_support,
        Rcpp::Named("abs_diff_prob_any_missing") = abs_diff_pam,
        Rcpp::Named("K") = static_cast<int>(K),
        Rcpp::Named("total_mass_S") = total_mass_S
    );
}


