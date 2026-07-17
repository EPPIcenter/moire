#pragma once

#include <Rcpp.h>

#include <algorithm>
#include <cmath>
#include <map>
#include <string>
#include <utility>
#include <vector>

/// Compact Initialization diagnostics accumulated across Ill-conditioned starts.
struct InitializationDiagnostics
{
    int max_initialization_tries = 0;
    int chains_attempted = 0;
    int total_ill_conditioned = 0;
    int unknown_genotyping_count = 0;

    // 1-based sample / locus counts of first non-finite genotyping terms
    std::map<int, int> locus_counts{};
    std::map<int, int> sample_counts{};

    // Secondary: non-finite prior term counts (does not drive the retry gate)
    std::map<std::string, int> prior_nonfinite_counts{};

    // Up to a few example (sample, locus) pairs, 1-based
    std::vector<std::pair<int, int>> examples{};

    // Per-chain retry attempts that were ill-conditioned
    std::vector<int> ill_conditioned_per_chain{};

    std::string classification{"hard_starting_set"};
    std::string concentration{""};  // "locus", "sample", or ""
    int dominant_locus = NA_INTEGER;
    int dominant_sample = NA_INTEGER;
    double dominant_share = 0.0;
    int secondary_sample = NA_INTEGER;

    void record_ill_conditioned(int sample_1based, int locus_1based)
    {
        ++total_ill_conditioned;
        if (sample_1based > 0 && locus_1based > 0)
        {
            ++sample_counts[sample_1based];
            ++locus_counts[locus_1based];
            if (examples.size() < 5)
            {
                examples.emplace_back(sample_1based, locus_1based);
            }
        }
        else
        {
            ++unknown_genotyping_count;
        }
    }

    void record_prior_nonfinite(const std::string &term)
    {
        ++prior_nonfinite_counts[term];
    }

    void classify(double majority_threshold = 0.8)
    {
        if (total_ill_conditioned <= 0)
        {
            classification = "hard_starting_set";
            return;
        }

        int best_locus = NA_INTEGER;
        int best_locus_count = 0;
        for (const auto &kv : locus_counts)
        {
            if (kv.second > best_locus_count)
            {
                best_locus = kv.first;
                best_locus_count = kv.second;
            }
        }

        int best_sample = NA_INTEGER;
        int best_sample_count = 0;
        for (const auto &kv : sample_counts)
        {
            if (kv.second > best_sample_count)
            {
                best_sample = kv.first;
                best_sample_count = kv.second;
            }
        }

        const double locus_share =
            static_cast<double>(best_locus_count) / total_ill_conditioned;
        const double sample_share =
            static_cast<double>(best_sample_count) / total_ill_conditioned;

        const bool locus_consistent = best_locus_count > 0 &&
                                      locus_share >= majority_threshold;
        const bool sample_consistent = best_sample_count > 0 &&
                                       sample_share >= majority_threshold;

        // Locus preferred when both clear the bar (rare-allele guidance).
        if (locus_consistent)
        {
            classification = "consistent_failure_cause";
            concentration = "locus";
            dominant_locus = best_locus;
            dominant_share = locus_share;
            if (sample_consistent)
            {
                secondary_sample = best_sample;
            }
        }
        else if (sample_consistent)
        {
            classification = "consistent_failure_cause";
            concentration = "sample";
            dominant_sample = best_sample;
            dominant_share = sample_share;
        }
        else
        {
            classification = "hard_starting_set";
            concentration = "";
            dominant_share = std::max(locus_share, sample_share);
            if (best_locus_count >= best_sample_count)
            {
                dominant_locus = best_locus;
            }
            else
            {
                dominant_sample = best_sample;
            }
        }
    }

    Rcpp::List to_list() const
    {
        Rcpp::IntegerVector locus_idx;
        Rcpp::IntegerVector locus_n;
        for (const auto &kv : locus_counts)
        {
            locus_idx.push_back(kv.first);
            locus_n.push_back(kv.second);
        }

        Rcpp::IntegerVector sample_idx;
        Rcpp::IntegerVector sample_n;
        for (const auto &kv : sample_counts)
        {
            sample_idx.push_back(kv.first);
            sample_n.push_back(kv.second);
        }

        Rcpp::CharacterVector prior_terms;
        Rcpp::IntegerVector prior_n;
        for (const auto &kv : prior_nonfinite_counts)
        {
            prior_terms.push_back(kv.first);
            prior_n.push_back(kv.second);
        }

        Rcpp::IntegerVector ex_sample;
        Rcpp::IntegerVector ex_locus;
        for (const auto &ex : examples)
        {
            ex_sample.push_back(ex.first);
            ex_locus.push_back(ex.second);
        }

        return Rcpp::List::create(
            Rcpp::Named("max_initialization_tries") = max_initialization_tries,
            Rcpp::Named("chains_attempted") = chains_attempted,
            Rcpp::Named("total_ill_conditioned") = total_ill_conditioned,
            Rcpp::Named("unknown_genotyping_count") = unknown_genotyping_count,
            Rcpp::Named("locus_counts") = Rcpp::List::create(
                Rcpp::Named("locus") = locus_idx,
                Rcpp::Named("n") = locus_n),
            Rcpp::Named("sample_counts") = Rcpp::List::create(
                Rcpp::Named("sample") = sample_idx,
                Rcpp::Named("n") = sample_n),
            Rcpp::Named("prior_nonfinite_counts") = Rcpp::List::create(
                Rcpp::Named("term") = prior_terms,
                Rcpp::Named("n") = prior_n),
            Rcpp::Named("examples") = Rcpp::List::create(
                Rcpp::Named("sample") = ex_sample,
                Rcpp::Named("locus") = ex_locus),
            Rcpp::Named("ill_conditioned_per_chain") = ill_conditioned_per_chain,
            Rcpp::Named("classification") = classification,
            Rcpp::Named("concentration") = concentration,
            Rcpp::Named("dominant_locus") = dominant_locus,
            Rcpp::Named("dominant_sample") = dominant_sample,
            Rcpp::Named("dominant_share") = dominant_share,
            Rcpp::Named("secondary_sample") = secondary_sample);
    }
};
