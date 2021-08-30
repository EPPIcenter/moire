#pragma once

#ifndef CHAIN_H_
#define CHAIN_H_

#include "combination_indices_generator.h"
#include "genotyping_data.h"
#include "lookup.h"
#include "parameters.h"
#include "prob_any_missing.h"
#include "sampler.h"

#include <Rcpp.h>

class Chain
{
   private:
    GenotypingData genotyping_data;
    Lookup lookup;
    Parameters params;
    Sampler sampler;
    probAnyMissingFunctor probAnyMissing_;
    std::vector<double> prVec_{};

    CombinationIndicesGenerator allele_index_generator_;

    void initialize_p();
    void initialize_m();
    void initialize_mean_coi();
    void initialize_eps_neg();
    void initialize_eps_pos();
    void initialize_likelihood();

    std::vector<double> reweight_allele_frequencies(
        std::vector<double> const &allele_frequencies,
        std::vector<int> const &observed_genotype, double epsilon_neg,
        double epsilon_pos, int coi);

    // std::vector<double> calc_genotype_log_pmf(
    //     std::vector<std::vector<int>> const &genotypes, int coi,
    // std::vector<double> const &allele_frequencies, int num_genotypes);

    double calc_transmission_process(
        std::vector<int> const &allele_index_vec,
        std::vector<double> const &allele_frequencies, int coi);

    double calc_observation_process(std::vector<int> const &allele_index_vec,
                                    std::vector<int> const &obs_genotype,
                                    int coi, double epsilon_neg,
                                    double epsilon_pos);

    double calc_genotype_log_pmf(std::vector<int> const &allele_index_vec,
                                 std::vector<int> const &obs_genotype,
                                 double epsilon_pos, double epsilon_neg,
                                 int coi,
                                 std::vector<double> const &allele_frequencies);

    std::vector<double> calc_obs_genotype_lliks(
        std::vector<int> const &obs_genotype,
        std::vector<std::vector<int>> const &true_genotypes, double epsilon_neg,
        double epsilon_pos, int num_genotypes);

    long double calc_genotype_marginal_llik(
        std::vector<int> const &obs_genotype, int coi,
        std::vector<double> const &allele_frequencies, double epsilon_neg,
        double epsilon_pos, bool importance_sample = false);

    long double calc_exact_genotype_marginal_llik(
        std::vector<int> const &obs_genotype, int coi,
        std::vector<double> const &allele_frequencies, double epsilon_neg,
        double epsilon_pos);

    long double importance_sample(std::vector<int> const &obs_genotype, int coi,
                                  double epsilon_neg, double epsilon_pos,
                                  std::vector<double> const &allele_frequencies,
                                  int sampling_depth);

    long double importance_sample2(
        std::vector<int> const &obs_genotype, int coi, double epsilon_neg,
        double epsilon_pos, std::vector<double> const &allele_frequencies,
        int sampling_depth
    );

    long double monte_carlo_sample(std::vector<int> const &obs_genotype,
                                   int coi, double epsilon_neg,
                                   double epsilon_pos,
                                   std::vector<double> const &true_distribution,
                                   int sampling_depth);

    long double calc_estimated_genotype_marginal_llik(
        std::vector<int> const &obs_genotype, int coi,
        std::vector<double> const &allele_frequencies, double epsilon_neg,
        double epsilon_pos, int sampling_depth, bool imp_sample);

   public:
    std::vector<std::vector<double>> llik_old{};
    std::vector<std::vector<double>> llik_new{};
    double llik;

    // Mean COI
    // TODO: Allow for other priors on complexity of infection
    // std::string prior;
    // double poisson_prior_lambda;
    double mean_coi;

    // COI
    std::vector<int> m{};
    int prop_m;
    std::vector<int> m_accept{};
    std::vector<double> m_prop_mean{};

    // Allele Frequencies
    std::vector<std::vector<double>> p{};
    std::vector<double> prop_p{};
    std::vector<double> p_prop_var{};
    std::vector<int> p_accept{};

    // Epsilon Positive
    // double eps_pos;
    std::vector<double> eps_pos{};
    double prop_eps_pos;
    std::vector<int> eps_pos_accept{};
    double eps_pos_var;

    // Epsilon Negative
    // double eps_neg;
    std::vector<double> eps_neg{};
    double prop_eps_neg;
    std::vector<int> eps_neg_accept{};
    double eps_neg_var;

    std::vector<int> individual_accept{};

    Chain(GenotypingData genotyping_data, Lookup lookup, Parameters params);
    void update_m(int iteration);
    void update_mean_coi(int iteration);
    void update_p(int iteration);
    void update_eps(int iteration);
    void update_eps_pos(int iteration);
    void update_eps_neg(int iteration);
    void update_individual_parameters(int iteration);
    void calculate_llik();
    double get_llik();
};

#endif  // CHAIN_H_
