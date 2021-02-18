#pragma once

#ifndef CHAIN_H_
#define CHAIN_H_

#include <Rcpp.h>
#include "genotyping_data.h"
#include "lookup.h"
#include "parameters.h"
#include "sampler.h"

class Chain {
private:
    GenotypingData genotyping_data;
    Lookup lookup;
    Parameters params;
    Sampler sampler;

    std::vector<std::vector<int>> importance_sample;
    std::tuple<int, int> genotype_cache_key;


    void initialize_p();
    void initialize_m();
    void initialize_mean_coi();
    void initialize_eps_neg();
    void initialize_eps_pos();
    void initialize_likelihood();

    std::vector<double> reweight_allele_frequencies(std::vector<double> const &allele_frequencies, std::vector<int> const &observed_genotype, double epsilon_neg, double epsilon_pos, int coi);
    std::vector<double> calc_genotype_log_pmf(std::vector<std::vector<int> > const &genotypes, int coi, std::vector<double> const &allele_frequencies, int num_genotypes);
    std::vector<double> calc_obs_genotype_lliks(std::vector<int> const &obs_genotype, std::vector<std::vector<int> > const &true_genotypes, double epsilon_neg, double epsilon_pos, int num_genotypes);
    long double calc_genotype_marginal_llik(std::vector<int> const &obs_genotype, int coi, std::vector<double> const &allele_frequencies, double epsilon_neg, double epsilon_pos);
    long double calc_exact_genotype_marginal_llik(std::vector<int> const &obs_genotype, int coi, std::vector<double> const &allele_frequencies, double epsilon_neg, double epsilon_pos);
    long double calc_estimated_genotype_marginal_llik(std::vector<int> const &obs_genotype, int coi, std::vector<double> const &allele_frequencies, double epsilon_neg, double epsilon_pos);
    void generate_possible_genotypes_helper(std::vector<int> &chosen, std::vector<int> &arr, int index, int r, int n, int start, int end);
    void generate_possible_genotypes(int coi, int total_alleles);

public:
    std::map<std::tuple<int, int>, std::vector<std::vector<int> > > true_genotypes_cache;
    std::vector<std::vector<double> > llik_old;
    std::vector<std::vector<double> > llik_new;
    double llik;

    // Mean COI
    // TODO: Allow for other priors on complexity of infection
    // std::string prior;
    // double poisson_prior_lambda;
    double mean_coi;

    // COI
    std::vector<int> m;
    int prop_m;
    std::vector<int> m_accept;
    std::vector<double> m_prop_mean;

    // Allele Frequencies
    std::vector<std::vector<double> > p;
    std::vector<double> prop_p;
    std::vector<double> p_prop_var;
    std::vector<int> p_accept;

    // Epsilon Positive
    // double eps_pos;
    std::vector<double> eps_pos;
    double prop_eps_pos;
    std::vector<int> eps_pos_accept;
    double eps_pos_var;

    // Epsilon Negative
    // double eps_neg;
    std::vector<double> eps_neg;
    double prop_eps_neg;
    std::vector<int> eps_neg_accept;
    double eps_neg_var;

    Chain() {};
    Chain(GenotypingData genotyping_data, Lookup lookup, Parameters params);

    void update_m(int iteration);
    void update_mean_coi(int iteration);
    void update_p(int iteration);
    void update_eps(int iteration);
    void update_eps_pos(int iteration);
    void update_eps_neg(int iteration);
    void calculate_llik();
    double get_llik();

};

#endif // CHAIN_H_
