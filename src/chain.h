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

    void initialize_p();
    void initialize_m();
    void initialize_eps_neg();
    void initialize_eps_pos();
    void initialize_likelihood();
    void initialize_sampler();

    std::vector<double> reweight_allele_frequencies(std::vector<double> const &allele_frequencies, std::vector<int> const &observed_genotype, double epsilon_neg, double epsilon_pos);
    std::vector<double> calc_genotype_log_pmf(std::vector<std::vector<int> > const &genotypes, int coi, std::vector<double> const &allele_frequencies);
    std::vector<double> calc_obs_genotype_lliks(std::vector<int> const &obs_genotype, std::vector<std::vector<int> > const &true_genotypes, double epsilon_neg, double epsilon_pos);
    long double calc_genotype_marginal_llik(std::vector<int> const &obs_genotype, int coi, std::vector<double> const &allele_frequencies, double epsilon_neg, double epsilon_pos, int importance_sampling_depth);

public:
    std::vector<std::vector<double> > llik_old;
    std::vector<std::vector<double> > llik_new;
    double llik;

    // COI
    std::vector<int> m;
    int prop_m;
    std::vector<int> m_accept;
    std::vector<double> m_prop_mean;

    // Allele Frequencies
    std::vector<std::vector<double> > p;
    std::vector<double> prop_p;
    std::vector<int> p_accept;

    // Epsilon Positive
    double eps_pos;
    // double prop_eps_pos;
    int eps_pos_accept;
    double eps_pos_var;

    // Epsilon Negative
    double eps_neg;
    // double prop_eps_neg;
    int eps_neg_accept;
    double eps_neg_var;

    Chain() {};
    Chain(GenotypingData genotyping_data, Lookup lookup, Parameters params);

    void update_m(int iteration);
    void update_p();
    void update_eps_pos(int iteration);
    void update_eps_neg(int iteration);
    void calculate_llik();
    double get_llik();

};

#endif // CHAIN_H_