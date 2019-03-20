#pragma once

#ifndef SAMPLER_H_
#define SAMPLER_H_

#include <random>


class Sampler {
private:
    static std::random_device rd; // obtain a random number from hardware
    static std::ranlux24_base eng; // seed the generator

    static std::uniform_int_distribution<int> unif_int_distr;
    static std::normal_distribution<double> norm_distr;
    static std::gamma_distribution<double> gamma_distr;
    static std::discrete_distribution<int> discrete_distr;
    static std::uniform_real_distribution<double> unif_distr;
    static std::bernoulli_distribution ber_distr;
    static std::geometric_distribution<int> geom_distr;

    double rgamma(double shape, double rate);
    std::vector<double> rdirichlet(std::vector<double> const &shape_vec);

    static std::map<int, std::vector<std::vector<int> > > genotype_samples;

public:

    double sample_epsilon(double curr_epsilon, double variance);
    double sample_epsilon_pos(double curr_epsilon_pos, double variance);
    double sample_epsilon_neg(double curr_epsilon_neg, double variance);
    int sample_coi(int curr_coi, int delta, int max_coi);
    int sample_coi_delta(double coi_prop_mean);
    std::vector<double> sample_allele_frequencies(std::vector<double> const &curr_allele_frequencies, double alpha);
    std::vector<std::vector<int> >& sample_genotype(int coi, std::vector<double> const &allele_frequencies, int num_samples);
    double sample_log_mh_acceptance();
    double runif_0_1();

    Sampler() { };
    Sampler(int genotype_sample_depth, std::vector<int> const &num_alleles);
    // Sampler(int seed);

};

#endif // SAMPLER_H_