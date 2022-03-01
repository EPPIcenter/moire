#ifndef SAMPLER_H_
#define SAMPLER_H_

#include "lookup.h"

#include <algorithm>

#include <boost/random.hpp>

#include <random>

struct LatentGenotype
{
    std::vector<int> value;
    double log_prob;
};

class Sampler
{
   private:
    std::uniform_int_distribution<int> unif_int_distr;
    std::normal_distribution<double> norm_distr;
    std::gamma_distribution<double> gamma_distr;
    std::discrete_distribution<int> discrete_distr;
    std::uniform_real_distribution<double> unif_distr;
    std::bernoulli_distribution ber_distr;
    std::geometric_distribution<int> geom_distr;

    double dbeta(double x, double alpha, double beta, bool return_log);
    double dpois(int x, double mean, bool return_log);
    double dztpois(int x, double mean);
    double dgamma(double x, double shape, double scale, bool return_log);
    double rgamma(double alpha, double beta);
    double rgamma2(double shape, double rate);

    std::vector<double> rdirichlet(std::vector<double> const &shape_vec);
    std::vector<double> rlogit_norm(std::vector<double> const &p,
                                    double variance);

   public:
    static std::random_device rd;
    std::ranlux24_base eng;
    boost::random::mt19937 r;

    double get_epsilon_log_prior(double x, double alpha, double beta);
    double get_coi_log_prob(int coi, double mean);
    double get_coi_mean_log_prior(double mean, double shape, double scale);

    double sample_epsilon(double curr_epsilon, double variance);
    std::tuple<double, double> sample_constrained(double curr, double var,
                                                  double lower, double upper);
    double sample_epsilon_pos(double curr_epsilon_pos, double variance);
    double sample_epsilon_neg(double curr_epsilon_neg, double variance);

    int sample_coi_delta(double coi_prop_mean);
    int sample_coi_delta();
    double sample_mean_coi(double coi_mean_shape, double coi_mean_rate);

    int sample_random_int(int lower, int upper);
    std::vector<double> sample_allele_frequencies(
        std::vector<double> const &curr_allele_frequencies, double alpha);
    std::vector<double> sample_allele_frequencies2(
        std::vector<double> const &curr_allele_frequencies, double variance);
    std::vector<std::vector<int>> &sample_genotype(
        int coi, std::vector<double> const &allele_frequencies,
        int num_samples);

    double sample_log_mh_acceptance();
    double sample_unif();
    void shuffle_vec(std::vector<int> &vec);

    LatentGenotype sample_latent_genotype(const std::vector<int> &obs_genotype,
                                          int coi, double epsilon_pos,
                                          double epsilon_neg);

    Sampler();
};

#endif  // SAMPLER_H_
