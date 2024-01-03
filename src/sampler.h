#ifndef SAMPLER_H_
#define SAMPLER_H_

#include <algorithm>

#include <boost/random.hpp>

#include <random>

struct LatentGenotype
{
    std::vector<int> value;
    float log_prob;
};

class Sampler
{
   private:
    std::uniform_int_distribution<int> unif_int_distr;
    std::normal_distribution<float> norm_distr;
    std::gamma_distribution<float> gamma_distr;
    std::discrete_distribution<int> discrete_distr;
    std::uniform_real_distribution<float> unif_distr;
    std::bernoulli_distribution ber_distr;
    std::geometric_distribution<int> geom_distr;
    std::array<float, 128> lgamma_lut;

    std::vector<float> rdirichlet(std::vector<float> const &shape_vec);
    std::vector<float> rlogit_norm(std::vector<float> const &p, float variance);

   public:
    static std::random_device rd;
    std::ranlux24_base eng;
    boost::random::mt19937 r;

    float get_epsilon_log_prior(float x, float alpha, float beta);
    float get_relatedness_log_prior(float x, float alpha, float beta);
    float get_coi_log_prior(int coi, float mean);
    float get_coi_mean_log_hyper_prior(float mean, float shape, float scale);

    float sample_epsilon(float curr_epsilon, float variance);
    std::tuple<float, float> sample_constrained(float curr, float var,
                                                float lower, float upper);
    float sample_epsilon_pos(float curr_epsilon_pos, float variance);
    float sample_epsilon_neg(float curr_epsilon_neg, float variance);

    int sample_coi_delta(float coi_prop_mean);
    int sample_coi_delta();
    float sample_mean_coi(float coi_mean_shape, float coi_mean_rate);

    int sample_random_int(int lower, int upper);
    std::vector<float> sample_allele_frequencies(
        std::vector<float> const &curr_allele_frequencies, float alpha);
    std::vector<float> sample_allele_frequencies2(
        std::vector<float> const &curr_allele_frequencies, float variance);
    std::vector<std::vector<int>> &sample_genotype(
        int coi, std::vector<float> const &allele_frequencies, int num_samples);

    float sample_log_mh_acceptance();
    float sample_unif();
    void shuffle_vec(std::vector<int> &vec);

    float dbeta(float x, float alpha, float beta, bool return_log);
    float dpois(int x, float mean, bool return_log);
    float dbinom(int x, int size, float prob, bool return_log);
    float dztpois(int x, float mean);
    float dgamma(float x, float shape, float scale, bool return_log);
    float rgamma(float alpha, float beta);
    float rgamma2(float shape, float rate);

    LatentGenotype sample_latent_genotype(const std::vector<int> &obs_genotype,
                                          int coi, float epsilon_pos,
                                          float epsilon_neg);

    Sampler();
};

#endif  // SAMPLER_H_
