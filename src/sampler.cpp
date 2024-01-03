
#include "sampler.h"

#include "mcmc_utils.h"

#include <Rcpp.h>
#include <Rmath.h>
#include <cmath>
#include <algorithm>
#include <random>
#include <tuple>

std::random_device Sampler::rd;

Sampler::Sampler()
{
    eng = std::ranlux24_base(rd());
    unif_distr = std::uniform_real_distribution<float>(0, 1);
    ber_distr = std::bernoulli_distribution(.5);

    for (int i = 0; i < 128; i++)
    {
        lgamma_lut[i] = std::lgamma(i + 1);
    }
}

float Sampler::dbeta(float x, float alpha, float beta, bool return_log)
{
    return R::dbeta(x, alpha, beta, return_log);
}

float Sampler::dpois(int x, float mean, bool return_log)
{
    return R::dpois(x, mean, return_log);
}

float Sampler::dztpois(int x, float lambda)
{
    return x * std::log(lambda) - std::log(std::exp(lambda) - 1) -
           lgamma_lut[x];
}

float Sampler::dbinom(int x, int size, float prob, bool return_log)
{
    float log_p = x * std::log(prob) + (size - x) * std::log(1 - prob);
    float log_coef = lgamma_lut[size] - lgamma_lut[x] - lgamma_lut[size - x];

    return return_log ? log_p + log_coef : std::exp(log_p + log_coef);
}

float Sampler::dgamma(float x, float shape, float scale, bool return_log)
{
    return R::dgamma(x, shape, scale, return_log);
}

float Sampler::rgamma(float alpha, float beta)
{
    gamma_distr.param(std::gamma_distribution<float>::param_type(alpha, beta));
    float x = gamma_distr(eng);

    if (x < UNDERFLO)
    {
        x = UNDERFLO;
    }
    else if (x > OVERFLO)
    {
        x = OVERFLO;
    }

    return x;
};

float Sampler::rgamma2(float shape, float rate)
{
    return R::rgamma(shape, 1 / rate);
}

std::vector<float> Sampler::rdirichlet(std::vector<float> const &shape_vec)
{
    int n = shape_vec.size();
    std::vector<float> res(n);

    float res_sum = 0;
    for (int i = 0; i < n; i++)
    {
        res[i] = rgamma(shape_vec[i], 1.0);
        res_sum += res[i];
    }

    float res_sum_inv = 1.0 / res_sum;
    for (size_t i = 0; i < res.size(); i++)
    {
        res[i] *= res_sum_inv;
    }

    return res;
};

std::vector<float> Sampler::rlogit_norm(std::vector<float> const &p,
                                        float variance)
{
    int n = p.size() - 1;

    std::vector<float> ret(n + 1);

    float tmp1 = 0;
    for (int i = 0; i < n; i++)
    {
        norm_distr.param(std::normal_distribution<float>::param_type(
            log(p[i] / p[n]), variance));
        ret[i] = exp(norm_distr(eng));
        tmp1 += ret[i];
    }

    float tmp2 = 1.0 / (1.0 + tmp1);
    for (int i = 0; i < n; i++)
    {
        ret[i] *= tmp2;
    }

    ret[n] = tmp2;

    return ret;
}

float Sampler::sample_mean_coi(float mean_shape, float mean_scale)
{
    return rgamma(mean_shape, mean_scale);
}

int Sampler::sample_random_int(int lower, int upper)
{
    unif_int_distr.param(
        std::uniform_int_distribution<>::param_type(lower, upper));
    return unif_int_distr(eng);
}

float Sampler::get_coi_log_prior(int coi, float mean)
{
    return dztpois(coi, mean);
}

float Sampler::get_coi_mean_log_hyper_prior(float mean, float shape,
                                            float scale)
{
    return dgamma(mean, shape, scale, true);
}

int Sampler::sample_coi_delta() { return (2 * ber_distr(eng) - 1); }

int Sampler::sample_coi_delta(float coi_prop_mean)
{
    geom_distr.param(std::geometric_distribution<int>::param_type(
        1.0 / (1.0 + coi_prop_mean)));
    // abs delta >= 1
    return (2 * ber_distr(eng) - 1) * (geom_distr(eng));
}

float Sampler::get_epsilon_log_prior(float x, float alpha, float beta)
{
    // return dgamma(x, shape, scale, true);
    return dbeta(x, alpha, beta, true);
}

float Sampler::get_relatedness_log_prior(float x, float alpha, float beta)
{
    return dbeta(x, alpha, beta, true);
}

float Sampler::sample_epsilon(float curr_epsilon, float variance)
{
    norm_distr.param(
        std::normal_distribution<float>::param_type(curr_epsilon, variance));
    float prop = norm_distr(eng);
    return prop;
};

std::tuple<float, float> Sampler::sample_constrained(float curr, float var,
                                                     float lower, float upper)
{
    norm_distr.param(std::normal_distribution<float>::param_type(0, var));
    float eps = norm_distr(eng);
    float unconstrained = std::log(curr - lower) - std::log(upper - curr);
    float exp_prop = std::exp(eps + unconstrained);
    float prop = (upper * exp_prop + lower) / (exp_prop + 1);
    prop = UtilFunctions::clamp(prop, lower, upper);

    float adj = std::log(prop - lower) + std::log(upper - prop) -
                std::log(curr - lower) - std::log(upper - curr);

    return std::make_tuple(prop, adj);
}

float Sampler::sample_epsilon_pos(float curr_epsilon_pos, float variance)
{
    return sample_epsilon(curr_epsilon_pos, variance);
};

float Sampler::sample_epsilon_neg(float curr_epsilon_neg, float variance)
{
    return sample_epsilon(curr_epsilon_neg, variance);
};

std::vector<float> Sampler::sample_allele_frequencies(
    std::vector<float> const &curr_allele_frequencies, float alpha)
{
    std::vector<float> shape_vec(curr_allele_frequencies.size());

    for (size_t i = 0; i < shape_vec.size(); i++)
    {
        shape_vec[i] = curr_allele_frequencies[i] * alpha;
    }

    return rdirichlet(shape_vec);
};

std::vector<float> Sampler::sample_allele_frequencies2(
    std::vector<float> const &curr_allele_frequencies, float variance)
{
    return rlogit_norm(curr_allele_frequencies, variance);
};

void Sampler::shuffle_vec(std::vector<int> &vec)
{
    std::shuffle(vec.begin(), vec.end(), eng);
}

float Sampler::sample_unif() { return unif_distr(eng); };

float Sampler::sample_log_mh_acceptance() { return std::log(unif_distr(eng)); };

LatentGenotype Sampler::sample_latent_genotype(
    const std::vector<int> &obs_genotype, int coi, float epsilon_pos,
    float epsilon_neg)
{
    int total_alleles = obs_genotype.size();
    int total_obs_positives = 0;
    int total_obs_negatives = 0;
    std::vector<int> obs_positive_indices{};
    std::vector<int> obs_negative_indices{};

    for (int ii = 0; ii < total_alleles; ++ii)
    {
        if (obs_genotype[ii] == 1)
        {
            total_obs_positives++;
            obs_positive_indices.push_back(ii);
        }
        else
        {
            total_obs_negatives++;
            obs_negative_indices.push_back(ii);
        }
    }

    // there must be at least one allele, so if all obs_positives are considered
    // false positives then there must be at least one false negative
    int min_false_negatives = std::max(0, 1 * (total_obs_positives == 0));

    int max_false_negatives =
        std::max(min_false_negatives, std::min(coi, total_obs_negatives));

    int total_false_negatives = min_false_negatives;
    for (int ii = total_false_negatives; ii < max_false_negatives; ++ii)
    {
        total_false_negatives +=
            (sample_unif() < (epsilon_neg / obs_genotype.size()));
    }
    int total_true_negatives = total_obs_negatives - total_false_negatives;

    float log_prob_total_false_negatives =
        std::log(boost::math::binomial_coefficient<float>(
            total_obs_negatives, total_false_negatives - min_false_negatives)) +
        (total_false_negatives - min_false_negatives) *
            std::log(epsilon_neg / obs_genotype.size()) +
        total_true_negatives *
            std::log(1 - (epsilon_neg / obs_genotype.size()));

    // if the observed number of positives exceeds the COI, then some number of
    // them must be false positives
    int min_false_positives =
        std::max(0, (total_obs_positives + total_false_negatives) - coi);

    int max_false_positives =
        std::min(total_obs_positives, total_false_negatives / 2);

    int total_false_positives = min_false_positives;
    for (int ii = total_false_positives; ii < max_false_positives; ++ii)
    {
        total_false_positives +=
            (sample_unif() < (epsilon_pos / obs_genotype.size()));
    }
    int total_true_positives = total_obs_positives - total_false_positives;

    float log_prob_total_false_positives =
        std::log(boost::math::binomial_coefficient<float>(
            total_obs_positives, total_false_positives - min_false_positives)) +
        (total_false_positives - min_false_positives) *
            std::log(epsilon_pos / obs_genotype.size()) +
        total_true_positives *
            std::log(1 - (epsilon_pos / obs_genotype.size()));

    shuffle_vec(obs_positive_indices);
    shuffle_vec(obs_negative_indices);

    std::vector<int> allele_index_vec{};
    allele_index_vec.insert(
        allele_index_vec.end(), obs_positive_indices.begin(),
        obs_positive_indices.begin() + total_true_positives);

    allele_index_vec.insert(
        allele_index_vec.end(), obs_negative_indices.begin(),
        obs_negative_indices.begin() + total_false_negatives);

    std::sort(allele_index_vec.begin(), allele_index_vec.end());

    assert(allele_index_vec.size() ==
           total_true_positives + total_false_negatives);

    float log_prob_positive_indices =
        -std::log(boost::math::binomial_coefficient<float>(
            total_obs_positives, total_true_positives));
    float log_prob_negative_indices =
        -std::log(boost::math::binomial_coefficient<float>(
            total_obs_negatives, total_false_negatives));

    float log_prob = log_prob_positive_indices + log_prob_negative_indices +
                     log_prob_total_false_positives +
                     log_prob_total_false_negatives;

    return LatentGenotype{allele_index_vec, log_prob};
}
