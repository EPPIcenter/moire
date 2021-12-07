
#include "sampler.h"

#include "mcmc_utils.h"

#include <Rcpp.h>
#include <Rmath.h>
#include <algorithm>
#include <random>
#include <tuple>

std::random_device Sampler::rd;

Sampler::Sampler(Lookup lookup) : lookup(lookup)
{
    eng = std::ranlux24_base(rd());
    unif_distr = std::uniform_real_distribution<double>(0, 1);
    ber_distr = std::bernoulli_distribution(.5);
}

double Sampler::dbeta(double x, double alpha, double beta, bool return_log)
{
    return R::dbeta(x, alpha, beta, return_log);
}

double Sampler::dpois(int x, double mean, bool return_log)
{
    return R::dpois(x, mean, return_log);
}

double Sampler::dztpois(int x, double lambda)
{
    return x * std::log(lambda) - std::log(std::exp(lambda) - 1) -
           std::lgamma(x + 1);
}

double Sampler::dgamma(double x, double shape, double scale, bool return_log)
{
    return R::dgamma(x, shape, scale, return_log);
}

double Sampler::rgamma(double alpha, double beta)
{
    gamma_distr.param(std::gamma_distribution<double>::param_type(alpha, beta));
    double x = gamma_distr(eng);

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

double Sampler::rgamma2(double shape, double rate)
{
    return R::rgamma(shape, 1 / rate);
}

std::vector<double> Sampler::rdirichlet(std::vector<double> const &shape_vec)
{
    int n = shape_vec.size();
    std::vector<double> res(n);

    double res_sum = 0;
    for (int i = 0; i < n; i++)
    {
        res[i] = rgamma(shape_vec[i], 1.0);
        res_sum += res[i];
    }

    double res_sum_inv = 1.0 / res_sum;
    for (size_t i = 0; i < res.size(); i++)
    {
        res[i] *= res_sum_inv;
    }

    return res;
};

std::vector<double> Sampler::rlogit_norm(std::vector<double> const &p,
                                         double variance)
{
    int n = p.size() - 1;

    std::vector<double> ret(n + 1);

    double tmp1 = 0;
    for (int i = 0; i < n; i++)
    {
        norm_distr.param(std::normal_distribution<double>::param_type(
            log(p[i] / p[n]), variance));
        ret[i] = exp(norm_distr(eng));
        tmp1 += ret[i];
    }

    double tmp2 = 1.0 / (1.0 + tmp1);
    for (int i = 0; i < n; i++)
    {
        ret[i] *= tmp2;
    }

    ret[n] = tmp2;

    return ret;
}

double Sampler::sample_mean_coi(double mean_shape, double mean_rate)
{
    return rgamma2(mean_shape, mean_rate) + 1;
}

int Sampler::sample_random_int(int lower, int upper)
{
    unif_int_distr.param(
        std::uniform_int_distribution<>::param_type(lower, upper));
    return unif_int_distr(eng);
}

double Sampler::get_coi_log_prob(int coi, double mean)
{
    return dztpois(coi, mean);
}

double Sampler::get_coi_mean_log_prior(double mean, double shape, double scale)
{
    return dgamma(mean, shape, scale, true);
}

int Sampler::sample_coi_delta() { return (2 * ber_distr(eng) - 1); }

int Sampler::sample_coi_delta(double coi_prop_mean)
{
    geom_distr.param(std::geometric_distribution<int>::param_type(
        1.0 / (1.0 + coi_prop_mean)));
    // abs delta >= 1
    return (2 * ber_distr(eng) - 1) * (geom_distr(eng));
}

double Sampler::get_epsilon_log_prior(double x, double shape, double scale)
{
    return dgamma(x, shape, scale, true);
}

double Sampler::sample_epsilon(double curr_epsilon, double variance)
{
    norm_distr.param(
        std::normal_distribution<double>::param_type(curr_epsilon, variance));
    double prop = norm_distr(eng);
    return prop;
};

std::tuple<double, double> Sampler::sample_constrained(double curr, double var,
                                                       double lower,
                                                       double upper)
{
    norm_distr.param(std::normal_distribution<double>::param_type(0, var));
    double eps = norm_distr(eng);
    double unconstrained = std::log(curr - lower) - std::log(upper - curr);
    double exp_prop = std::exp(eps + unconstrained);
    double prop = (upper * exp_prop + lower) / (exp_prop + 1);
    prop = UtilFunctions::clamp(prop, lower, upper);

    double adj = std::log(prop - lower) + std::log(upper - prop) -
                 std::log(curr - lower) - std::log(upper - curr);
    return std::make_tuple(prop, adj);
}

double Sampler::sample_epsilon_pos(double curr_epsilon_pos, double variance)
{
    return sample_epsilon(curr_epsilon_pos, variance);
};

double Sampler::sample_epsilon_neg(double curr_epsilon_neg, double variance)
{
    return sample_epsilon(curr_epsilon_neg, variance);
};

std::vector<double> Sampler::sample_allele_frequencies(
    std::vector<double> const &curr_allele_frequencies, double alpha)
{
    std::vector<double> shape_vec(curr_allele_frequencies.size());

    for (size_t i = 0; i < shape_vec.size(); i++)
    {
        shape_vec[i] = curr_allele_frequencies[i] * alpha;
    }

    return rdirichlet(shape_vec);
};

std::vector<double> Sampler::sample_allele_frequencies2(
    std::vector<double> const &curr_allele_frequencies, double variance)
{
    return rlogit_norm(curr_allele_frequencies, variance);
};

void Sampler::shuffle_vec(std::vector<int> &vec)
{
    std::shuffle(vec.begin(), vec.end(), eng);
}

double Sampler::sample_unif() { return unif_distr(eng); };

double Sampler::sample_log_mh_acceptance()
{
    return std::log(unif_distr(eng));
};

LatentGenotype Sampler::sample_latent_genotype(
    const std::vector<int> &obs_genotype, int coi, double epsilon_pos,
    double epsilon_neg)
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

    double log_prob_total_false_negatives =
        std::log(boost::math::binomial_coefficient<double>(
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

    double log_prob_total_false_positives =
        std::log(boost::math::binomial_coefficient<double>(
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

    double log_prob_positive_indices =
        -std::log(boost::math::binomial_coefficient<double>(
            total_obs_positives, total_true_positives));
    double log_prob_negative_indices =
        -std::log(boost::math::binomial_coefficient<double>(
            total_obs_negatives, total_false_negatives));

    double log_prob = log_prob_positive_indices + log_prob_negative_indices +
                      log_prob_total_false_positives +
                      log_prob_total_false_negatives;

    return LatentGenotype{allele_index_vec, log_prob};
}
