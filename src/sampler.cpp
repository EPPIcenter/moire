#include "sampler.h"

#include "mcmc_utils.h"

#include <Rcpp.h>
#include <Rmath.h>
#include <cmath>
#include <algorithm>
#include <random>
#include <tuple>
#include <span>

std::random_device Sampler::rd;

Sampler::Sampler()
{
    eng = std::ranlux24_base(rd());
    unif_distr = std::uniform_real_distribution<float>(0, 1);
    ber_distr = std::bernoulli_distribution(.5);

    for (std::size_t i = 0; i < std::size(lgamma_lut); i++)
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

float Sampler::dztnbinom(int x, float p, float r)
{
    return lgamma_lut[x + r] + r * std::log(p) + x * std::log(1 - p) -
           lgamma_lut[r] - lgamma_lut[x + 1] - std::log(1 - std::pow(p, r));
}

float Sampler::dbinom(int x, int size, float prob)
{
    float log_p = x * std::log(prob) + (size - x) * std::log(1 - prob);
    float log_coef = lgamma_lut[size] - lgamma_lut[x] - lgamma_lut[size - x];

    return log_p + log_coef;
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

float Sampler::sample_gamma(float mean_shape, float mean_scale)
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

float Sampler::get_coi_log_prior(int coi, float p, float r)
{
    return dztnbinom(coi, p, r);
}

float Sampler::get_coi_mean_log_hyper_prior(float mean, float shape,
                                            float scale)
{
    return dgamma(mean, shape, scale, true);
}

float Sampler::get_gamma_log_prior(float variance, float shape,
                                                float scale)
{
    return dgamma(variance, shape, scale, true);
}

int Sampler::sample_coi_delta() { return (2 * ber_distr(eng) - 1); }

int Sampler::sample_coi_delta(float coi_prop_mean)
{
    geom_distr.param(std::geometric_distribution<int>::param_type(1.0 / (1.0 + coi_prop_mean)));
    return (2 * ber_distr(eng) - 1) * (geom_distr(eng));
}

float Sampler::get_beta_log_prior(float x, float alpha, float beta)
{
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

std::vector<float> Sampler::sample_dirichlet(
    std::vector<float> const &curr_allele_frequencies, float alpha)
{
    std::vector<float> shape_vec;
    shape_vec.reserve(curr_allele_frequencies.size());

    std::transform(curr_allele_frequencies.begin(), curr_allele_frequencies.end(), 
                   std::back_inserter(shape_vec), 
                   [alpha](float freq) { return freq * alpha; });

    return rdirichlet(shape_vec);
};

float Sampler::unnormalized_dirichlet_log_prior(std::span<float> x, std::span<float> alpha)
{
    float res = 0.0f;
    const size_t n = x.size();
    
    switch (n) {
        case 1:
            return (alpha[0] - 1.0f) * std::log(x[0]);
        case 2:
            return (alpha[0] - 1.0f) * std::log(x[0]) +
                   (alpha[1] - 1.0f) * std::log(x[1]);
        case 3:
            return (alpha[0] - 1.0f) * std::log(x[0]) +
                   (alpha[1] - 1.0f) * std::log(x[1]) +
                   (alpha[2] - 1.0f) * std::log(x[2]);
    }
    
    // For sizes 4 and above, use the unrolled loop
    size_t i = 0;
    for (; i + 3 < n; i += 4) {
        res += (alpha[i] - 1.0f) * std::log(x[i]);
        res += (alpha[i+1] - 1.0f) * std::log(x[i+1]);
        res += (alpha[i+2] - 1.0f) * std::log(x[i+2]);
        res += (alpha[i+3] - 1.0f) * std::log(x[i+3]);
    }
    
    // Handle remaining elements
    for (; i < n; ++i) {
        res += (alpha[i] - 1.0f) * std::log(x[i]);
    }
    
    return res;
}

float Sampler::dirichlet_log_prior(std::span<float> x, std::span<float> alpha)
{
    return unnormalized_dirichlet_log_prior(x, alpha) - 
        std::lgamma(std::reduce(alpha.begin(), alpha.end())) + 
        std::transform_reduce(alpha.begin(), alpha.end(), 0.0f, std::plus<float>(), [](float a) { return std::lgamma(a); });
}

std::vector<float> Sampler::sample_logit_norm(
    std::vector<float> const &curr_allele_frequencies, float variance)
{
    return rlogit_norm(curr_allele_frequencies, variance);
};

float Sampler::sample_unif() { return unif_distr(eng); };

float Sampler::sample_log_mh_acceptance() { return std::log(unif_distr(eng)); };

std::vector<int> Sampler::sample_random_sequence(int min, int max)
{
    std::vector<int> indices(max - min);
    std::iota(std::begin(indices), std::end(indices), min);
    std::shuffle(indices.begin(), indices.end(), eng);
    return indices;
}

