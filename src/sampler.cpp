
#include <random>
#include <algorithm>
#include <Rcpp.h>
#include "sampler.h"
#include "mcmc_utils.h"

std::random_device Sampler::rd;

Sampler::Sampler(int genotype_sampling_depth, std::vector<int> const &num_alleles) {
    eng = std::ranlux24_base(rd());
    unif_distr = std::uniform_real_distribution<double>(0, 1);
    ber_distr = std::bernoulli_distribution(.5);

    for(size_t i = 0; i < num_alleles.size(); i++)
    {
        int allele_key = num_alleles[i];
        auto search = genotype_samples.find(allele_key);
        if (search == genotype_samples.end()) {
            genotype_samples[allele_key] = std::vector<std::vector<int > >(genotype_sampling_depth, std::vector<int>(allele_key));
        }
    }
}

double Sampler::dbeta(double x, double alpha, double beta, bool return_log) {
    return R::dbeta(x, alpha, beta, return_log);
}

double Sampler::dpois(int x, double mean, bool return_log) {
    return R::dpois(x, mean, return_log);
}

double Sampler::rgamma(double alpha, double beta) {
    gamma_distr.param(std::gamma_distribution<double>::param_type(alpha, beta));
    double x = gamma_distr(eng);

    if (x < UNDERFLO) {
        x = UNDERFLO;
    } else if (x > OVERFLO) {
        x = OVERFLO;
    }

    return x;
};

double Sampler::rgamma2(double shape, double rate) {
    return R::rgamma(shape, 1 / rate);
}

std::vector<double> Sampler::rdirichlet(std::vector<double> const &shape_vec) {
    int n = shape_vec.size();
    std::vector<double> res(n);

    double res_sum = 0;
    for(int i = 0; i < n; i++) {
        res[i] = rgamma(shape_vec[i], 1.0);
        res_sum += res[i];
    }

    double res_sum_inv = 1.0 / res_sum;
    for(size_t i = 0; i < res.size(); i++) {
        res[i] *= res_sum_inv;
    }

    return res;
    
};

std::vector<double> Sampler::rlogit_norm(std::vector<double> const &p, double variance) {
    int n = p.size() - 1;
    
    std::vector<double> ret(n+1);
    
    double tmp1 = 0;  
    for (int i = 0; i < n; i++) {
        norm_distr.param(std::normal_distribution<double>::param_type(log(p[i] / p[n]), variance));
        ret[i] = exp(norm_distr(eng));
        tmp1 += ret[i];
    }
    
    double tmp2 = 1.0 / (1.0 + tmp1);
    for (int i = 0; i < n; i++) {
        ret[i] *= tmp2;
    }
    
    ret[n] = tmp2;
    
    return ret;
}

// int Sampler::sample_coi(int curr_coi, int delta, int max_coi) {
//     unif_int_distr.param(std::uniform_int_distribution<int>::param_type(curr_coi - delta, std::min(curr_coi + delta, max_coi)));
//     return unif_int_distr(eng);
// };


double Sampler::sample_mean_coi(double mean_shape, double mean_rate) {
    return rgamma2(mean_shape, mean_rate) + 1;
}

double Sampler::get_coi_log_prior(int coi, double mean) {
    return dpois(coi, mean, true);
}

int Sampler::sample_coi_delta(double coi_prop_mean) {
    geom_distr.param(std::geometric_distribution<int>::param_type(1.0 / (1.0 + coi_prop_mean)));
    return (2 * ber_distr(eng) - 1) * geom_distr(eng);
}

double Sampler::get_epsilon_log_prior(double x, double alpha, double beta) {
    return dbeta(x, alpha, beta, true);
}

// double Sampler::sample_epsilon(double curr_epsilon, double variance) {
//     norm_distr.param(std::normal_distribution<double>::param_type(log(curr_epsilon / (1 - curr_epsilon)), variance));
//     double prop = norm_distr(eng);
//     return exp(prop) / (1 + exp(prop));
// };

double Sampler::sample_epsilon(double curr_epsilon, double variance) {
    norm_distr.param(std::normal_distribution<double>::param_type(curr_epsilon, variance));
    double prop = norm_distr(eng);
    return prop;
};

double Sampler::sample_epsilon_pos(double curr_epsilon_pos, double variance) {
    return sample_epsilon(curr_epsilon_pos, variance);
};

double Sampler::sample_epsilon_neg(double curr_epsilon_neg, double variance) {
    return sample_epsilon(curr_epsilon_neg, variance);
};

std::vector<double> Sampler::sample_allele_frequencies(std::vector<double> const &curr_allele_frequencies, double alpha) {
    std::vector<double> shape_vec(curr_allele_frequencies.size());

    for(size_t i = 0; i < shape_vec.size(); i++) {
        shape_vec[i] = curr_allele_frequencies[i] * alpha;
    }

    return rdirichlet(shape_vec);
};

std::vector<double> Sampler::sample_allele_frequencies2(std::vector<double> const &curr_allele_frequencies, double variance) {
    return rlogit_norm(curr_allele_frequencies, variance);
};

std::vector<std::vector<int > >& Sampler::sample_genotype(int coi, std::vector<double> const &allele_frequencies, int num_samples) {
    discrete_distr.param(std::discrete_distribution<int>::param_type(allele_frequencies.begin(), allele_frequencies.end()));
    for(size_t i = 0; i < num_samples; i++)
    {
        std::fill(genotype_samples[allele_frequencies.size()][i].begin(), genotype_samples[allele_frequencies.size()][i].end(), 0);
        for(int j = 0; j < coi; j++)
        {
            genotype_samples[allele_frequencies.size()][i][discrete_distr(eng)] += 1;
        }
    }
    return genotype_samples[allele_frequencies.size()];
}

double Sampler::sample_log_mh_acceptance() {
    return log(unif_distr(eng));
};