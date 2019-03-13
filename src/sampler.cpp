
#include <random>
#include <algorithm>
#include <Rcpp.h>
#include "sampler.h"
#include "mcmc_utils.h"

std::random_device Sampler::rd;
std::ranlux24_base Sampler::eng(rd());

std::uniform_int_distribution<int> Sampler::unif_int_distr;
std::normal_distribution<double> Sampler::norm_distr;
std::gamma_distribution<double> Sampler::gamma_distr;
std::discrete_distribution<int> Sampler::discrete_distr;
std::uniform_real_distribution<double> Sampler::unif_distr(0, 1);

Sampler::Sampler(int seed) {
    eng.seed(seed);
};

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

int Sampler::sample_coi(int curr_coi, int delta, int max_coi) {
    unif_int_distr.param(std::uniform_int_distribution<int>::param_type(curr_coi - delta, std::min(curr_coi + delta, max_coi)));
    return unif_int_distr(eng);
};

double Sampler::sample_epsilon(double curr_epsilon, double variance) {
    norm_distr.param(std::normal_distribution<double>::param_type(curr_epsilon, variance));
    return norm_distr(eng);
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

std::vector<std::vector<int> > Sampler::sample_genotype(int coi, std::vector<double> const &allele_frequencies, int num_samples) {
    discrete_distr.param(std::discrete_distribution<int>::param_type(allele_frequencies.begin(), allele_frequencies.end()));
    std::vector<std::vector<int> > res(num_samples);
    for(size_t i = 0; i < num_samples; i++)
    {
        res[i] = std::vector<int>(allele_frequencies.size(), 0);
        for(int j = 0; j < coi; j++)
        {
            res[i][discrete_distr(eng)] += 1;
        }
    }
    return res;
    
}

double Sampler::sample_log_mh_acceptance() {
    return UtilFunctions::fastlog(unif_distr(eng));
};