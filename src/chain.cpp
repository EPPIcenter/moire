#include <algorithm>

#include "chain.h"
#include "sampler.h"
#include "mcmc_utils.h"


// Initialize P with empirical allele frequencies
void Chain::initialize_p() {
    p_accept.resize(genotyping_data.num_loci, 0);

    std::vector<std::vector<int> > total_locus_alleles(genotyping_data.num_loci);
    std::vector<int> total_alleles(genotyping_data.num_loci);

    for (size_t i = 0; i < genotyping_data.num_loci; i++) {
        total_locus_alleles.push_back(std::vector<int>(genotyping_data.num_alleles[i]));
        total_alleles.push_back(0);
        p.push_back(std::vector<double>(genotyping_data.num_alleles[i]));

        auto obs_locus_genotypes = genotyping_data.observed_alleles[i];

        for (size_t j = 0; j < obs_locus_genotypes.size(); j++) {
            auto sample_genotype = obs_locus_genotypes[j];
            for (size_t k = 0; k < sample_genotype.size(); k++) {

                if(j == 0) {
                    total_locus_alleles[i].push_back(0);
                }

                total_locus_alleles[i][k] += sample_genotype[k];
                total_alleles[i] += sample_genotype[k];
            }
        }
        for(size_t j = 0; j < genotyping_data.num_alleles[i]; j++) {
            p[i][j] = (total_locus_alleles[i][j] + 1) / ((double) total_alleles[i] + genotyping_data.num_alleles[i]); // Make sure at least 1 allele everywhere
        }
    }
};

void Chain::initialize_m() {
    m = genotyping_data.observed_coi;
    m_accept.resize(genotyping_data.num_samples, 0);
}

void Chain::initialize_eps_neg() {
    eps_neg = params.eps_neg_0;
    eps_neg_accept = 0;
}

void Chain::initialize_eps_pos() {
    eps_pos = params.eps_pos_0;
    eps_pos_accept = 0;
}

void Chain::initialize_mean_coi() {
    mean_coi = 0;
    for(size_t i = 0; i < m.size(); i++) {
        mean_coi += m[i];
    };
    mean_coi = mean_coi / m.size();
}

// void Chain::initialize_sampler() {
//     sampler = Sampler(params.importance_sampling_depth, genotyping_data.num_alleles);
// }

void Chain::update_m(int iteration) {
    for(int i = 0; i < genotyping_data.num_samples; i++) {
        int prop_m =  m[i] + sampler.sample_coi_delta(m_prop_mean[i]);

        // Accept automatically if COI is unchanged
        if(prop_m == m[i]) {
            m_prop_mean[i] += (1 - 0.23) / sqrt(double(iteration));
            m_accept[i] += 1;
            continue;
        }

        if(params.max_coi >= prop_m && prop_m > 0) {
            double sum_can = 0;
            double sum_orig = 0;
            for(size_t j = 0; j < genotyping_data.num_loci; j++) {
                llik_new[j][i] = calc_genotype_marginal_llik(genotyping_data.observed_alleles[j][i], prop_m, p[j], eps_neg, eps_pos);
                sum_can += llik_new[j][i];
                sum_orig += llik_old[j][i];
            }
            sum_can += sampler.get_coi_log_prior(prop_m, mean_coi);
            sum_orig += sampler.get_coi_log_prior(m[i], mean_coi);

            // Accept
            if(sampler.sample_log_mh_acceptance() <= (sum_can - sum_orig)) {
                // UtilFunctions::print("Accept M", m[i], prop_m, m_prop_mean[i]);
                m[i] = prop_m;
                m_prop_mean[i] += (1 - 0.23) / sqrt(double(iteration));
                m_accept[i] += 1;
                for(size_t j = 0; j < genotyping_data.num_loci; j++) {
                    llik_old[j][i] = llik_new[j][i];
                }
            } else {
                m_prop_mean[i] -= 0.23 / sqrt(double(iteration));
                m_prop_mean[i] = (m_prop_mean[i] < 0) ? 0 : m_prop_mean[i];
            }
        }

    }

}

// Gibbs update
void Chain::update_mean_coi(int iteration) {
    double coi_mean_shape = .25;
    double coi_mean_rate = .25;

    for(size_t i = 0; i < m.size(); i++)
    {
        coi_mean_shape += m[i] - 1;
        coi_mean_rate += 1;
    }

    mean_coi = sampler.sample_mean_coi(coi_mean_shape, coi_mean_rate);
}

void Chain::update_p(int iteration) {
    for(size_t j = 0; j < genotyping_data.num_loci; j++) {
        // prop_p = sampler.sample_allele_frequencies(p[j], params.alpha);
        prop_p = sampler.sample_allele_frequencies2(p[j], p_prop_var[j]);
        // UtilFunctions::print("P", j, p_prop_var[j]);
        // UtilFunctions::print_vector(p[j]);
        double sum_can = 0;
        double sum_orig = 0;
        for(size_t i = 0; i < genotyping_data.num_samples; i++) {
            llik_new[j][i] = calc_genotype_marginal_llik(genotyping_data.observed_alleles[j][i], m[i], prop_p, eps_neg, eps_pos);
            sum_can += llik_new[j][i];
            sum_orig += llik_old[j][i];
        }

        // Accept
        if(sampler.sample_log_mh_acceptance() <= (sum_can - sum_orig)) {
            // UtilFunctions::print("Updating P", j, p_accept[j] / iteration);
            // UtilFunctions::print_vector(p[j]);
            p[j] = prop_p;
            p_accept[j] += 1;
            p_prop_var[j] = exp(log(p_prop_var[j]) + (1 - 0.23) / sqrt(iteration));
            for(size_t i = 0; i < genotyping_data.num_samples; i++) {
                llik_old[j][i] = llik_new[j][i];
            }
        } else {
            p_prop_var[j] = exp(log(p_prop_var[j]) - 0.23 / sqrt(iteration));
        }
    }

}

void Chain::update_eps(int iteration) {
    double prop_eps_pos = sampler.sample_epsilon_pos(eps_pos, eps_pos_var);
    double prop_eps_neg = sampler.sample_epsilon_neg(eps_neg, eps_neg_var);
    // UtilFunctions::print("Epsilon Pos:", eps_pos, prop_eps_pos, eps_pos_var);
    // UtilFunctions::print("Epsilon Neg:", eps_neg, prop_eps_neg, eps_neg_var);

    if (prop_eps_pos < params.max_eps_pos && prop_eps_pos > 0 && prop_eps_neg < params.max_eps_neg && prop_eps_neg > 0) {
        double sum_can = 0;
        double sum_orig = 0;
        for(size_t j = 0; j < genotyping_data.num_loci; j++) {
            for(size_t i = 0; i < genotyping_data.num_samples; i++) {
                llik_new[j][i] = calc_genotype_marginal_llik(genotyping_data.observed_alleles[j][i], m[i], p[j], prop_eps_neg, prop_eps_pos);
                sum_can += llik_new[j][i];
                sum_orig += llik_old[j][i];
            }
        }

        // Incorporate prior    
        sum_can += sampler.get_epsilon_log_prior(prop_eps_neg, params.eps_neg_alpha, params.eps_neg_beta);
        sum_can += sampler.get_epsilon_log_prior(prop_eps_pos, params.eps_pos_alpha, params.eps_pos_beta);
        sum_orig += sampler.get_epsilon_log_prior(eps_neg, params.eps_neg_alpha, params.eps_neg_beta);
        sum_orig += sampler.get_epsilon_log_prior(eps_pos, params.eps_pos_alpha, params.eps_pos_beta);

        // Accept
        if(sampler.sample_log_mh_acceptance() <= (sum_can - sum_orig)) {
            // UtilFunctions::print("Updating Eps", prop_eps_pos, prop_eps_neg);
            eps_pos = prop_eps_pos;
            // eps_pos_var += (1-0.23)/sqrt(double(iteration));
            eps_pos_accept += 1;

            eps_neg = prop_eps_neg;
            // eps_neg_var += (1-0.23)/sqrt(double(iteration));
            eps_neg_accept += 1;

            for(size_t j = 0; j < genotyping_data.num_loci; j++) {
                for(size_t i = 0; i < genotyping_data.num_samples; i++) {
                    llik_old[j][i] = llik_new[j][i];
                }
            }
        } else {
            // eps_pos_var -= 0.23/sqrt(double(iteration));
            // if (eps_pos_var < UNDERFLO) {
            //     eps_pos_var = UNDERFLO;
            // }

            // eps_neg_var -= 0.23/sqrt(double(iteration));
            // if (eps_neg_var < UNDERFLO) {
            //     eps_neg_var = UNDERFLO;
            // }
        }
    }
}

void Chain::update_eps_pos(int iteration) {
    double prop_eps_pos = sampler.sample_epsilon_pos(eps_pos, eps_pos_var);
    UtilFunctions::print("Epsilon Pos:", eps_pos, prop_eps_pos, eps_pos_var);

    if (prop_eps_pos < params.max_eps_pos && prop_eps_pos > 0) {
        double sum_can = 0;
        double sum_orig = 0;
        for(size_t j = 0; j < genotyping_data.num_loci; j++) {
            for(size_t i = 0; i < genotyping_data.num_samples; i++) {
                llik_new[j][i] = calc_genotype_marginal_llik(genotyping_data.observed_alleles[j][i], m[i], p[j], eps_neg, prop_eps_pos);
                sum_can += llik_new[j][i];
                sum_orig += llik_old[j][i];
            }
        }

        // Incorporate prior    
        sum_can += sampler.get_epsilon_log_prior(prop_eps_pos, params.eps_pos_alpha, params.eps_pos_beta);
        sum_orig += sampler.get_epsilon_log_prior(eps_pos, params.eps_pos_alpha, params.eps_pos_beta);

        // Accept
        if(sampler.sample_log_mh_acceptance() <= (sum_can - sum_orig)) {
            UtilFunctions::print("Updating Eps Pos", prop_eps_pos);
            eps_pos = prop_eps_pos;
            eps_pos_var += (1-0.23)/sqrt(double(iteration));
            eps_pos_accept += 1;
            for(size_t j = 0; j < genotyping_data.num_loci; j++) {
                for(size_t i = 0; i < genotyping_data.num_samples; i++) {
                    llik_old[j][i] = llik_new[j][i];
                }
            }
        } else {
            eps_pos_var -= 0.23/sqrt(double(iteration));
            if (eps_pos_var < UNDERFLO) {
                eps_pos_var = UNDERFLO;
            }
        }
    }
    
    

}

void Chain::update_eps_neg(int iteration) {
    double prop_eps_neg = sampler.sample_epsilon_neg(eps_neg, eps_neg_var);
    if (prop_eps_neg < params.max_eps_neg && prop_eps_neg > 0) {
        double sum_can = 0;
        double sum_orig = 0;
        for(size_t j = 0; j < genotyping_data.num_loci; j++) {
            for(size_t i = 0; i < genotyping_data.num_samples; i++) {
                llik_new[j][i] = calc_genotype_marginal_llik(genotyping_data.observed_alleles[j][i], m[i], p[j], prop_eps_neg, eps_pos);
                sum_can += llik_new[j][i];
                sum_orig += llik_old[j][i];
            }

        }

        // // Incorporate prior    
        sum_can += sampler.get_epsilon_log_prior(prop_eps_neg, params.eps_neg_alpha, params.eps_neg_beta);
        sum_orig += sampler.get_epsilon_log_prior(eps_neg, params.eps_neg_alpha, params.eps_neg_beta);

        // Accept
        if(sampler.sample_log_mh_acceptance() <= (sum_can - sum_orig)) {
            UtilFunctions::print("Updating Eps Neg", prop_eps_neg);
            eps_neg = prop_eps_neg;
            eps_neg_var += (1-0.23)/sqrt(double(iteration));
            eps_neg_accept += 1;
            for(size_t j = 0; j < genotyping_data.num_loci; j++) {
                for(size_t i = 0; i < genotyping_data.num_samples; i++) {
                    llik_old[j][i] = llik_new[j][i];
                }
            }
        } else {
            eps_neg_var -= 0.23/sqrt(double(iteration));
            if (eps_neg_var < UNDERFLO) {
                eps_neg_var = UNDERFLO;
            }
        }
    }
    
}

std::vector<double> Chain::reweight_allele_frequencies(std::vector<double> const &allele_frequencies, std::vector<int> const &observed_genotype, double epsilon_neg, double epsilon_pos, int coi) {
    // std::vector<double> res(allele_frequencies.size(), 0);
    // double res_sum = 0;
    // int total_obs_alleles = std::accumulate(observed_genotype.begin(), observed_genotype.end(), 0);
    // for(size_t i = 0; i < allele_frequencies.size(); i++) {
    //     res[i] = allele_frequencies[i] * ((observed_genotype[i] * (1 - epsilon_neg)) + ((1 - observed_genotype[i]) * epsilon_neg)) + (epsilon_pos * coi) / allele_frequencies.size() ;
    //     res_sum += res[i];
    // }

    // double res_sum_inv = 1.0 / res_sum;
    // // UtilFunctions::print("Original Allele Frequencies");
    // // UtilFunctions::print_vector(allele_frequencies);
    // // UtilFunctions::print("Res Sum Inv:", res_sum_inv, res_sum);
    // for(size_t i = 0; i < res.size(); i++)
    // {
    //     res[i] *= res_sum_inv;
    // }
    // return res;

    std::vector<double> res(allele_frequencies.size(), 0);
    std::vector<double> tp(allele_frequencies.size(), 0);
    std::vector<double> fn(allele_frequencies.size(), 0);
    double tp_sum = 0;
    double fn_sum = 0;

    for(size_t i = 0; i < allele_frequencies.size(); i++)
    {
        if (observed_genotype[i]) {
            tp[i] = allele_frequencies[i];
            tp_sum += allele_frequencies[i];
        } else {
            fn[i] = allele_frequencies[i];
            fn_sum += allele_frequencies[i];
        }

    }

    double inv_tp_sum = 1.0 / tp_sum;
    double inv_fn_sum = 1.0 / fn_sum;

    for(size_t i = 0; i < allele_frequencies.size(); i++)
    {
        if (observed_genotype[i]) {
            res[i] += tp[i] * inv_tp_sum * (1 - eps_neg);
        } else {
            res[i] += fn[i] * inv_fn_sum * (eps_neg);
        }
    }
    return res;

}

std::vector<double> Chain::calc_genotype_log_pmf(std::vector<std::vector<int> > const &genotypes, int coi, std::vector<double> const &allele_frequencies, int num_genotypes) {
    std::vector<double> res(num_genotypes, lookup.lookup_lgamma[coi + 1]);
    for(size_t i = 0; i < allele_frequencies.size(); i++) {
        double log_prob = log(allele_frequencies[i] + 1e-12);
        for(size_t j = 0; j < num_genotypes; j++) {
            if(genotypes[j][i] > 0) {
                res[j] += (genotypes[j][i] * log_prob) - lookup.lookup_lgamma[genotypes[j][i] + 1];
            }
        }
    }

    return res;

};

std::vector<double> Chain::calc_obs_genotype_lliks(std::vector<int> const &obs_genotype, std::vector<std::vector<int> > const &true_genotypes, double epsilon_neg, double epsilon_pos, int num_genotypes) {
    std::vector<double> lliks(true_genotypes.size(), 0);

    double tp = log(1 - epsilon_neg);
    double tn = log(1 - epsilon_pos);
    double fp = log(epsilon_pos);
    double fn = log(epsilon_neg);

    for(size_t i = 0; i < num_genotypes; i++) // Iterate over genotypes
    {
        for(size_t j = 0; j < true_genotypes[i].size(); j++) // Iterate over alleles
        {
            if(obs_genotype[j]) {
                if(true_genotypes[i][j]) {
                    lliks[i] += true_genotypes[i][j] * tp;
                } else {
                    lliks[i] += fp;
                }
            } else {
                if(true_genotypes[i][j]) {
                    lliks[i] += true_genotypes[i][j] * fn;
                } else {
                    lliks[i] += tn;
                }
            }
        }
    }
    return lliks;
}

void Chain::generate_possible_genotypes_helper(std::vector<int> &chosen, std::vector<int> &arr, int index, int r, int n, int start, int end) {
    if(index == r) {
        std::vector<int> genotype(n);
        for(size_t i = 0; i < r; i++)
        {
            genotype[arr[chosen[i]]] += 1;
        }
        Rcpp::checkUserInterrupt();
        // UtilFunctions::print("Genotype: ");
        // UtilFunctions::print_vector(genotype);
        true_genotypes_cache[std::make_tuple(r, n)].push_back(genotype);
        return;
    }

    for(int i = start; i <= end; i++)
    {
        chosen[index] = i;
        generate_possible_genotypes_helper(chosen, arr, index + 1, r, n, i, end);
    }
    
}

void Chain::generate_possible_genotypes(int coi, int total_alleles) {
    std::vector<int> chosen(coi + 1);
    std::vector<int> arr(total_alleles);
    std::iota(arr.begin(), arr.end(), 0);
    UtilFunctions::print("Generating Possible Genotypes:", coi, total_alleles);
    UtilFunctions::print("Total Genotypes:", lookup.lookup_lgamma[coi + total_alleles], lookup.lookup_lgamma[coi + 1], lookup.lookup_lgamma[total_alleles]);
    generate_possible_genotypes_helper(chosen, arr, 0, coi, total_alleles, 0, total_alleles - 1);
}

long double Chain::calc_exact_genotype_marginal_llik(std::vector<int> const &obs_genotype, int coi, std::vector<double> const &allele_frequencies, double epsilon_neg, double epsilon_pos) {
    // UtilFunctions::print("Calculating Exact Genotype");
    // std::tuple<int, int> cache_key = std::make_tuple(coi, obs_genotype.size());
    std::get<0>(genotype_cache_key) = coi;
    std::get<1>(genotype_cache_key) = obs_genotype.size();
    
    true_genotypes_iter = true_genotypes_cache.find(genotype_cache_key);
    
    if(true_genotypes_iter == true_genotypes_cache.end()) {

        generate_possible_genotypes(coi, obs_genotype.size());
    }

    genotypes = true_genotypes_cache[genotype_cache_key];
    g_given_g_star_probabilities = calc_obs_genotype_lliks(obs_genotype, genotypes, epsilon_neg, epsilon_pos, genotypes.size());
    g_star_probabilities = calc_genotype_log_pmf(genotypes, coi, allele_frequencies, genotypes.size());

    long double res = 0;
    for(size_t i = 0; i < genotypes.size(); i++) {
        res += exp(g_given_g_star_probabilities[i] + g_star_probabilities[i]);
    }

    return log(res);
}

long double Chain::calc_estimated_genotype_marginal_llik(std::vector<int> const &obs_genotype, int coi, std::vector<double> const &allele_frequencies, double epsilon_neg, double epsilon_pos) {
    // int importance_sampling_depth = std::min(params.importance_sampling_depth, lookup.lookup_sampling_depth[coi][allele_frequencies.size()]);
    // auto importance_reweighted_frequencies = allele_frequencies;
    int importance_sampling_depth = params.importance_sampling_depth;
    importance_reweighted_frequencies = reweight_allele_frequencies(allele_frequencies, obs_genotype, epsilon_neg, epsilon_pos, coi);

    sample_true_genotypes = sampler.sample_genotype(coi, importance_reweighted_frequencies, importance_sampling_depth);
    importance_probabilities = calc_genotype_log_pmf(sample_true_genotypes, coi, importance_reweighted_frequencies, importance_sampling_depth);
    g_star_probabilities = calc_genotype_log_pmf(sample_true_genotypes, coi, allele_frequencies, importance_sampling_depth);
    g_given_g_star_probabilities = calc_obs_genotype_lliks(obs_genotype, sample_true_genotypes, epsilon_neg, epsilon_pos, importance_sampling_depth);

    long double res = 0;
    for(size_t i = 0; i < importance_sampling_depth; i++) {
        res += exp(g_given_g_star_probabilities[i] + g_star_probabilities[i] - importance_probabilities[i]);
        // UtilFunctions::print("Sample I:", g_given_g_star_probabilities[i], g_star_probabilities[i], importance_probabilities[i]);
        // UtilFunctions::print("Curr Res:", res);
    }

    // UtilFunctions::print("Pre-normed Marginal llik:", UtilFunctions::fastlog(res));
    res = res / importance_sampling_depth;

    // double variance = 0;
    // for(size_t i = 0; i < importance_sampling_depth; i++) {
    //     variance += pow(exp(g_given_g_star_probabilities[i] + g_star_probabilities[i] - importance_probabilities[i]) - res, 2);
    // }
    // variance = variance / importance_sampling_depth;
    // UtilFunctions::print("Importance Samping Depth Variance:", log(variance));
    // UtilFunctions::print("Normed Marginal llik:", log(res));
    // UtilFunctions::print("---------------------------------------------------");
    return log(res);
}

long double Chain::calc_genotype_marginal_llik(std::vector<int> const &obs_genotype, int coi, std::vector<double> const &allele_frequencies, double epsilon_neg, double epsilon_pos) {
    double log_total_combinations = lookup.lookup_lgamma[coi + allele_frequencies.size()] - lookup.lookup_lgamma[coi + 1] - lookup.lookup_lgamma[allele_frequencies.size()];
    double log_importance_sampling = log(params.importance_sampling_depth * 10);
    // UtilFunctions::print("Exact or Estimate:", log_total_combinations, log_importance_sampling);
    if(log_total_combinations <= log_importance_sampling) { // Fix this condition
        return calc_exact_genotype_marginal_llik(obs_genotype, coi, allele_frequencies, epsilon_neg, epsilon_pos);
    } else {
        return calc_estimated_genotype_marginal_llik(obs_genotype, coi, allele_frequencies, epsilon_neg, epsilon_pos);
    }
}

void Chain::initialize_likelihood() {
    // UtilFunctions::print("Initializing Likelihood");
    for(size_t j = 0; j < genotyping_data.num_loci; j++) {
        llik_old.push_back(std::vector<double>(genotyping_data.num_samples));
        llik_new.push_back(std::vector<double>(genotyping_data.num_samples));
        for(size_t i = 0; i < genotyping_data.num_samples; i++) {
            auto observed_alleles = genotyping_data.observed_alleles[j][i];
            double marginal_llik = calc_genotype_marginal_llik(observed_alleles, m[i], p[j], eps_neg, eps_pos);
            llik_old[j][i] = marginal_llik;
            llik_new[j][i] = marginal_llik;
            // UtilFunctions::print("LogLik:", marginal_llik, i, j, llik_old[j][i]);
        }
    }
    // UtilFunctions::print("Initializing Likelihood Finished");
};

void Chain::calculate_llik() {
    llik = 0;
    for(size_t j = 0; j < genotyping_data.num_loci; j++) {
        for(size_t i = 0; i < genotyping_data.num_samples; i++) {
            llik += llik_old[j][i];
        }
    }
}

double Chain::get_llik() {
    calculate_llik();
    return llik;
}

Chain::Chain(GenotypingData genotyping_data, Lookup lookup, Parameters params) :
    genotyping_data(genotyping_data),
    lookup(lookup),
    params(params),
    sampler(params.importance_sampling_depth, genotyping_data.num_alleles)
    {
        UtilFunctions::print("Starting Sampler...");
        genotype_cache_key = std::make_tuple(0, 0);

        llik = 0;
        eps_pos_var = params.eps_pos_var;
        eps_neg_var = params.eps_neg_var;
        m_prop_mean = std::vector<double>(genotyping_data.num_samples, 1);
        p_prop_var = std::vector<double>(genotyping_data.num_loci, 1);
        
        // initialize_sampler();
        initialize_p();
        initialize_m();
        initialize_eps_neg();
        initialize_eps_pos();

        UtilFunctions::print("M List Length:", m.size());
        UtilFunctions::print("P List Length:", p.size());

        initialize_likelihood();

    };
