
#include "chain.h"

#include "mcmc_utils.h"
#include "sampler.h"

#include <algorithm>

// Initialize P with empirical allele frequencies
void Chain::initialize_p()
{
    p_accept.resize(genotyping_data.num_loci, 0);

    std::vector<std::vector<int>> total_locus_alleles(genotyping_data.num_loci);
    std::vector<int> total_alleles(genotyping_data.num_loci);

    for (size_t i = 0; i < genotyping_data.num_loci; i++)
    {
        total_locus_alleles.push_back(
            std::vector<int>(genotyping_data.num_alleles[i]));
        total_alleles.push_back(0);
        p.push_back(std::vector<double>(genotyping_data.num_alleles[i]));

        auto obs_locus_genotypes = genotyping_data.observed_alleles[i];

        for (size_t j = 0; j < obs_locus_genotypes.size(); j++)
        {
            auto sample_genotype = obs_locus_genotypes[j];
            for (size_t k = 0; k < sample_genotype.size(); k++)
            {
                if (j == 0)
                {
                    total_locus_alleles[i].push_back(0);
                }

                total_locus_alleles[i][k] += sample_genotype[k];
                total_alleles[i] += sample_genotype[k];
            }
        }
        for (size_t j = 0; j < genotyping_data.num_alleles[i]; j++)
        {
            p[i][j] = (total_locus_alleles[i][j] + 1) /
                      ((double)total_alleles[i] +
                       genotyping_data.num_alleles[i]);  // Make sure at least 1
                                                         // allele everywhere
        }
    }
};

void Chain::initialize_m()
{
    m = genotyping_data.observed_coi;
    m_accept.resize(genotyping_data.num_samples, 0);
    individual_accept.resize(genotyping_data.num_samples, 0);
}

void Chain::initialize_eps_neg()
{
    // eps_neg = params.eps_neg_0;
    // eps_neg_accept = 0;
    eps_neg.resize(genotyping_data.num_samples, params.eps_neg_0);
    eps_neg_accept.resize(genotyping_data.num_samples, 0);
}

void Chain::initialize_eps_pos()
{
    // eps_pos = params.eps_pos_0;
    // eps_pos_accept = 0;
    eps_pos.resize(genotyping_data.num_samples, params.eps_pos_0);
    eps_pos_accept.resize(genotyping_data.num_samples, 0);
}

void Chain::initialize_mean_coi()
{
    mean_coi = 0;
    for (size_t i = 0; i < m.size(); i++)
    {
        mean_coi += m[i];
    };
    mean_coi = mean_coi / m.size();
}

// void Chain::initialize_sampler() {
//     sampler = Sampler(params.importance_sampling_depth,
//     genotyping_data.num_alleles);
// }

void Chain::update_m(int iteration)
{
    for (int i = 0; i < genotyping_data.num_samples; i++)
    {
        int prop_m = m[i] + sampler.sample_coi_delta(m_prop_mean[i]);
        // int prop_m = m[i] + sampler.sample_coi_delta(2);
        // UtilFunctions::print("Update M: ", m[i], prop_m);
        // UtilFunctions::print("M mean: ", m_prop_mean[i]);
        // Accept automatically if COI is unchanged
        if (prop_m == m[i])
        {
            m_prop_mean[i] += (1 - 0.23) / sqrt(double(iteration));
            m_accept[i] += 1;
            continue;
        }

        if (params.max_coi >= prop_m && prop_m > 0)
        {
            double sum_can = 0;
            double sum_orig = 0;
            for (size_t j = 0; j < genotyping_data.num_loci; j++)
            {
                llik_new[j][i] = calc_genotype_marginal_llik(
                    genotyping_data.observed_alleles[j][i], prop_m, p[j],
                    eps_neg[i], eps_pos[i]);
                sum_can += llik_new[j][i];
                sum_orig += llik_old[j][i];
            }

            // Poisson prior on COI
            // sum_can += sampler.get_coi_log_prior(prop_m, mean_coi);
            // sum_orig += sampler.get_coi_log_prior(m[i], mean_coi);

            // if (m[i] < prop_m)
            // {
            //     UtilFunctions::print("Going from ", m[i], "to", prop_m, "with
            //     probability", exp(sum_can - sum_orig));
            // }

            // Accept
            if (sampler.sample_log_mh_acceptance() <= (sum_can - sum_orig))
            {
                m[i] = prop_m;
                m_prop_mean[i] += (1 - 0.23) / sqrt(double(iteration));
                m_accept[i] += 1;
                for (size_t j = 0; j < genotyping_data.num_loci; j++)
                {
                    llik_old[j][i] = llik_new[j][i];
                }
            }
            else
            {
                m_prop_mean[i] -= 0.23 / sqrt(double(iteration));
                m_prop_mean[i] = (m_prop_mean[i] < 0) ? 0 : m_prop_mean[i];
            }
        }
    }
}

// Gibbs update
void Chain::update_mean_coi(int iteration)
{
    double coi_mean_shape = .25;
    double coi_mean_rate = .25;

    for (size_t i = 0; i < m.size(); i++)
    {
        coi_mean_shape += m[i] - 1;
        coi_mean_rate += 1;
    }

    mean_coi = sampler.sample_mean_coi(coi_mean_shape, coi_mean_rate);
}

void Chain::update_p(int iteration)
{
    for (size_t j = 0; j < genotyping_data.num_loci; j++)
    {
        // prop_p = sampler.sample_allele_frequencies2(p[j], p_prop_var[j]);
        prop_p = sampler.sample_allele_frequencies2(p[j], .1);
        double sum_can = 0;
        double sum_orig = 0;
        for (size_t i = 0; i < genotyping_data.num_samples; i++)
        {
            llik_new[j][i] = calc_genotype_marginal_llik(
                genotyping_data.observed_alleles[j][i], m[i], prop_p,
                eps_neg[i], eps_pos[i]);
            sum_can += llik_new[j][i];
            sum_orig += llik_old[j][i];
        }

        // Accept
        if (sampler.sample_log_mh_acceptance() <= (sum_can - sum_orig))
        {
            p[j] = prop_p;
            p_accept[j] += 1;
            // p_prop_var[j] = exp(log(p_prop_var[j]) + (1 - 0.23) /
            // sqrt(iteration));
            for (size_t i = 0; i < genotyping_data.num_samples; i++)
            {
                llik_old[j][i] = llik_new[j][i];
            }
            // Reject
        }
        // else
        // {
        // p_prop_var[j] = exp(log(p_prop_var[j]) - 0.23 / sqrt(iteration));
        // }
    }
}

// unused at the moment, updating eps_pos/eps_neg independently
void Chain::update_eps(int iteration)
{
    for (size_t i = 0; i < m.size(); i++)
    {
        double prop_eps_pos =
            sampler.sample_epsilon_pos(eps_pos[i], eps_pos_var);
        double prop_eps_neg =
            sampler.sample_epsilon_neg(eps_neg[i], eps_neg_var);

        if (prop_eps_pos < params.max_eps_pos && prop_eps_pos > 0 &&
            prop_eps_neg < params.max_eps_neg && prop_eps_neg > 0)
        {
            double sum_can = 0;
            double sum_orig = 0;
            for (size_t j = 0; j < genotyping_data.num_loci; j++)
            {
                llik_new[j][i] = calc_genotype_marginal_llik(
                    genotyping_data.observed_alleles[j][i], m[i], p[j],
                    prop_eps_neg, prop_eps_pos);
                sum_can += llik_new[j][i];
                sum_orig += llik_old[j][i];
            }

            // Incorporate prior
            sum_can += sampler.get_epsilon_log_prior(
                prop_eps_neg, params.eps_neg_alpha, params.eps_neg_beta);
            sum_can += sampler.get_epsilon_log_prior(
                prop_eps_pos, params.eps_pos_alpha, params.eps_pos_beta);
            sum_orig += sampler.get_epsilon_log_prior(
                eps_neg[i], params.eps_neg_alpha, params.eps_neg_beta);
            sum_orig += sampler.get_epsilon_log_prior(
                eps_pos[i], params.eps_pos_alpha, params.eps_pos_beta);

            // Accept
            if (sampler.sample_log_mh_acceptance() <= (sum_can - sum_orig))
            {
                eps_pos[i] = prop_eps_pos;
                // eps_pos_var += (1-0.23)/sqrt(double(iteration));
                eps_pos_accept[i] += 1;

                eps_neg[i] = prop_eps_neg;
                // eps_neg_var += (1-0.23)/sqrt(double(iteration));
                eps_neg_accept[i] += 1;

                for (size_t j = 0; j < genotyping_data.num_loci; j++)
                {
                    llik_old[j][i] = llik_new[j][i];
                }
            }
        }
    }
}

void Chain::update_eps_pos(int iteration)
{
    for (size_t i = 0; i < m.size(); i++)
    {
        double prop_eps_pos =
            sampler.sample_epsilon_pos(eps_pos[i], eps_pos_var);
        // UtilFunctions::print("Eps Pos prop: ", prop_eps_pos, prop_eps_pos -
        // eps_pos[i]);

        if (prop_eps_pos < params.max_eps_pos && prop_eps_pos > 0)
        {
            double sum_can = 0;
            double sum_orig = 0;
            for (size_t j = 0; j < genotyping_data.num_loci; j++)
            {
                llik_new[j][i] = calc_genotype_marginal_llik(
                    genotyping_data.observed_alleles[j][i], m[i], p[j],
                    eps_neg[i], prop_eps_pos);
                sum_can += llik_new[j][i];
                sum_orig += llik_old[j][i];
            }

            // Incorporate prior
            sum_can += sampler.get_epsilon_log_prior(
                prop_eps_pos, params.eps_pos_alpha, params.eps_pos_beta);
            sum_orig += sampler.get_epsilon_log_prior(
                eps_pos[i], params.eps_pos_alpha, params.eps_pos_beta);

            // UtilFunctions::print("Eps pos:", prop_eps_pos, exp(sum_can -
            // sum_orig));

            // Accept
            if (sampler.sample_log_mh_acceptance() <= (sum_can - sum_orig))
            {
                eps_pos[i] = prop_eps_pos;
                // eps_pos_var += (1-0.23) / sqrt(double(iteration));
                eps_pos_accept[i] += 1;
                for (size_t j = 0; j < genotyping_data.num_loci; j++)
                {
                    llik_old[j][i] = llik_new[j][i];
                }
            }
            else
            {
                // eps_pos_var -= 0.23/sqrt(double(iteration));
                // if (eps_pos_var < .005) { // Minimum variance, should make
                // input parameter instead of hardcoding eps_pos_var = .005;
                // }
            }
        }
    }
}

void Chain::update_eps_neg(int iteration)
{
    for (size_t i = 0; i < m.size(); i++)
    {
        double prop_eps_neg =
            sampler.sample_epsilon_neg(eps_neg[i], eps_neg_var);
        // UtilFunctions::print("Eps Neg prop: ", prop_eps_neg, prop_eps_neg -
        // eps_neg[i]);

        if (prop_eps_neg < params.max_eps_neg && prop_eps_neg > 0)
        {
            double sum_can = 0;
            double sum_orig = 0;
            for (size_t j = 0; j < genotyping_data.num_loci; j++)
            {
                llik_new[j][i] = calc_genotype_marginal_llik(
                    genotyping_data.observed_alleles[j][i], m[i], p[j],
                    prop_eps_neg, eps_pos[i]);
                sum_can += llik_new[j][i];
                sum_orig += llik_old[j][i];
            }

            // // Incorporate prior
            sum_can += sampler.get_epsilon_log_prior(
                prop_eps_neg, params.eps_neg_alpha, params.eps_neg_beta);
            sum_orig += sampler.get_epsilon_log_prior(
                eps_neg[i], params.eps_neg_alpha, params.eps_neg_beta);

            // UtilFunctions::print("Eps neg:", prop_eps_neg, exp(sum_can -
            // sum_orig));

            // Accept
            if (sampler.sample_log_mh_acceptance() <= (sum_can - sum_orig))
            {
                eps_neg[i] = prop_eps_neg;
                // eps_neg_var += (1-0.23) / sqrt(double(iteration));
                eps_neg_accept[i] += 1;
                for (size_t j = 0; j < genotyping_data.num_loci; j++)
                {
                    llik_old[j][i] = llik_new[j][i];
                }
            }
            else
            {
                // eps_neg_var -= 0.23/sqrt(double(iteration));
                // if (eps_neg_var < .005) { // Minimum variance, should make
                // input parameter instead of hardcoding
                //     eps_neg_var = .005;
                // }
            }
        }
    }
}

void Chain::update_individual_parameters(int iteration)
{
    for (size_t i = 0; i < m.size(); i++)
    {
        int prop_m = m[i] + sampler.sample_coi_delta(2);
        double prop_eps_neg =
            sampler.sample_epsilon_neg(eps_neg[i], eps_neg_var);
        double prop_eps_pos =
            sampler.sample_epsilon_neg(eps_pos[i], eps_pos_var);

        if (prop_eps_neg < params.max_eps_neg && prop_eps_neg > 0 &&
            prop_eps_pos < params.max_eps_pos && prop_eps_pos > 0 &&
            params.max_coi >= prop_m && prop_m > 0)
        {
            double sum_can = 0;
            double sum_orig = 0;
            for (size_t j = 0; j < genotyping_data.num_loci; j++)
            {
                llik_new[j][i] = calc_genotype_marginal_llik(
                    genotyping_data.observed_alleles[j][i], prop_m, p[j],
                    prop_eps_neg, prop_eps_pos);
                sum_can += llik_new[j][i];
                sum_orig += llik_old[j][i];
            }

            // Incorporate priors
            sum_can += sampler.get_epsilon_log_prior(
                prop_eps_neg, params.eps_neg_alpha, params.eps_neg_beta);
            sum_orig += sampler.get_epsilon_log_prior(
                eps_neg[i], params.eps_neg_alpha, params.eps_neg_beta);
            sum_can += sampler.get_epsilon_log_prior(
                prop_eps_pos, params.eps_pos_alpha, params.eps_pos_beta);
            sum_orig += sampler.get_epsilon_log_prior(
                eps_pos[i], params.eps_pos_alpha, params.eps_pos_beta);
            // sum_can += sampler.get_coi_log_prior(prop_m, mean_coi);
            // sum_orig += sampler.get_coi_log_prior(m[i], mean_coi);

            // Accept
            if (sampler.sample_log_mh_acceptance() <= (sum_can - sum_orig))
            {
                m[i] = prop_m;
                eps_neg[i] = prop_eps_neg;
                eps_pos[i] = prop_eps_pos;
                individual_accept[i] += 1;
                // m_prop_mean[i] += (1 - 0.23) / sqrt(double(iteration));
                // m_accept[i] += 1;
                for (size_t j = 0; j < genotyping_data.num_loci; j++)
                {
                    llik_old[j][i] = llik_new[j][i];
                }
            }
        }
    }
}

std::vector<double> Chain::reweight_allele_frequencies(
    std::vector<double> const &allele_frequencies,
    std::vector<int> const &observed_genotype, double epsilon_neg,
    double epsilon_pos, int coi)
{
    std::vector<double> res(allele_frequencies.size(), 0);
    std::vector<double> tp(allele_frequencies.size(), 0);
    std::vector<double> fn(allele_frequencies.size(), 0);
    double tp_sum = 0;
    double fn_sum = 0;

    for (size_t i = 0; i < allele_frequencies.size(); i++)
    {
        if (observed_genotype[i])
        {
            tp[i] = allele_frequencies[i] * coi;
            tp_sum += allele_frequencies[i] * coi;
        }
        else
        {
            fn[i] = allele_frequencies[i];
            fn_sum += allele_frequencies[i];
        }
    }

    double inv_tp_sum = 1.0 / tp_sum;
    double inv_fn_sum = 1.0 / fn_sum;
    double total_ = 0;

    for (size_t i = 0; i < allele_frequencies.size(); i++)
    {
        if (observed_genotype[i])
        {
            res[i] += tp[i] * inv_tp_sum * (1 - epsilon_neg);
        }
        else
        {
            res[i] += fn[i] * inv_fn_sum * (epsilon_neg);
        }
        total_ += res[i];
    }

    for (size_t i = 0; i < res.size(); i++)
    {
        res[i] = res[i] / total_;
    }

    return res;
}

std::vector<double> Chain::calc_genotype_log_pmf(
    std::vector<std::vector<int>> const &genotypes, int coi,
    std::vector<double> const &allele_frequencies, int num_genotypes)
{
    std::vector<double> res(num_genotypes, lookup.lookup_lgamma[coi + 1]);
    for (size_t i = 0; i < allele_frequencies.size(); i++)
    {
        double log_prob = log(allele_frequencies[i] + 1e-12);
        for (size_t j = 0; j < num_genotypes; j++)
        {
            if (genotypes[j][i] > 0)
            {
                res[j] += (genotypes[j][i] * log_prob) -
                          lookup.lookup_lgamma[genotypes[j][i] + 1];
            }
        }
    }

    return res;
};

// std::vector<double>
// Chain::calc_genotype_log_pmf(std::vector<std::vector<int>> const &genotypes,
// int coi, std::vector<double> const &allele_frequencies, int num_genotypes)
// {
//     std::vector<double> res(num_genotypes);
//     for (size_t i = 0; i < num_genotypes; i++)
//     {
//         res[i] = gsl_ran_multinomial_lnpdf(allele_frequencies.size(),
//         allele_frequencies.data(), (unsigned int *)genotypes[i].data());
//     }
//     return res;
// }

std::vector<double> Chain::calc_obs_genotype_lliks(
    std::vector<int> const &obs_genotype,
    std::vector<std::vector<int>> const &true_genotypes, double epsilon_neg,
    double epsilon_pos, int num_genotypes)
{
    std::vector<double> lliks(true_genotypes.size(), 0);

    double tp = log(1 - epsilon_neg);
    double tn = log(1 - epsilon_pos);
    double fp = log(epsilon_pos);
    double fn = log(epsilon_neg);

    for (size_t i = 0; i < num_genotypes; i++)  // Iterate over genotypes
    {
        for (size_t j = 0; j < true_genotypes[i].size();
             j++)  // Iterate over alleles
        {
            if (obs_genotype[j])
            {
                if (true_genotypes[i][j])
                {
                    lliks[i] += true_genotypes[i][j] * tp;
                }
                else
                {
                    lliks[i] += fp;
                }
            }
            else
            {
                if (true_genotypes[i][j])
                {
                    lliks[i] += true_genotypes[i][j] * fn;
                }
                else
                {
                    lliks[i] += tn;
                }
            }
        }
    }
    return lliks;
}

void Chain::generate_possible_genotypes_helper(std::vector<int> &chosen,
                                               std::vector<int> &arr, int index,
                                               int r, int n, int start, int end)
{
    if (index == r)
    {
        std::vector<int> genotype(n);
        for (size_t i = 0; i < r; i++)
        {
            genotype[arr[chosen[i]]] += 1;
        }
        Rcpp::checkUserInterrupt();
        true_genotypes_cache[std::make_tuple(r, n)].push_back(genotype);
        return;
    }

    for (int i = start; i <= end; i++)
    {
        chosen[index] = i;
        generate_possible_genotypes_helper(chosen, arr, index + 1, r, n, i,
                                           end);
    }
}

void Chain::generate_possible_genotypes(int coi, int total_alleles)
{
    std::vector<int> chosen(coi + 1);
    std::vector<int> arr(total_alleles);
    std::iota(arr.begin(), arr.end(), 0);
    generate_possible_genotypes_helper(chosen, arr, 0, coi, total_alleles, 0,
                                       total_alleles - 1);
}

long double Chain::calc_exact_genotype_marginal_llik(
    std::vector<int> const &obs_genotype, int coi,
    std::vector<double> const &allele_frequencies, double epsilon_neg,
    double epsilon_pos)
{
    // UtilFunctions::print("----------------");
    // UtilFunctions::print("Calculating Exact LLik");
    // UtilFunctions::print("----------------");

    std::get<0>(genotype_cache_key) = coi;
    std::get<1>(genotype_cache_key) = obs_genotype.size();

    const auto &true_genotypes_iter =
        true_genotypes_cache.find(genotype_cache_key);

    if (true_genotypes_iter == true_genotypes_cache.end())
    {
        generate_possible_genotypes(coi, obs_genotype.size());
    }

    const auto &genotypes = true_genotypes_cache[genotype_cache_key];
    // UtilFunctions::print("Cached Genotypes --", coi);
    // for (size_t i = 0; i < genotypes.size(); i++)
    // {
    //     UtilFunctions::print_vector(genotypes[i]);
    // }
    // UtilFunctions::print("-------------");

    const auto &g_given_g_star_probabilities = calc_obs_genotype_lliks(
        obs_genotype, genotypes, epsilon_neg, epsilon_pos, genotypes.size());
    const auto &g_star_probabilities = calc_genotype_log_pmf(
        genotypes, coi, allele_frequencies, genotypes.size());

    long double res = 0;
    for (size_t i = 0; i < genotypes.size(); i++)
    {
        res += exp(g_given_g_star_probabilities[i] + g_star_probabilities[i]);
    }
    // UtilFunctions::print("Exact GT Llik: ", res);
    return log(res);
}

long double Chain::calc_estimated_genotype_marginal_llik(
    std::vector<int> const &obs_genotype, int coi,
    std::vector<double> const &allele_frequencies, double epsilon_neg,
    double epsilon_pos)
{
    int importance_sampling_depth = params.importance_sampling_depth * coi;
    const auto &importance_reweighted_frequencies = reweight_allele_frequencies(
        allele_frequencies, obs_genotype, epsilon_neg, epsilon_pos, coi);

    const auto &sample_true_genotypes = sampler.sample_genotype(
        coi, importance_reweighted_frequencies, importance_sampling_depth);
    const auto &importance_probabilities = calc_genotype_log_pmf(
        sample_true_genotypes, coi, importance_reweighted_frequencies,
        importance_sampling_depth);
    const auto &g_star_probabilities =
        calc_genotype_log_pmf(sample_true_genotypes, coi, allele_frequencies,
                              importance_sampling_depth);
    const auto &g_given_g_star_probabilities = calc_obs_genotype_lliks(
        obs_genotype, sample_true_genotypes, epsilon_neg, epsilon_pos,
        importance_sampling_depth);

    long double res = 0;
    // long double importance_prob = 0;
    // long double target_prob = 0;
    // long double divergence = 0;
    for (size_t i = 0; i < importance_sampling_depth; i++)
    {
        // divergence += exp(g_star_probabilities[i] -
        // importance_probabilities[i]); importance_prob +=
        // exp(importance_probabilities[i]); target_prob +=
        // exp(g_given_g_star_probabilities[i] + g_star_probabilities[i]);
        // UtilFunctions::print("importance reweighted: ",
        // exp(g_given_g_star_probabilities[i] + g_star_probabilities[i] -
        // importance_probabilities[i]));
        res += exp(g_given_g_star_probabilities[i] + g_star_probabilities[i] -
                   importance_probabilities[i]);
    }

    // UtilFunctions::print("---------------");
    // UtilFunctions::print("Eps: ", epsilon_neg, epsilon_pos);
    // UtilFunctions::print("Importance Depth:", importance_prob);
    // UtilFunctions::print("True Depth:", target_prob);
    // UtilFunctions::print("Importance Prob:", importance_prob);
    // UtilFunctions::print("Divergence:", divergence /
    // importance_sampling_depth); UtilFunctions::print_vector(obs_genotype);
    // UtilFunctions::print("COI:", coi);

    // UtilFunctions::print("Reweighted Freqs: ");
    // UtilFunctions::print_vector(importance_reweighted_frequencies);
    // UtilFunctions::print("---------------");

    res = res / importance_sampling_depth;
    // UtilFunctions::print("Approx GT Llik: ", res);
    return log(res);
}

long double Chain::calc_genotype_marginal_llik(
    std::vector<int> const &obs_genotype, int coi,
    std::vector<double> const &allele_frequencies, double epsilon_neg,
    double epsilon_pos)
{
    // Implement a way to deal with missing data here --> could either
    // completely ignore missing data, or bias towards low COI
    double log_total_combinations =
        lookup.lookup_lgamma[coi + allele_frequencies.size()] -
        lookup.lookup_lgamma[coi + 1] -
        lookup.lookup_lgamma[allele_frequencies.size()];
    double log_importance_sampling =
        log(params.importance_sampling_depth * coi * 10);
    if (log_total_combinations <= log_importance_sampling)
    {
        return calc_exact_genotype_marginal_llik(
            obs_genotype, coi, allele_frequencies, epsilon_neg, epsilon_pos);
    }
    else
    {
        return calc_estimated_genotype_marginal_llik(
            obs_genotype, coi, allele_frequencies, epsilon_neg, epsilon_pos);
    }
}

void Chain::initialize_likelihood()
{
    for (size_t j = 0; j < genotyping_data.num_loci; j++)
    {
        llik_old.push_back(std::vector<double>(genotyping_data.num_samples));
        llik_new.push_back(std::vector<double>(genotyping_data.num_samples));
        for (size_t i = 0; i < genotyping_data.num_samples; i++)
        {
            auto observed_alleles = genotyping_data.observed_alleles[j][i];
            double marginal_llik = calc_genotype_marginal_llik(
                observed_alleles, m[i], p[j], eps_neg[i], eps_pos[i]);
            llik_old[j][i] = marginal_llik;
            llik_new[j][i] = marginal_llik;
        }
    }
};

void Chain::calculate_llik()
{
    llik = 0;
    for (size_t j = 0; j < genotyping_data.num_loci; j++)
    {
        for (size_t i = 0; i < genotyping_data.num_samples; i++)
        {
            llik += llik_old[j][i];
        }
    }
}

double Chain::get_llik()
{
    calculate_llik();
    return llik;
}

Chain::Chain(GenotypingData genotyping_data, Lookup lookup, Parameters params)
    : genotyping_data(genotyping_data),
      lookup(lookup),
      params(params),
      sampler(params.importance_sampling_depth * params.max_coi,
              genotyping_data.num_alleles)
{
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
    initialize_likelihood();
};
