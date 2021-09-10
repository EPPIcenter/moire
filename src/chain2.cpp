#include "chain.h"
#include "mcmc_utils.h"
#include "sampler.h"

#include <cmath>
#include <algorithm>

#include <boost/math/special_functions/binomial.hpp>
#include <boost/math/special_functions/math_fwd.hpp>

#include <limits>
#include <map>
#include <numeric>

// Initialize P with empirical allele frequencies
void Chain::initialize_p()
{
    std::vector<std::vector<int>> total_locus_alleles(genotyping_data.num_loci);
    std::vector<int> total_alleles(genotyping_data.num_loci);

    for (size_t i = 0; i < genotyping_data.num_loci; i++)
    {
        total_locus_alleles.push_back(
            std::vector<int>(genotyping_data.num_alleles[i]));
        total_alleles.push_back(0);
        p.push_back(std::vector<double>(genotyping_data.num_alleles[i]));
        for (size_t j = 0; j < genotyping_data.num_samples; j++)
        {
            const auto &sample_genotype =
                genotyping_data.get_observed_alleles(i, j);
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

    p_prop_var.resize(genotyping_data.num_loci);
    p_accept.resize(genotyping_data.num_loci);
    p_attempt.resize(genotyping_data.num_loci);

    for (int jj = 0; jj < genotyping_data.num_loci; ++jj)
    {
        int total_alleles = p[jj].size();
        p_prop_var[jj] =
            std::vector<double>(total_alleles, params.allele_freq_vars[jj]);
        p_accept[jj] = std::vector<int>(total_alleles, 0);
        p_attempt[jj] = std::vector<int>(total_alleles, 0);
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
    eps_neg.resize(genotyping_data.num_samples);

    std::fill(eps_neg.begin(), eps_neg.end(), params.eps_pos_0);
    eps_neg_accept.resize(genotyping_data.num_samples, 0);
}

void Chain::initialize_eps_pos()
{
    eps_pos.resize(genotyping_data.num_samples, params.eps_pos_0);
    eps_pos_accept.resize(genotyping_data.num_samples, 0);
}

void Chain::initialize_mean_coi()
{
    mean_coi = params.mean_coi_prior_shape * params.mean_coi_prior_scale;
}

void Chain::update_mean_coi(int iteration)
{
    double prop_mean_coi =
        sampler.sample_epsilon(mean_coi, params.mean_coi_var);

    if (prop_mean_coi > 0)
    {
        double prev_mean_coi = mean_coi;
        mean_coi = prop_mean_coi;
        calculate_mean_coi_likelihood();

        double new_llik = calc_new_likelihood();

        if (std::isnan(new_llik))
        {
            UtilFunctions::print("Encountered NaN -- mean_coi", prev_mean_coi,
                                 prop_mean_coi);
        }

        double alpha = sampler.sample_log_mh_acceptance();
        if (alpha > (new_llik - llik))
        {
            mean_coi = prev_mean_coi;
            restore_mean_coi_likelihood();
        }
        else
        {
            if (new_llik < llik)
            {
                UtilFunctions::print("Accepted (mean COI):", new_llik, llik,
                                     alpha);
            }

            llik = new_llik;
            save_mean_coi_likelihood();
        }
    }
}

void Chain::update_m(int iteration)
{
    for (size_t ii = 0; ii < genotyping_data.num_samples; ii++)
    {
        int prop_m = m[ii] + sampler.sample_coi_delta(1);
        if (prop_m > 0)
        {
            double prev_m = m[ii];
            m[ii] = prop_m;

            calculate_coi_likelihood(ii);

            for (int i = 0; i < genotyping_data.num_samples; ++i)
            {
                for (int jj = 0; jj < genotyping_data.num_loci; ++jj)
                {
                    calculate_genotype_likelihood(i, jj);
                }
            }

            double new_llik = calc_new_likelihood();
            if (std::isnan(new_llik))
            {
                UtilFunctions::print("Encountered NaN -- coi", prev_m, prop_m);
            }

            if (prev_m == prop_m)
            {
                UtilFunctions::print(llik, new_llik);
            }

            double alpha = sampler.sample_log_mh_acceptance();
            if (!std::isnan(new_llik) and alpha <= (new_llik - llik))
            {
                // if (new_llik < llik)
                // {
                //     UtilFunctions::print("Accepted (COI):", new_llik, llik,
                //                          alpha);
                // }

                llik = new_llik;
                save_coi_likelihood(ii);
                for (int jj = 0; jj < genotyping_data.num_loci; ++jj)
                {
                    save_genotype_likelihood(ii, jj);
                }
            }
            else
            {
                m[ii] = prev_m;
                restore_coi_likelihood(ii);
                for (int jj = 0; jj < genotyping_data.num_loci; ++jj)
                {
                    restore_genotype_likelihood(ii, jj);
                }
            }
        }
    }
}

/*
 * SALT Sampler approach.
 * https://doi.org/10.1080/00949655.2017.1376063
 */
void Chain::update_p(int iteration)
{
    for (size_t j = 0; j < genotyping_data.num_loci; j++)
    {
        int k = p[j].size();
        int rep = k;
        while (--rep >= 0)
        {
            const int idx = sampler.sample_random_int(0, k - 1);

            auto logitPropP = UtilFunctions::logitVec(p[j]);

            double logitCurr = logitPropP[idx];
            double logitProp =
                sampler.sample_epsilon(logitCurr, p_prop_var[j][idx]);

            auto currLogPQ = UtilFunctions::log_pq(logitCurr);
            auto propLogPQ = UtilFunctions::log_pq(logitProp);

            logitPropP.erase(logitPropP.begin() + idx);

            double ls = propLogPQ.second - UtilFunctions::logitSum(logitPropP);
            logitPropP = UtilFunctions::logitScale(logitPropP, ls);
            logitPropP.insert(logitPropP.begin() + idx, logitProp);

            double logAdj = (currLogPQ.first - propLogPQ.first) +
                            (k - 1) * (currLogPQ.second - propLogPQ.second);

            auto prop_p = UtilFunctions::expitVec(logitPropP);
            // check to make sure the proposed simplex is within a bounded range
            bool sub_threshold_flag = false;
            for (const auto &el : prop_p)
            {
                if (el < 1e-3)
                {
                    sub_threshold_flag = true;
                    break;
                }
            }

            if (sub_threshold_flag)
            {
                break;
            }

            auto prev_p = p[j];

            for (int ii = 0; ii < genotyping_data.num_samples; ++ii)
            {
                calculate_genotype_likelihood(ii, j);
            }
            double new_llik2 = calc_new_likelihood();

            p[j] = prop_p;

            for (int ii = 0; ii < genotyping_data.num_samples; ++ii)
            {
                calculate_genotype_likelihood(ii, j);
            }

            double new_llik = calc_new_likelihood();

            UtilFunctions::print("Prop:", llik, new_llik, new_llik2);
            UtilFunctions::print_vector(prev_p);
            UtilFunctions::print_vector(prop_p);

            if (std::isnan(new_llik))
            {
                UtilFunctions::print("Encountered NaN -- allele_freqs");
            }
            p_attempt[j][idx] += 1;

            double acceptanceRatio = new_llik - llik + logAdj;

            UtilFunctions::print("Freq:", acceptanceRatio, new_llik, llik,
                                 logAdj);
            if (std::isnan(new_llik) or
                sampler.sample_log_mh_acceptance() > acceptanceRatio)
            {
                p[j] = prev_p;
                for (int ii = 0; ii < genotyping_data.num_samples; ++ii)
                {
                    restore_genotype_likelihood(ii, j);
                }
            }
            else
            {
                llik = new_llik;
                for (int ii = 0; ii < genotyping_data.num_samples; ++ii)
                {
                    save_genotype_likelihood(ii, j);
                }
                p_accept[j][idx] += 1;
            }

            // UtilFunctions::print("Freq:", llik, p_accept[j][idx],
            //                      p_prop_var[j][idx]);

            if (params.adapt_allele_freq_vars and iteration < params.burnin and
                p_attempt[j][idx] > 15)  // don't start adapting until there are
                                         // at least a few samples
            {
                double acceptanceRate = std::min(
                    (p_accept[j][idx] + 1) / (double(p_attempt[j][idx]) + 1),
                    1.0);
                p_prop_var[j][idx] += (acceptanceRate - .43) /
                                      std::pow(p_attempt[j][idx] + 1, .5);
                // (acceptanceRate - .23) / (iteration + 1);
                p_prop_var[j][idx] = std::max(p_prop_var[j][idx], .01);
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

        if (prop_eps_pos < params.max_eps_pos && prop_eps_pos > 1e-5)
        {
            double prev_eps_pos = eps_pos[i];
            eps_pos[i] = prop_eps_pos;
            calculate_eps_pos_likelihood(i);

            for (int jj = 0; jj < genotyping_data.num_loci; ++jj)
            {
                calculate_genotype_likelihood(i, jj);
            }

            double new_llik = calc_new_likelihood();
            if (std::isnan(new_llik))
            {
                UtilFunctions::print("Encountered NaN -- eps_pos");
            }

            // Reject
            if (std::isnan(new_llik) or
                sampler.sample_log_mh_acceptance() > (new_llik - llik))
            {
                eps_pos[i] = prev_eps_pos;
                restore_eps_pos_likelihood(i);
                for (int jj = 0; jj < genotyping_data.num_loci; ++jj)
                {
                    restore_genotype_likelihood(i, jj);
                }
            }
            else
            {
                llik = new_llik;
                save_eps_pos_likelihood(i);
                for (int jj = 0; jj < genotyping_data.num_loci; ++jj)
                {
                    save_genotype_likelihood(i, jj);
                }
                eps_pos_accept[i] += 1;
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

        if (prop_eps_neg < params.max_eps_neg && prop_eps_neg > 1e-5)
        {
            double prev_eps_neg = eps_neg[i];
            eps_neg[i] = prop_eps_neg;
            calculate_eps_neg_likelihood(i);

            for (int jj = 0; jj < genotyping_data.num_loci; ++jj)
            {
                calculate_genotype_likelihood(i, jj);
            }

            double new_llik = calc_new_likelihood();
            if (std::isnan(new_llik))
            {
                UtilFunctions::print("Encountered NaN -- eps_neg");
            }

            // Reject
            if (std::isnan(new_llik) or
                sampler.sample_log_mh_acceptance() > (new_llik - llik))
            {
                eps_neg[i] = prev_eps_neg;
                restore_eps_neg_likelihood(i);
                for (int jj = 0; jj < genotyping_data.num_loci; ++jj)
                {
                    restore_genotype_likelihood(i, jj);
                }
            }
            else
            {
                llik = new_llik;
                save_eps_neg_likelihood(i);
                for (int jj = 0; jj < genotyping_data.num_loci; ++jj)
                {
                    save_genotype_likelihood(i, jj);
                }
                eps_neg_accept[i] += 1;
            }
        }
    }
}

void Chain::update_individual_parameters(int iteration)
{
    for (size_t i = 0; i < m.size(); i++)
    {
        int prop_m = m[i] + sampler.sample_coi_delta(1);
        double prop_eps_neg =
            sampler.sample_epsilon_neg(eps_neg[i], eps_neg_var);
        double prop_eps_pos =
            sampler.sample_epsilon_neg(eps_pos[i], eps_pos_var);

        if (prop_eps_neg < params.max_eps_neg and prop_eps_neg > 1e-5 and
            prop_eps_pos < params.max_eps_pos and prop_eps_pos > 1e-5 and
            (prop_m > 0))
        {
            double prev_m = m[i];
            double prev_eps_pos = eps_pos[i];
            double prev_eps_neg = eps_neg[i];

            m[i] = prop_m;
            eps_neg[i] = prop_eps_neg;
            eps_pos[i] = prop_eps_pos;

            calculate_coi_likelihood(i);
            calculate_eps_neg_likelihood(i);
            calculate_eps_pos_likelihood(i);

            for (int jj = 0; jj < genotyping_data.num_loci; ++jj)
            {
                calculate_genotype_likelihood(i, jj);
            }

            double new_llik = calc_new_likelihood();
            if (std::isnan(new_llik))
            {
                UtilFunctions::print("Encountered NaN -- indiv_params");
            }

            // Reject
            if (std::isnan(new_llik) or
                sampler.sample_log_mh_acceptance() > (new_llik - llik))
            {
                m[i] = prev_m;
                eps_neg[i] = prev_eps_neg;
                eps_pos[i] = prev_eps_pos;

                restore_coi_likelihood(i);
                restore_eps_neg_likelihood(i);
                restore_eps_pos_likelihood(i);
                for (int jj = 0; jj < genotyping_data.num_loci; ++jj)
                {
                    restore_genotype_likelihood(i, jj);
                }
            }
            else
            {
                llik = new_llik;
                save_coi_likelihood(i);
                save_eps_neg_likelihood(i);
                save_eps_pos_likelihood(i);
                for (int jj = 0; jj < genotyping_data.num_loci; ++jj)
                {
                    restore_genotype_likelihood(i, jj);
                }
                individual_accept[i] += 1;
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
            tp[i] = allele_frequencies[i];
            tp_sum += allele_frequencies[i];
        }
        else
        {
            fn[i] = allele_frequencies[i];
            fn_sum += allele_frequencies[i];
        }
    }

    double inv_tp_sum = 1.0 / tp_sum;
    double inv_fn_sum = 1.0 / fn_sum;

    // double prop_pos = ((1 - epsilon_pos) / (1 - epsilon_pos +
    // epsilon_neg)); double prop_neg = (epsilon_neg / (1 - epsilon_pos +
    // epsilon_neg));

    // double obs_pos_mass = inv_tp_sum * prop_pos;
    // double obs_neg_mass = inv_fn_sum * prop_neg;

    double obs_pos_mass = inv_tp_sum * .999;
    double obs_neg_mass = inv_fn_sum * .001;

    for (size_t i = 0; i < allele_frequencies.size(); i++)
    {
        if (observed_genotype[i])
        {
            res[i] += tp[i] * obs_pos_mass;
        }
        else
        {
            res[i] += fn[i] * obs_neg_mass;
        }
    }

    return res;
}

double Chain::calc_transmission_process(
    std::vector<int> const &allele_index_vec,
    std::vector<double> const &allele_frequencies, int coi)
{
    // transmission process - prob that after "coi" number of draws, all
    // alleles are drawn at least once conditional on all draws come
    // from the constrained set, where the constrained set is the set of
    // positive alleles in the latent genotype

    double constrained_set_total_prob = 0;
    double res = 0;
    prVec_.clear();
    prVec_.reserve(allele_index_vec.size());

    for (size_t j = 0; j < allele_index_vec.size(); j++)
    {
        prVec_.push_back(allele_frequencies[allele_index_vec[j]]);
        constrained_set_total_prob += prVec_.back();
    }

    // normalize the vector
    for (double &k : prVec_)
    {
        k = k / constrained_set_total_prob;
    }

    res = std::log(1 - probAnyMissing_(prVec_, coi)) +
          log(constrained_set_total_prob) * coi;

    return res;
}

double Chain::calc_observation_process(std::vector<int> const &allele_index_vec,
                                       std::vector<int> const &obs_genotype,
                                       int coi, double epsilon_neg,
                                       double epsilon_pos)
{
    double res = 0;
    int fp = 0;
    int tp = 0;
    int fn = 0;
    int tn = 0;

    int vec_pointer = 0;
    int next_allele_index = allele_index_vec[vec_pointer];
    int total_alleles = allele_index_vec.size();

    int j = 0;
    bool is_allele;
    for (const auto &e : obs_genotype)
    {
        is_allele = (j == next_allele_index);
        fp += e == 1 and !is_allele;
        tp += e == 1 and is_allele;
        fn += e == 0 and is_allele;
        tn += e == 0 and !is_allele;
        vec_pointer += is_allele;

        if (vec_pointer < total_alleles)
        {
            next_allele_index = allele_index_vec[vec_pointer];
        }
        else
        {
            next_allele_index = -1;
        }
        ++j;
    }

    res += std::log(epsilon_neg) * fn;
    res += std::log(1 - epsilon_neg) * tn;
    res += std::log(epsilon_pos) * fp;
    res += std::log(1 - epsilon_pos) * tp;

    return res;
};

double Chain::calc_genotype_log_pmf(
    std::vector<int> const &allele_index_vec,
    std::vector<int> const &obs_genotype, double epsilon_pos,
    double epsilon_neg, int coi, std::vector<double> const &allele_frequencies)
{
    double res = 0.0;
    res += calc_transmission_process(allele_index_vec, allele_frequencies, coi);

    res += calc_observation_process(allele_index_vec, obs_genotype, coi,
                                    epsilon_neg, epsilon_pos);

    return res;
}

long double Chain::importance_sample(
    std::vector<int> const &obs_genotype, int coi, double epsilon_neg,
    double epsilon_pos, std::vector<double> const &allele_frequencies,
    int sampling_depth)
{
    int i = sampling_depth;
    double true_prob;
    double importance_prob;
    double obs_prob;
    double accumulated_weight = 0.0;
    double accumulated_prob = 0.0;
    double est = 0.0;

    std::map<std::vector<int>, std::tuple<double, double, double>> pqf_memo{};

    std::vector<int> allele_index_vec;

    std::vector<double> importance_dist = reweight_allele_frequencies(
        allele_frequencies, obs_genotype, epsilon_neg, epsilon_pos, coi);

    while (--i >= 0)
    {
        allele_index_vec = sampler.sample_latent_genotype(coi, importance_dist);
        auto pqf_iter = pqf_memo.find(allele_index_vec);
        if (pqf_iter != pqf_memo.end())
        {
            auto pqf = pqf_iter->second;
            accumulated_prob +=
                std::get<0>(pqf) / std::get<1>(pqf) * std::get<2>(pqf);
            accumulated_weight += std::get<0>(pqf) / std::get<1>(pqf);
        }
        else
        {
            true_prob = std::exp(calc_transmission_process(
                allele_index_vec, allele_frequencies, coi));
            importance_prob = std::exp(calc_transmission_process(
                allele_index_vec, importance_dist, coi));
            obs_prob = std::exp(calc_observation_process(
                allele_index_vec, obs_genotype, coi, epsilon_neg, epsilon_pos));
            pqf_memo[allele_index_vec] = {true_prob, importance_prob, obs_prob};
            accumulated_prob += true_prob / importance_prob * obs_prob;
            accumulated_weight += true_prob / importance_prob;
        }
    }
    est = std::log(accumulated_prob / accumulated_weight);
    return est;
}

long double Chain::monte_carlo_sample(
    std::vector<int> const &obs_genotype, int coi, double epsilon_neg,
    double epsilon_pos, std::vector<double> const &true_distribution,
    int sampling_depth)
{
    int i = sampling_depth;
    double true_prob = 0.0;
    double obs_prob = 0.0;

    std::vector<int> allele_index_vec;

    double est = 0.0;
    std::map<std::vector<int>, double> pf_memo{};

    double accumulated_prob = 0.0;
    while (--i >= 0)
    {
        allele_index_vec =
            sampler.sample_latent_genotype(coi, true_distribution);
        auto pf_iter = pf_memo.find(allele_index_vec);
        if (pf_iter != pf_memo.end())
        {
            accumulated_prob += pf_iter->second;
        }
        else
        {
            true_prob = std::exp(calc_transmission_process(
                allele_index_vec, true_distribution, coi));
            obs_prob = std::exp(calc_observation_process(
                allele_index_vec, obs_genotype, coi, epsilon_pos, epsilon_neg));

            pf_memo[allele_index_vec] = true_prob * obs_prob;
            accumulated_prob += true_prob * obs_prob;
        }
    }

    est = std::log(accumulated_prob / sampling_depth);
    return est;
}

long double Chain::calc_exact_genotype_marginal_llik(
    std::vector<int> const &obs_genotype, int coi,
    std::vector<double> const &allele_frequencies, double epsilon_neg,
    double epsilon_pos)
{
    UtilFunctions::print("Calculating Exact");
    long double res = 0.0;
    long double prob = 0.0;
    for (int i = 1; i <= coi; i++)
    {
        Rcpp::checkUserInterrupt();
        allele_index_generator_.reset(allele_frequencies.size(), i);
        while (!allele_index_generator_.completed)
        {
            prob = std::exp(calc_genotype_log_pmf(
                allele_index_generator_.curr, obs_genotype, epsilon_pos,
                epsilon_neg, coi, allele_frequencies));

            res += prob;
            allele_index_generator_.next();
        }
    }
    return log(res);
}

double Chain::calc_estimated_genotype_marginal_llik(
    std::vector<int> const &obs_genotype, int coi,
    std::vector<double> const &allele_frequencies, double epsilon_neg,
    double epsilon_pos, bool verbose)
{
    std::vector<double> importance_dist = reweight_allele_frequencies(
        allele_frequencies, obs_genotype, epsilon_neg, epsilon_pos, coi);
    std::vector<int> allele_index_vec =
        sampler.sample_latent_genotype(coi, importance_dist);

    if (verbose)
    {
        UtilFunctions::print("Genotypes:");
        UtilFunctions::print_vector(obs_genotype);
        UtilFunctions::print_vector(allele_index_vec);
    }

    double true_prob =
        calc_transmission_process(allele_index_vec, allele_frequencies, coi);
    double importance_prob =
        calc_transmission_process(allele_index_vec, importance_dist, coi);
    double obs_prob = calc_observation_process(allele_index_vec, obs_genotype,
                                               coi, epsilon_neg, epsilon_pos);

    double est = true_prob + obs_prob - importance_prob;

    return est;
}

double Chain::calc_estimated_genotype_marginal_llik2(
    std::vector<int> const &obs_genotype, int coi,
    std::vector<double> const &allele_frequencies, double epsilon_neg,
    double epsilon_pos, bool verbose)
{
    // double epsilon_pos = .1;
    // double epsilon_neg = .2;

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

    // if the observed number of positives exceeds the COI, then some number of
    // them must be false positives
    int min_false_positives = std::max(0, total_obs_positives - coi);
    int max_false_positives = total_obs_positives - (total_obs_negatives == 0);

    int total_false_positives = min_false_positives;
    for (int ii = total_false_positives; ii < max_false_positives; ++ii)
    {
        total_false_positives += (sampler.sample_unif() < epsilon_pos);
    }
    int total_true_positives = total_obs_positives - total_false_positives;

    double log_prob_total_false_positives =
        std::log(boost::math::binomial_coefficient<double>(
            total_obs_positives, total_false_positives - min_false_positives)) +
        (total_false_positives - min_false_positives) * std::log(epsilon_pos) +
        total_true_positives * std::log(1 - epsilon_pos);

    // there must be at least one allele, so if all obs_positives are considered
    // false positives then there must be at least one false negative
    int min_false_negatives =
        std::min(total_obs_negatives,
                 1 * (total_false_positives == total_obs_positives));
    // also, there can't be more than min((coi - total_true_positives),
    // total_obs_negatives) false negatives
    int max_false_negatives =
        std::min(coi - total_true_positives, total_obs_negatives);

    int total_false_negatives = min_false_negatives;
    for (int ii = total_false_negatives; ii < max_false_negatives; ++ii)
    {
        total_false_negatives += (sampler.sample_unif() < epsilon_neg);
    }
    int total_true_negatives = total_obs_negatives - total_false_negatives;

    double log_prob_total_false_negatives =
        std::log(boost::math::binomial_coefficient<double>(
            total_obs_negatives, total_false_negatives - min_false_negatives)) +
        (total_false_negatives - min_false_negatives) * std::log(epsilon_neg) +
        total_true_negatives * std::log(1 - epsilon_neg);

    sampler.shuffle_vec(obs_positive_indices);
    sampler.shuffle_vec(obs_negative_indices);

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

    double true_prob =
        calc_transmission_process(allele_index_vec, allele_frequencies, coi);
    double importance_prob =
        log_prob_positive_indices + log_prob_negative_indices +
        log_prob_total_false_positives + log_prob_total_false_negatives;

    double obs_prob = calc_observation_process(allele_index_vec, obs_genotype,
                                               coi, epsilon_neg, epsilon_pos);

    double est = true_prob + obs_prob - importance_prob;

    return est;
}

long double Chain::calc_genotype_marginal_llik(
    std::vector<int> const &obs_genotype, int coi,
    std::vector<double> const &allele_frequencies, double epsilon_neg,
    double epsilon_pos)
{
    double log_total_combinations =
        lookup.get_sampling_depth(coi, allele_frequencies.size());
    double llik_calc;

    if (log_total_combinations <= std::log(params.complexity_limit))
    {
        llik_calc = calc_exact_genotype_marginal_llik(
            obs_genotype, coi, allele_frequencies, epsilon_neg, epsilon_pos);
    }
    else
    {
        UtilFunctions::print("estimating Llik");
        llik_calc = calc_estimated_genotype_marginal_llik(
            obs_genotype, coi, allele_frequencies, epsilon_neg, epsilon_pos);
        double exact_llik = calc_exact_genotype_marginal_llik(
            obs_genotype, coi, allele_frequencies, epsilon_neg, epsilon_pos);
        UtilFunctions::print("Est:", llik_calc, exact_llik);
    }

    return llik;
}

double Chain::calc_old_likelihood()
{
    double base_llik = 0;

    for (int ii = 0; ii < genotyping_data.num_samples; ++ii)
    {
        base_llik += coi_prior_old[ii];
        base_llik += eps_neg_prior_old[ii];
        base_llik += eps_pos_prior_old[ii];
    }

    std::vector<double> llik_samples(params.importance_sampling_depth,
                                     base_llik);

    for (int kk = 0; kk < params.importance_sampling_depth; ++kk)
    {
        const auto &llik_vec = genotyping_llik_old[kk];
        for (const double llik : llik_vec)
        {
            llik_samples[kk] += llik;
        }
    }

    return UtilFunctions::logSumExp(llik_samples) -
           std::log(params.importance_sampling_depth);
}

double Chain::calc_new_likelihood()
{
    double base_llik = 0;

    for (int ii = 0; ii < genotyping_data.num_samples; ++ii)
    {
        base_llik += coi_prior_new[ii];
        base_llik += eps_neg_prior_new[ii];
        base_llik += eps_pos_prior_new[ii];
    }

    std::vector<double> llik_samples(params.importance_sampling_depth,
                                     base_llik);

    for (int kk = 0; kk < params.importance_sampling_depth; ++kk)
    {
        const auto &llik_vec = genotyping_llik_new[kk];
        for (const double llik : llik_vec)
        {
            llik_samples[kk] += llik;
        }
    }

    double est_llik = UtilFunctions::logSumExp(llik_samples) -
                      std::log(params.importance_sampling_depth);
    // UtilFunctions::print("New Llik:", est_llik, llik);
    return est_llik;
}

// double Chain::calculate_llik(int num_samples)
// {
//     llik_store.resize(num_samples, 0);
//     std::fill(llik_store.begin(), llik_store.end(), 0);

//     for (int kk = 0; kk < num_samples; ++kk)
//     {
//         for (size_t j = 0; j < genotyping_data.num_loci; j++)
//         {
//             for (size_t i = 0; i < genotyping_data.num_samples; i++)
//             {
//                 double marginal_llik = 0;
//                 if (!genotyping_data.is_missing(j, i))
//                 {
//                     double log_total_combinations =
//                         lookup.get_sampling_depth(m[i], p[j].size());

//                     // cache exact genotype calculations
//                     if (log_total_combinations <=
//                         std::log(params.complexity_limit))
//                     {
//                         if (kk == 0)
//                         {
//                             for (double &d : llik_store)
//                             {
//                                 d += calc_exact_genotype_marginal_llik(
//                                     genotyping_data.get_observed_alleles(j,
//                                     i), m[i], p[j], eps_neg[i], eps_pos[i]);
//                             }
//                         }
//                     }
//                     else
//                     {
//                         llik_store[kk] +=
//                         calc_estimated_genotype_marginal_llik(
//                             genotyping_data.get_observed_alleles(j, i), m[i],
//                             p[j], eps_neg[i], eps_pos[i]);
//                     }
//                 }
//             }
//         }

//         for (size_t i = 0; i < genotyping_data.num_samples; i++)
//         {
//             auto eps_neg_prior = sampler.get_epsilon_log_prior(
//                 eps_neg[i], params.eps_neg_alpha, params.eps_neg_beta);
//             auto eps_pos_prior = sampler.get_epsilon_log_prior(
//                 eps_pos[i], params.eps_pos_alpha, params.eps_pos_beta);
//             auto coi_prior = sampler.get_coi_log_prob(m[i], mean_coi);

//             llik_store[kk] += eps_neg_prior + eps_pos_prior + coi_prior;
//         }

//         llik_store[kk] += sampler.get_coi_mean_log_prior(
//             mean_coi, params.mean_coi_prior_shape,
//             params.mean_coi_prior_scale);
//     }

//     double out = UtilFunctions::logSumExp(llik_store);
//     out = out - std::log(num_samples);
//     return out;
// }

double Chain::get_llik() { return llik; }

void Chain::calculate_genotype_likelihood(int sample_idx, int locus_idx)
{
    int sampling_depth = params.importance_sampling_depth;
    int idx = sample_idx * genotyping_data.num_loci + locus_idx;
    double log_total_combinations =
        lookup.get_sampling_depth(m[sample_idx], p[locus_idx].size());

    if (genotyping_data.is_missing(locus_idx, sample_idx))
    {
        for (int kk = 0; kk < sampling_depth; ++kk)
        {
            genotyping_llik_new[kk][idx] = 0;
        }
    }
    else if (log_total_combinations <= std::log(params.complexity_limit))
    {
        double exact_llik = calc_exact_genotype_marginal_llik(
            genotyping_data.get_observed_alleles(locus_idx, sample_idx),
            m[sample_idx], p[locus_idx], eps_neg[sample_idx],
            eps_pos[sample_idx]);
        for (int kk = 0; kk < sampling_depth; ++kk)
        {
            genotyping_llik_new[kk][idx] = exact_llik;
        }
    }
    else
    {
        // double exact_llik = calc_exact_genotype_marginal_llik(
        //     genotyping_data.get_observed_alleles(locus_idx, sample_idx),
        //     m[sample_idx], p[locus_idx], eps_neg[sample_idx],
        //     eps_pos[sample_idx]);

        for (int kk = 0; kk < sampling_depth; ++kk)
        {
            double est_llik = calc_estimated_genotype_marginal_llik2(
                genotyping_data.get_observed_alleles(locus_idx, sample_idx),
                m[sample_idx], p[locus_idx], eps_neg[sample_idx],
                eps_pos[sample_idx]);

            // UtilFunctions::print("Llik Calc:", exact_llik, est_llik);

            genotyping_llik_new[kk][idx] = est_llik;
        }
    }
}

void Chain::calculate_coi_likelihood(int sample_idx)
{
    coi_prior_new[sample_idx] =
        sampler.get_coi_log_prob(m[sample_idx], mean_coi);
}

void Chain::calculate_mean_coi_likelihood()
{
    mean_coi_prior_new = sampler.get_coi_mean_log_prior(
        mean_coi, params.mean_coi_prior_shape, params.mean_coi_prior_scale);
}

void Chain::calculate_eps_neg_likelihood(int sample_idx)
{
    eps_neg_prior_new[sample_idx] = sampler.get_epsilon_log_prior(
        eps_neg[sample_idx], params.eps_neg_alpha, params.eps_neg_beta);
}

void Chain::calculate_eps_pos_likelihood(int sample_idx)
{
    eps_pos_prior_new[sample_idx] = sampler.get_epsilon_log_prior(
        eps_pos[sample_idx], params.eps_pos_alpha, params.eps_pos_beta);
}

void Chain::initialize_likelihood()
{
    int num_samples = genotyping_data.num_samples;
    int num_loci = genotyping_data.num_loci;
    int sampling_depth = params.importance_sampling_depth;
    int row_idx, idx;

    genotyping_llik_new.resize(sampling_depth);
    genotyping_llik_old.resize(sampling_depth);

    coi_prior_old.resize(num_samples);
    coi_prior_new.resize(num_samples);
    eps_neg_prior_new.resize(num_samples);
    eps_neg_prior_old.resize(num_samples);
    eps_pos_prior_new.resize(num_samples);
    eps_pos_prior_old.resize(num_samples);
    for (int kk = 0; kk < sampling_depth; ++kk)
    {
        genotyping_llik_new[kk].resize(num_samples * num_loci);
        genotyping_llik_old[kk].resize(num_samples * num_loci);

        for (int ii = 0; ii < genotyping_data.num_samples; ++ii)
        {
            row_idx = ii * num_loci;
            for (int jj = 0; jj < genotyping_data.num_loci; ++jj)
            {
                idx = row_idx + jj;
                double log_total_combinations =
                    lookup.get_sampling_depth(m[ii], p[jj].size());

                if (genotyping_data.is_missing(jj, ii))
                {
                    genotyping_llik_new[kk][idx] = 0;
                }
                else if (log_total_combinations <=
                         std::log(params.complexity_limit))
                {
                    if (kk == 0)
                    {
                        genotyping_llik_new[kk][idx] =
                            calc_exact_genotype_marginal_llik(
                                genotyping_data.get_observed_alleles(jj, ii),
                                m[ii], p[jj], eps_neg[ii], eps_pos[ii]);
                    }
                    else
                    {
                        genotyping_llik_new[kk][idx] =
                            genotyping_llik_new[0][idx];
                    }
                }
                else
                {
                    genotyping_llik_new[kk][idx] =
                        calc_estimated_genotype_marginal_llik2(
                            genotyping_data.get_observed_alleles(jj, ii), m[ii],
                            p[jj], eps_neg[ii], eps_pos[ii]);
                }
                genotyping_llik_old[kk][idx] = genotyping_llik_new[kk][idx];
            }
        }
    }

    for (size_t ii = 0; ii < genotyping_data.num_samples; ii++)
    {
        double eps_neg_prior = sampler.get_epsilon_log_prior(
            eps_neg[ii], params.eps_neg_alpha, params.eps_neg_beta);
        double eps_pos_prior = sampler.get_epsilon_log_prior(
            eps_pos[ii], params.eps_pos_alpha, params.eps_pos_beta);
        double coi_prior = sampler.get_coi_log_prob(m[ii], mean_coi);

        eps_neg_prior_old[ii] = eps_neg_prior;
        eps_neg_prior_new[ii] = eps_neg_prior;
        eps_pos_prior_old[ii] = eps_pos_prior;
        eps_pos_prior_new[ii] = eps_pos_prior;
        coi_prior_old[ii] = coi_prior;
        coi_prior_new[ii] = coi_prior;
    }

    mean_coi_prior_new = sampler.get_coi_mean_log_prior(
        mean_coi, params.mean_coi_prior_shape, params.mean_coi_prior_scale);
    mean_coi_prior_old = mean_coi_prior_new;
    llik = calc_new_likelihood();
}

void Chain::save_genotype_likelihood(int sample_idx, int locus_idx)
{
    int idx = sample_idx * genotyping_data.num_loci + locus_idx;
    for (int kk = 0; kk < params.importance_sampling_depth; ++kk)
    {
        genotyping_llik_old.at(kk).at(idx) = genotyping_llik_new.at(kk).at(idx);
    }
}

void Chain::save_mean_coi_likelihood()
{
    mean_coi_prior_old = mean_coi_prior_new;
}

void Chain::save_coi_likelihood(int sample_idx)
{
    coi_prior_old[sample_idx] = coi_prior_new[sample_idx];
}

void Chain::save_eps_neg_likelihood(int sample_idx)
{
    eps_neg_prior_old[sample_idx] = eps_neg_prior_new[sample_idx];
}

void Chain::save_eps_pos_likelihood(int sample_idx)
{
    eps_pos_prior_old[sample_idx] = eps_pos_prior_new[sample_idx];
}

void Chain::restore_genotype_likelihood(int sample_idx, int locus_idx)
{
    int idx = sample_idx * genotyping_data.num_loci + locus_idx;
    for (int kk = 0; kk < params.importance_sampling_depth; ++kk)
    {
        genotyping_llik_new.at(kk).at(idx) = genotyping_llik_old.at(kk).at(idx);
    }
}

void Chain::restore_mean_coi_likelihood()
{
    mean_coi_prior_new = mean_coi_prior_old;
}

void Chain::restore_coi_likelihood(int sample_idx)
{
    coi_prior_new[sample_idx] = coi_prior_old[sample_idx];
}

void Chain::restore_eps_neg_likelihood(int sample_idx)
{
    eps_neg_prior_new[sample_idx] = eps_neg_prior_old[sample_idx];
}

void Chain::restore_eps_pos_likelihood(int sample_idx)
{
    eps_pos_prior_new[sample_idx] = eps_pos_prior_old[sample_idx];
}

Chain::Chain(GenotypingData genotyping_data, Lookup lookup, Parameters params)
    : genotyping_data(genotyping_data),
      lookup(lookup),
      params(params),
      sampler(lookup)

{
    eps_pos_var = params.eps_pos_var;
    eps_neg_var = params.eps_neg_var;
    m_prop_mean = std::vector<double>(genotyping_data.num_samples, 1);
    p_prop_var = std::vector<std::vector<double>>(genotyping_data.num_loci);
    p_accept = std::vector<std::vector<int>>(genotyping_data.num_loci);

    llik = std::numeric_limits<double>::lowest();

    initialize_p();
    initialize_m();
    initialize_eps_neg();
    initialize_eps_pos();
    initialize_mean_coi();
    initialize_likelihood();
}
