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
};

void Chain::initialize_m()
{
    m = genotyping_data.observed_coi;
    m_accept.resize(genotyping_data.num_samples, 0);
    individual_accept.resize(genotyping_data.num_samples, 0);
}

void Chain::initialize_eps_neg()
{
    eps_neg.resize(genotyping_data.num_samples, params.eps_neg_0);
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
        double sum_can = 0;
        double sum_orig = 0;
        for (size_t ii = 0; ii < genotyping_data.num_samples; ii++)
        {
            sum_can += sampler.get_coi_log_prob(m[ii], prop_mean_coi);
            sum_orig += sampler.get_coi_log_prob(m[ii], mean_coi);
        }

        sum_can += sampler.get_coi_mean_log_prior(prop_mean_coi,
                                                  params.mean_coi_prior_shape,
                                                  params.mean_coi_prior_scale);

        sum_orig += sampler.get_coi_mean_log_prior(
            mean_coi, params.mean_coi_prior_shape, params.mean_coi_prior_scale);

        if (sampler.sample_log_mh_acceptance() <= (sum_can - sum_orig))
        {
            mean_coi = prop_mean_coi;
        }
    }
}

void Chain::update_m(int iteration)
{
    for (size_t i = 0; i < genotyping_data.num_samples; i++)
    {
        int prop_m = m[i] + sampler.sample_coi_delta(3);

        if (prop_m > 0)
        {
            double sum_can = 0;
            double sum_orig = 0;

            for (size_t j = 0; j < genotyping_data.num_loci; j++)
            {
                if (!genotyping_data.is_missing(j, i))
                {
                    llik_new[j][i] = calc_genotype_marginal_llik(
                        genotyping_data.get_observed_alleles(j, i), prop_m,
                        p[j], eps_neg[i], eps_pos[i]);
                    sum_can += llik_new[j][i];
                    sum_orig += llik_old[j][i];
                }
            }

            // ZTPoisson prior on COI
            sum_can += sampler.get_coi_log_prob(prop_m, mean_coi);
            sum_orig += sampler.get_coi_log_prob(m[i], mean_coi);

            // Accept
            if (sampler.sample_log_mh_acceptance() <= (sum_can - sum_orig))
            {
                for (size_t j = 0; j < genotyping_data.num_loci; j++)
                {
                    llik_old[j][i] = llik_new[j][i];
                }
                m_accept[i] += 1;
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
        int rep = 1;
        while (--rep >= 0)
        {
            int k = p[j].size();
            const int idx = sampler.sample_random_int(0, k - 1);

            auto logitPropP = UtilFunctions::logitVec(p[j]);

            double logitCurr = logitPropP[idx];
            double logitProp =
                sampler.sample_epsilon(logitCurr, params.allele_freq_var);

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
            for (const auto &el : prop_p)
            {
                if (el < 1e-12)
                {
                    return;
                }
            }

            double sum_can = 0;
            double sum_orig = 0;
            for (size_t i = 0; i < genotyping_data.num_samples; i++)
            {
                if (!genotyping_data.is_missing(j, i))
                {
                    const auto observed_alleles =
                        genotyping_data.get_observed_alleles(j, i);
                    auto emphasized_alleles = observed_alleles;
                    emphasized_alleles[idx] = 1;

                    llik_new[j][i] = calc_genotype_marginal_llik(
                        observed_alleles, emphasized_alleles, m[i], prop_p,
                        eps_neg[i], eps_pos[i]);
                    sum_can += llik_new[j][i];
                    sum_orig += llik_old[j][i];
                }
            }

            double acceptanceRatio = sum_can - sum_orig + logAdj;
            if (sampler.sample_log_mh_acceptance() <= acceptanceRatio)
            {
                p[j] = prop_p;
                p_accept[j] += 1;
                for (size_t i = 0; i < genotyping_data.num_samples; i++)
                {
                    llik_old[j][i] = llik_new[j][i];
                }
            }
        }
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
                if (!genotyping_data.is_missing(j, i))
                {
                    llik_new[j][i] = calc_genotype_marginal_llik(
                        genotyping_data.get_observed_alleles(j, i), m[i], p[j],
                        prop_eps_neg, prop_eps_pos);
                    sum_can += llik_new[j][i];
                    sum_orig += llik_old[j][i];
                }
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
                eps_pos_accept[i] += 1;

                eps_neg[i] = prop_eps_neg;
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

        if (prop_eps_pos < params.max_eps_pos && prop_eps_pos > 0)
        {
            double sum_can = 0;
            double sum_orig = 0;

            for (size_t j = 0; j < genotyping_data.num_loci; j++)
            {
                if (!genotyping_data.is_missing(j, i))
                {
                    llik_new[j][i] = calc_genotype_marginal_llik(
                        genotyping_data.get_observed_alleles(j, i), m[i], p[j],
                        eps_neg[i], prop_eps_pos);
                    sum_can += llik_new[j][i];
                    sum_orig += llik_old[j][i];
                }
            }

            // Incorporate prior
            sum_can += sampler.get_epsilon_log_prior(
                prop_eps_pos, params.eps_pos_alpha, params.eps_pos_beta);
            sum_orig += sampler.get_epsilon_log_prior(
                eps_pos[i], params.eps_pos_alpha, params.eps_pos_beta);

            // Accept
            if (sampler.sample_log_mh_acceptance() <= (sum_can - sum_orig))
            {
                eps_pos[i] = prop_eps_pos;
                eps_pos_accept[i] += 1;
                for (size_t j = 0; j < genotyping_data.num_loci; j++)
                {
                    llik_old[j][i] = llik_new[j][i];
                }
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

        if (prop_eps_neg < params.max_eps_neg && prop_eps_neg > 0)
        {
            double sum_can = 0;
            double sum_orig = 0;

            for (size_t j = 0; j < genotyping_data.num_loci; j++)
            {
                if (!genotyping_data.is_missing(j, i))
                {
                    llik_new[j][i] = calc_genotype_marginal_llik(
                        genotyping_data.get_observed_alleles(j, i), m[i], p[j],
                        prop_eps_neg, eps_pos[i]);
                    sum_can += llik_new[j][i];
                    sum_orig += llik_old[j][i];
                }
            }

            // // Incorporate prior
            sum_can += sampler.get_epsilon_log_prior(
                prop_eps_neg, params.eps_neg_alpha, params.eps_neg_beta);
            sum_orig += sampler.get_epsilon_log_prior(
                eps_neg[i], params.eps_neg_alpha, params.eps_neg_beta);

            // Accept
            if (sampler.sample_log_mh_acceptance() <= (sum_can - sum_orig))
            {
                eps_neg[i] = prop_eps_neg;
                eps_neg_accept[i] += 1;
                for (size_t j = 0; j < genotyping_data.num_loci; j++)
                {
                    llik_old[j][i] = llik_new[j][i];
                }
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
            prop_eps_pos < params.max_eps_pos && prop_eps_pos > 0 && prop_m &&
            prop_m > 0)
        {
            double sum_can = 0;
            double sum_orig = 0;

            for (size_t j = 0; j < genotyping_data.num_loci; j++)
            {
                if (!genotyping_data.is_missing(j, i))
                {
                    llik_new[j][i] = calc_genotype_marginal_llik(
                        genotyping_data.get_observed_alleles(j, i), prop_m,
                        p[j], prop_eps_neg, prop_eps_pos);
                    sum_can += llik_new[j][i];
                    sum_orig += llik_old[j][i];
                }
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
            sum_can += sampler.get_coi_log_prob(prop_m, mean_coi);
            sum_orig += sampler.get_coi_log_prob(m[i], mean_coi);

            // Accept
            if (sampler.sample_log_mh_acceptance() <= (sum_can - sum_orig))
            {
                m[i] = prop_m;
                eps_neg[i] = prop_eps_neg;
                eps_pos[i] = prop_eps_pos;
                individual_accept[i] += 1;
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

    // double obs_pos_mass =
    //     inv_tp_sum * ((1 - epsilon_pos) / (1 - epsilon_pos + epsilon_neg));
    // double obs_neg_mass =
    //     inv_fn_sum * (epsilon_neg / (1 - epsilon_pos + epsilon_neg));
    //
    double obs_pos_mass = inv_tp_sum * .9;
    double obs_neg_mass = inv_fn_sum * .1;

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

    // for (size_t j = 0; j < total_alleles; j++)
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

        if (vec_pointer < coi)
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

long double Chain::calc_exact_genotype_marginal_llik(
    std::vector<int> const &obs_genotype, int coi,
    std::vector<double> const &allele_frequencies, double epsilon_neg,
    double epsilon_pos)
{
    long double res = 0;
    for (int i = 1; i <= coi; i++)
    {
        Rcpp::checkUserInterrupt();
        allele_index_generator_.reset(allele_frequencies.size(), i);
        while (!allele_index_generator_.completed)
        {
            res += std::exp(calc_genotype_log_pmf(
                allele_index_generator_.curr, obs_genotype, epsilon_pos,
                epsilon_neg, coi, allele_frequencies));
            allele_index_generator_.next();
        }
    }
    return log(res);
}

long double Chain::calc_estimated_genotype_marginal_llik(
    std::vector<int> const &obs_genotype,
    std::vector<int> const &emphasized_alleles, int coi,
    std::vector<double> const &allele_frequencies, double epsilon_neg,
    double epsilon_pos, int sampling_depth)
{
    int i = sampling_depth;
    double importance_weight = 0;

    std::vector<int> allele_index_vec;
    std::vector<double> reweighted_allele_frequencies =
        reweight_allele_frequencies(allele_frequencies, emphasized_alleles,
                                    epsilon_neg, epsilon_pos, coi);

    double est = 0.0;
    double val = 0.0;
    std::map<std::vector<int>, double> memo{};

    while (--i >= 0)
    {
        allele_index_vec =
            sampler.sample_latent_genotype(coi, reweighted_allele_frequencies);
        auto s = memo.find(allele_index_vec);
        if (s != memo.end())
        {
            val = s->second;
        }
        else
        {
            importance_weight = calc_transmission_process(
                allele_index_vec, reweighted_allele_frequencies, coi);
            val = std::exp(calc_genotype_log_pmf(allele_index_vec, obs_genotype,
                                                 epsilon_pos, epsilon_neg, coi,
                                                 allele_frequencies) -
                           importance_weight);
            memo[allele_index_vec] = val;
        }
        est += val;
    }

    est = std::log(est / sampling_depth);
    return est;
}

long double Chain::calc_genotype_marginal_llik(
    std::vector<int> const &obs_genotype,
    std::vector<int> const &emphasized_alleles, int coi,
    std::vector<double> const &allele_frequencies, double epsilon_neg,
    double epsilon_pos)
{
    double log_total_combinations =
        lookup.get_sampling_depth(coi, allele_frequencies.size());

    if (log_total_combinations <= std::log(params.complexity_limit))
    {
        return calc_exact_genotype_marginal_llik(
            obs_genotype, coi, allele_frequencies, epsilon_neg, epsilon_pos);
    }
    else
    {
        double approx = calc_estimated_genotype_marginal_llik(
            obs_genotype, emphasized_alleles, coi, allele_frequencies,
            epsilon_neg, epsilon_pos,
            params.importance_sampling_depth +
                coi * params.importance_sampling_scaling_factor);

        return approx;
    }
}

long double Chain::calc_genotype_marginal_llik(
    std::vector<int> const &obs_genotype, int coi,
    std::vector<double> const &allele_frequencies, double epsilon_neg,
    double epsilon_pos)
{
    return calc_genotype_marginal_llik(obs_genotype, obs_genotype, coi,
                                       allele_frequencies, epsilon_neg,
                                       epsilon_pos);
}

void Chain::initialize_likelihood()
{
    for (size_t j = 0; j < genotyping_data.num_loci; j++)
    {
        llik_old.push_back(std::vector<double>(genotyping_data.num_samples));
        llik_new.push_back(std::vector<double>(genotyping_data.num_samples));
        for (size_t i = 0; i < genotyping_data.num_samples; i++)
        {
            double marginal_llik = 0;
            if (!genotyping_data.is_missing(j, i))
            {
                marginal_llik = calc_genotype_marginal_llik(
                    genotyping_data.get_observed_alleles(j, i), m[i], p[j],
                    eps_neg[i], eps_pos[i]);
            }
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

    for (size_t i = 0; i < genotyping_data.num_samples; i++)
    {
        llik += sampler.get_epsilon_log_prior(eps_neg[i], params.eps_neg_alpha,
                                              params.eps_neg_beta);
        llik += sampler.get_epsilon_log_prior(eps_pos[i], params.eps_pos_alpha,
                                              params.eps_pos_beta);
        llik += sampler.get_coi_log_prob(m[i], mean_coi);
    }

    llik += sampler.get_coi_mean_log_prior(
        mean_coi, params.mean_coi_prior_shape, params.mean_coi_prior_scale);
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
      sampler(lookup)

{
    llik = 0;
    eps_pos_var = params.eps_pos_var;
    eps_neg_var = params.eps_neg_var;
    m_prop_mean = std::vector<double>(genotyping_data.num_samples, 1);
    p_prop_var = std::vector<double>(genotyping_data.num_loci, 1);

    initialize_p();
    initialize_m();
    initialize_eps_neg();
    initialize_eps_pos();
    initialize_mean_coi();
    initialize_likelihood();
};
