#include "chain.h"

#include "mcmc_utils.h"
#include "sampler.h"

#include <algorithm>
#include <map>

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
        int prop_m = m[i] + sampler.sample_coi_delta(2);

        if (prop_m > 0)
        {
            double sum_can = 0;
            double sum_orig = 0;

            // std::vector<double> coi_range(12, 0);
            for (size_t j = 0; j < genotyping_data.num_loci; j++)
            {
                // std::vector<double> coi_llik(12, 0);
                if (!genotyping_data.is_missing(j, i))
                {
                    // UtilFunctions::print_vector(
                        // genotyping_data.get_observed_alleles(j, i));
                    // for (int ii = 1; ii <= 12; ++ii)
                    // {
                    //     coi_llik[ii - 1] = calc_genotype_marginal_llik(
                    //         genotyping_data.get_observed_alleles(j, i), ii,
                    //         p[j], eps_neg[i], eps_pos[i], true);
                    //     coi_range[ii - 1] += coi_llik[ii - 1];
                    // }
                    // UtilFunctions::print_vector(coi_llik);
                    llik_new[j][i] = calc_genotype_marginal_llik(
                        genotyping_data.get_observed_alleles(j, i), prop_m,
                        p[j], eps_neg[i], eps_pos[i], true);
                    sum_can += llik_new[j][i];
                    sum_orig += llik_old[j][i];
                }
            }
            // UtilFunctions::print("COI Range:", m[i], prop_m);
            // UtilFunctions::print_vector(coi_range);

            // ZTPoisson prior on COI
            sum_can += sampler.get_coi_log_prob(prop_m, mean_coi);
            sum_orig += sampler.get_coi_log_prob(m[i], mean_coi);

            // UtilFunctions::print("M:", prop_m, m[i], sum_can, sum_orig);

            // Accept
            if ((sampler.sample_log_mh_acceptance() <= (sum_can - sum_orig)) or
                prop_m == m[i])
            {
                m[i] = prop_m;
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
                if (el < 1e-6)
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

                    llik_new[j][i] = calc_genotype_marginal_llik(
                        observed_alleles, m[i], prop_p, eps_neg[i], eps_pos[i],
                        true);

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
            UtilFunctions::print("Accept AF:", p_accept[j] / double(iteration));
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
                        prop_eps_neg, prop_eps_pos, false);
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
                        eps_neg[i], prop_eps_pos, false);
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
                        prop_eps_neg, eps_pos[i], false);
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
                        p[j], prop_eps_neg, prop_eps_pos, false);
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

    // double prop_pos = ((1 - epsilon_pos) / (1 - epsilon_pos + epsilon_neg));
    // double prop_neg = (epsilon_neg / (1 - epsilon_pos + epsilon_neg));

    // double obs_pos_mass = inv_tp_sum * prop_pos;
    // double obs_neg_mass = inv_fn_sum * prop_neg;

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

long double Chain::calc_exact_genotype_marginal_llik(
    std::vector<int> const &obs_genotype, int coi,
    std::vector<double> const &allele_frequencies, double epsilon_neg,
    double epsilon_pos)
{
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

long double Chain::importance_sample2(
    std::vector<int> const &obs_genotype, int coi, double epsilon_neg,
    double epsilon_pos, std::vector<double> const &allele_frequencies,
    int sampling_depth)
{
    int i = sampling_depth;
    double importance_prob;
    double true_prob;
    double obs_prob;
    double accumulated_prob = 0;
    double accumulated_weight = 0;

    std::vector<int> obs_pos_indices{};
    std::vector<int> obs_neg_indices{};
    std::vector<int> allele_index_vec{};

    // count the total number of positive alleles
    int observed_positives = 0;
    int observed_negatives = 0;
    int idx = 0;
    for (int allele : obs_genotype)
    {
        if (allele == 1)
        {
            observed_positives++;
            obs_pos_indices.push_back(idx);
        }
        else
        {
            observed_negatives++;
            obs_neg_indices.push_back(idx);
        }
        ++idx;
    }

    double fp_adjustment = .5;
    double fn_adjustment = .5;
    // Calculate the distribution over how many false positives there might be.
    // There must be a minimum of observed_positives - coi false positives
    double total = 1 - fp_adjustment;
    std::vector<double> fp_dist{1 - fp_adjustment};

    if (observed_positives > 0)
    {
        fp_dist.push_back(fp_adjustment);
        total = 1;
    }

    for (int ii = 2; ii <= observed_positives; ++ii)
    {
        fp_dist.push_back(fp_dist[ii - 1] * fp_adjustment);
        total += fp_dist.back();
    }

    // zero out impossible numbers of false positives, i.e. if the COI is
    // less than the number of observed alleles, then some number of these
    // alleles must be false positives.
    if (observed_positives > coi)
    {
        for (int ii = 0; ii < (observed_positives - coi); ++ii)
        {
            total -= fp_dist[ii];
            fp_dist[ii] = 0;
        }
    }

    for (int ii = 0; ii < fp_dist.size(); ++ii)
    {
        fp_dist[ii] = fp_dist[ii] / total;
    }

    // importance sampling
    std::vector<double> fn_dist{};
    while (--i >= 0)
    {
        double p = sampler.sample_unif();
        double cumsum = fp_dist[0];
        int total_false_positives = 0;
        while (cumsum <= p)
        {
            total_false_positives++;
            cumsum += fp_dist[total_false_positives];
        }
        double false_positive_prob = fp_dist[total_false_positives];
        int total_false_positive_sequences =
            std::tgamma(total_false_positives + 1);

        // Upper limit on the number of possible false negatives.
        // I.e., COI is 5, we observe 6, but sample 2 false positives -> total
        // possible false negatives is 1
        int possible_false_negatives =
            std::min(coi - observed_positives + total_false_positives,
                     observed_negatives);

        total = 1 - fn_adjustment;
        fn_dist.clear();
        fn_dist.push_back(1 - fn_adjustment);
        if (possible_false_negatives > 0)
        {
            fn_dist.push_back(fn_adjustment);
            total = 1;
        }

        for (int ii = 2; ii <= possible_false_negatives; ++ii)
        {
            fn_dist.push_back(fn_dist[ii - 1] * fn_adjustment);
            total += fn_dist.back();
        }

        for (int ii = 0; ii < fn_dist.size(); ++ii)
        {
            fn_dist[ii] = fn_dist[ii] / total;
        }

        p = sampler.sample_unif();
        cumsum = fn_dist[0];
        int total_false_negatives = 0;

        while (cumsum <= p)
        {
            total_false_negatives++;
            cumsum += fn_dist[total_false_negatives];
        }

        double false_negative_prob = fn_dist[total_false_negatives];
        int total_false_negative_sequences =
            std::tgamma(total_false_negatives + 1);

        std::vector<int> latent_genotype(obs_genotype);

        if (obs_pos_indices.size() > 0)
        {
            sampler.shuffle_vec(obs_pos_indices);
        }

        if (obs_neg_indices.size() > 0)
        {
            sampler.shuffle_vec(obs_neg_indices);
        }

        for (int ii = 0; ii < total_false_positives; ++ii)
        {
            latent_genotype[obs_pos_indices.at(ii)] = 0;
        }
        for (int ii = 0; ii < total_false_negatives; ++ii)
        {
            latent_genotype[obs_neg_indices[ii]] = 1;
        }

        allele_index_vec.clear();
        for (int ii = 0; ii < latent_genotype.size(); ++ii)
        {
            if (latent_genotype[ii] == 1)
            {
                allele_index_vec.push_back(ii);
            }
        }

        // require at least one allele to be observed. discard samples that do
        // not satisfy requirement.
        if (allele_index_vec.size() > 0)
        {
            importance_prob =
                total_false_positive_sequences * false_positive_prob *
                total_false_negative_sequences * false_negative_prob;

            true_prob = std::exp(calc_transmission_process(
                allele_index_vec, allele_frequencies, coi));

            obs_prob = std::exp(calc_observation_process(
                allele_index_vec, obs_genotype, coi, epsilon_neg, epsilon_pos));
            // UtilFunctions::print("Prob:", true_prob, obs_prob,
            // importance_prob, total_false_positives, total_false_negatives);
            double weight = true_prob / importance_prob;
            if (!std::isnan(weight) and !std::isnan(obs_prob))
            {
                accumulated_prob += (weight * obs_prob);
                accumulated_weight += weight;
            }
        }
    }
    double est = std::log(accumulated_prob / accumulated_weight);
    if (std::isnan(est))
    {
        UtilFunctions::print("NaN Encountered:", accumulated_prob,
                             accumulated_weight);
        UtilFunctions::print_vector(obs_genotype);
        UtilFunctions::print_vector(fp_dist);
        UtilFunctions::print_vector(allele_index_vec);
        UtilFunctions::print_vector(allele_frequencies);
    }
    return est;
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

long double Chain::calc_estimated_genotype_marginal_llik(
    std::vector<int> const &obs_genotype, int coi,
    std::vector<double> const &allele_frequencies, double epsilon_neg,
    double epsilon_pos, int sampling_depth, bool imp_sample)
{
    double est;
    // if (imp_sample)
    // {
    est = importance_sample(obs_genotype, coi, epsilon_neg, epsilon_pos,
                             allele_frequencies, sampling_depth);
    // }
    // else
    // {
    // est = monte_carlo_sample(obs_genotype, coi, epsilon_neg, epsilon_pos,
    // allele_frequencies, sampling_depth);
    // }

    return est;
}

long double Chain::calc_genotype_marginal_llik(
    std::vector<int> const &obs_genotype, int coi,
    std::vector<double> const &allele_frequencies, double epsilon_neg,
    double epsilon_pos, bool importance_sample)
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
            obs_genotype, coi, allele_frequencies, epsilon_neg, epsilon_pos,
            params.importance_sampling_depth +
                coi * params.importance_sampling_scaling_factor,
            importance_sample);

        return approx;
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
        auto eps_neg_prior = sampler.get_epsilon_log_prior(
            eps_neg[i], params.eps_neg_alpha, params.eps_neg_beta);
        auto eps_pos_prior = sampler.get_epsilon_log_prior(
            eps_pos[i], params.eps_pos_alpha, params.eps_pos_beta);
        auto coi_prior = sampler.get_coi_log_prob(m[i], mean_coi);

        llik += eps_neg_prior + eps_pos_prior + coi_prior;
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
