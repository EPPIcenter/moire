#include "chain.h"

#include "mcmc_utils.h"
#include "sampler.h"

#include <cmath>

#include <Rcpp/stats/binom.h>

#include <algorithm>
#include <limits>
#include <map>
#include <numeric>

void Chain::initialize_latent_genotypes()
{
    for (int jj = 0; jj < genotyping_data.num_loci; ++jj)
    {
        latent_genotypes_new.push_back(std::vector<std::vector<int>>{});
        latent_genotypes_old.push_back(std::vector<std::vector<int>>{});
        lg_adj_old.push_back(std::vector<double>{});
        lg_adj_new.push_back(std::vector<double>{});

        for (int ii = 0; ii < genotyping_data.num_samples; ++ii)
        {
            auto obs_alleles = genotyping_data.get_observed_alleles(jj, ii);
            std::vector<int> allele_index_vec{};
            auto lg =
                sampler.sample_latent_genotype(obs_alleles, m[ii], .2, .2);

            latent_genotypes_new.back().push_back(lg.value);
            lg_adj_new.back().push_back(lg.log_prob);

            latent_genotypes_old.back().push_back(lg.value);
            lg_adj_old.back().push_back(lg.log_prob);
        }
    }
}

// Initialize P with random allele frequencies
void Chain::initialize_p()
{
    for (size_t i = 0; i < genotyping_data.num_loci; i++)
    {
        p.push_back(sampler.sample_allele_frequencies(
            std::vector<double>(genotyping_data.num_alleles[i], 1), 1));
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
    m = std::vector<int>();
    for (const auto coi : genotyping_data.observed_coi)
    {
        int m_coi = coi + sampler.sample_coi_delta(3);
        m_coi = std::min(m_coi, params.max_coi);
        m_coi = std::max(1, m_coi);
        m.push_back(m_coi);
    }
    m_accept.resize(genotyping_data.num_samples, 0);
    sample_accept.resize(genotyping_data.num_samples, 0);
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

void Chain::initialize_r()
{
    if (params.allow_relatedness)
    {
        r.resize(genotyping_data.num_samples, 0.5);
    }
    else
    {
        r.resize(genotyping_data.num_samples, 0.0);
    }
    r_accept.resize(genotyping_data.num_samples, 0);
    r_var.resize(genotyping_data.num_samples, 1);
}

void Chain::update_m(int iteration)
{
    for (size_t ii = 0; ii < genotyping_data.num_samples; ii++)
    {
        int prop_m = m[ii] + sampler.sample_coi_delta(2);
        if (prop_m > 0 and prop_m <= params.max_coi)
        {
            double prev_m = m[ii];
            m[ii] = prop_m;

            double adj_ratio = 0;

            for (int jj = 0; jj < genotyping_data.num_loci; ++jj)
            {
                auto lg = sampler.sample_latent_genotype(
                    genotyping_data.get_observed_alleles(jj, ii), m[ii],
                    eps_pos[ii], eps_neg[ii]);
                latent_genotypes_new[jj][ii] = lg.value;
                lg_adj_new[jj][ii] = lg.log_prob;
                adj_ratio = adj_ratio + lg_adj_new[jj][ii] - lg_adj_old[jj][ii];
                calculate_genotype_likelihood(ii, jj);
            }

            double new_llik = calc_new_likelihood();
            double alpha = sampler.sample_log_mh_acceptance();
            if (!std::isnan(new_llik) and
                alpha <= (new_llik - llik + adj_ratio))
            {
                llik = new_llik;
                for (int jj = 0; jj < genotyping_data.num_loci; ++jj)
                {
                    save_genotype_likelihood(ii, jj);
                    lg_adj_old[jj][ii] = lg_adj_new[jj][ii];
                    latent_genotypes_old[jj][ii] = latent_genotypes_new[jj][ii];
                }
                m_accept[ii] += 1;
            }
            else
            {
                m[ii] = prev_m;
                for (int jj = 0; jj < genotyping_data.num_loci; ++jj)
                {
                    restore_genotype_likelihood(ii, jj);
                    lg_adj_new[jj][ii] = lg_adj_old[jj][ii];
                    latent_genotypes_new[jj][ii] = latent_genotypes_old[jj][ii];
                }
            }
        }
    }
}

void Chain::update_r(int iteration)
{
    for (size_t i = 0; i < m.size(); i++)
    {
        auto prop_adj = sampler.sample_constrained(r[i], r_var[i], 0, 1);
        double prop_r = std::get<0>(prop_adj);
        double adj = std::get<1>(prop_adj);

        if (prop_r < 1 && prop_r > 1e-32)
        {
            double prev_r = r[i];
            r[i] = prop_r;

            for (int jj = 0; jj < genotyping_data.num_loci; ++jj)
            {
                calculate_genotype_likelihood(i, jj);
            }

            double new_llik = calc_new_likelihood();

            // Reject
            if (std::isnan(new_llik) or
                sampler.sample_log_mh_acceptance() > (new_llik - llik + adj))
            {
                r[i] = prev_r;
                for (int jj = 0; jj < genotyping_data.num_loci; ++jj)
                {
                    restore_genotype_likelihood(i, jj);
                }
            }
            else
            {
                llik = new_llik;
                for (int jj = 0; jj < genotyping_data.num_loci; ++jj)
                {
                    save_genotype_likelihood(i, jj);
                }
                r_accept[i] += 1;
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
            p_attempt[j][idx] += 1;

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
                if (el < params.allele_freq_threshold)
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
            p[j] = prop_p;

            for (int ii = 0; ii < genotyping_data.num_samples; ++ii)
            {
                calculate_genotype_likelihood(ii, j);
            }

            double new_llik = calc_new_likelihood();

            double acceptanceRatio = new_llik - llik + logAdj;

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

            if (params.adapt_allele_freq_vars and iteration < params.burnin and
                p_attempt[j][idx] > 15)  // don't start adapting until there are
                                         // at least a few samples
            {
                double acceptanceRate = std::min(
                    (p_accept[j][idx] + 1) / (double(p_attempt[j][idx]) + 1),
                    1.0);
                p_prop_var[j][idx] += (acceptanceRate - .43) /
                                      std::pow(p_attempt[j][idx] + 1, .5);
                p_prop_var[j][idx] = std::max(p_prop_var[j][idx], .01);
            }
        }
    }
}

void Chain::update_eps_pos(int iteration)
{
    for (size_t i = 0; i < m.size(); i++)
    {
        auto prop_adj =
            sampler.sample_constrained(eps_pos[i], eps_pos_var, 0, 1);
        double prop_eps_pos = std::get<0>(prop_adj);
        double adj = std::get<1>(prop_adj);

        if (prop_eps_pos < 1 && prop_eps_pos > 1e-32)
        {
            double prev_eps_pos = eps_pos[i];
            eps_pos[i] = prop_eps_pos;
            calculate_eps_pos_likelihood(i);

            for (int jj = 0; jj < genotyping_data.num_loci; ++jj)
            {
                calculate_genotype_likelihood(i, jj);
            }

            double new_llik = calc_new_likelihood();

            // Reject
            if (std::isnan(new_llik) or
                sampler.sample_log_mh_acceptance() > (new_llik - llik + adj))
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
        auto prop_adj =
            sampler.sample_constrained(eps_neg[i], eps_neg_var, 0, 1);
        double prop_eps_neg = std::get<0>(prop_adj);
        double adj = std::get<1>(prop_adj);

        if (prop_eps_neg < 1 && prop_eps_neg > 1e-32)
        {
            double prev_eps_neg = eps_neg[i];
            eps_neg[i] = prop_eps_neg;
            calculate_eps_neg_likelihood(i);

            for (int jj = 0; jj < genotyping_data.num_loci; ++jj)
            {
                calculate_genotype_likelihood(i, jj);
            }

            double new_llik = calc_new_likelihood();

            // Reject
            if (std::isnan(new_llik) or
                sampler.sample_log_mh_acceptance() > (new_llik - llik + adj))
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

void Chain::update_samples(int iteration)
{
    for (size_t ii = 0; ii < m.size(); ii++)
    {
        auto eps_neg_prop_adj =
            sampler.sample_constrained(eps_neg[ii], eps_neg_var, 0, 1);
        double prop_eps_neg = std::get<0>(eps_neg_prop_adj);
        double eps_neg_adj = std::get<1>(eps_neg_prop_adj);
        bool valid_prop_eps_neg = prop_eps_neg < 1 && prop_eps_neg > 1e-32;

        auto eps_pos_prop_adj =
            sampler.sample_constrained(eps_pos[ii], eps_pos_var, 0, 1);
        double prop_eps_pos = std::get<0>(eps_pos_prop_adj);
        double eps_pos_adj = std::get<1>(eps_pos_prop_adj);
        bool valid_prop_eps_pos = prop_eps_pos < 1 && prop_eps_pos > 1e-32;

        double prop_r = 0;
        double r_adj = 1;
        bool valid_prop_r = true;
        if (params.allow_relatedness)
        {
            auto r_prop_adj =
                sampler.sample_constrained(r[ii], r_var[ii], 0, 1);
            double prop_r = std::get<0>(r_prop_adj);
            double r_adj = std::get<1>(r_prop_adj);
            bool valid_prop_r = prop_r < 1 && prop_r > 1e-32;
        }

        int prop_m = m[ii] + sampler.sample_coi_delta(2);
        bool valid_prop_m = prop_m > 0 and prop_m <= params.max_coi;

        if (valid_prop_eps_neg && valid_prop_eps_pos && valid_prop_r &&
            valid_prop_m)
        {
            double prev_eps_pos = eps_pos[ii];
            eps_pos[ii] = prop_eps_pos;
            calculate_eps_pos_likelihood(ii);

            double prev_eps_neg = eps_neg[ii];
            eps_neg[ii] = prop_eps_neg;
            calculate_eps_neg_likelihood(ii);

            double prev_r = r[ii];
            r[ii] = prop_r;

            double prev_m = m[ii];
            m[ii] = prop_m;

            double adj_ratio = eps_neg_adj + eps_pos_adj + r_adj;

            for (int jj = 0; jj < genotyping_data.num_loci; ++jj)
            {
                auto lg = sampler.sample_latent_genotype(
                    genotyping_data.get_observed_alleles(jj, ii), m[ii],
                    eps_pos[ii], eps_neg[ii]);
                latent_genotypes_new[jj][ii] = lg.value;
                lg_adj_new[jj][ii] = lg.log_prob;
                adj_ratio = adj_ratio + lg_adj_new[jj][ii] - lg_adj_old[jj][ii];
                calculate_genotype_likelihood(ii, jj);
            }

            double new_llik = calc_new_likelihood();
            double alpha = sampler.sample_log_mh_acceptance();

            if (!std::isnan(new_llik) and
                alpha <= (new_llik - llik + adj_ratio))
            {
                llik = new_llik;
                save_eps_neg_likelihood(ii);
                save_eps_pos_likelihood(ii);
                for (int jj = 0; jj < genotyping_data.num_loci; ++jj)
                {
                    save_genotype_likelihood(ii, jj);
                    lg_adj_old[jj][ii] = lg_adj_new[jj][ii];
                    latent_genotypes_old[jj][ii] = latent_genotypes_new[jj][ii];
                }
                sample_accept[ii] += 1;
            }
            else
            {
                m[ii] = prev_m;
                eps_pos[ii] = prev_eps_pos;
                eps_neg[ii] = prev_eps_neg;
                r[ii] = prev_r;
                restore_eps_neg_likelihood(ii);
                restore_eps_pos_likelihood(ii);
                for (int jj = 0; jj < genotyping_data.num_loci; ++jj)
                {
                    restore_genotype_likelihood(ii, jj);
                    lg_adj_new[jj][ii] = lg_adj_old[jj][ii];
                    latent_genotypes_new[jj][ii] = latent_genotypes_old[jj][ii];
                }
            }
        }
    }
}

double Chain::calc_transmission_process(
    std::vector<int> const &allele_index_vec,
    std::vector<double> const &allele_frequencies, int coi, double relatedness)
{
    // transmission process - prob that after "coi" number of draws, all
    // alleles are drawn at least once conditional on all draws come
    // from the constrained set, where the constrained set is the set of
    // positive alleles in the latent genotype

    double constrained_set_total_prob = 0;
    std::vector<double> res{};
    int total_alleles = allele_index_vec.size();

    prVec_.clear();
    prVec_.reserve(total_alleles);
    res.reserve(coi - total_alleles + 1);
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

    if (!params.allow_relatedness)
    {
        return std::log(1 - probAnyMissing_(prVec_, coi)) +
               std::log(constrained_set_total_prob) * coi;
    }

    for (size_t i = total_alleles; i <= coi; ++i)
    {
        double pr = R::dbinom(i, coi, 1 - relatedness, true);
        double i_res = std::log(1 - probAnyMissing_(prVec_, i)) +
                       std::log(constrained_set_total_prob) * i;
        res.push_back(pr + i_res);
    }

    double out = UtilFunctions::logSumExp(res);
    return out;
}

double Chain::calc_observation_process(std::vector<int> const &allele_index_vec,
                                       std::vector<int> const &obs_genotype,
                                       int coi, double epsilon_neg,
                                       double epsilon_pos)
{
    double res = 0;
    unsigned int fp = 0;
    unsigned int tp = 0;
    unsigned int fn = 0;
    unsigned int tn = 0;

    unsigned int vec_pointer = 0;
    unsigned int next_allele_index = allele_index_vec[vec_pointer];
    unsigned int total_alleles = allele_index_vec.size();

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

    int total_possible_alleles = obs_genotype.size();
    res += std::log(1 - (epsilon_neg * params.max_eps_neg /
                         total_possible_alleles)) *
           tp;
    res += std::log(epsilon_neg * params.max_eps_neg / total_possible_alleles) *
           fn;
    res += std::log(1 - (epsilon_pos * params.max_eps_pos /
                         total_possible_alleles)) *
           tn;
    res += std::log(epsilon_pos * params.max_eps_pos / total_possible_alleles) *
           fp;

    return res;
};

double Chain::calc_genotype_log_pmf(
    std::vector<int> const &allele_index_vec,
    std::vector<int> const &obs_genotype, double epsilon_pos,
    double epsilon_neg, int coi, double relatedness,
    std::vector<double> const &allele_frequencies)
{
    double res = 0.0;
    res += calc_transmission_process(allele_index_vec, allele_frequencies, coi,
                                     relatedness);
    res += calc_observation_process(allele_index_vec, obs_genotype, coi,
                                    epsilon_neg, epsilon_pos);

    return res;
}

double Chain::calc_old_likelihood()
{
    double old_llik = 0;

    for (int ii = 0; ii < genotyping_data.num_samples; ++ii)
    {
        old_llik += eps_neg_prior_old[ii];
        old_llik += eps_pos_prior_old[ii];
    }

    for (const double ll : genotyping_llik_old)
    {
        old_llik += ll;
    }

    return old_llik;
}

double Chain::calc_new_likelihood()
{
    double new_llik = 0;

    for (int ii = 0; ii < genotyping_data.num_samples; ++ii)
    {
        new_llik += eps_neg_prior_new[ii];
        new_llik += eps_pos_prior_new[ii];
    }

    for (const double ll : genotyping_llik_new)
    {
        new_llik += ll;
    }

    return new_llik;
}

double Chain::get_llik() { return llik; }

void Chain::calculate_genotype_likelihood(int sample_idx, int locus_idx)
{
    int idx = sample_idx * genotyping_data.num_loci + locus_idx;
    if (genotyping_data.is_missing(locus_idx, sample_idx))
    {
        genotyping_llik_new[idx] = 0;
    }
    else
    {
        double res;
        double transmission_prob = calc_transmission_process(
            latent_genotypes_new[locus_idx][sample_idx], p[locus_idx],
            m[sample_idx], r[sample_idx]);
        double obs_prob = calc_observation_process(
            latent_genotypes_new[locus_idx][sample_idx],
            genotyping_data.get_observed_alleles(locus_idx, sample_idx),
            m[sample_idx], eps_neg[sample_idx], eps_pos[sample_idx]);
        res = transmission_prob + obs_prob;
        genotyping_llik_new[idx] = res;
    }
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
    int row_idx, idx;

    genotyping_llik_new.resize(num_samples * num_loci);
    genotyping_llik_old.resize(num_samples * num_loci);

    eps_neg_prior_new.resize(num_samples);
    eps_neg_prior_old.resize(num_samples);
    eps_pos_prior_new.resize(num_samples);
    eps_pos_prior_old.resize(num_samples);
    for (int ii = 0; ii < genotyping_data.num_samples; ++ii)
    {
        row_idx = ii * num_loci;
        for (int jj = 0; jj < genotyping_data.num_loci; ++jj)
        {
            idx = row_idx + jj;
            calculate_genotype_likelihood(ii, jj);
            genotyping_llik_old[idx] = genotyping_llik_new[idx];
        }
    }

    for (size_t ii = 0; ii < genotyping_data.num_samples; ii++)
    {
        double eps_neg_prior = sampler.get_epsilon_log_prior(
            eps_neg[ii], params.eps_neg_alpha, params.eps_neg_beta);
        double eps_pos_prior = sampler.get_epsilon_log_prior(
            eps_pos[ii], params.eps_pos_alpha, params.eps_pos_beta);

        eps_neg_prior_old[ii] = eps_neg_prior;
        eps_neg_prior_new[ii] = eps_neg_prior;
        eps_pos_prior_old[ii] = eps_pos_prior;
        eps_pos_prior_new[ii] = eps_pos_prior;
    }

    llik = calc_new_likelihood();
}

void Chain::save_genotype_likelihood(int sample_idx, int locus_idx)
{
    int idx = sample_idx * genotyping_data.num_loci + locus_idx;

    genotyping_llik_old.at(idx) = genotyping_llik_new.at(idx);
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
    genotyping_llik_new.at(idx) = genotyping_llik_old.at(idx);
}

void Chain::restore_eps_neg_likelihood(int sample_idx)
{
    eps_neg_prior_new[sample_idx] = eps_neg_prior_old[sample_idx];
}

void Chain::restore_eps_pos_likelihood(int sample_idx)
{
    eps_pos_prior_new[sample_idx] = eps_pos_prior_old[sample_idx];
}

Chain::Chain(GenotypingData genotyping_data, Parameters params)
    : genotyping_data(genotyping_data), params(params), sampler()

{
    eps_pos_var = params.eps_pos_var;
    eps_neg_var = params.eps_neg_var;
    m_prop_mean = std::vector<double>(genotyping_data.num_samples, 1);
    p_prop_var = std::vector<std::vector<double>>(genotyping_data.num_loci);
    p_accept = std::vector<std::vector<int>>(genotyping_data.num_loci);

    llik = std::numeric_limits<double>::lowest();

    initialize_m();
    initialize_eps_neg();
    initialize_eps_pos();
    initialize_r();
    initialize_latent_genotypes();
    initialize_p();
    initialize_likelihood();
}
