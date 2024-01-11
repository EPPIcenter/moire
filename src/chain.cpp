#include "chain.h"

#include "mcmc_utils.h"
#include "sampler.h"

#include <cmath>
#include <algorithm>

#if defined(__cpp_lib_execution) && __cpp_lib_execution >= 201603
#include <execution>
#define HAS_EXECUTION 1
#endif

#include <limits>
#include <map>
#include <numeric>

constexpr float min_sampled = std::numeric_limits<float>::min();

void Chain::initialize_latent_genotypes()
{
    for (std::size_t jj = 0; jj < genotyping_data.num_loci; ++jj)
    {
        latent_genotypes_new.push_back(std::vector<std::vector<int>>{});
        latent_genotypes_old.push_back(std::vector<std::vector<int>>{});
        lg_adj_old.push_back(std::vector<float>{});
        lg_adj_new.push_back(std::vector<float>{});

        for (std::size_t ii = 0; ii < genotyping_data.num_samples; ++ii)
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
            std::vector<float>(genotyping_data.num_alleles[i], 1), 1));
    }

    p_prop_var.resize(genotyping_data.num_loci);
    p_accept.resize(genotyping_data.num_loci);
    p_attempt.resize(genotyping_data.num_loci);

    for (std::size_t jj = 0; jj < genotyping_data.num_loci; ++jj)
    {
        const int total_alleles = p[jj].size();
        p_prop_var[jj] = std::vector<float>(total_alleles, 1);
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
    eps_neg.resize(genotyping_data.num_samples, sampler.sample_unif() * .1);
    eps_neg_accept.resize(genotyping_data.num_samples, 0);
    eps_neg_var.resize(genotyping_data.num_samples, 1);
}

void Chain::initialize_eps_pos()
{
    eps_pos.resize(genotyping_data.num_samples, sampler.sample_unif() * .1);
    eps_pos_accept.resize(genotyping_data.num_samples, 0);
    eps_pos_var.resize(genotyping_data.num_samples, 1);
}

void Chain::initialize_r()
{
    if (params.allow_relatedness)
    {
        for (size_t i = 0; i < genotyping_data.num_samples; ++i)
        {
            r.push_back(sampler.sample_unif() * .99);
        }
    }
    else
    {
        r.resize(genotyping_data.num_samples, 0.0);
    }
    r_accept.resize(genotyping_data.num_samples, 0);
    r_var.resize(genotyping_data.num_samples, 1);
    m_r_accept.resize(genotyping_data.num_samples, 0);
    m_r_var.resize(genotyping_data.num_samples, 1);
}

void Chain::update_m(int iteration)
{
    auto indices = std::vector<int>(genotyping_data.num_samples);
    std::iota(indices.begin(), indices.end(), 0);
    sampler.shuffle_vec(indices);

    for (const auto ii : indices)
    {
        const int prop_m = m[ii] + sampler.sample_coi_delta(2);
        if (prop_m > 0 and prop_m <= params.max_coi)
        {
            const float prev_m = m[ii];
            m[ii] = prop_m;
            calculate_coi_likelihood(ii);

            float adj_ratio = 0;

            for (std::size_t jj = 0; jj < genotyping_data.num_loci; ++jj)
            {
                auto lg = sampler.sample_latent_genotype(
                    genotyping_data.get_observed_alleles(jj, ii), m[ii],
                    eps_pos[ii], eps_neg[ii]);
                latent_genotypes_new[jj][ii] = lg.value;
                lg_adj_new[jj][ii] = lg.log_prob;
                adj_ratio = adj_ratio + lg_adj_new[jj][ii] - lg_adj_old[jj][ii];
                calculate_genotype_likelihood(ii, jj);
            }

            const float new_llik = calc_new_likelihood();
            const float new_prior = calc_new_prior();
            const float new_post = new_llik * temp + new_prior;
            const float alpha = sampler.sample_log_mh_acceptance();
            if (!std::isnan(new_post) and
                alpha <= (new_post - get_posterior() + adj_ratio))
            {
                llik = new_llik;
                prior = new_prior;
                save_coi_likelihood(ii);
                for (std::size_t jj = 0; jj < genotyping_data.num_loci; ++jj)
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
                restore_coi_likelihood(ii);
                for (std::size_t jj = 0; jj < genotyping_data.num_loci; ++jj)
                {
                    restore_genotype_likelihood(ii, jj);
                    lg_adj_new[jj][ii] = lg_adj_old[jj][ii];
                    latent_genotypes_new[jj][ii] = latent_genotypes_old[jj][ii];
                }
            }
        }
    }
}

void Chain::update_eff_coi(int iteration)
{
    auto indices = std::vector<int>(genotyping_data.num_samples);
    std::iota(indices.begin(), indices.end(), 0);
    sampler.shuffle_vec(indices);

    for (const int i : indices)
    {
        const float curr_eff_coi = (m[i] - 1) * (1.0f - r[i]) + 1.0f;
        const auto prop_adj = sampler.sample_constrained(
            curr_eff_coi, m_r_var[i], 1, params.max_coi);

        const float prop_eff_coi = std::get<0>(prop_adj);
        const float adj = std::get<1>(prop_adj);

        const int prop_m = m[i] + sampler.sample_coi_delta(2);
        const float prop_r = 1.0f - (prop_eff_coi - 1.0f) / (prop_m - 1.0f);

        if (prop_m <= 0 or prop_m > params.max_coi or prop_r > 1.0f or
            prop_r < 1e-32)
            continue;

        const int prev_m = m[i];
        const float prev_r = r[i];
        m[i] = prop_m;
        r[i] = prop_r;
        calculate_coi_likelihood(i);
        calculate_relatedness_likelihood(i);

        float adj_ratio = adj;

        for (std::size_t jj = 0; jj < genotyping_data.num_loci; ++jj)
        {
            const auto &lg = sampler.sample_latent_genotype(
                genotyping_data.get_observed_alleles(jj, i), m[i], eps_pos[i],
                eps_neg[i]);
            latent_genotypes_new[jj][i] = lg.value;
            lg_adj_new[jj][i] = lg.log_prob;
            adj_ratio = adj_ratio + lg_adj_new[jj][i] - lg_adj_old[jj][i];
            calculate_genotype_likelihood(i, jj);
        }

        const float new_llik = calc_new_likelihood();
        const float new_prior = calc_new_prior();
        const float new_post = new_llik * temp + new_prior;
        const double mh_ratio = new_post - get_posterior() + adj_ratio;
        const double alpha = sampler.sample_log_mh_acceptance();

        // Reject
        if (std::isnan(new_post) or std::isnan(mh_ratio) or alpha > mh_ratio)
        {
            m[i] = prev_m;
            r[i] = prev_r;
            restore_coi_likelihood(i);
            restore_relatedness_likelihood(i);
            for (std::size_t jj = 0; jj < genotyping_data.num_loci; ++jj)
            {
                restore_genotype_likelihood(i, jj);
                lg_adj_new[jj][i] = lg_adj_old[jj][i];
                latent_genotypes_new[jj][i] = latent_genotypes_old[jj][i];
            }
        }
        else
        {
            llik = new_llik;
            prior = new_prior;
            save_relatedness_likelihood(i);
            save_coi_likelihood(i);
            for (std::size_t jj = 0; jj < genotyping_data.num_loci; ++jj)
            {
                save_genotype_likelihood(i, jj);
                lg_adj_old[jj][i] = lg_adj_new[jj][i];
                latent_genotypes_old[jj][i] = latent_genotypes_new[jj][i];
            }
            m_r_accept[i] += 1;
        }

        if (iteration < params.burnin and
            iteration > 15)  // don't start adapting until there are
                             // at least a few samples
        {
            const float acceptanceRate = m_r_accept[i] / float(iteration);
            const float update =
                (acceptanceRate - .23) / std::pow(iteration + 1, .5);
            m_r_var[i] = std::max(m_r_var[i] + update, .01f);
        }
    }
}

void Chain::update_r(int iteration)
{
    auto indices = std::vector<int>(genotyping_data.num_samples);
    std::iota(indices.begin(), indices.end(), 0);
    sampler.shuffle_vec(indices);

    for (const auto i : indices)
    {
        const auto prop_adj =
            sampler.sample_constrained(r[i], r_var[i], min_sampled, .99);
        const float prop_r = std::get<0>(prop_adj);
        const float adj = std::get<1>(prop_adj);

        const float prev_r = r[i];
        r[i] = prop_r;
        calculate_relatedness_likelihood(i);

        for (std::size_t jj = 0; jj < genotyping_data.num_loci; ++jj)
        {
            calculate_genotype_likelihood(i, jj);
        }

        const float new_llik = calc_new_likelihood();
        const float new_prior = calc_new_prior();
        const float new_post = new_llik * temp + new_prior;

        // Reject
        if (std::isnan(new_post) or sampler.sample_log_mh_acceptance() >
                                        (new_post - get_posterior() + adj))
        {
            r[i] = prev_r;
            restore_relatedness_likelihood(i);
            for (std::size_t jj = 0; jj < genotyping_data.num_loci; ++jj)
            {
                restore_genotype_likelihood(i, jj);
            }
        }
        else
        {
            llik = new_llik;
            prior = new_prior;
            save_relatedness_likelihood(i);
            for (std::size_t jj = 0; jj < genotyping_data.num_loci; ++jj)
            {
                save_genotype_likelihood(i, jj);
            }
            r_accept[i] += 1;
        }

        if (iteration < params.burnin and
            iteration > 15)  // don't start adapting until there are
                             // at least a few samples
        {
            const float acceptanceRate = r_accept[i] / float(iteration);
            const float update =
                (acceptanceRate - .23) / std::pow(iteration + 1, .5);
            r_var[i] = std::max(r_var[i] + update, .01f);
        }
    }
}

void Chain::update_m_r(int iteration)
{
    auto indices = std::vector<int>(genotyping_data.num_samples);
    std::iota(indices.begin(), indices.end(), 0);
    sampler.shuffle_vec(indices);

    for (const auto i : indices)
    {
        const auto prop_adj =
            sampler.sample_constrained(r[i], m_r_var[i], min_sampled, .99);
        const float prop_r = std::get<0>(prop_adj);
        const float adj = std::get<1>(prop_adj);

        const int prop_m = m[i] + sampler.sample_coi_delta(1);

        if (prop_m <= 0 or prop_m > params.max_coi) return;

        const float prev_m = m[i];
        const float prev_r = r[i];
        r[i] = prop_r;
        calculate_relatedness_likelihood(i);
        m[i] = prop_m;
        calculate_coi_likelihood(i);

        float adj_ratio = adj;

        for (std::size_t jj = 0; jj < genotyping_data.num_loci; ++jj)
        {
            const auto &lg = sampler.sample_latent_genotype(
                genotyping_data.get_observed_alleles(jj, i), m[i], eps_pos[i],
                eps_neg[i]);
            // const auto &lg = sampler.sample_latent_genotype(
            //     latent_genotypes_old[jj][i], m[i], eps_pos[i], eps_neg[i]);
            latent_genotypes_new[jj][i] = lg.value;
            lg_adj_new[jj][i] = lg.log_prob;
            adj_ratio = adj_ratio + lg_adj_new[jj][i] - lg_adj_old[jj][i];
            calculate_genotype_likelihood(i, jj);
        }

        const float new_llik = calc_new_likelihood();
        const float new_prior = calc_new_prior();
        const float new_post = new_llik * temp + new_prior;

        // Reject
        if (std::isnan(new_post) or
            sampler.sample_log_mh_acceptance() >
                (new_post - get_posterior() + adj_ratio))
        {
            r[i] = prev_r;
            restore_relatedness_likelihood(i);
            m[i] = prev_m;
            restore_coi_likelihood(i);
            for (std::size_t jj = 0; jj < genotyping_data.num_loci; ++jj)
            {
                restore_genotype_likelihood(i, jj);
                lg_adj_new[jj][i] = lg_adj_old[jj][i];
                latent_genotypes_new[jj][i] = latent_genotypes_old[jj][i];
            }
        }
        else
        {
            llik = new_llik;
            prior = new_prior;
            save_relatedness_likelihood(i);
            save_coi_likelihood(i);
            for (std::size_t jj = 0; jj < genotyping_data.num_loci; ++jj)
            {
                save_genotype_likelihood(i, jj);
                lg_adj_old[jj][i] = lg_adj_new[jj][i];
                latent_genotypes_old[jj][i] = latent_genotypes_new[jj][i];
            }
            m_r_accept[i] += 1;
        }

        if (iteration < params.burnin and
            iteration > 15)  // don't start adapting until there are
                             // at least a few samples
        {
            const float acceptanceRate = m_r_accept[i] / float(iteration);
            const float update =
                (acceptanceRate - .23) / std::pow(iteration + 1, .5);
            m_r_var[i] = std::max(m_r_var[i] + update, .01f);
        }
    }
}

/*
 * SALT Sampler approach.
 * https://doi.org/10.1080/00949655.2017.1376063
 */
void Chain::update_p(int iteration)
{
    auto indices = std::vector<int>(genotyping_data.num_loci);
    std::iota(indices.begin(), indices.end(), 0);
    sampler.shuffle_vec(indices);

    for (const auto j : indices)
    {
        const int k = p[j].size();
        int rep = k;
        while (--rep >= 0)
        {
            const int idx = sampler.sample_random_int(0, k - 1);
            p_attempt[j][idx] += 1;

            auto logitPropP = UtilFunctions::logitVec(p[j]);

            const float logitCurr = logitPropP[idx];
            const float logitProp =
                sampler.sample_epsilon(logitCurr, p_prop_var[j][idx]);

            const auto currLogPQ = UtilFunctions::log_pq(logitCurr);
            const auto propLogPQ = UtilFunctions::log_pq(logitProp);

            logitPropP.erase(logitPropP.begin() + idx);

            const float ls =
                propLogPQ.second - UtilFunctions::logitSum(logitPropP);
            logitPropP = UtilFunctions::logitScale(logitPropP, ls);
            logitPropP.insert(logitPropP.begin() + idx, logitProp);

            const float logAdj =
                (currLogPQ.first - propLogPQ.first) +
                (k - 1) * (currLogPQ.second - propLogPQ.second);

            const auto prop_p = UtilFunctions::expitVec(logitPropP);
            // check to make sure the proposed simplex is within a bounded range
            bool sub_threshold_flag = false;
            for (const auto el : prop_p)
            {
                // below a very small threshold that can cause numerical
                // instability
                if (el < 1e-12)
                {
                    sub_threshold_flag = true;
                    break;
                }
            }

            if (sub_threshold_flag)
            {
                break;
            }

            const auto prev_p = p[j];
            p[j] = prop_p;

            for (std::size_t ii = 0; ii < genotyping_data.num_samples; ++ii)
            {
                calculate_genotype_likelihood(ii, j);
            }

            const float new_llik = calc_new_likelihood();
            const float new_prior = calc_new_prior();
            const float new_post = new_llik * temp + new_prior;

            const float acceptanceRatio = new_post - get_posterior() + logAdj;

            if (std::isnan(new_post) or
                sampler.sample_log_mh_acceptance() > acceptanceRatio)
            {
                p[j] = prev_p;
                for (std::size_t ii = 0; ii < genotyping_data.num_samples; ++ii)
                {
                    restore_genotype_likelihood(ii, j);
                }
            }
            else
            {
                llik = new_llik;
                prior = new_prior;
                for (std::size_t ii = 0; ii < genotyping_data.num_samples; ++ii)
                {
                    save_genotype_likelihood(ii, j);
                }
                p_accept[j][idx] += 1;
            }

            // don't start adapting until there are at least a few samples
            if (iteration < params.burnin and p_attempt[j][idx] > 15)
            {
                const float acceptanceRate =
                    (p_accept[j][idx] + 1) / (float(p_attempt[j][idx]) + 1);
                p_prop_var[j][idx] += (acceptanceRate - .23) /
                                      std::pow(p_attempt[j][idx] + 1, .5);
                p_prop_var[j][idx] = std::max(p_prop_var[j][idx], .01f);
            }
        }
    }
}

void Chain::update_eps_pos(int iteration)
{
    auto indices = std::vector<int>(genotyping_data.num_samples);
    std::iota(indices.begin(), indices.end(), 0);
    sampler.shuffle_vec(indices);

    for (const auto i : indices)
    {
        const auto prop_adj = sampler.sample_constrained(
            eps_pos[i], eps_pos_var[i], min_sampled, 1);
        const float prop_eps_pos = std::get<0>(prop_adj);
        const float adj = std::get<1>(prop_adj);

        if (prop_eps_pos < 1 && prop_eps_pos > 1e-32)
        {
            const float prev_eps_pos = eps_pos[i];
            eps_pos[i] = prop_eps_pos;
            calculate_eps_pos_likelihood(i);

            for (std::size_t jj = 0; jj < genotyping_data.num_loci; ++jj)
            {
                calculate_genotype_likelihood(i, jj);
            }

            const float new_llik = calc_new_likelihood();
            const float new_prior = calc_new_prior();
            const float new_post = new_llik * temp + new_prior;

            // Reject
            if (std::isnan(new_post) or sampler.sample_log_mh_acceptance() >
                                            (new_post - get_posterior() + adj))
            {
                eps_pos[i] = prev_eps_pos;
                restore_eps_pos_likelihood(i);
                for (std::size_t jj = 0; jj < genotyping_data.num_loci; ++jj)
                {
                    restore_genotype_likelihood(i, jj);
                }
            }
            else
            {
                llik = new_llik;
                prior = new_prior;
                save_eps_pos_likelihood(i);
                for (std::size_t jj = 0; jj < genotyping_data.num_loci; ++jj)
                {
                    save_genotype_likelihood(i, jj);
                }
                eps_pos_accept[i] += 1;
            }

            if (iteration < params.burnin and
                iteration > 15)  // don't start adapting until there are
                                 // at least a few samples
            {
                const float acceptanceRate =
                    eps_pos_accept[i] / float(iteration);
                const float update =
                    (acceptanceRate - .23) / std::pow(iteration + 1, .5);
                eps_pos_var[i] = std::max(eps_pos_var[i] + update, .01f);
            }
        }
    }
}

void Chain::update_eps_neg(int iteration)
{
    auto indices = std::vector<int>(genotyping_data.num_samples);
    std::iota(indices.begin(), indices.end(), 0);
    sampler.shuffle_vec(indices);

    for (const auto i : indices)
    {
        const auto prop_adj = sampler.sample_constrained(
            eps_neg[i], eps_neg_var[i], min_sampled, 1);
        const float prop_eps_neg = std::get<0>(prop_adj);
        const float adj = std::get<1>(prop_adj);

        if (prop_eps_neg < 1 && prop_eps_neg > 1e-32)
        {
            const float prev_eps_neg = eps_neg[i];
            eps_neg[i] = prop_eps_neg;
            calculate_eps_neg_likelihood(i);

            for (std::size_t jj = 0; jj < genotyping_data.num_loci; ++jj)
            {
                calculate_genotype_likelihood(i, jj);
            }

            const float new_llik = calc_new_likelihood();
            const float new_prior = calc_new_prior();
            const float new_post = new_llik * temp + new_prior;

            // Reject
            if (std::isnan(new_post) or sampler.sample_log_mh_acceptance() >
                                            (new_post - get_posterior() + adj))
            {
                eps_neg[i] = prev_eps_neg;
                restore_eps_neg_likelihood(i);
                for (std::size_t jj = 0; jj < genotyping_data.num_loci; ++jj)
                {
                    restore_genotype_likelihood(i, jj);
                }
            }
            else
            {
                llik = new_llik;
                prior = new_prior;
                save_eps_neg_likelihood(i);
                for (std::size_t jj = 0; jj < genotyping_data.num_loci; ++jj)
                {
                    save_genotype_likelihood(i, jj);
                }
                eps_neg_accept[i] += 1;
            }

            if (iteration < params.burnin and
                iteration > 15)  // don't start adapting until there are
                                 // at least a few samples
            {
                const float acceptanceRate =
                    eps_neg_accept[i] / float(iteration);
                const float update =
                    (acceptanceRate - .23) / std::pow(iteration + 1, .5);
                eps_neg_var[i] = std::max(eps_neg_var[i] + update, .01f);
            }
        }
    }
}

void Chain::update_samples(int iteration)
{
    auto indices = std::vector<int>(genotyping_data.num_samples);
    std::iota(indices.begin(), indices.end(), 0);
    sampler.shuffle_vec(indices);

    for (const auto ii : indices)
    {
        const auto eps_neg_prop_adj = sampler.sample_constrained(
            eps_neg[ii], eps_neg_var[ii], min_sampled, 1);
        const float prop_eps_neg = std::get<0>(eps_neg_prop_adj);
        const float eps_neg_adj = std::get<1>(eps_neg_prop_adj);
        const bool valid_prop_eps_neg =
            prop_eps_neg < 1 && prop_eps_neg > 1e-32;

        const auto eps_pos_prop_adj = sampler.sample_constrained(
            eps_pos[ii], eps_pos_var[ii], min_sampled, 1);
        const float prop_eps_pos = std::get<0>(eps_pos_prop_adj);
        const float eps_pos_adj = std::get<1>(eps_pos_prop_adj);
        const bool valid_prop_eps_pos =
            prop_eps_pos < 1 && prop_eps_pos > 1e-32;

        float prop_r = 0;
        float r_adj = 0;
        bool valid_prop_r = true;
        if (params.allow_relatedness)
        {
            auto r_prop_adj =
                sampler.sample_constrained(r[ii], r_var[ii], min_sampled, .99);
            prop_r = std::get<0>(r_prop_adj);
            r_adj = std::get<1>(r_prop_adj);
            valid_prop_r = prop_r < 1 && prop_r > 1e-32;
        }

        const int prop_m = m[ii] + sampler.sample_coi_delta(2);
        const bool valid_prop_m = prop_m > 0 and prop_m <= params.max_coi;

        if (valid_prop_eps_neg && valid_prop_eps_pos && valid_prop_r &&
            valid_prop_m)
        {
            const float prev_eps_pos = eps_pos[ii];
            eps_pos[ii] = prop_eps_pos;
            calculate_eps_pos_likelihood(ii);

            const float prev_eps_neg = eps_neg[ii];
            eps_neg[ii] = prop_eps_neg;
            calculate_eps_neg_likelihood(ii);

            const float prev_r = r[ii];
            r[ii] = prop_r;
            calculate_relatedness_likelihood(ii);

            const float prev_m = m[ii];
            m[ii] = prop_m;
            calculate_coi_likelihood(ii);

            float adj_ratio = eps_neg_adj + eps_pos_adj + r_adj;

            for (std::size_t jj = 0; jj < genotyping_data.num_loci; ++jj)
            {
                auto lg = sampler.sample_latent_genotype(
                    genotyping_data.get_observed_alleles(jj, ii), m[ii],
                    eps_pos[ii], eps_neg[ii]);
                latent_genotypes_new[jj][ii] = lg.value;
                lg_adj_new[jj][ii] = lg.log_prob;
                adj_ratio = adj_ratio + lg_adj_new[jj][ii] - lg_adj_old[jj][ii];
                calculate_genotype_likelihood(ii, jj);
            }

            const float new_llik = calc_new_likelihood();
            const float new_prior = calc_new_prior();
            const float new_post = new_llik * temp + new_prior;
            const float alpha = sampler.sample_log_mh_acceptance();

            if (!std::isnan(new_post) and
                alpha <= (new_post - get_posterior() + adj_ratio))
            {
                llik = new_llik;
                prior = new_prior;
                save_eps_neg_likelihood(ii);
                save_eps_pos_likelihood(ii);
                save_relatedness_likelihood(ii);
                save_coi_likelihood(ii);
                for (std::size_t jj = 0; jj < genotyping_data.num_loci; ++jj)
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
                restore_relatedness_likelihood(ii);
                restore_coi_likelihood(ii);
                for (std::size_t jj = 0; jj < genotyping_data.num_loci; ++jj)
                {
                    restore_genotype_likelihood(ii, jj);
                    lg_adj_new[jj][ii] = lg_adj_old[jj][ii];
                    latent_genotypes_new[jj][ii] = latent_genotypes_old[jj][ii];
                }
            }
        }
    }
}

void Chain::update_mean_coi(int iteration)
{
    const auto prop_adj = sampler.sample_constrained(mean_coi, mean_coi_var,
                                                     1e-5, params.max_coi);
    const float prop_mean_coi = std::get<0>(prop_adj);
    const float adj = std::get<1>(prop_adj);

    const float prev_mean_coi = mean_coi;
    mean_coi = prop_mean_coi;
    calculate_mean_coi_likelihood();
    for (std::size_t ii = 0; ii < genotyping_data.num_samples; ++ii)
    {
        calculate_coi_likelihood(ii);
    }
    const float new_llik = calc_new_likelihood();
    const float new_prior = calc_new_prior();
    const float new_post = new_llik * temp + new_prior;

    const float alpha = sampler.sample_log_mh_acceptance();
    if (!std::isnan(new_post) and alpha <= (new_post - get_posterior() + adj))
    {
        llik = new_llik;
        prior = new_prior;
        save_mean_coi_likelihood();
        for (std::size_t ii = 0; ii < genotyping_data.num_samples; ++ii)
        {
            save_coi_likelihood(ii);
        }
        mean_coi_accept += 1;
    }
    else
    {
        mean_coi = prev_mean_coi;
        restore_mean_coi_likelihood();
        for (std::size_t ii = 0; ii < genotyping_data.num_samples; ++ii)
        {
            restore_coi_likelihood(ii);
        }
    }

    if (iteration < params.burnin and iteration > 15)
    {
        const float acceptance_rate = mean_coi_accept / (iteration + 1);
        const float update =
            (acceptance_rate - .23) / std::pow(iteration + 1, .5);
        mean_coi_var = std::max(mean_coi_var + update, .1f);
    }
}

float Chain::calc_transmission_process(
    std::vector<int> const &allele_index_vec,
    std::vector<float> const &allele_frequencies, int coi, float relatedness)
{
    // transmission process - prob that after "coi" number of draws, all
    // alleles are drawn at least once conditional on all draws come
    // from the constrained set, where the constrained set is the set of
    // positive alleles in the latent genotype

    float constrained_set_total_prob = 0;
    std::vector<float> res{};
    const int total_alleles = allele_index_vec.size();

    prVec_.clear();
    prVec_.reserve(total_alleles);
    res.reserve(coi - total_alleles + 1);
    for (size_t j = 0; j < allele_index_vec.size(); j++)
    {
        prVec_.push_back(allele_frequencies[allele_index_vec[j]]);
        constrained_set_total_prob += prVec_.back();
    }

    float log_constrained_set_total_prob = std::log(constrained_set_total_prob);

    // normalize the vector
    for (float &k : prVec_)
    {
        k = k / constrained_set_total_prob;
    }

    if (!params.allow_relatedness)
    {
        return std::log(1 - probAnyMissing_(prVec_, coi)) +
               log_constrained_set_total_prob * coi;
    }

    // Only go up to coi - total_alleles because there must be at least
    // total_alleles unrelated strains at this locus
    std::vector<float> pamVec = probAnyMissing_.vectorized(prVec_, coi);
    for (int i = 0; i <= coi - total_alleles; ++i)
    {
        // Calculate the probability of `i` related strains being present, where
        // there are at most coi - total_alleles related strains

        // prob of i related strains
        // const float pr = R::dbinom(i, coi - 1, relatedness, true);
        const float pr = sampler.dbinom(i, coi - 1, relatedness, true);

        const float i_res = std::log(1 - pamVec[coi - i - 1]) +
                            log_constrained_set_total_prob * (coi - i);

        res.push_back(pr + i_res);
    }

    const float out = UtilFunctions::logSumExp(res);
    return out;
}

float Chain::calc_observation_process(std::vector<int> const &allele_index_vec,
                                      std::vector<int> const &obs_genotype,
                                      float epsilon_neg, float epsilon_pos)
{
    float res = 0;
    unsigned int fp = 0;
    unsigned int tp = 0;
    unsigned int fn = 0;
    unsigned int tn = 0;

    unsigned int vec_pointer = 0;
    unsigned int next_allele_index = allele_index_vec[vec_pointer];
    unsigned int total_alleles = allele_index_vec.size();

    unsigned int j = 0;
    for (const auto &e : obs_genotype)
    {
        fp += (e & 1) & !(j == next_allele_index);
        tp += (e & 1) & (j == next_allele_index);
        fn += !(e & 1) & (j == next_allele_index);
        tn += !(e & 1) & !(j == next_allele_index);
        vec_pointer += (j == next_allele_index);

        next_allele_index = -1 + (vec_pointer < total_alleles) *
                                     (allele_index_vec[vec_pointer] + 1);
        ++j;
    }

    const float norm_factor = 1.0f / obs_genotype.size();
    res += std::log(1 - (epsilon_neg * params.max_eps_neg * norm_factor)) * tp;
    res += std::log(epsilon_neg * params.max_eps_neg * norm_factor) * fn;
    res += std::log(1 - (epsilon_pos * params.max_eps_pos * norm_factor)) * tn;
    res += std::log(epsilon_pos * params.max_eps_pos * norm_factor) * fp;

    return res;
};

float Chain::calc_genotype_log_pmf(const std::vector<int> &allele_index_vec,
                                   const std::vector<int> &obs_genotype,
                                   float epsilon_pos, float epsilon_neg,
                                   int coi, float relatedness,
                                   const std::vector<float> &allele_frequencies)
{
    float res = 0.0;
    res += calc_transmission_process(allele_index_vec, allele_frequencies, coi,
                                     relatedness);
    res += calc_observation_process(allele_index_vec, obs_genotype, epsilon_neg,
                                    epsilon_pos);

    return res;
}

float Chain::calc_new_likelihood()
{
    #ifdef HAS_EXECUTION
        return std::reduce(std::execution::unseq, genotyping_llik_new.begin(),
                            genotyping_llik_new.end(), 0.0f);
    #else
        float sum = 0.0f;
        #pragma omp simd reduction(+:sum)
        for (int i = 0; i < genotyping_llik_new.size(); ++i) {
            sum += genotyping_llik_new[i];
        }
        return sum;
    #endif
}

float Chain::calc_new_prior()
{
    #ifdef HAS_EXECUTION
        return std::reduce(std::execution::unseq, eps_neg_prior_new.begin(),
                           eps_neg_prior_new.end(), 0.0f) +
               std::reduce(std::execution::unseq, eps_pos_prior_new.begin(),
                           eps_pos_prior_new.end(), 0.0f) +
               std::reduce(std::execution::unseq, relatedness_prior_new.begin(),
                           relatedness_prior_new.end(), 0.0f) +
               std::reduce(std::execution::unseq, coi_prior_new.begin(),
                           coi_prior_new.end(), 0.0f) +
               mean_coi_hyper_prior_new;
    #else
        float sum = 0.0f;
        #pragma omp simd reduction(+:sum)
        for (int i = 0; i < eps_neg_prior_new.size(); ++i) {
            sum += eps_neg_prior_new[i];
        }
        #pragma omp simd reduction(+:sum)
        for (int i = 0; i < eps_pos_prior_new.size(); ++i) {
            sum += eps_pos_prior_new[i];
        }
        #pragma omp simd reduction(+:sum)
        for (int i = 0; i < relatedness_prior_new.size(); ++i) {
            sum += relatedness_prior_new[i];
        }
        #pragma omp simd reduction(+:sum)
        for (int i = 0; i < coi_prior_new.size(); ++i) {
            sum += coi_prior_new[i];
        }
        sum += mean_coi_hyper_prior_new;
        return sum;
    #endif
}

float Chain::get_llik() { return llik; }
float Chain::get_prior() { return prior; }
float Chain::get_posterior() { return llik * temp + prior; }

void Chain::calculate_genotype_likelihood(int sample_idx, int locus_idx)
{
    int idx = sample_idx * genotyping_data.num_loci + locus_idx;
    if (genotyping_data.is_missing(locus_idx, sample_idx))
    {
        genotyping_llik_new[idx] = 0;
    }
    else
    {
        float res;
        const float transmission_prob = calc_transmission_process(
            latent_genotypes_new[locus_idx][sample_idx], p[locus_idx],
            m[sample_idx], r[sample_idx]);
        const float obs_prob = calc_observation_process(
            latent_genotypes_new[locus_idx][sample_idx],
            genotyping_data.get_observed_alleles(locus_idx, sample_idx),
            eps_neg[sample_idx], eps_pos[sample_idx]);
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

void Chain::calculate_relatedness_likelihood(int sample_idx)
{
    relatedness_prior_new[sample_idx] = sampler.get_relatedness_log_prior(
        r[sample_idx], params.r_alpha, params.r_beta);
}

void Chain::calculate_coi_likelihood(int sample_idx)
{
    coi_prior_new[sample_idx] =
        sampler.get_coi_log_prior(m[sample_idx], mean_coi);
}

void Chain::calculate_mean_coi_likelihood()
{
    mean_coi_hyper_prior_new = sampler.get_coi_mean_log_hyper_prior(
        mean_coi, params.mean_coi_shape, params.mean_coi_scale);
}

void Chain::initialize_likelihood()
{
    const int num_samples = genotyping_data.num_samples;
    const int num_loci = genotyping_data.num_loci;

    genotyping_llik_new.resize(num_samples * num_loci);
    genotyping_llik_old.resize(num_samples * num_loci);

    eps_neg_prior_new.resize(num_samples);
    eps_neg_prior_old.resize(num_samples);
    eps_pos_prior_new.resize(num_samples);
    eps_pos_prior_old.resize(num_samples);
    relatedness_prior_new.resize(num_samples);
    relatedness_prior_old.resize(num_samples);
    coi_prior_new.resize(num_samples);
    coi_prior_old.resize(num_samples);

    for (std::size_t ii = 0; ii < genotyping_data.num_samples; ++ii)
    {
        const int row_idx = ii * num_loci;
        for (std::size_t jj = 0; jj < genotyping_data.num_loci; ++jj)
        {
            const int idx = row_idx + jj;
            calculate_genotype_likelihood(ii, jj);
            genotyping_llik_old[idx] = genotyping_llik_new[idx];
        }
    }

    for (size_t ii = 0; ii < genotyping_data.num_samples; ii++)
    {
        const float eps_neg_prior = sampler.get_epsilon_log_prior(
            eps_neg[ii], params.eps_neg_alpha, params.eps_neg_beta);
        const float eps_pos_prior = sampler.get_epsilon_log_prior(
            eps_pos[ii], params.eps_pos_alpha, params.eps_pos_beta);

        const float relatedness_prior = sampler.get_relatedness_log_prior(
            r[ii], params.r_alpha, params.r_beta);

        const float coi_prior = sampler.get_coi_log_prior(m[ii], mean_coi);

        eps_neg_prior_old[ii] = eps_neg_prior;
        eps_neg_prior_new[ii] = eps_neg_prior;
        eps_pos_prior_old[ii] = eps_pos_prior;
        eps_pos_prior_new[ii] = eps_pos_prior;
        relatedness_prior_old[ii] = relatedness_prior;
        relatedness_prior_new[ii] = relatedness_prior;
        coi_prior_old[ii] = coi_prior;
        coi_prior_new[ii] = coi_prior;
    }

    mean_coi_hyper_prior_old = sampler.get_coi_mean_log_hyper_prior(
        mean_coi, params.mean_coi_shape, params.mean_coi_scale);
    mean_coi_hyper_prior_new = sampler.get_coi_mean_log_hyper_prior(
        mean_coi, params.mean_coi_shape, params.mean_coi_scale);

    llik = calc_new_likelihood();
    prior = calc_new_prior();
}

void Chain::save_genotype_likelihood(int sample_idx, int locus_idx)
{
    const int idx = sample_idx * genotyping_data.num_loci + locus_idx;

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

void Chain::save_relatedness_likelihood(int sample_idx)
{
    relatedness_prior_old[sample_idx] = relatedness_prior_new[sample_idx];
}

void Chain::save_coi_likelihood(int sample_idx)
{
    coi_prior_old[sample_idx] = coi_prior_new[sample_idx];
}

void Chain::save_mean_coi_likelihood()
{
    mean_coi_hyper_prior_old = mean_coi_hyper_prior_new;
}

void Chain::restore_genotype_likelihood(int sample_idx, int locus_idx)
{
    const int idx = sample_idx * genotyping_data.num_loci + locus_idx;
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

void Chain::restore_relatedness_likelihood(int sample_idx)
{
    relatedness_prior_new[sample_idx] = relatedness_prior_old[sample_idx];
}

void Chain::restore_coi_likelihood(int sample_idx)
{
    coi_prior_new[sample_idx] = coi_prior_old[sample_idx];
}

void Chain::restore_mean_coi_likelihood()
{
    mean_coi_hyper_prior_new = mean_coi_hyper_prior_old;
}

void Chain::set_llik(float llik) { this->llik = llik; }

void Chain::set_temp(float temp) { this->temp = temp; }

float Chain::get_temp() { return this->temp; }

Chain::Chain(GenotypingData genotyping_data, Parameters params, float temp)
    : genotyping_data(genotyping_data), params(params), sampler()

{
    p_prop_var = std::vector<std::vector<float>>(genotyping_data.num_loci);
    p_accept = std::vector<std::vector<int>>(genotyping_data.num_loci);
    mean_coi_accept = 0;
    mean_coi_var = 1;
    mean_coi =
        sampler.sample_mean_coi(params.mean_coi_shape, params.mean_coi_scale) +
        1;

    llik = std::numeric_limits<float>::lowest();
    this->temp = temp;

    initialize_m();
    initialize_eps_neg();
    initialize_eps_pos();
    initialize_r();
    initialize_latent_genotypes();
    initialize_p();
    initialize_likelihood();
}
