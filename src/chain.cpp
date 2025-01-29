#include "chain.h"

#include "mcmc_utils.h"
#include "sampler.h"

#include <cmath>
#include <algorithm>

#if defined(__cpp_lib_execution) && __cpp_lib_execution >= 201603
#include <execution>
#include "RcppParallel.h"
#define HAS_EXECUTION 1
#endif

#include <limits>
#include <map>
#include <numeric>
#include <span>

constexpr float min_sampled = std::numeric_limits<float>::min();

void Chain::initialize_latent_genotypes()
{
    latent_genotypes_new.clear();
    latent_genotypes_new.resize(std::array{genotyping_data.num_samples, genotyping_data.num_loci}, genotyping_data.num_alleles, -1);
    latent_genotypes_old.clear();
    latent_genotypes_old.resize(std::array{genotyping_data.num_samples, genotyping_data.num_loci}, genotyping_data.num_alleles, -1);
    lg_adj_old.clear();
    lg_adj_old.resize(std::array{genotyping_data.num_samples, genotyping_data.num_loci});
    lg_adj_new.clear();
    lg_adj_new.resize(std::array{genotyping_data.num_samples, genotyping_data.num_loci});
    for (std::size_t locus_idx = 0; locus_idx < genotyping_data.num_loci; ++locus_idx)
    {
        for (std::size_t sample_idx = 0; sample_idx < genotyping_data.num_samples; ++sample_idx)
        {
            const auto obs_alleles = genotyping_data.get_observed_alleles(locus_idx, sample_idx);
            const auto lg = sampler.sample_latent_genotype(obs_alleles, m[sample_idx], .2, .2);

            latent_genotypes_new.inner_fill({sample_idx, locus_idx}, lg.value);
            lg_adj_new.at({sample_idx, locus_idx}) = lg.log_prob;

            latent_genotypes_old.inner_fill({sample_idx, locus_idx}, lg.value);
            lg_adj_old.at({sample_idx, locus_idx}) = lg.log_prob;
        }
    }
}

// Initialize P with random allele frequencies
void Chain::initialize_p()
{
    p.clear();
    p.resize({genotyping_data.num_loci}, genotyping_data.num_alleles);
    p_prop_var.clear();
    p_accept.clear();
    p_attempt.clear();

    for (size_t i = 0; i < genotyping_data.num_loci; i++)
    {
        p.inner_fill({i}, sampler.sample_allele_frequencies(
            std::vector<float>(genotyping_data.num_alleles[i], 1), 10));
    }

    p_prop_var.resize({genotyping_data.num_loci}, genotyping_data.num_alleles);
    p_accept.resize({genotyping_data.num_loci}, genotyping_data.num_alleles);
    p_attempt.resize({genotyping_data.num_loci}, genotyping_data.num_alleles);
};

void Chain::initialize_m()
{
    m.clear();
    m_accept.clear();
    sample_accept.clear();

    m.reserve(genotyping_data.num_samples);
    for (const auto coi : genotyping_data.observed_coi)
    {
        std::size_t m_coi = std::clamp(coi + sampler.sample_coi_delta(3), std::size_t(1), params.max_coi);
        m.push_back(m_coi);
    }
    m_accept.resize(genotyping_data.num_samples, 0);
    sample_accept.resize(genotyping_data.num_samples, 0);
}

void Chain::initialize_eps_neg()
{
    eps_neg.clear();
    eps_neg_accept.clear();
    eps_neg_var.clear();

    eps_neg.reserve(genotyping_data.num_samples);
    for (size_t i = 0; i < genotyping_data.num_samples; ++i)
    {
        eps_neg.push_back(sampler.sample_unif() * .1);
    }
    eps_neg_accept.resize(genotyping_data.num_samples, 0);
    eps_neg_var.resize(genotyping_data.num_samples, 1);
}

void Chain::initialize_eps_pos()
{
    eps_pos.clear();
    eps_pos_accept.clear();
    eps_pos_var.clear();

    eps_pos.reserve(genotyping_data.num_samples);
    for (size_t i = 0; i < genotyping_data.num_samples; ++i)
    {
        eps_pos.push_back(sampler.sample_unif() * .1);
    }
    eps_pos_accept.resize(genotyping_data.num_samples, 0);
    eps_pos_var.resize(genotyping_data.num_samples, 1);
}

void Chain::initialize_r()
{
    r.clear();
    r_accept.clear();
    r_var.clear();
    m_r_accept.clear();
    m_r_var.clear();

    r.reserve(genotyping_data.num_samples);
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

void Chain::initialize_parameters()
{
    initialize_m();
    initialize_eps_neg();
    initialize_eps_pos();
    initialize_r();
    initialize_latent_genotypes();
    initialize_p();
    initialize_likelihood();
    initialize_likelihood();
}

void Chain::update_m(int iteration)
{
    auto indices = std::vector<std::size_t>(genotyping_data.num_samples);
    std::iota(indices.begin(), indices.end(), 0);
    sampler.shuffle_vec(indices);

    for (const auto sample_idx : indices)
    {
        const std::size_t prop_m = m[sample_idx] + sampler.sample_coi_delta(2);
        if (prop_m > 0 and prop_m <= params.max_coi)
        {
            const float prev_m = m[sample_idx];
            m[sample_idx] = prop_m;
            calculate_coi_likelihood(sample_idx);

            float adj_ratio = 0;

            for (std::size_t locus_idx = 0; locus_idx < genotyping_data.num_loci; ++locus_idx)
            {
                const auto& lg = sampler.sample_latent_genotype(
                    genotyping_data.get_observed_alleles(locus_idx, sample_idx), m[sample_idx],
                    eps_pos[sample_idx], eps_neg[sample_idx]);
                latent_genotypes_new.inner_fill({sample_idx, locus_idx}, -1);
                latent_genotypes_new.inner_fill({sample_idx, locus_idx}, lg.value);
                lg_adj_new.at({sample_idx, locus_idx}) = lg.log_prob;
                adj_ratio = adj_ratio + lg_adj_new.at({sample_idx, locus_idx}) - lg_adj_old.at({sample_idx, locus_idx});
                calculate_genotype_likelihood(sample_idx, locus_idx);
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
                save_coi_likelihood(sample_idx);
                for (std::size_t locus_idx = 0; locus_idx < genotyping_data.num_loci; ++locus_idx)
                {
                    save_genotype_likelihood(sample_idx, locus_idx);
                    lg_adj_old.at({sample_idx, locus_idx}) = lg_adj_new.at({sample_idx, locus_idx});
                    auto [begin, end] = latent_genotypes_new.inner_iterators({sample_idx, locus_idx});
                    std::copy(begin, end, latent_genotypes_old.inner_begin({sample_idx, locus_idx}));
                }
                ++m_accept[sample_idx];
            }
            else
            {
                m[sample_idx] = prev_m;
                restore_coi_likelihood(sample_idx);
                for (std::size_t locus_idx = 0; locus_idx < genotyping_data.num_loci; ++locus_idx)
                {
                    restore_genotype_likelihood(sample_idx, locus_idx);
                    lg_adj_new.at({sample_idx, locus_idx}) = lg_adj_old.at({sample_idx, locus_idx});
                    auto [begin, end] = latent_genotypes_old.inner_iterators({sample_idx, locus_idx});
                    std::copy(begin, end, latent_genotypes_new.inner_begin({sample_idx, locus_idx}));
                }
            }
        }
    }
}

void Chain::update_eff_coi(int iteration)
{
    auto indices = std::vector<std::size_t>(genotyping_data.num_samples);
    std::iota(indices.begin(), indices.end(), 0);
    sampler.shuffle_vec(indices);

    for (const auto sample_idx : indices)
    {
        const float curr_eff_coi = (m[sample_idx] - 1) * (1.0f - r[sample_idx]) + 1.0f;
        const auto prop_adj = sampler.sample_constrained(
            curr_eff_coi, m_r_var[sample_idx], 1, params.max_coi);

        const float prop_eff_coi = std::get<0>(prop_adj);
        const float adj = std::get<1>(prop_adj);

        const std::size_t prop_m = m[sample_idx] + sampler.sample_coi_delta(2);
        const float prop_r = 1.0f - (prop_eff_coi - 1.0f) / (prop_m - 1.0f);

        if (prop_m <= 0 || prop_m > params.max_coi || prop_r > 1.0f || prop_r < 1e-32) {
            continue;
        }

        const int prev_m = m[sample_idx];
        const float prev_r = r[sample_idx];
        m[sample_idx] = prop_m;
        r[sample_idx] = prop_r;
        calculate_coi_likelihood(sample_idx);
        calculate_relatedness_likelihood(sample_idx);

        float adj_ratio = adj;

        for (std::size_t locus_idx = 0; locus_idx < genotyping_data.num_loci; ++locus_idx)
        {
            const auto &lg = sampler.sample_latent_genotype(
                genotyping_data.get_observed_alleles(locus_idx, sample_idx), m[sample_idx], eps_pos[sample_idx],
                eps_neg[sample_idx]);
            latent_genotypes_new.inner_fill({sample_idx, locus_idx}, -1);
            latent_genotypes_new.inner_fill({sample_idx, locus_idx}, lg.value);
            lg_adj_new.at({sample_idx, locus_idx}) = lg.log_prob;
            adj_ratio = adj_ratio + lg_adj_new.at({sample_idx, locus_idx}) - lg_adj_old.at({sample_idx, locus_idx});
            calculate_genotype_likelihood(sample_idx, locus_idx);
        }

        const float new_llik = calc_new_likelihood();
        const float new_prior = calc_new_prior();
        const float new_post = new_llik * temp + new_prior;
        const double mh_ratio = new_post - get_posterior() + adj_ratio;
        const double alpha = sampler.sample_log_mh_acceptance();

        // Reject
        if (std::isnan(new_post) or std::isnan(mh_ratio) or alpha > mh_ratio)
        {
            m[sample_idx] = prev_m;
            r[sample_idx] = prev_r;
            restore_coi_likelihood(sample_idx);
            restore_relatedness_likelihood(sample_idx);
            for (std::size_t locus_idx = 0; locus_idx < genotyping_data.num_loci; ++locus_idx)
            {
                restore_genotype_likelihood(sample_idx, locus_idx);
                lg_adj_new.at({sample_idx, locus_idx}) = lg_adj_old.at({sample_idx, locus_idx});
                const auto [begin, end] = latent_genotypes_old.inner_iterators({sample_idx, locus_idx});
                std::copy(begin, end, latent_genotypes_new.inner_begin({sample_idx, locus_idx}));
            }
        }
        else
        {
            llik = new_llik;
            prior = new_prior;
            save_relatedness_likelihood(sample_idx);
            save_coi_likelihood(sample_idx);
            for (std::size_t locus_idx = 0; locus_idx < genotyping_data.num_loci; ++locus_idx)
            {
                save_genotype_likelihood(sample_idx, locus_idx);
                lg_adj_old.at({sample_idx, locus_idx}) = lg_adj_new.at({sample_idx, locus_idx});
                const auto [begin, end] = latent_genotypes_new.inner_iterators({sample_idx, locus_idx});
                std::copy(begin, end, latent_genotypes_old.inner_begin({sample_idx, locus_idx}));
            }
            ++m_r_accept[sample_idx];
        }

        if (iteration < params.burnin and
            iteration > 15)  // don't start adapting until there are
                             // at least a few samples
        {
            const float acceptanceRate = m_r_accept[sample_idx] / float(iteration);
            const float update =
                (acceptanceRate - .23) / std::pow(iteration + 1, .5);
            m_r_var[sample_idx] = std::max(m_r_var[sample_idx] + update, .01f);
        }
    }
}

void Chain::update_r(int iteration)
{
    auto indices = std::vector<int>(genotyping_data.num_samples);
    std::iota(indices.begin(), indices.end(), 0);
    sampler.shuffle_vec(indices);

    for (const auto sample_idx : indices)
    {
        const auto prop_adj =
            sampler.sample_constrained(r[sample_idx], r_var[sample_idx], min_sampled, .99);
        const float prop_r = std::get<0>(prop_adj);
        const float adj = std::get<1>(prop_adj);

        const float prev_r = r[sample_idx];
        r[sample_idx] = prop_r;
        calculate_relatedness_likelihood(sample_idx);

        for (std::size_t locus_idx = 0; locus_idx < genotyping_data.num_loci; ++locus_idx)
        {
            calculate_genotype_likelihood(sample_idx, locus_idx);
        }

        const float new_llik = calc_new_likelihood();
        const float new_prior = calc_new_prior();
        const float new_post = new_llik * temp + new_prior;

        // Reject
        if (std::isnan(new_post) or sampler.sample_log_mh_acceptance() >
                                        (new_post - get_posterior() + adj))
        {
            r[sample_idx] = prev_r;
            restore_relatedness_likelihood(sample_idx);
            for (std::size_t locus_idx = 0; locus_idx < genotyping_data.num_loci; ++locus_idx)
            {
                restore_genotype_likelihood(sample_idx, locus_idx);
            }
        }
        else
        {
            llik = new_llik;
            prior = new_prior;
            save_relatedness_likelihood(sample_idx);
            for (std::size_t locus_idx = 0; locus_idx < genotyping_data.num_loci; ++locus_idx)
            {
                save_genotype_likelihood(sample_idx, locus_idx);
            }
            ++r_accept[sample_idx];
        }

        if (iteration < params.burnin and
            iteration > 15)  // don't start adapting until there are
                             // at least a few samples
        {
            const float acceptanceRate = r_accept[sample_idx] / float(iteration);
            const float update =
                (acceptanceRate - .23) / std::pow(iteration + 1, .5);
            r_var[sample_idx] = std::max(r_var[sample_idx] + update, .0001f);
        }
    }
}

void Chain::update_m_r(int iteration)
{
    auto indices = std::vector<std::size_t>(genotyping_data.num_samples);
    std::iota(indices.begin(), indices.end(), 0);
    sampler.shuffle_vec(indices);

    for (const auto sample_idx : indices)
    {
        const auto prop_adj =
            sampler.sample_constrained(r[sample_idx], m_r_var[sample_idx], min_sampled, .99);
        const float prop_r = std::get<0>(prop_adj);
        const float adj = std::get<1>(prop_adj);

        const std::size_t prop_m = m[sample_idx] + sampler.sample_coi_delta(1);

        if (prop_m <= 0 or prop_m > params.max_coi) return;

        const float prev_m = m[sample_idx];
        const float prev_r = r[sample_idx];
        r[sample_idx] = prop_r;
        calculate_relatedness_likelihood(sample_idx);
        m[sample_idx] = prop_m;
        calculate_coi_likelihood(sample_idx);

        float adj_ratio = adj;

        for (std::size_t locus_idx = 0; locus_idx < genotyping_data.num_loci; ++locus_idx)
        {
            const auto &lg = sampler.sample_latent_genotype(
                genotyping_data.get_observed_alleles(locus_idx, sample_idx), m[sample_idx], eps_pos[sample_idx],
                eps_neg[sample_idx]);
            latent_genotypes_new.inner_fill({sample_idx, locus_idx}, -1);
            latent_genotypes_new.inner_fill({sample_idx, locus_idx}, lg.value);
            lg_adj_new.at({sample_idx, locus_idx}) = lg.log_prob;
            adj_ratio = adj_ratio + lg_adj_new.at({sample_idx, locus_idx}) - lg_adj_old.at({sample_idx, locus_idx});
            calculate_genotype_likelihood(sample_idx, locus_idx);
        }

        const float new_llik = calc_new_likelihood();
        const float new_prior = calc_new_prior();
        const float new_post = new_llik * temp + new_prior;

        // Reject
        if (std::isnan(new_post) or
            sampler.sample_log_mh_acceptance() >
                (new_post - get_posterior() + adj_ratio))
        {
            r[sample_idx] = prev_r;
            restore_relatedness_likelihood(sample_idx);
            m[sample_idx] = prev_m;
            restore_coi_likelihood(sample_idx);
            for (std::size_t locus_idx = 0; locus_idx < genotyping_data.num_loci; ++locus_idx)
            {
                restore_genotype_likelihood(sample_idx, locus_idx);
                lg_adj_new.at({sample_idx, locus_idx}) = lg_adj_old.at({sample_idx, locus_idx});
                auto [begin, end] = latent_genotypes_old.inner_iterators({sample_idx, locus_idx});
                std::copy(begin, end, latent_genotypes_new.inner_begin({sample_idx, locus_idx}));
            }
        }
        else
        {
            llik = new_llik;
            prior = new_prior;
            save_relatedness_likelihood(sample_idx);
            save_coi_likelihood(sample_idx);
            for (std::size_t locus_idx = 0; locus_idx < genotyping_data.num_loci; ++locus_idx)
            {
                save_genotype_likelihood(sample_idx, locus_idx);
                lg_adj_old.at({sample_idx, locus_idx}) = lg_adj_new.at({sample_idx, locus_idx});
                auto [begin, end] = latent_genotypes_new.inner_iterators({sample_idx, locus_idx});
                std::copy(begin, end, latent_genotypes_old.inner_begin({sample_idx, locus_idx}));
            }
            ++m_r_accept[sample_idx];
        }

        if (iteration < params.burnin and
            iteration > 15)  // don't start adapting until there are
                             // at least a few samples
        {
            const float acceptanceRate = m_r_accept[sample_idx] / float(iteration);
            const float update =
                (acceptanceRate - .23) / std::pow(iteration + 1, .5);
            m_r_var[sample_idx] = std::max(m_r_var[sample_idx] + update, .0001f);
        }
    }
}

/*
 * SALT Sampler approach.
 * https://doi.org/10.1080/00949655.2017.1376063
 */
void Chain::update_p(int iteration)
{
    auto indices = std::vector<std::size_t>(genotyping_data.num_loci);
    std::iota(indices.begin(), indices.end(), 0);
    sampler.shuffle_vec(indices);

    for (const std::size_t locus_idx : indices)
    {
        const std::size_t num_alleles = p.ragged_dimensions(locus_idx);
        std::size_t rep = num_alleles; // +1 for post order decrement, total of attem

        const auto [begin, end] = p.inner_iterators({locus_idx});
        while (rep-- > 0)
        {
            const size_t allele_idx = sampler.sample_random_int(0, num_alleles - 1);

            ++p_attempt.at({locus_idx, allele_idx});

            auto logitPropP = UtilFunctions::logitVec(begin, end);

            const float logitCurr = logitPropP[allele_idx];
            const float logitProp =
                sampler.sample_epsilon(logitCurr, p_prop_var.at({locus_idx, allele_idx}));

            const auto currLogPQ = UtilFunctions::log_pq(logitCurr);
            const auto propLogPQ = UtilFunctions::log_pq(logitProp);

            logitPropP.erase(logitPropP.begin() + allele_idx);

            const float ls =
                propLogPQ.second - UtilFunctions::logitSum(logitPropP);
            logitPropP = UtilFunctions::logitScale(logitPropP, ls);
            logitPropP.insert(logitPropP.begin() + allele_idx, logitProp);

            const float logAdj =
                (currLogPQ.first - propLogPQ.first) +
                (num_alleles - 1) * (currLogPQ.second - propLogPQ.second);

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

            const auto prev_p = std::vector<float>(begin, end);
            p.inner_fill({locus_idx}, prop_p);

            for (std::size_t sample_idx = 0; sample_idx < genotyping_data.num_samples; ++sample_idx)
            {
                calculate_genotype_likelihood(sample_idx, locus_idx);
            }

            const float new_llik = calc_new_likelihood();
            const float new_prior = calc_new_prior();
            const float new_post = new_llik * temp + new_prior;

            const float acceptanceRatio = new_post - get_posterior() + logAdj;

            if (std::isnan(new_post) or
                sampler.sample_log_mh_acceptance() > acceptanceRatio)
            {
                p.inner_fill({locus_idx}, prev_p);
                for (std::size_t sample_idx = 0; sample_idx < genotyping_data.num_samples; ++sample_idx)
                {
                    restore_genotype_likelihood(sample_idx, locus_idx);
                }
            }
            else
            {
                llik = new_llik;
                prior = new_prior;
                for (std::size_t sample_idx = 0; sample_idx < genotyping_data.num_samples; ++sample_idx)
                {
                    save_genotype_likelihood(sample_idx, locus_idx);
                }
                ++p_accept.at({locus_idx, allele_idx});
            }

            // don't start adapting until there are at least a few samples
            if (iteration < params.burnin and p_attempt.at({locus_idx, allele_idx}) > 15)
            {
                const float acceptanceRate =
                    (p_accept.at({locus_idx, allele_idx}) + 1) / (float(p_attempt.at({locus_idx, allele_idx})) + 1);
                p_prop_var.at({locus_idx, allele_idx}) += (acceptanceRate - .23) /
                                      std::pow(p_attempt.at({locus_idx, allele_idx}) + 1, .5);
                p_prop_var.at({locus_idx, allele_idx}) = std::max(p_prop_var.at({locus_idx, allele_idx}), .01f);
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
                ++eps_pos_accept[i];
            }

            if (iteration < params.burnin and
                iteration > 15)  // don't start adapting until there are
                                 // at least a few samples
            {
                const float acceptanceRate =
                    eps_pos_accept[i] / float(iteration);
                const float update =
                    (acceptanceRate - .23) / std::pow(iteration + 1, .5);
                eps_pos_var[i] = std::max(eps_pos_var[i] + update, .0001f);
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
                ++eps_neg_accept[i];
            }

            if (iteration < params.burnin and
                iteration > 15)  // don't start adapting until there are
                                 // at least a few samples
            {
                const float acceptanceRate =
                    eps_neg_accept[i] / float(iteration);
                const float update =
                    (acceptanceRate - .23) / std::pow(iteration + 1, .5);
                eps_neg_var[i] = std::max(eps_neg_var[i] + update, .0001f);
            }
        }
    }
}

void Chain::update_samples(int iteration)
{
    auto indices = std::vector<std::size_t>(genotyping_data.num_samples);
    std::iota(indices.begin(), indices.end(), 0);
    sampler.shuffle_vec(indices);

    for (const auto sample_idx : indices)
    {
        const auto eps_neg_prop_adj = sampler.sample_constrained(
            eps_neg[sample_idx], eps_neg_var[sample_idx], min_sampled, 1);
        const float prop_eps_neg = std::get<0>(eps_neg_prop_adj);
        const float eps_neg_adj = std::get<1>(eps_neg_prop_adj);
        const bool valid_prop_eps_neg =
            prop_eps_neg < 1 && prop_eps_neg > 1e-32;

        const auto eps_pos_prop_adj = sampler.sample_constrained(
            eps_pos[sample_idx], eps_pos_var[sample_idx], min_sampled, 1);
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
                sampler.sample_constrained(r[sample_idx], r_var[sample_idx], min_sampled, .99);
            prop_r = std::get<0>(r_prop_adj);
            r_adj = std::get<1>(r_prop_adj);
            valid_prop_r = prop_r < 1 && prop_r > 1e-32;
        }

        const std::size_t prop_m = m[sample_idx] + sampler.sample_coi_delta(2);
        const bool valid_prop_m = prop_m > 0 && prop_m <= params.max_coi;

        if (valid_prop_eps_neg && valid_prop_eps_pos && valid_prop_r &&
            valid_prop_m)
        {
            const float prev_eps_pos = eps_pos[sample_idx];
            eps_pos[sample_idx] = prop_eps_pos;
            calculate_eps_pos_likelihood(sample_idx);

            const float prev_eps_neg = eps_neg[sample_idx];
            eps_neg[sample_idx] = prop_eps_neg;
            calculate_eps_neg_likelihood(sample_idx);

            const float prev_r = r[sample_idx];
            r[sample_idx] = prop_r;
            calculate_relatedness_likelihood(sample_idx);

            const float prev_m = m[sample_idx];
            m[sample_idx] = prop_m;
            calculate_coi_likelihood(sample_idx);

            float adj_ratio = eps_neg_adj + eps_pos_adj + r_adj;

            for (std::size_t locus_idx = 0; locus_idx < genotyping_data.num_loci; ++locus_idx)
            {
                auto lg = sampler.sample_latent_genotype(
                    genotyping_data.get_observed_alleles(locus_idx, sample_idx), m[sample_idx],
                    eps_pos[sample_idx], eps_neg[sample_idx]);
                latent_genotypes_new.inner_fill({sample_idx, locus_idx}, -1);
                latent_genotypes_new.inner_fill({sample_idx, locus_idx}, lg.value);
                lg_adj_new.at({sample_idx, locus_idx}) = lg.log_prob;
                adj_ratio = adj_ratio + lg_adj_new.at({sample_idx, locus_idx}) - lg_adj_old.at({sample_idx, locus_idx});
                calculate_genotype_likelihood(sample_idx, locus_idx);
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
                save_eps_neg_likelihood(sample_idx);
                save_eps_pos_likelihood(sample_idx);
                save_relatedness_likelihood(sample_idx);
                save_coi_likelihood(sample_idx);
                for (std::size_t locus_idx = 0; locus_idx < genotyping_data.num_loci; ++locus_idx)
                {
                    save_genotype_likelihood(sample_idx, locus_idx);
                    lg_adj_old.at({sample_idx, locus_idx}) = lg_adj_new.at({sample_idx, locus_idx});
                    auto [begin, end] = latent_genotypes_new.inner_iterators({sample_idx, locus_idx});
                    std::copy(begin, end, latent_genotypes_old.inner_begin({sample_idx, locus_idx}));
                }
                ++sample_accept[sample_idx];
            }
            else
            {
                m[sample_idx] = prev_m;
                eps_pos[sample_idx] = prev_eps_pos;
                eps_neg[sample_idx] = prev_eps_neg;
                r[sample_idx] = prev_r;
                restore_eps_neg_likelihood(sample_idx);
                restore_eps_pos_likelihood(sample_idx);
                restore_relatedness_likelihood(sample_idx);
                restore_coi_likelihood(sample_idx);
                for (std::size_t locus_idx = 0; locus_idx < genotyping_data.num_loci; ++locus_idx)
                {
                    restore_genotype_likelihood(sample_idx, locus_idx);
                    lg_adj_new.at({sample_idx, locus_idx}) = lg_adj_old.at({sample_idx, locus_idx});
                    auto [begin, end] = latent_genotypes_old.inner_iterators({sample_idx, locus_idx});
                    std::copy(begin, end, latent_genotypes_new.inner_begin({sample_idx, locus_idx}));
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
        ++mean_coi_accept;
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
        mean_coi_var = std::max(mean_coi_var + update, .0001f);
    }
}

float Chain::calc_transmission_process(
    std::span<int const> full_allele_index_vec,
    std::span<float const> allele_frequencies, int coi, float relatedness)
{
    // transmission process - prob that after "coi" number of draws, all
    // alleles are drawn at least once conditional on all draws come
    // from the constrained set, where the constrained set is the set of
    // positive alleles in the latent genotype

    // get span up to first instance of -1
    auto first_instance = std::find(full_allele_index_vec.begin(), full_allele_index_vec.end(), -1);
    std::span<int const> allele_index_vec = full_allele_index_vec.subspan(0, first_instance - full_allele_index_vec.begin());

    float constrained_set_total_prob = 0;
    std::vector<float> res{};
    const size_t total_alleles = allele_index_vec.size();
    std::vector<float> prVec_{};

    res.reserve(coi - total_alleles + 1);
    prVec_.reserve(total_alleles);
    for (const auto &e : allele_index_vec)
    {
        prVec_.push_back(allele_frequencies[e]);
        constrained_set_total_prob += prVec_.back();
    }

    const float log_constrained_set_total_prob = std::log(constrained_set_total_prob);

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
    const std::vector<float> pamVec = probAnyMissing_.vectorized(prVec_, coi);
    for (size_t i = 0; i <= coi - total_alleles; ++i)
    {
        // Calculate the probability of `i` related strains being present, where
        // there are at most coi - total_alleles related strains

        // prob of i related strains
        const float pr = sampler.dbinom(i, coi - 1, relatedness);

        const float i_res = std::log(1 - pamVec[coi - i - 1]) +
                            log_constrained_set_total_prob * (coi - i);

        res.push_back(pr + i_res);
    }

    return UtilFunctions::logSumExp(res);
}

float Chain::calc_observation_process(std::span<int const> full_allele_index_vec,
                                      std::span<int const> obs_genotype,
                                      float epsilon_neg, float epsilon_pos)
{
    float res = 0;
    unsigned int fp = 0;
    unsigned int tp = 0;
    unsigned int fn = 0;
    unsigned int tn = 0;

    auto first_instance = std::find(full_allele_index_vec.begin(), full_allele_index_vec.end(), -1);
    std::span<int const> allele_index_vec = full_allele_index_vec.subspan(0, first_instance - full_allele_index_vec.begin());

    unsigned int total_alleles = allele_index_vec.size(); 
    unsigned int vec_pointer = 0;
    int j = 0;
    for (const auto &e : obs_genotype)
    {
        auto next_allele_index = -1;
        if (vec_pointer < total_alleles)
        {
            next_allele_index = allele_index_vec[vec_pointer];
        }

        const auto mask = (j == next_allele_index);
        fp += (e & 1) & !mask;
        tp += (e & 1) & mask;
        fn += !(e & 1) & mask;
        tn += !(e & 1) & !mask;
        vec_pointer += mask;
        
        ++j;
    }

    const float norm_factor = 1.0f / obs_genotype.size();
    res += std::log(1 - (epsilon_neg * params.max_eps_neg * norm_factor)) * tp;
    res += std::log(epsilon_neg * params.max_eps_neg * norm_factor) * fn;
    res += std::log(1 - (epsilon_pos * params.max_eps_pos * norm_factor)) * tn;
    res += std::log(epsilon_pos * params.max_eps_pos * norm_factor) * fp;

    return res;
};

float Chain::calc_genotype_log_pmf(std::span<int const> allele_index_vec,
                                   std::span<int const> obs_genotype,
                                   float epsilon_pos, float epsilon_neg,
                                   int coi, float relatedness,
                                   std::span<float const> allele_frequencies)
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

float Chain::get_llik(int sample) { 
    int idx = sample * genotyping_data.num_loci;

    #ifdef HAS_EXECUTION
        return std::reduce(std::execution::unseq, genotyping_llik_new.begin() + idx,
                           genotyping_llik_new.begin() + idx + genotyping_data.num_loci);
    #else
        float llik = 0.0f;
        #pragma omp simd reduction(+:llik)
        for (int i = 0; i < genotyping_data.num_loci; ++i) {
            llik += genotyping_llik_new[idx + i];
        }
        return llik;
    #endif
}
float Chain::get_prior(int sample) { 
    float prior = 0.0f;
    if (params.allow_relatedness) {
        prior += relatedness_prior_new[sample];
    }
    prior += coi_prior_new[sample];
    prior += eps_neg_prior_new[sample];
    prior += eps_pos_prior_new[sample];
    return prior;
}

float Chain::get_posterior(int sample) { 
    float posterior = get_llik(sample) * temp + get_prior(sample);
    return posterior;
}

void Chain::calculate_genotype_likelihood(std::size_t sample_idx, std::size_t locus_idx)
{
    const std::size_t idx = sample_idx * genotyping_data.num_loci + locus_idx;
    if (genotyping_data.is_missing(locus_idx, sample_idx))
    {
        genotyping_llik_new[idx] = 0;
    }
    else
    {
        auto [latent_genotypes_begin, latent_genotypes_end] = latent_genotypes_new.inner_iterators({sample_idx, locus_idx});
        auto [p_begin, p_end] = p.inner_iterators({locus_idx});
        const float transmission_prob = calc_transmission_process(
            std::span(latent_genotypes_begin, latent_genotypes_end), 
            std::span(p_begin, p_end), m[sample_idx], r[sample_idx]
        );
        const float obs_prob = calc_observation_process(
            std::span(latent_genotypes_begin, latent_genotypes_end),
            genotyping_data.get_observed_alleles(locus_idx, sample_idx),
            eps_neg[sample_idx], eps_pos[sample_idx]
        );
        genotyping_llik_new[idx] = transmission_prob + obs_prob;
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
    const std::size_t num_samples = genotyping_data.num_samples;
    const std::size_t num_loci = genotyping_data.num_loci;

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

    for (std::size_t sample_idx = 0; sample_idx < genotyping_data.num_samples; ++sample_idx)
    {
        const int row_idx = sample_idx * num_loci;
        for (std::size_t locus_idx = 0; locus_idx < genotyping_data.num_loci; ++locus_idx)
        {
            const int idx = row_idx + locus_idx;
            calculate_genotype_likelihood(sample_idx, locus_idx);
            genotyping_llik_old[idx] = genotyping_llik_new[idx];
        }
    }

    for (size_t sample_idx = 0; sample_idx < genotyping_data.num_samples; sample_idx++)
    {
        const float eps_neg_prior = sampler.get_epsilon_log_prior(
            eps_neg[sample_idx], params.eps_neg_alpha, params.eps_neg_beta);
        const float eps_pos_prior = sampler.get_epsilon_log_prior(
            eps_pos[sample_idx], params.eps_pos_alpha, params.eps_pos_beta);

        const float relatedness_prior = sampler.get_relatedness_log_prior(
            r[sample_idx], params.r_alpha, params.r_beta);

        const float coi_prior = sampler.get_coi_log_prior(m[sample_idx], mean_coi);

        eps_neg_prior_old[sample_idx] = eps_neg_prior;
        eps_neg_prior_new[sample_idx] = eps_neg_prior;
        eps_pos_prior_old[sample_idx] = eps_pos_prior;
        eps_pos_prior_new[sample_idx] = eps_pos_prior;
        relatedness_prior_old[sample_idx] = relatedness_prior;
        relatedness_prior_new[sample_idx] = relatedness_prior;
        coi_prior_old[sample_idx] = coi_prior;
        coi_prior_new[sample_idx] = coi_prior;
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
    mean_coi_accept = 0;
    mean_coi_var = 1;
    mean_coi =
        sampler.sample_mean_coi(params.mean_coi_shape, params.mean_coi_scale) +
        1;

    llik = std::numeric_limits<float>::lowest();
    this->temp = temp;

    initialize_parameters();
}
