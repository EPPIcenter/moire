#include "chain.h"

#include "mcmc_utils.h"
#include "sampler.h"

#include <cmath>
#include <algorithm>

#include "env_defs.h"

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
    for (std::size_t sample_idx = 0; sample_idx < genotyping_data.num_samples; ++sample_idx)
    {
        for (std::size_t locus_idx = 0; locus_idx < genotyping_data.num_loci; ++locus_idx)
        {
            const auto obs_alleles = genotyping_data.get_observed_alleles(sample_idx, locus_idx);

            const auto lg = sampler.sample_latent_genotype(obs_alleles, m.at({sample_idx}), .1, .1);

            latent_genotypes_new.inner_fill({sample_idx, locus_idx}, -1);
            latent_genotypes_new.inner_fill({sample_idx, locus_idx}, lg.value);
            lg_adj_new.at({sample_idx, locus_idx}) = lg.log_prob;

            latent_genotypes_old.inner_fill({sample_idx, locus_idx}, -1);
            latent_genotypes_old.inner_fill({sample_idx, locus_idx}, lg.value);
            lg_adj_old.at({sample_idx, locus_idx}) = lg.log_prob;
        }
    }
}


void Chain::initialize_population_responsibility()
{
    population_responsibility_vector.clear();
    population_responsibility_vector.resize({params.num_populations});

    // sample from dirichlet distribution and sort in descending order to ensure there is no label switching
    auto vec = sampler.sample_dirichlet(params.population_responsibility_vector_alpha, 1);
    std::sort(vec.begin(), vec.end(), std::greater<float>());

    population_responsibility_vector.inner_fill(std::span<float const>(vec));

    population_responsibility_vector_prop_var.clear();
    population_responsibility_vector_prop_var.resize({params.num_populations}, 1);
    population_responsibility_vector_accept.clear();
    population_responsibility_vector_accept.resize({params.num_populations}, 0);
    population_responsibility_vector_attempt.clear();
    population_responsibility_vector_attempt.resize({params.num_populations}, 0);
}

// Initialize P with random allele frequencies
void Chain::initialize_p()
{
    p.clear();
    p.resize({params.num_populations, genotyping_data.num_loci}, genotyping_data.num_alleles);
    p_prop_var.clear();
    p_accept.clear();
    p_attempt.clear();

    p = UtilFunctions::calculate_clustered_allele_frequencies<float>(genotyping_data, params.num_populations, sampler);

    // for (size_t pop_idx = 0; pop_idx < params.num_populations; ++pop_idx)
    // {
    //     for (size_t locus_idx = 0; locus_idx < genotyping_data.num_loci; ++locus_idx)
    //     {
    //         p.inner_fill({pop_idx, locus_idx}, sampler.sample_dirichlet(
    //             std::vector<float>(genotyping_data.num_alleles[locus_idx], 1), 10));
    //     }
    // }

    p_prop_var.resize({params.num_populations, genotyping_data.num_loci}, genotyping_data.num_alleles, 1);
    p_accept.resize({params.num_populations, genotyping_data.num_loci}, genotyping_data.num_alleles, 0);
    p_attempt.resize({params.num_populations, genotyping_data.num_loci}, genotyping_data.num_alleles, 0);
};

void Chain::initialize_m()
{
    m.clear();
    m_accept.clear();
    sample_accept.clear();

    m.resize({genotyping_data.num_samples});

    for (std::size_t sample_idx = 0; sample_idx < genotyping_data.num_samples; ++sample_idx)
    {
        std::size_t m_coi = std::clamp(genotyping_data.observed_coi[sample_idx] + sampler.sample_coi_delta(3), std::size_t(1), params.max_coi);
        m.at({sample_idx}) = m_coi;
    }
    m_accept.resize({genotyping_data.num_samples}, 0);
    sample_accept.resize({genotyping_data.num_samples}, 0);
}

void Chain::initialize_eps_neg()
{
    eps_neg.clear();
    eps_neg_accept.clear();
    eps_neg_var.clear();

    eps_neg.resize({genotyping_data.num_samples});
    for (size_t i = 0; i < genotyping_data.num_samples; ++i)
    {
        eps_neg.at({i}) = sampler.sample_unif() * .1;
    }
    eps_neg_accept.resize({genotyping_data.num_samples}, 0);
    eps_neg_var.resize({genotyping_data.num_samples}, 1);
}

void Chain::initialize_eps_pos()
{
    eps_pos.clear();
    eps_pos_accept.clear();
    eps_pos_var.clear();

    eps_pos.resize({genotyping_data.num_samples});
    for (size_t i = 0; i < genotyping_data.num_samples; ++i)
    {
        eps_pos.at({i}) = sampler.sample_unif() * .1;
    }
    eps_pos_accept.resize({genotyping_data.num_samples}, 0);
    eps_pos_var.resize({genotyping_data.num_samples}, 1);
}

void Chain::initialize_r()
{
    r.clear();
    r_accept.clear();
    r_var.clear();
    m_r_accept.clear();
    m_r_var.clear();

    r.resize({genotyping_data.num_samples});
    if (params.allow_relatedness)
    {
        for (size_t i = 0; i < genotyping_data.num_samples; ++i)
        {
            r.at({i}) = sampler.sample_unif() * .99;
        }
    }
    else
    {
        r.resize({genotyping_data.num_samples}, 0.0);
    }
    r_accept.resize({genotyping_data.num_samples}, 0);
    r_var.resize({genotyping_data.num_samples}, 1);
    m_r_accept.resize({genotyping_data.num_samples}, 0);
    m_r_var.resize({genotyping_data.num_samples}, 1);
}

void Chain::initialize_population_coi()
{
    // population_mean_coi.clear();
    // population_mean_coi.resize({params.num_populations});
    // for (size_t pop_idx = 0; pop_idx < params.num_populations; ++pop_idx)
    // {
    //     population_mean_coi.at({pop_idx}) = sampler.sample_unif() * 10;
    // }

    // population_mean_coi_var.resize({params.num_populations}, 1);
    // population_mean_coi_accept.resize({params.num_populations}, 0);

    population_coi_p = .5;
    population_coi_r = 2;
    population_coi_p_accept = 0;
    population_coi_r_accept = 0;
    population_coi_p_sampling_variance = 1;
    population_coi_r_sampling_variance = 1;
}



void Chain::initialize_parameters()
{
    initialize_m();
    initialize_eps_neg();
    initialize_eps_pos();
    initialize_r();
    initialize_latent_genotypes();
    initialize_population_responsibility();
    initialize_p();
    initialize_population_coi();
    initialize_likelihood();
}

void Chain::update_m(int iteration)
{
    auto indices = std::vector<std::size_t>(genotyping_data.num_samples);
    std::iota(indices.begin(), indices.end(), 0);
    sampler.shuffle_vec(indices);

    for (const auto sample_idx : indices)
    {
        calculate_coi_likelihood(sample_idx);
        const std::size_t prop_m = m.at({sample_idx}) + sampler.sample_coi_delta(2);
        if (prop_m > 0 and prop_m <= params.max_coi)
        {
            const float prev_m = m.at({sample_idx});
            m.at({sample_idx}) = prop_m;

            float adj_ratio = 0;

            for (std::size_t pop_idx = 0; pop_idx < params.num_populations; ++pop_idx)
            {
                for (std::size_t locus_idx = 0; locus_idx < genotyping_data.num_loci; ++locus_idx)
                {
                    if (pop_idx == 0)
                    {
                        const auto& lg = sampler.sample_latent_genotype(
                            genotyping_data.get_observed_alleles(sample_idx, locus_idx), 
                            m.at({sample_idx}), 
                            eps_pos.at({sample_idx}), 
                            eps_neg.at({sample_idx})
                        );
                        latent_genotypes_new.inner_fill({sample_idx, locus_idx}, -1);
                        latent_genotypes_new.inner_fill({sample_idx, locus_idx}, lg.value);
                        lg_adj_new.at({sample_idx, locus_idx}) = lg.log_prob;
                        adj_ratio = adj_ratio + lg_adj_new.at({sample_idx, locus_idx}) - lg_adj_old.at({sample_idx, locus_idx});
                        calculate_observation_likelihood(sample_idx, locus_idx);
                    }
                    calculate_transmission_likelihood(sample_idx, pop_idx, locus_idx);
                }
            }

            const float new_llik = calc_new_likelihood();
            const float new_prior = calc_new_prior();
            const float new_post = new_llik * temp + new_prior;
            const float alpha = sampler.sample_log_mh_acceptance();
            const float mh_ratio = new_post - get_posterior() + adj_ratio;
            if (!std::isnan(new_post) && alpha <= mh_ratio)
            // Accept
            {
                llik = new_llik;
                prior = new_prior;

                save_coi_likelihood(sample_idx);
                for (std::size_t pop_idx = 0; pop_idx < params.num_populations; ++pop_idx)
                {
                    for (std::size_t locus_idx = 0; locus_idx < genotyping_data.num_loci; ++locus_idx)
                    {
                        if (pop_idx == 0)
                        {
                            save_observation_likelihood(sample_idx, locus_idx);
                            lg_adj_old.at({sample_idx, locus_idx}) = lg_adj_new.at({sample_idx, locus_idx});
                            auto [begin, end] = latent_genotypes_new.inner_iterators({sample_idx, locus_idx});
                            std::copy(begin, end, latent_genotypes_old.inner_begin({sample_idx, locus_idx}));
                        }
                        save_transmission_likelihood(sample_idx, pop_idx, locus_idx);
                    }
                }
                m_accept.at({sample_idx}) = m_accept.at({sample_idx}) + 1;
            }
            else
            // Reject
            {
                m.at({sample_idx}) = prev_m;
                restore_coi_likelihood(sample_idx);
                for (std::size_t pop_idx = 0; pop_idx < params.num_populations; ++pop_idx)
                {
                    for (std::size_t locus_idx = 0; locus_idx < genotyping_data.num_loci; ++locus_idx)
                    {
                        if (pop_idx == 0)
                        {
                            restore_observation_likelihood(sample_idx, locus_idx);
                            lg_adj_new.at({sample_idx, locus_idx}) = lg_adj_old.at({sample_idx, locus_idx});
                            auto [begin, end] = latent_genotypes_old.inner_iterators({sample_idx, locus_idx});
                            std::copy(begin, end, latent_genotypes_new.inner_begin({sample_idx, locus_idx}));
                        }
                        restore_transmission_likelihood(sample_idx, pop_idx, locus_idx);
                    }
                }
            }
        }
    }
}

void Chain::update_eff_coi(int iteration)
{
    auto sample_indices = std::vector<std::size_t>(genotyping_data.num_samples);
    std::iota(sample_indices.begin(), sample_indices.end(), 0);
    sampler.shuffle_vec(sample_indices);

    for (const auto sample_idx : sample_indices)
    {
        const float curr_eff_coi = (m.at({sample_idx}) - 1) * (1.0f - r.at({sample_idx})) + 1.0f;
        const auto prop_adj = sampler.sample_constrained(
            curr_eff_coi, m_r_var.at({sample_idx}), 1, params.max_coi);

        const float prop_eff_coi = std::get<0>(prop_adj);
        const float adj = std::get<1>(prop_adj);

        const std::size_t prop_m = m.at({sample_idx}) + sampler.sample_coi_delta(2);
        const float prop_r = 1.0f - (prop_eff_coi - 1.0f) / (prop_m - 1.0f);

        if (prop_m <= 0 || prop_m > params.max_coi || prop_r > .99999 || prop_r < .00001 || std::isnan(prop_r)) {
            continue;
        }

        const int prev_m = m.at({sample_idx});
        const float prev_r = r.at({sample_idx});
        m.at({sample_idx}) = prop_m;
        r.at({sample_idx}) = prop_r;


        float adj_ratio = adj;
        calculate_relatedness_likelihood(sample_idx);
        calculate_coi_likelihood(sample_idx);

        for (std::size_t pop_idx = 0; pop_idx < params.num_populations; ++pop_idx)
        {

            for (std::size_t locus_idx = 0; locus_idx < genotyping_data.num_loci; ++locus_idx)
            {
                if (pop_idx == 0)
                {
                    const auto &lg = sampler.sample_latent_genotype(
                        genotyping_data.get_observed_alleles(sample_idx, locus_idx), 
                        m.at({sample_idx}), 
                        eps_pos.at({sample_idx}),
                        eps_neg.at({sample_idx})
                    );
                    latent_genotypes_new.inner_fill({sample_idx, locus_idx}, -1);
                    latent_genotypes_new.inner_fill({sample_idx, locus_idx}, lg.value);
                    lg_adj_new.at({sample_idx, locus_idx}) = lg.log_prob;
                    adj_ratio = adj_ratio + lg_adj_new.at({sample_idx, locus_idx}) - lg_adj_old.at({sample_idx, locus_idx});
                    calculate_observation_likelihood(sample_idx, locus_idx);
                }
                calculate_transmission_likelihood(sample_idx, pop_idx, locus_idx);
            }

        }

        const float new_llik = calc_new_likelihood();
        const float new_prior = calc_new_prior();
        const float new_post = new_llik * temp + new_prior;
        const double mh_ratio = new_post - get_posterior() + adj_ratio;
        const double alpha = sampler.sample_log_mh_acceptance();

        // Reject
        if (std::isnan(new_post) || std::isnan(mh_ratio) || alpha > mh_ratio)
        {
            m.at({sample_idx}) = prev_m;
            r.at({sample_idx}) = prev_r;
            restore_relatedness_likelihood(sample_idx);
            restore_coi_likelihood(sample_idx);
            for (std::size_t pop_idx = 0; pop_idx < params.num_populations; ++pop_idx)
            {
                for (std::size_t locus_idx = 0; locus_idx < genotyping_data.num_loci; ++locus_idx)
                {
                    restore_transmission_likelihood(sample_idx, pop_idx, locus_idx);
                    if (pop_idx == 0)
                    {
                        lg_adj_new.at({sample_idx, locus_idx}) = lg_adj_old.at({sample_idx, locus_idx});
                        const auto [begin, end] = latent_genotypes_old.inner_iterators({sample_idx, locus_idx});
                        std::copy(begin, end, latent_genotypes_new.inner_begin({sample_idx, locus_idx}));
                        restore_observation_likelihood(sample_idx, locus_idx);
                    }
                }
            }
        }
        // Accept
        else
        {
            llik = new_llik;
            prior = new_prior;
            save_relatedness_likelihood(sample_idx);
            save_coi_likelihood(sample_idx);
            for (std::size_t pop_idx = 0; pop_idx < params.num_populations; ++pop_idx)
            {
                for (std::size_t locus_idx = 0; locus_idx < genotyping_data.num_loci; ++locus_idx)
                {
                    save_transmission_likelihood(sample_idx, pop_idx, locus_idx);
                    if (pop_idx == 0)
                    {
                        lg_adj_old.at({sample_idx, locus_idx}) = lg_adj_new.at({sample_idx, locus_idx});
                        const auto [begin, end] = latent_genotypes_new.inner_iterators({sample_idx, locus_idx});
                        std::copy(begin, end, latent_genotypes_old.inner_begin({sample_idx, locus_idx}));
                        save_observation_likelihood(sample_idx, locus_idx);
                    }

                }
            }
            ++m_r_accept.at({sample_idx});
        }

        if (iteration < params.burnin and
            iteration > 15)  // don't start adapting until there are
                             // at least a few samples
        {
            const float acceptanceRate = m_r_accept.at({sample_idx}) / float(iteration);
            const float update =
                (acceptanceRate - .23) / std::pow(iteration + 1, .5);
            m_r_var.at({sample_idx}) = std::max(m_r_var.at({sample_idx}) + update, .01f);
        }
    }
}

void Chain::update_r(int iteration)
{
    auto sample_indices = std::vector<size_t>(genotyping_data.num_samples);
    std::iota(sample_indices.begin(), sample_indices.end(), 0);
    sampler.shuffle_vec(sample_indices);

    for (const auto sample_idx : sample_indices)
    {
        const auto prop_adj =
            sampler.sample_constrained(r.at({sample_idx}), r_var.at({sample_idx}), .00001, .99999);
        const float prop_r = std::get<0>(prop_adj);
        const float adj = std::get<1>(prop_adj);

        const float prev_r = r.at({sample_idx});
        r.at({sample_idx}) = prop_r;
        calculate_relatedness_likelihood(sample_idx);

        for (std::size_t locus_idx = 0; locus_idx < genotyping_data.num_loci; ++locus_idx)
        {
            for (std::size_t pop_idx = 0; pop_idx < params.num_populations; ++pop_idx)
            {
                calculate_transmission_likelihood(sample_idx, pop_idx, locus_idx);
            }
        }

        const float new_llik = calc_new_likelihood();
        const float new_prior = calc_new_prior();
        const float new_post = new_llik * temp + new_prior;

        // Reject
        if (std::isnan(new_post) or sampler.sample_log_mh_acceptance() >
                                        (new_post - get_posterior() + adj))
        {
            r.at({sample_idx}) = prev_r;
            restore_relatedness_likelihood(sample_idx);
            for (std::size_t locus_idx = 0; locus_idx < genotyping_data.num_loci; ++locus_idx)
            {
                for (std::size_t pop_idx = 0; pop_idx < params.num_populations; ++pop_idx)
                {
                    restore_transmission_likelihood(sample_idx, pop_idx, locus_idx);
                }
            }
        }
        else
        {
            llik = new_llik;
            prior = new_prior;
            save_relatedness_likelihood(sample_idx);
            for (std::size_t locus_idx = 0; locus_idx < genotyping_data.num_loci; ++locus_idx)
            {
                for (std::size_t pop_idx = 0; pop_idx < params.num_populations; ++pop_idx)
                {
                    save_transmission_likelihood(sample_idx, pop_idx, locus_idx);
                }
            }
            ++r_accept.at({sample_idx});
        }

        if (iteration < params.burnin and
            iteration > 15)  // don't start adapting until there are
                             // at least a few samples
        {
            const float acceptanceRate = r_accept.at({sample_idx}) / float(iteration);
            const float update =
                (acceptanceRate - .23) / std::pow(iteration + 1, .5);
            r_var.at({sample_idx}) = std::max(r_var.at({sample_idx}) + update, .0001f);
        }
    }
}

void Chain::update_m_r(int iteration)
{
    auto sample_indices = std::vector<std::size_t>(genotyping_data.num_samples);
    std::iota(sample_indices.begin(), sample_indices.end(), 0);
    sampler.shuffle_vec(sample_indices);

    for (const auto sample_idx : sample_indices)
    {
        const auto prop_adj =
            sampler.sample_constrained(r.at({sample_idx}), m_r_var.at({sample_idx}), .00001, .99999);
        const float prop_r = std::get<0>(prop_adj);
        const float adj = std::get<1>(prop_adj);

        const std::size_t prop_m = m.at({sample_idx}) + sampler.sample_coi_delta(1);

        if (prop_m <= 0 or prop_m > params.max_coi) {
            continue;
        }

        const float prev_m = m.at({sample_idx});
        const float prev_r = r.at({sample_idx});
        r.at({sample_idx}) = prop_r;
        calculate_relatedness_likelihood(sample_idx);
        m.at({sample_idx}) = prop_m;

        float adj_ratio = adj;

        calculate_coi_likelihood(sample_idx);
        for (std::size_t pop_idx = 0; pop_idx < params.num_populations; ++pop_idx)
        {
            for (std::size_t locus_idx = 0; locus_idx < genotyping_data.num_loci; ++locus_idx)
            {
                if (pop_idx == 0)
                {
                    const auto &lg = sampler.sample_latent_genotype(
                        genotyping_data.get_observed_alleles(sample_idx, locus_idx), 
                        m.at({sample_idx}), 
                        eps_pos.at({sample_idx}),
                        eps_neg.at({sample_idx})
                    );
                    latent_genotypes_new.inner_fill({sample_idx, locus_idx}, -1);
                    latent_genotypes_new.inner_fill({sample_idx, locus_idx}, lg.value);
                    lg_adj_new.at({sample_idx, locus_idx}) = lg.log_prob;
                    adj_ratio = adj_ratio + lg_adj_new.at({sample_idx, locus_idx}) - lg_adj_old.at({sample_idx, locus_idx});
                    calculate_observation_likelihood(sample_idx, locus_idx);
                }
                calculate_transmission_likelihood(sample_idx, pop_idx, locus_idx);
            }
        }

        const float new_llik = calc_new_likelihood();
        const float new_prior = calc_new_prior();
        const float new_post = new_llik * temp + new_prior;
        const float alpha = sampler.sample_log_mh_acceptance();
        const float mh_ratio = new_post - get_posterior() + adj_ratio;

        // Reject
        if (std::isnan(new_post) or alpha > mh_ratio)
        {
            r.at({sample_idx}) = prev_r;
            restore_relatedness_likelihood(sample_idx);
            m.at({sample_idx}) = prev_m;
            restore_coi_likelihood(sample_idx);
            for (std::size_t pop_idx = 0; pop_idx < params.num_populations; ++pop_idx)
            {
                for (std::size_t locus_idx = 0; locus_idx < genotyping_data.num_loci; ++locus_idx)
                {
                    if (pop_idx == 0)
                    {
                        restore_observation_likelihood(sample_idx, locus_idx);
                        lg_adj_new.at({sample_idx, locus_idx}) = lg_adj_old.at({sample_idx, locus_idx});
                        const auto [begin, end] = latent_genotypes_old.inner_iterators({sample_idx, locus_idx});
                        std::copy(begin, end, latent_genotypes_new.inner_begin({sample_idx, locus_idx}));
                    }
                    restore_transmission_likelihood(sample_idx, pop_idx, locus_idx);
                }
            }
        }
        else
        {
            llik = new_llik;
            prior = new_prior;
            save_relatedness_likelihood(sample_idx);
            save_coi_likelihood(sample_idx);
            for (std::size_t pop_idx = 0; pop_idx < params.num_populations; ++pop_idx)
            {
                for (std::size_t locus_idx = 0; locus_idx < genotyping_data.num_loci; ++locus_idx)
                {
                    if (pop_idx == 0)
                    {
                        save_observation_likelihood(sample_idx, locus_idx);
                        lg_adj_old.at({sample_idx, locus_idx}) = lg_adj_new.at({sample_idx, locus_idx});
                        const auto [begin, end] = latent_genotypes_new.inner_iterators({sample_idx, locus_idx});
                        std::copy(begin, end, latent_genotypes_old.inner_begin({sample_idx, locus_idx}));
                    }
                    save_transmission_likelihood(sample_idx, pop_idx, locus_idx);
                }
            }
            ++m_r_accept.at({sample_idx});
        }

        if (iteration < params.burnin and
            iteration > 15)  // don't start adapting until there are
                             // at least a few samples
        {
            const float acceptanceRate = m_r_accept.at({sample_idx}) / float(iteration);
            const float update =
                (acceptanceRate - .23) / std::pow(iteration + 1, .5);
            m_r_var.at({sample_idx}) = std::max(m_r_var.at({sample_idx}) + update, .0001f);
        }
    }
}

/*
 * SALT Sampler approach.
 * https://doi.org/10.1080/00949655.2017.1376063
 */
void Chain::update_p(int iteration)
{
    auto locus_indices = std::vector<std::size_t>(genotyping_data.num_loci);
    auto pop_indices = std::vector<std::size_t>(params.num_populations);
    std::iota(locus_indices.begin(), locus_indices.end(), 0);
    std::iota(pop_indices.begin(), pop_indices.end(), 0);
    sampler.shuffle_vec(locus_indices);
    sampler.shuffle_vec(pop_indices);

    for (const std::size_t pop_idx : pop_indices) {
        for (const std::size_t locus_idx : locus_indices)
        {
            const std::size_t num_alleles = p.ragged_dimensions(locus_idx);
            // std::size_t rep = num_alleles; 
            std::size_t rep = 3;
            const auto [begin, end] = p.inner_iterators({pop_idx, locus_idx});
            while (rep-- > 0)
            {
                const size_t allele_idx = sampler.sample_random_int(0, num_alleles - 1);

                ++p_attempt.at({pop_idx, locus_idx, allele_idx});

                auto logitPropP = UtilFunctions::logitVec(begin, end);

                const float logitCurr = logitPropP[allele_idx];
                const float logitProp =
                    sampler.sample_epsilon(logitCurr, p_prop_var.at({pop_idx, locus_idx, allele_idx}));

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
                p.inner_fill({pop_idx, locus_idx}, prop_p);

                for (std::size_t sample_idx = 0; sample_idx < genotyping_data.num_samples; ++sample_idx)
                {
                    calculate_transmission_likelihood(sample_idx, pop_idx, locus_idx);
                }

                const float new_llik = calc_new_likelihood();
                const float new_prior = calc_new_prior();
                const float new_post = new_llik * temp + new_prior;

                const float acceptanceRatio = new_post - get_posterior() + logAdj;

                if (std::isnan(new_post) or
                    sampler.sample_log_mh_acceptance() > acceptanceRatio)
                {
                    p.inner_fill({pop_idx, locus_idx}, prev_p);
                    for (std::size_t sample_idx = 0; sample_idx < genotyping_data.num_samples; ++sample_idx)
                    {
                        for (std::size_t pop_idx = 0; pop_idx < params.num_populations; ++pop_idx)
                        {
                            restore_transmission_likelihood(sample_idx, pop_idx, locus_idx);
                        }
                    }
                }
                else
                {
                    llik = new_llik;
                    prior = new_prior;
                    for (std::size_t sample_idx = 0; sample_idx < genotyping_data.num_samples; ++sample_idx)
                    {
                        for (std::size_t pop_idx = 0; pop_idx < params.num_populations; ++pop_idx)
                        {
                            save_transmission_likelihood(sample_idx, pop_idx, locus_idx);
                        }
                    }
                    ++p_accept.at({pop_idx, locus_idx, allele_idx});
                }

                // don't start adapting until there are at least a few samples
                if (iteration < params.burnin and p_attempt.at({pop_idx, locus_idx, allele_idx}) > 15)
                {
                    const float acceptanceRate =
                        (p_accept.at({pop_idx, locus_idx, allele_idx}) + 1) / (float(p_attempt.at({pop_idx, locus_idx, allele_idx})) + 1);
                    p_prop_var.at({pop_idx, locus_idx, allele_idx}) += (acceptanceRate - .23) /
                                        std::pow(p_attempt.at({pop_idx, locus_idx, allele_idx}) + 1, .5);
                    p_prop_var.at({pop_idx, locus_idx, allele_idx}) = std::max(p_prop_var.at({pop_idx, locus_idx, allele_idx}), .01f);
                }
            }
        }
    }
}

void Chain::update_eps_pos(int iteration)
{
    auto sample_indices = std::vector<std::size_t>(genotyping_data.num_samples);
    std::iota(sample_indices.begin(), sample_indices.end(), 0);
    sampler.shuffle_vec(sample_indices);

    for (const auto sample_idx : sample_indices)
    {
        const auto prop_adj = sampler.sample_constrained(
            eps_pos.at({sample_idx}), eps_pos_var.at({sample_idx}), min_sampled, 1);
        const float prop_eps_pos = std::get<0>(prop_adj);
        const float adj = std::get<1>(prop_adj);

        if (prop_eps_pos < 1 && prop_eps_pos > 1e-32)
        {
            const float prev_eps_pos = eps_pos.at({sample_idx});
            eps_pos.at({sample_idx}) = prop_eps_pos;
            calculate_eps_pos_likelihood(sample_idx);

            for (std::size_t locus_idx = 0; locus_idx < genotyping_data.num_loci; ++locus_idx)
            {
                calculate_observation_likelihood(sample_idx, locus_idx);
            }

            const float new_llik = calc_new_likelihood();
            const float new_prior = calc_new_prior();
            const float new_post = new_llik * temp + new_prior;

            // Reject
            if (std::isnan(new_post) or sampler.sample_log_mh_acceptance() >
                                            (new_post - get_posterior() + adj))
            {
                eps_pos.at({sample_idx}) = prev_eps_pos;
                restore_eps_pos_likelihood(sample_idx);
                for (std::size_t locus_idx = 0; locus_idx < genotyping_data.num_loci; ++locus_idx)
                {
                    restore_observation_likelihood(sample_idx, locus_idx);
                }
            }
            else
            {
                llik = new_llik;
                prior = new_prior;
                save_eps_pos_likelihood(sample_idx);
                for (std::size_t locus_idx = 0; locus_idx < genotyping_data.num_loci; ++locus_idx)
                {
                    save_observation_likelihood(sample_idx, locus_idx);
                }
                ++eps_pos_accept.at({sample_idx});
            }

            if (iteration < params.burnin and
                iteration > 15)  // don't start adapting until there are
                                 // at least a few samples
            {
                const float acceptanceRate =
                    eps_pos_accept.at({sample_idx}) / float(iteration);
                const float update =
                    (acceptanceRate - .23) / std::pow(iteration + 1, .5);
                eps_pos_var.at({sample_idx}) = std::max(eps_pos_var.at({sample_idx}) + update, .0001f);
            }
        }
    }
}

void Chain::update_eps_neg(int iteration)
{
    auto sample_indices = std::vector<std::size_t>(genotyping_data.num_samples);
    std::iota(sample_indices.begin(), sample_indices.end(), 0);
    sampler.shuffle_vec(sample_indices);

    for (const auto sample_idx : sample_indices)
    {
        const auto prop_adj = sampler.sample_constrained(
            eps_neg.at({sample_idx}), eps_neg_var.at({sample_idx}), min_sampled, 1);
        const float prop_eps_neg = std::get<0>(prop_adj);
        const float adj = std::get<1>(prop_adj);

        if (prop_eps_neg < 1 && prop_eps_neg > 1e-32)
        {
            const float prev_eps_neg = eps_neg.at({sample_idx});
            eps_neg.at({sample_idx}) = prop_eps_neg;
            calculate_eps_neg_likelihood(sample_idx);

            for (std::size_t locus_idx = 0; locus_idx < genotyping_data.num_loci; ++locus_idx)
            {
                calculate_observation_likelihood(sample_idx, locus_idx);
            }

            const float new_llik = calc_new_likelihood();
            const float new_prior = calc_new_prior();
            const float new_post = new_llik * temp + new_prior;

            // Reject
            if (std::isnan(new_post) or sampler.sample_log_mh_acceptance() >
                                            (new_post - get_posterior() + adj))
            {
                eps_neg.at({sample_idx}) = prev_eps_neg;
                restore_eps_neg_likelihood(sample_idx);
                for (std::size_t locus_idx = 0; locus_idx < genotyping_data.num_loci; ++locus_idx)
                {
                    restore_observation_likelihood(sample_idx, locus_idx);
                }
            }
            else
            {
                llik = new_llik;
                prior = new_prior;
                save_eps_neg_likelihood(sample_idx);
                for (std::size_t locus_idx = 0; locus_idx < genotyping_data.num_loci; ++locus_idx)
                {
                    save_observation_likelihood(sample_idx, locus_idx);
                }
                ++eps_neg_accept.at({sample_idx});
            }

            if (iteration < params.burnin and
                iteration > 15)  // don't start adapting until there are
                                 // at least a few samples
            {
                const float acceptanceRate =
                    eps_neg_accept.at({sample_idx}) / float(iteration);
                const float update =
                    (acceptanceRate - .23) / std::pow(iteration + 1, .5);
                eps_neg_var.at({sample_idx}) = std::max(eps_neg_var.at({sample_idx}) + update, .0001f);
            }
        }
    }
}

void Chain::update_samples(int iteration)
{
    auto sample_indices = std::vector<std::size_t>(genotyping_data.num_samples);
    std::iota(sample_indices.begin(), sample_indices.end(), 0);
    sampler.shuffle_vec(sample_indices);

    for (const auto sample_idx : sample_indices)
    {
        const auto eps_neg_prop_adj = sampler.sample_constrained(
            eps_neg.at({sample_idx}), eps_neg_var.at({sample_idx}), min_sampled, 1);
        const float prop_eps_neg = std::get<0>(eps_neg_prop_adj);
        const float eps_neg_adj = std::get<1>(eps_neg_prop_adj);
        const bool valid_prop_eps_neg =
            prop_eps_neg < 1 && prop_eps_neg > 1e-32;

        const auto eps_pos_prop_adj = sampler.sample_constrained(
            eps_pos.at({sample_idx}), eps_pos_var.at({sample_idx}), min_sampled, 1);
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
                sampler.sample_constrained(r.at({sample_idx}), r_var.at({sample_idx}), min_sampled, .99);
            prop_r = std::get<0>(r_prop_adj);
            r_adj = std::get<1>(r_prop_adj);
            valid_prop_r = prop_r < 1 && prop_r > 1e-32;
        }

        const std::size_t prop_m = m.at({sample_idx}) + sampler.sample_coi_delta(2);
        const bool valid_prop_m = prop_m > 0 && prop_m <= params.max_coi;

        if (valid_prop_eps_neg && valid_prop_eps_pos && valid_prop_r &&
            valid_prop_m)
        {
            const float prev_eps_pos = eps_pos.at({sample_idx});
            eps_pos.at({sample_idx}) = prop_eps_pos;
            calculate_eps_pos_likelihood(sample_idx);

            const float prev_eps_neg = eps_neg.at({sample_idx});
            eps_neg.at({sample_idx}) = prop_eps_neg;
            calculate_eps_neg_likelihood(sample_idx);

            const float prev_r = r.at({sample_idx});
            r.at({sample_idx}) = prop_r;
            calculate_relatedness_likelihood(sample_idx);

            const float prev_m = m.at({sample_idx});
            m.at({sample_idx}) = prop_m;

            float adj_ratio = eps_neg_adj + eps_pos_adj + r_adj;

            calculate_coi_likelihood(sample_idx);
            for (std::size_t pop_idx = 0; pop_idx < params.num_populations; ++pop_idx)
            {
                for (std::size_t locus_idx = 0; locus_idx < genotyping_data.num_loci; ++locus_idx)
                {
                    if (pop_idx == 0)
                    {
                        auto lg = sampler.sample_latent_genotype(
                            genotyping_data.get_observed_alleles(sample_idx, locus_idx), m.at({sample_idx}),
                            eps_pos.at({sample_idx}), eps_neg.at({sample_idx}));
                        latent_genotypes_new.inner_fill({sample_idx, locus_idx}, -1);
                        latent_genotypes_new.inner_fill({sample_idx, locus_idx}, lg.value);
                        lg_adj_new.at({sample_idx, locus_idx}) = lg.log_prob;
                        adj_ratio = adj_ratio + lg_adj_new.at({sample_idx, locus_idx}) - lg_adj_old.at({sample_idx, locus_idx});
                        calculate_observation_likelihood(sample_idx, locus_idx);
                    }
                    calculate_transmission_likelihood(sample_idx, pop_idx, locus_idx);
                }
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
                for (std::size_t pop_idx = 0; pop_idx < params.num_populations; ++pop_idx)
                {
                    for (std::size_t locus_idx = 0; locus_idx < genotyping_data.num_loci; ++locus_idx)
                    {
                        if (pop_idx == 0)
                        {
                            save_observation_likelihood(sample_idx, locus_idx);
                            lg_adj_old.at({sample_idx, locus_idx}) = lg_adj_new.at({sample_idx, locus_idx});
                            auto [begin, end] = latent_genotypes_new.inner_iterators({sample_idx, locus_idx});
                            std::copy(begin, end, latent_genotypes_old.inner_begin({sample_idx, locus_idx}));
                        }
                        save_transmission_likelihood(sample_idx, pop_idx, locus_idx);
                    }
                }
                ++sample_accept.at({sample_idx});
            }
            else
            {
                m.at({sample_idx}) = prev_m;
                eps_pos.at({sample_idx}) = prev_eps_pos;
                eps_neg.at({sample_idx}) = prev_eps_neg;
                r.at({sample_idx}) = prev_r;
                restore_eps_neg_likelihood(sample_idx);
                restore_eps_pos_likelihood(sample_idx);
                restore_relatedness_likelihood(sample_idx);
                restore_coi_likelihood(sample_idx);
                for (std::size_t pop_idx = 0; pop_idx < params.num_populations; ++pop_idx)
                {
                    for (std::size_t locus_idx = 0; locus_idx < genotyping_data.num_loci; ++locus_idx)
                    {
                        if (pop_idx == 0)
                        {
                            restore_observation_likelihood(sample_idx, locus_idx);
                            lg_adj_new.at({sample_idx, locus_idx}) = lg_adj_old.at({sample_idx, locus_idx});
                            auto [begin, end] = latent_genotypes_old.inner_iterators({sample_idx, locus_idx});
                            std::copy(begin, end, latent_genotypes_new.inner_begin({sample_idx, locus_idx}));
                        }
                        restore_transmission_likelihood(sample_idx, pop_idx, locus_idx);
                    }
                }
            }
        }
    }
}

void Chain::update_population_coi_p(int iteration)
{
    const auto prop_adj = sampler.sample_constrained(population_coi_p, population_coi_p_sampling_variance,
                                                    0.01, 0.99);
    const float prop_p = std::get<0>(prop_adj);
    const float adj = std::get<1>(prop_adj);

    const float prev_p = population_coi_p;
    population_coi_p = prop_p;
    calculate_population_coi_p_likelihood();
    for (std::size_t sample_idx = 0; sample_idx < genotyping_data.num_samples; ++sample_idx)
    {
        calculate_coi_likelihood(sample_idx);
    }
    const float new_llik = calc_new_likelihood();
    const float new_prior = calc_new_prior();
    const float new_post = new_llik * temp + new_prior;

    const float alpha = sampler.sample_log_mh_acceptance();
    if (!std::isnan(new_post) and alpha <= (new_post - get_posterior() + adj))
    {
        llik = new_llik;
        prior = new_prior;
        save_population_coi_p_likelihood();
        for (std::size_t sample_idx = 0; sample_idx < genotyping_data.num_samples; ++sample_idx)
        {
            save_coi_likelihood(sample_idx);
        }
        ++population_coi_p_accept;
    }
    else
    {
        population_coi_p = prev_p;
        restore_population_coi_p_likelihood();
        for (std::size_t sample_idx = 0; sample_idx < genotyping_data.num_samples; ++sample_idx)
        {
            restore_coi_likelihood(sample_idx);
        }
    }

    if (iteration < params.burnin and iteration > 15)
    {
        const float acceptance_rate = population_coi_p_accept / (iteration + 1);
        const float update =
            (acceptance_rate - .23) / std::pow(iteration + 1, .5);
        population_coi_p_sampling_variance = std::max(population_coi_p_sampling_variance + update, .0001f);
    }
}

void Chain::update_population_coi_r(int iteration)
{
    const auto prop_adj = sampler.sample_constrained(population_coi_r, population_coi_r_sampling_variance,
                                                    .01, 100);
    const float prop_r = std::get<0>(prop_adj);
    const float adj = std::get<1>(prop_adj);

    const float prev_r = population_coi_r;
    population_coi_r = prop_r;
    calculate_population_coi_r_likelihood();
    for (std::size_t sample_idx = 0; sample_idx < genotyping_data.num_samples; ++sample_idx)
    {
        calculate_coi_likelihood(sample_idx);
    }
    const float new_llik = calc_new_likelihood();
    const float new_prior = calc_new_prior();
    const float new_post = new_llik * temp + new_prior;

    const float alpha = sampler.sample_log_mh_acceptance();

    if (!std::isnan(new_post) and alpha <= (new_post - get_posterior() + adj))
    {
        llik = new_llik;
        prior = new_prior;
        save_population_coi_r_likelihood();
        for (std::size_t sample_idx = 0; sample_idx < genotyping_data.num_samples; ++sample_idx)
        {
            save_coi_likelihood(sample_idx);
        }
        ++population_coi_r_accept;
    }
    else
    {
        population_coi_r = prev_r;
        restore_population_coi_r_likelihood();
        for (std::size_t sample_idx = 0; sample_idx < genotyping_data.num_samples; ++sample_idx)
        {
            restore_coi_likelihood(sample_idx);
        }
    }

    if (iteration < params.burnin and iteration > 15)
    {
        const float acceptance_rate = population_coi_r_accept / (iteration + 1);
        const float update = (acceptance_rate - .23) / std::pow(iteration + 1, .5);
        population_coi_r_sampling_variance = std::max(population_coi_r_sampling_variance + update, .0001f);
    }
}

void Chain::update_population_responsibility_vector(int iteration) 
{
    const std::size_t num_populations = params.num_populations;
    if (num_populations == 1) {
        return;
    }
    std::size_t rep = 2;
    const auto [begin, end] = population_responsibility_vector.inner_iterators();
    while (rep-- > 0) 
    {
        const size_t pop_idx = sampler.sample_random_int(0, num_populations - 1);
        ++population_responsibility_vector_attempt.at({pop_idx});

        auto logitPropP = UtilFunctions::logitVec(begin, end);

        const float logitCurr = logitPropP[pop_idx];
        const float logitProp = sampler.sample_epsilon(logitCurr, population_responsibility_vector_prop_var.at({pop_idx}));

        const auto currLogPQ = UtilFunctions::log_pq(logitCurr);
        const auto propLogPQ = UtilFunctions::log_pq(logitProp);

        logitPropP.erase(logitPropP.begin() + pop_idx);

        const float ls = propLogPQ.second - UtilFunctions::logitSum(logitPropP);
        logitPropP = UtilFunctions::logitScale(logitPropP, ls);
        logitPropP.insert(logitPropP.begin() + pop_idx, logitProp);

        const float logAdj = (currLogPQ.first - propLogPQ.first) + (num_populations - 1) * (currLogPQ.second - propLogPQ.second);
        auto prop_p = UtilFunctions::expitVec(logitPropP);

        bool sub_threshold_flag = false;
        float sum = 0.0f;
        for (auto &el : prop_p)
        {
            el = std::max(el, 1e-6f);
            sum += el;
        }
        for (auto &el : prop_p) {
            el = el / sum;
        }

        // sort the proposed population responsibility vector 
        auto prop_p_begin = prop_p.begin();
        auto prop_p_end = prop_p.end();
        auto [var_begin, var_end] = population_responsibility_vector_prop_var.inner_iterators();
        // // Create a vector of indices
        std::vector<size_t> indices(params.num_populations);
        std::iota(indices.begin(), indices.end(), 0);
        
        // // Sort indices based on population responsibility values
        std::sort(indices.begin(), indices.end(), 
            [prop_p_begin](size_t i1, size_t i2) { return prop_p_begin[i1] > prop_p_begin[i2]; });
        
        std::vector<float> sorted_resp(params.num_populations);
        // std::vector<float> sorted_var(params.num_populations);

        for (size_t i = 0; i < params.num_populations; ++i) {
            sorted_resp[i] = prop_p_begin[indices[i]];
        }

        const auto prev_p = std::vector<float>(begin, end);
        // population_responsibility_vector.inner_fill(prop_p);
        population_responsibility_vector.inner_fill(sorted_resp);
        calculate_population_responsibility_vector_likelihood();

        const float new_llik = calc_new_likelihood();
        const float new_prior = calc_new_prior();
        const float new_post = new_llik * temp + new_prior;

        const float acceptance_ratio = new_post - get_posterior() + logAdj;

        bool reject = sampler.sample_log_mh_acceptance() > acceptance_ratio;

        if (std::isnan(new_post) or reject)
        {
            population_responsibility_vector.inner_fill(prev_p);
            restore_population_responsibility_vector_likelihood();
        }
        else
        {
            llik = new_llik;
            prior = new_prior;
            save_population_responsibility_vector_likelihood();
            ++population_responsibility_vector_accept.at({pop_idx});
            // for (size_t i = 0; i < params.num_populations; ++i) {
            //     sorted_var[i] = var_begin[indices[i]];
            // }
            // population_responsibility_vector_prop_var.inner_fill(sorted_var);
        }

        if (iteration < params.burnin and iteration > 15)
        {
            const float acceptance_rate = population_responsibility_vector_accept.at({pop_idx}) / (float(population_responsibility_vector_attempt.at({pop_idx})) + 1);
            const float update = (acceptance_rate - .23) / std::pow(float(population_responsibility_vector_attempt.at({pop_idx})) + 1, .5);
            population_responsibility_vector_prop_var.at({pop_idx}) = std::max(population_responsibility_vector_prop_var.at({pop_idx}) + update, .0001f);
        }
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

    if (allele_index_vec.size() > coi) 
    {
        throw std::invalid_argument("Allele index vec size is greater than COI");
    }


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

    // break early if relatedness is not allowed
    if (!params.allow_relatedness)
    {
        return std::log(1 - probAnyMissing_(prVec_, coi)) +
               log_constrained_set_total_prob * coi;
    }

    // Only go up to coi - total_alleles because there must be at least
    // total_alleles unrelated strains at this locus
    const std::vector<double> pamVec = probAnyMissing_.vectorized(prVec_, coi);
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

    // if (std::isnan(UtilFunctions::logSumExp(res))) {
    //     std::cout << "Transmission probability is NaN" << std::endl;
    //     std::cout << "Allele frequencies: ";
    //     for (const auto &e : allele_frequencies) {
    //         std::cout << e << " ";
    //     }
    //     std::cout << std::endl;
    //     std::cout << "COI: " << coi << std::endl;
    //     std::cout << "Relatedness: " << relatedness << std::endl;
    //     std::cout << "Allele index vec: ";
    //     for (const auto &e : allele_index_vec) {
    //         std::cout << e << " ";
    //     }
    //     std::cout << std::endl;
    //     std::cout << "PrVec_: ";
    //     for (const auto &e : prVec_) {
    //         std::cout << e << " ";
    //     }
    //     std::cout << std::endl;
    //     std::cout << "Res: ";
    //     for (const auto &e : res) {
    //         std::cout << e << " ";
    //     }
    //     std::cout << std::endl;
    //     std::cout << "PAMVec: ";
    //     for (const auto &e : pamVec) {
    //         std::cout << e << " ";
    //     }
    //     std::cout << std::endl;
    //     std::cout << "UtilFunctions::logSumExp(res): " << UtilFunctions::logSumExp(res) << std::endl;
    //     throw std::runtime_error("Transmission probability is NaN");
    // }

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


float Chain::calc_new_likelihood() {
    // All observation likelihoods are independent of each other
    auto observation_llik = observation_llik_new.full_sum();

    // std::cout << "Observation likelihood: " << observation_llik << std::endl;

    // Transmission likelihoods are independent across loci
    auto transmission_llik = (transmission_llik_new.parallel_sum() + coi_prior_new)
        // weight by population responsibility
        .parallel_element_add(population_responsibility_vector.parallel_log())
        // logsumexp over populations
        .parallel_logsumexp()
        // sum over samples
        .parallel_full_sum();

    return observation_llik + transmission_llik;
}

float Chain::calc_new_prior() {
    return eps_neg_prior_new.full_sum() + 
                eps_pos_prior_new.full_sum() + 
                relatedness_prior_new.full_sum() + 
                population_coi_p_hyper_prior_new + 
                population_coi_r_hyper_prior_new + 
                population_responsibility_vector_prior_new;
}

float Chain::get_llik() { return llik; }
float Chain::get_prior() { return prior; }
float Chain::get_posterior() { return llik * temp + prior; }

// float Chain::get_llik(int sample) { 
//     int idx = sample * genotyping_data.num_loci;

//     #ifdef HAS_EXECUTION
//         return std::reduce(std::execution::unseq, genotyping_llik_new.begin() + idx,
//                            genotyping_llik_new.begin() + idx + genotyping_data.num_loci);
//     #else
//         float llik = 0.0f;
//         #pragma omp simd reduction(+:llik)
//         for (int i = 0; i < genotyping_data.num_loci; ++i) {
//             llik += genotyping_llik_new[idx + i];
//         }
//         return llik;
//     #endif
// }
// float Chain::get_prior(int sample) { 
//     float prior = 0.0f;
//     if (params.allow_relatedness) {
//         prior += relatedness_prior_new[sample];
//     }
//     prior += coi_prior_new[sample];
//     prior += eps_neg_prior_new[sample];
//     prior += eps_pos_prior_new[sample];
//     return prior;
// }

// float Chain::get_posterior(int sample) { 
//     float posterior = get_llik(sample) * temp + get_prior(sample);
//     return posterior;
// }

// void Chain::calculate_genotype_likelihood(std::size_t sample_idx, std::size_t locus_idx)
// {
//     // const std::size_t idx = sample_idx * genotyping_data.num_loci + locus_idx;
//     if (genotyping_data.is_missing(locus_idx, sample_idx))
//     {
//         genotyping_llik_new[idx] = 0;
//     }
//     else
//     {
//         auto [latent_genotypes_begin, latent_genotypes_end] = latent_genotypes_new.inner_iterators({sample_idx, locus_idx});
//         auto [p_begin, p_end] = p.inner_iterators({locus_idx});
//         const float transmission_prob = calc_transmission_process(
//             std::span(latent_genotypes_begin, latent_genotypes_end), 
//             std::span(p_begin, p_end), m[sample_idx], r[sample_idx]
//         );
//         const float obs_prob = calc_observation_process(
//             std::span(latent_genotypes_begin, latent_genotypes_end),
//             genotyping_data.get_observed_alleles(locus_idx, sample_idx),
//             eps_neg[sample_idx], eps_pos[sample_idx]
//         );
//         genotyping_llik_new[idx] = transmission_prob + obs_prob;
//     }
// }

void Chain::calculate_observation_likelihood(std::size_t sample_idx, std::size_t locus_idx) {
    if (genotyping_data.is_missing(sample_idx, locus_idx)) {
        observation_llik_new.at({sample_idx, locus_idx}) = 0;
    } else {
        auto [latent_genotypes_begin, latent_genotypes_end] = latent_genotypes_new.inner_iterators({sample_idx, locus_idx});
        const float obs_prob = calc_observation_process(
            std::span(latent_genotypes_begin, latent_genotypes_end),
            genotyping_data.get_observed_alleles(sample_idx, locus_idx),
            eps_neg.at({sample_idx}), eps_pos.at({sample_idx})
        );
        observation_llik_new.at({sample_idx, locus_idx}) = obs_prob;
    }
}

void Chain::calculate_transmission_likelihood(std::size_t sample_idx, std::size_t population_idx, std::size_t locus_idx) {
    if (genotyping_data.is_missing(sample_idx, locus_idx)) {
        transmission_llik_new.at({sample_idx, population_idx, locus_idx}) = 0;
    } else {
        const auto [latent_genotypes_begin, latent_genotypes_end] = latent_genotypes_new.inner_iterators({sample_idx, locus_idx});
        const auto [p_begin, p_end] = p.inner_iterators({population_idx, locus_idx});
        const auto latent_genotypes = std::span(latent_genotypes_begin, latent_genotypes_end);
        const auto p = std::span(p_begin, p_end);
        float transmission_prob = calc_transmission_process(
            latent_genotypes,
            p, m.at({sample_idx}), r.at({sample_idx})
        );

        // if (std::isnan(transmission_prob)) {
        //     std::cout << "Transmission probability is NaN" << std::endl;
        //     std::cout << "Latent genotypes: ";
        //     for (const auto &e : latent_genotypes) {
        //         std::cout << e << " ";
        //     }
        //     std::cout << std::endl;
        //     std::cout << "P: ";
        //     for (const auto &e : p) {
        //         std::cout << e << " ";
        //     }
        //     std::cout << std::endl;
        //     std::cout << "M: " << m.at({sample_idx}) << std::endl;
        //     std::cout << "R: " << r.at({sample_idx}) << std::endl;
        //     throw std::runtime_error("Transmission probability is NaN");
        // }

        transmission_llik_new.at({sample_idx, population_idx, locus_idx}) = transmission_prob;
    }
}

void Chain::calculate_eps_neg_likelihood(std::size_t sample_idx)
{
    eps_neg_prior_new.at({sample_idx}) = sampler.get_beta_log_prior(
        eps_neg.at({sample_idx}), params.eps_neg_alpha, params.eps_neg_beta);
}

void Chain::calculate_eps_pos_likelihood(std::size_t sample_idx)
{
    eps_pos_prior_new.at({sample_idx}) = sampler.get_beta_log_prior(
        eps_pos.at({sample_idx}), params.eps_pos_alpha, params.eps_pos_beta);
}

void Chain::calculate_relatedness_likelihood(std::size_t sample_idx)
{
    relatedness_prior_new.at({sample_idx}) = sampler.get_relatedness_log_prior(
        r.at({sample_idx}), params.r_alpha, params.r_beta);
}

void Chain::calculate_coi_likelihood(std::size_t sample_idx)
{
    coi_prior_new.at({sample_idx}) =
        sampler.get_coi_log_prior(m.at({sample_idx}), population_coi_p, population_coi_r);
}

void Chain::calculate_population_coi_p_likelihood()
{
    population_coi_p_hyper_prior_new = sampler.get_beta_log_prior(
        population_coi_p, params.population_coi_p_alpha, params.population_coi_p_beta);
}

void Chain::calculate_population_coi_r_likelihood()
{
    population_coi_r_hyper_prior_new = sampler.get_gamma_log_prior(
        population_coi_r, params.population_coi_r_shape, params.population_coi_r_rate);
}

void Chain::calculate_population_responsibility_vector_likelihood()
{
    population_responsibility_vector_prior_new = sampler.unnormalized_dirichlet_log_prior(
        population_responsibility_vector, params.population_responsibility_vector_alpha
    );
}

void Chain::save_population_responsibility_vector_likelihood()
{
    population_responsibility_vector_prior_old = population_responsibility_vector_prior_new;
}

void Chain::initialize_likelihood()
{
    const std::size_t num_samples = genotyping_data.num_samples;
    const std::size_t num_loci = genotyping_data.num_loci;
    const std::size_t num_populations = params.num_populations;

    observation_llik_new.resize({num_samples, num_loci});
    observation_llik_old.resize({num_samples, num_loci});

    transmission_llik_new.resize({num_samples, num_populations, num_loci});
    transmission_llik_old.resize({num_samples, num_populations, num_loci});

    eps_neg_prior_new.resize({num_samples});
    eps_neg_prior_old.resize({num_samples});

    eps_pos_prior_new.resize({num_samples});
    eps_pos_prior_old.resize({num_samples});

    relatedness_prior_new.resize({num_samples});
    relatedness_prior_old.resize({num_samples});
    
    coi_prior_new.resize({num_samples, num_populations});
    coi_prior_old.resize({num_samples, num_populations});

    // population_coi_mean_hyper_prior_new.resize({num_populations});
    // population_coi_mean_hyper_prior_old.resize({num_populations});

    // Calculate likelihoods of sample level parameters
    for (std::size_t sample_idx = 0; sample_idx < genotyping_data.num_samples; ++sample_idx) {
        calculate_eps_neg_likelihood(sample_idx);
        save_eps_neg_likelihood(sample_idx);
        
        calculate_eps_pos_likelihood(sample_idx);
        save_eps_pos_likelihood(sample_idx);
        
        calculate_relatedness_likelihood(sample_idx);
        save_relatedness_likelihood(sample_idx);
        
        for (std::size_t locus_idx = 0; locus_idx < genotyping_data.num_loci; ++locus_idx) {
            calculate_observation_likelihood(sample_idx, locus_idx);
            save_observation_likelihood(sample_idx, locus_idx);
        }
    }

    // Calculate likelihoods of sample specific population varying parameters
    for (std::size_t sample_idx = 0; sample_idx < genotyping_data.num_samples; ++sample_idx)
    {
        calculate_coi_likelihood(sample_idx);
        save_coi_likelihood(sample_idx);
        for (std::size_t pop_idx = 0; pop_idx < params.num_populations; ++pop_idx) {
            for (std::size_t locus_idx = 0; locus_idx < genotyping_data.num_loci; ++locus_idx)
            {
                calculate_transmission_likelihood(sample_idx, pop_idx, locus_idx);
                save_transmission_likelihood(sample_idx, pop_idx, locus_idx);
            }
        }
    }

    // Calculate likelihoods of population mean coi
    calculate_population_coi_p_likelihood();
    save_population_coi_p_likelihood();
    calculate_population_coi_r_likelihood();
    save_population_coi_r_likelihood();

    // Calculate likelihood of population responsibility vector
    calculate_population_responsibility_vector_likelihood();
    save_population_responsibility_vector_likelihood();
    llik = calc_new_likelihood();
    prior = calc_new_prior();
}

void Chain::save_observation_likelihood(std::size_t sample_idx, std::size_t locus_idx)
{
    observation_llik_old.at({sample_idx, locus_idx}) = observation_llik_new.at({sample_idx, locus_idx});
}

void Chain::save_transmission_likelihood(std::size_t sample_idx, std::size_t population_idx, std::size_t locus_idx)
{
    transmission_llik_old.at({sample_idx, population_idx, locus_idx}) = transmission_llik_new.at({sample_idx, population_idx, locus_idx});
}

void Chain::save_eps_neg_likelihood(std::size_t sample_idx)
{
    eps_neg_prior_old.at({sample_idx}) = eps_neg_prior_new.at({sample_idx});
}

void Chain::save_eps_pos_likelihood(std::size_t sample_idx)
{
    eps_pos_prior_old.at({sample_idx}) = eps_pos_prior_new.at({sample_idx});
}

void Chain::save_relatedness_likelihood(std::size_t sample_idx)
{
    relatedness_prior_old.at({sample_idx}) = relatedness_prior_new.at({sample_idx});
}

void Chain::save_coi_likelihood(std::size_t sample_idx)
{
    coi_prior_old.at({sample_idx}) = coi_prior_new.at({sample_idx});
}

void Chain::save_population_coi_p_likelihood()
{
    population_coi_p_hyper_prior_old = population_coi_p_hyper_prior_new;
}

void Chain::save_population_coi_r_likelihood()
{
    population_coi_r_hyper_prior_old = population_coi_r_hyper_prior_new;
}

void Chain::restore_observation_likelihood(std::size_t sample_idx, std::size_t locus_idx)
{
    observation_llik_new.at({sample_idx, locus_idx}) = observation_llik_old.at({sample_idx, locus_idx});
}

void Chain::restore_transmission_likelihood(std::size_t sample_idx, std::size_t population_idx, std::size_t locus_idx)
{
    transmission_llik_new.at({sample_idx, population_idx, locus_idx}) = transmission_llik_old.at({sample_idx, population_idx, locus_idx});
}

void Chain::restore_eps_neg_likelihood(std::size_t sample_idx)
{
    eps_neg_prior_new.at({sample_idx}) = eps_neg_prior_old.at({sample_idx});
}

void Chain::restore_eps_pos_likelihood(std::size_t sample_idx)
{
    eps_pos_prior_new.at({sample_idx}) = eps_pos_prior_old.at({sample_idx});
}

void Chain::restore_relatedness_likelihood(std::size_t sample_idx)
{
    relatedness_prior_new.at({sample_idx}) = relatedness_prior_old.at({sample_idx});
}

void Chain::restore_coi_likelihood(std::size_t sample_idx)
{
    coi_prior_new.at({sample_idx}) = coi_prior_old.at({sample_idx});
}

void Chain::restore_population_coi_p_likelihood()
{
    population_coi_p_hyper_prior_new = population_coi_p_hyper_prior_old;
}

void Chain::restore_population_coi_r_likelihood()
{
    population_coi_r_hyper_prior_new = population_coi_r_hyper_prior_old;
}

void Chain::restore_population_responsibility_vector_likelihood()
{
    population_responsibility_vector_prior_new = population_responsibility_vector_prior_old;
}

void Chain::set_llik(float llik) { this->llik = llik; }

void Chain::set_temp(float temp) { this->temp = temp; }

float Chain::get_temp() { return this->temp; }

Chain::Chain(GenotypingData genotyping_data, Parameters params, float temp)
    : genotyping_data(genotyping_data), params(params), temp(temp), llik(std::numeric_limits<float>::lowest()), sampler()

{
    initialize_parameters();
}
