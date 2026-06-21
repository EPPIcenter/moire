#include "chain.h"

#include "prob_any_missing.h"
#include "multivector_fused.h"
#include "mcmc_utils.h"
#include "sampler.h"
#include "profiler.h"

#include <cmath>
#include <algorithm>

#include "env_defs.h"

#include <cstdlib>
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
    latent_support_k_.clear();
    latent_support_k_.resize(std::array{genotyping_data.num_samples, genotyping_data.num_loci});
    for (std::size_t sample_idx = 0; sample_idx < genotyping_data.num_samples; ++sample_idx)
    {
        for (std::size_t locus_idx = 0; locus_idx < genotyping_data.num_loci; ++locus_idx)
        {
            const auto obs_alleles = genotyping_data.get_observed_alleles(sample_idx, locus_idx);

            const auto lg = observation_model_->sample_latent_genotype(
                sampler, obs_alleles, m.at({sample_idx}), .1f, .1f);

            assign_latent_genotype_new(sample_idx, locus_idx, lg.value);
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
    // Initialize cache as invalid (will be computed on first calc_new_likelihood call)
    population_responsibility_vector_log_valid_ = false;
    tx_llik_cache_valid_ = false;

    population_responsibility_vector_prop_var.clear();
    population_responsibility_vector_prop_var.resize({params.num_populations}, 1);
    population_responsibility_vector_accept.clear();
    population_responsibility_vector_accept.resize({params.num_populations}, 0);
    population_responsibility_vector_attempt.clear();
    population_responsibility_vector_attempt.resize({params.num_populations}, 0);
}

// Initialize P with allele frequencies
void Chain::initialize_p()
{
    p.clear();
    p.resize({params.num_populations, genotyping_data.num_loci}, genotyping_data.num_alleles);
    p_prop_var.clear();
    p_accept.clear();
    p_attempt.clear();

    if (params.use_initial_allele_frequencies) {
        // Use pre-computed allele frequencies from R
        UtilFunctions::print("=== Using Initial Allele Frequencies from R ===");
        UtilFunctions::print("Number of populations:", params.num_populations);
        UtilFunctions::print("Number of loci:", genotyping_data.num_loci);
        
        for (size_t pop_idx = 0; pop_idx < params.num_populations; ++pop_idx) {
            Rcpp::List pop_frequencies = Rcpp::as<Rcpp::List>(params.initial_allele_frequencies[pop_idx]);
            UtilFunctions::print("Population", pop_idx + 1, "has", pop_frequencies.size(), "loci");
            if (pop_frequencies.size() != genotyping_data.num_loci) {
                Rcpp::stop(
                    "initial_allele_frequencies: population %d has %d loci but data has %d loci.",
                    static_cast<int>(pop_idx + 1),
                    static_cast<int>(pop_frequencies.size()),
                    static_cast<int>(genotyping_data.num_loci)
                );
            }
            for (size_t locus_idx = 0; locus_idx < genotyping_data.num_loci; ++locus_idx) {
                Rcpp::NumericVector locus_frequencies = Rcpp::as<Rcpp::NumericVector>(pop_frequencies[locus_idx]);
                std::vector<float> freq_vec(locus_frequencies.begin(), locus_frequencies.end());
                p.inner_fill({pop_idx, locus_idx}, freq_vec);
                
                // Log first few loci for verification
                if (locus_idx < 3) {
                    std::string freq_str = "[";
                    float sum_freq = 0.0f;
                    for (size_t i = 0; i < freq_vec.size(); ++i) {
                        if (i > 0) freq_str += ", ";
                        freq_str += std::to_string(freq_vec[i]);
                        sum_freq += freq_vec[i];
                    }
                    freq_str += "]";
                    UtilFunctions::print("  Locus", locus_idx + 1, ":", freq_str, "(sum:", sum_freq, ")");
                }
            }
        }
        UtilFunctions::print("=== Initial Allele Frequencies Loaded Successfully ===");
    } else {
        // Fall back to clustering-based initialization
        UtilFunctions::print("=== Using Clustering-Based Initialization ===");
        p = UtilFunctions::calculate_clustered_allele_frequencies<float>(genotyping_data, params.num_populations, sampler);
    }

    p_prop_var.resize({params.num_populations, genotyping_data.num_loci}, genotyping_data.num_alleles, 1);
    p_accept.resize({params.num_populations, genotyping_data.num_loci}, genotyping_data.num_alleles, 0);
    p_attempt.resize({params.num_populations, genotyping_data.num_loci}, genotyping_data.num_alleles, 0);

    std::size_t max_alleles_per_locus = 0;
    for (std::size_t locus_idx = 0; locus_idx < genotyping_data.num_loci; ++locus_idx) {
        max_alleles_per_locus =
            std::max(max_alleles_per_locus, p.ragged_dimensions(locus_idx));
    }
    update_p_prev_p_ws_.resize(max_alleles_per_locus);
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

    // eCOI mode pools the false-positive rate to one value per locus.
    if (params.marginal_ecoi)
    {
        const std::size_t L = genotyping_data.num_loci;
        eps_pos_locus.clear();
        eps_pos_locus_accept.clear();
        eps_pos_locus_var.clear();
        eps_pos_locus.resize({L});
        // Start each locus at the (anchored) prior mean so the chain begins at the
        // assay rate rather than an arbitrary high value it would have to walk down.
        const float prior_mean =
            params.eps_pos_locus_alpha /
            std::max(1e-6f, params.eps_pos_locus_alpha + params.eps_pos_locus_beta);
        for (std::size_t l = 0; l < L; ++l)
        {
            eps_pos_locus.at({l}) = prior_mean;
        }
        eps_pos_locus_accept.resize({L}, 0);
        eps_pos_locus_var.resize({L}, 1);
    }
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

    // Effective COI, initialized from the (m, r) starting point.
    eff_coi.clear();
    eff_coi_accept.clear();
    eff_coi_var.clear();
    eff_coi.resize({genotyping_data.num_samples});
    for (std::size_t i = 0; i < genotyping_data.num_samples; ++i)
    {
        eff_coi.at({i}) =
            (m.at({i}) - 1) * (1.0f - r.at({i})) + 1.0f;  // = eCOI(m, r)
    }
    eff_coi_accept.resize({genotyping_data.num_samples}, 0);
    eff_coi_var.resize({genotyping_data.num_samples}, 1);
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



void Chain::initialize_population_e()
{
    ecoi_support_cache_init();

    const std::size_t n = genotyping_data.num_samples;
    const std::size_t P = params.num_populations;

    // Data-driven start: place each sample's effective COI at the argmax of its
    // own transmission-only marginal profile (uniform relatedness, m integrated),
    // evaluated on the initial allele frequencies and observed supports. The
    // per-sample and per-allele-frequency conditionals are individually well
    // identified, but the joint (e, p, G) posterior is multimodal; from a cold
    // clustering start the chain otherwise settles into a partially-collapsed
    // mode (some high-eCOI samples explained as e~1) that the local moves cannot
    // escape. Seeding e at its data-supported MAP places the chain in the correct
    // basin. NOTE: empirically this alone does NOT fix the k underestimate (the
    // collapse is a stationary phenomenon, not a cold-start trap), so it is
    // opt-in via MOIRE_ECOI_SMART_INIT and the (m, r) start remains the default.
    const bool smart_init = std::getenv("MOIRE_ECOI_SMART_INIT") != nullptr;
    if (smart_init)
    {
        const double e_hi = std::max(1.1, static_cast<double>(params.max_coi) - 1e-2);
        for (std::size_t i = 0; i < n; ++i)
        {
            double best_e = static_cast<double>(eff_coi.at({i}));
            double best_val = -std::numeric_limits<double>::infinity();
            for (double e = 1.05; e <= e_hi; e += 0.25)
            {
                const double v = ecoi_sample_rel_transmission_at_e(i, e);
                if (std::isfinite(v) && v > best_val)
                {
                    best_val = v;
                    best_e = e;
                }
            }
            eff_coi.at({i}) = static_cast<float>(best_e);
        }
    }

    // Every sample's effective COI is continuous (e > 1); there is no atom. Seed
    // any sample initialized at e <= 1 (e.g. a starting COI of 1) just above 1 so
    // the continuous prior and the m-marginalization are well defined.
    double sum_excess = 0.0;
    for (std::size_t i = 0; i < n; ++i)
    {
        if (!(eff_coi.at({i}) > 1.0f))
        {
            eff_coi.at({i}) = 1.0f + 1e-2f;
        }
        sum_excess += static_cast<double>(eff_coi.at({i})) - 1.0;
    }

    // Seed the per-population e-prior parameters from the starting configuration,
    // clamped to sane interiors so f_p(e) and its hyperpriors are finite.
    const float mu0 = static_cast<float>(
        std::max(0.5, n > 0 ? sum_excess / double(n) : 3.0));
    ecoi_mu_plus.assign(P, mu0);
    ecoi_k.assign(P, 2.0f);
    ecoi_mu_log_sd = 0.15f;
    ecoi_k_log_sd = 0.2f;
    ecoi_pop_accept = 0;
    ecoi_pop_attempt = 0;

    recompute_ecoi_prior();
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
    if (params.marginal_ecoi)
    {
        initialize_population_e();
    }
    initialize_likelihood();
}
