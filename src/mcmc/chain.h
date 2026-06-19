#pragma once

#include "observation_model_factory.h"

#include <cstdint>
#include <memory>
#include <span>
#include <string>
#include "combination_indices_generator.h"
#include "genotyping_data.h"
#include "parameters.h"
#include "sampler.h"
#include "multivector.h"

#include <Rcpp.h>


class Chain
{
   private:
    GenotypingData genotyping_data;
    Parameters params;
    std::unique_ptr<ObservationModel> observation_model_;
    Sampler sampler;
    bool hot = false;
    float temp;
    float llik;
    float prior;

    CombinationIndicesGenerator allele_index_generator_;

    void initialize_latent_genotypes();
    void initialize_population_responsibility();
    void initialize_p();
    void initialize_m();
    void initialize_eps_neg();
    void initialize_eps_pos();
    void initialize_r();
    void initialize_population_coi();
    void initialize_likelihood();

    /// Number of non-padding alleles in latent_genotypes_new[sample, locus].
    void refresh_latent_support_k(std::size_t sample_idx, std::size_t locus_idx);
    void assign_latent_genotype_new(std::size_t sample_idx, std::size_t locus_idx,
                                    std::span<const int> value);
    void restore_latent_genotype_new(std::size_t sample_idx, std::size_t locus_idx);
    std::span<const int> latent_allele_support(std::size_t sample_idx,
                                               std::size_t locus_idx) const;
    std::span<const int> latent_allele_support_old(std::size_t sample_idx,
                                                   std::size_t locus_idx) const;
    bool locus_tx_inputs_unchanged(std::size_t sample_idx,
                                   std::size_t locus_idx,
                                   int prev_coi,
                                   float prev_r) const;

    float calc_transmission_process(
        std::span<int const> allele_index_vec,
        std::span<float const> allele_frequencies, int coi,
        float relatedness);

    float calc_transmission_process_after_r_change(
        std::span<int const> allele_index_vec,
        std::span<float const> allele_frequencies,
        int coi,
        float relatedness);

    float calc_transmission_process_after_p_change(
        std::span<int const> allele_index_vec,
        std::span<float const> allele_frequencies,
        int coi,
        float relatedness);

    float calc_new_likelihood();
    float calc_new_prior();
    float calc_transmission_llik_sum();
    void invalidate_transmission_llik_cache();
    void rebuild_transmission_llik_cache();
    float apply_transmission_cell_change(
        std::size_t sample_idx, std::size_t pop_idx, float old_val, float new_val);
    /// Recompute tx_sample_logsumexp[sample] after coi_prior_new changed (loci sums unchanged).
    void refresh_sample_tx_after_coi_change(std::size_t sample_idx);
    /// Recompute all per-sample tx logsumexp terms (coi prior and/or population log weights changed).
    void refresh_all_samples_tx_logsumexp();
    /// Incrementally update transmission cells for one sample (all pop x loci).
    void recalculate_transmission_for_sample_incremental(std::size_t sample_idx,
                                                         int prev_coi,
                                                         float prev_r);
    void restore_transmission_for_sample_incremental(std::size_t sample_idx);

    void calculate_observation_likelihood(std::size_t sample_idx, std::size_t locus_idx);
    void sync_obs_sum_for_sample(std::size_t sample_idx);
    void calculate_transmission_likelihood(std::size_t population_idx, std::size_t sample_idx, std::size_t locus_idx);
    void calculate_transmission_likelihood_after_r_change(
        std::size_t population_idx,
        std::size_t sample_idx,
        std::size_t locus_idx);
    /// Recompute transmission_llik_new[:, pop, locus] after a p proposal (grouped PAM).
    void recalculate_transmission_at_locus_after_p_change(
        std::size_t population_idx,
        std::size_t locus_idx,
        std::span<const float> p_old);
    void calculate_eps_neg_likelihood(std::size_t sample_idx);
    void calculate_eps_pos_likelihood(std::size_t sample_idx);
    void calculate_coi_likelihood(std::size_t sample_idx);
    void calculate_relatedness_likelihood(std::size_t sample_idx);
    void calculate_population_coi_p_likelihood();
    void calculate_population_coi_r_likelihood();
    void calculate_population_responsibility_vector_likelihood();

    void save_observation_likelihood(std::size_t sample_idx, std::size_t locus_idx);
    void save_transmission_likelihood(std::size_t population_idx, std::size_t sample_idx, std::size_t locus_idx);
    void save_eps_neg_likelihood(std::size_t sample_idx);
    void save_eps_pos_likelihood(std::size_t sample_idx);
    void save_coi_likelihood(std::size_t sample_idx);
    void save_relatedness_likelihood(std::size_t sample_idx);
    void save_population_coi_p_likelihood();
    void save_population_coi_r_likelihood();
    void save_population_responsibility_vector_likelihood();

    void restore_observation_likelihood(std::size_t sample_idx, std::size_t locus_idx);
    void restore_transmission_likelihood(std::size_t sample_idx, std::size_t population_idx, std::size_t locus_idx);
    void restore_eps_neg_likelihood(std::size_t sample_idx);
    void restore_eps_pos_likelihood(std::size_t sample_idx);
    void restore_coi_likelihood(std::size_t sample_idx);
    void restore_relatedness_likelihood(std::size_t sample_idx);
    void restore_population_coi_p_likelihood();
    void restore_population_coi_r_likelihood();
    void restore_population_responsibility_vector_likelihood();

   public:
   
    // Transmission likelihood per sample
    // indexed by sample, population, locus
    MultiVector<float, 3> transmission_llik_old{};
    // Transmission likelihood per sample
    // indexed by sample, population, locus
    MultiVector<float, 3> transmission_llik_new{};

    // Observation likelihood per sample
    // indexed by sample, locus
    MultiVector<float, 2> observation_llik_old{};
    // Observation likelihood per sample
    // indexed by sample, locus
    MultiVector<float, 2> observation_llik_new{};

    // Epsilon negative prior
    // indexed by sample
    MultiVector<float, 1> eps_neg_prior_old{};
    // Epsilon negative prior
    // indexed by sample
    MultiVector<float, 1> eps_neg_prior_new{};

    // Epsilon positive prior
    // indexed by sample
    MultiVector<float, 1> eps_pos_prior_old{};
    // Epsilon positive prior
    // indexed by sample
    MultiVector<float, 1> eps_pos_prior_new{};

    // COI prior
    // indexed by sample, population
    MultiVector<float, 2> coi_prior_new{};
    // COI prior
    // indexed by sample, population
    MultiVector<float, 2> coi_prior_old{};

    // Relatedness prior
    // indexed by sample
    MultiVector<float, 1> relatedness_prior_new{};
    // Relatedness prior
    // indexed by sample
    MultiVector<float, 1> relatedness_prior_old{};



    // Latent Genotypes
    // indexed by sample, locus, allele
    RaggedMultiVector<int, 3> latent_genotypes_old{};
    // Latent Genotypes
    // indexed by sample, locus, allele
    RaggedMultiVector<int, 3> latent_genotypes_new{};
    // Latent Genotype Adjustment   
    // indexed by sample, locus
    MultiVector<float, 2> lg_adj_old{};
    // Latent Genotype Adjustment
    // indexed by sample, locus
    MultiVector<float, 2> lg_adj_new{};

    // Cached support size per (sample, locus); avoids std::find on every transmission eval.
    MultiVector<std::uint16_t, 2> latent_support_k_{};

    // COI ~ ZTNB(population_coi_mean, population_coi_variance)
    // indexed by sample
    MultiVector<int, 1> m{};
    // COI acceptance
    // indexed by sample
    MultiVector<int, 1> m_accept{};


    // Population COI parameters
    float population_coi_p;
    float population_coi_r;
    // Hyper parameters for population COI mean ~ Gamma(shape, rate)
    float population_coi_p_alpha;
    float population_coi_p_beta;
    // Hyper parameters for population COI variance ~ Gamma(shape, rate)
    float population_coi_r_shape;
    float population_coi_r_rate;

    float population_coi_p_sampling_variance;
    float population_coi_r_sampling_variance;
    float population_coi_p_accept;
    float population_coi_r_accept;
    float population_coi_p_hyper_prior_old;
    float population_coi_p_hyper_prior_new;
    float population_coi_r_hyper_prior_old;
    float population_coi_r_hyper_prior_new;

    // Population responsibility vector
    // indexed by population
    MultiVector<float, 1> population_responsibility_vector{};
    
    // Cached log of population responsibility vector (invalidated when vector changes)
    // This avoids recomputing log() on every calc_new_likelihood() call
    mutable MultiVector<float, 1> population_responsibility_vector_log_{};
    mutable bool population_responsibility_vector_log_valid_ = false;

    // Cached transmission log-likelihood decomposition for incremental updates.
    // tx_loci_sum_new[sample, pop] = sum_loci transmission_llik_new[sample, pop, :]
    // tx_sample_logsumexp_new[sample] = logsumexp_pop(loci_sum + coi_prior + pop_log)
    // tx_llik_sum_new = sum_sample tx_sample_logsumexp_new
    MultiVector<float, 2> tx_loci_sum_new{};
    MultiVector<float, 1> tx_sample_logsumexp_new{};
    float tx_llik_sum_new{0.f};
    bool tx_llik_cache_valid_{false};

    // Running scalar sums maintained incrementally to avoid O(N) / O(N*L)
    // reductions on every Metropolis proposal. Each tracks the sum of the
    // corresponding *_new array and is rebuilt from scratch at the end of
    // initialize_likelihood(). Observation updates happen in parallel over loci,
    // so its sum is resynced per-sample (sequentially) via
    // sync_obs_sum_for_sample(); the per-sample priors are updated in place
    // inside calculate_*/restore_* since those run sequentially.
    float obs_llik_sum_new_{0.f};
    std::vector<float> obs_row_sum_new_{};
    float eps_neg_prior_sum_new_{0.f};
    float eps_pos_prior_sum_new_{0.f};
    float relatedness_prior_sum_new_{0.f};

    // Reused buffers for update_p SALT proposals (sized to max alleles per locus).
    std::vector<float> update_p_prev_p_ws_{};

    // Population responsibility vector proposal variance
    // indexed by population
    MultiVector<float, 1> population_responsibility_vector_prop_var{};

    // Population responsibility vector acceptance
    // indexed by population
    MultiVector<float, 1> population_responsibility_vector_accept{};

    // Population responsibility vector attempt
    // indexed by population
    MultiVector<float, 1> population_responsibility_vector_attempt{};

    // Population responsibility vector prior
    float population_responsibility_vector_prior_old{};
    
    // Population responsibility vector prior
    float population_responsibility_vector_prior_new{};

    // Relatedness parameter
    // indexed by sample
    MultiVector<float, 1> r{};
    // Relatedness proposal acceptance
    // indexed by sample
    MultiVector<int, 1> r_accept{};
    // Relatedness proposal variance
    // indexed by sample
    MultiVector<float, 1> r_var{};

    // COI and relatedness acceptance
    // indexed by sample
    MultiVector<int, 1> m_r_accept{};
    // COI and relatedness proposal variance
    // indexed by sample
    MultiVector<float, 1> m_r_var{};

    // Allele Frequencies Parameter
    // indexed by population, locus, allele
    RaggedMultiVector<float, 3> p{};
    // Allele frequency proposal variance
    // indexed by population, locus, allele 
    RaggedMultiVector<float, 3> p_prop_var{};
    // Allele frequency acceptance
    // indexed by population, locus, allele
    RaggedMultiVector<int, 3> p_accept{};
    // Allele frequency attempt
    // indexed by population, locus, allele
    RaggedMultiVector<int, 3> p_attempt{};

    // Epsilon Positive Parameter
    // indexed by sample
    MultiVector<float, 1> eps_pos{};
    // Epsilon Positive acceptance
    // indexed by sample
    MultiVector<int, 1> eps_pos_accept{};
    // Epsilon Positive proposal variance
    // indexed by sample
    MultiVector<float, 1> eps_pos_var{};

    // Epsilon Negative Parameter
    // indexed by sample
    MultiVector<float, 1> eps_neg{};
    // Epsilon Negative acceptance
    // indexed by sample
    MultiVector<int, 1> eps_neg_accept{};
    // Epsilon Negative proposal variance
    // indexed by sample
    MultiVector<float, 1> eps_neg_var{};

    // Sample update acceptance
    // indexed by sample
    MultiVector<int, 1> sample_accept{};

    Chain(GenotypingData genotyping_data, Parameters params, float temp = 1.0);
    Chain()
        : observation_model_(make_observation_model(ObservationModelKind::Binary, 2.0f, 2.0f))
    {
    }
    void update_m(int iteration);
    void update_r(int iteration);
    void update_m_r(int iteration);
    void update_eff_coi(int iteration);
    void update_p(int iteration);
    void update_eps(int iteration);
    void update_eps_pos(int iteration);
    void update_eps_neg(int iteration);
    void update_samples(int iteration);
    void update_population_coi_p(int iteration);
    void update_population_coi_r(int iteration);
    void update_population_responsibility_vector(int iteration);
    void initialize_parameters();
    float get_llik();
    float get_prior();
    float get_posterior();

    void set_llik(float llik);
    void set_temp(float temp);
    float get_temp();
    void set_hot(bool hot) { this->hot = hot; };
    bool is_hot() { return hot; };
};
