#pragma once

#ifndef CHAIN_H_
#define CHAIN_H_

#include "combination_indices_generator.h"
#include "genotyping_data.h"
#include "parameters.h"
#include "prob_any_missing.h"
#include "sampler.h"
#include "multivector.h"

#include <Rcpp.h>


class Chain
{
   private:
    GenotypingData genotyping_data;
    Parameters params;
    Sampler sampler;
    probAnyMissingFunctor probAnyMissing_;
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

    float calc_transmission_process(
        std::span<int const> allele_index_vec,
        std::span<float const> allele_frequencies, int coi,
        float relatedness);

    float calc_observation_process(std::span<int const> allele_index_vec,
                                   std::span<int const> obs_genotype,
                                   float epsilon_neg, float epsilon_pos);

    float calculate_llik(int num_samples);
    float calc_new_likelihood();
    float calc_new_prior();
    float calc_new_posterior();

    void calculate_observation_likelihood(std::size_t sample_idx, std::size_t locus_idx);
    void calculate_transmission_likelihood(std::size_t sample_idx, std::size_t population_idx, std::size_t locus_idx);
    void calculate_eps_neg_likelihood(std::size_t sample_idx);
    void calculate_eps_pos_likelihood(std::size_t sample_idx);
    void calculate_coi_likelihood(std::size_t sample_idx);
    void calculate_relatedness_likelihood(std::size_t sample_idx);
    void calculate_population_coi_p_likelihood();
    void calculate_population_coi_r_likelihood();
    void calculate_population_responsibility_vector_likelihood();

    void save_observation_likelihood(std::size_t sample_idx, std::size_t locus_idx);
    void save_transmission_likelihood(std::size_t sample_idx, std::size_t population_idx, std::size_t locus_idx);
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


    // // Population COI
    // // indexed by population
    // MultiVector<float, 1> population_mean_coi{};
    // // Population COI variance
    // // indexed by population
    // MultiVector<float, 1> population_mean_coi_var{};
    // // Population COI acceptance
    // // indexed by population
    // MultiVector<float, 1> population_mean_coi_accept{};
    // // Population COI hyper prior
    // // indexed by population
    // MultiVector<float, 1> population_mean_coi_hyper_prior_old{};
    // // Population COI hyper prior
    // // indexed by population
    // MultiVector<float, 1> population_mean_coi_hyper_prior_new{};
    // // Population COI variance hyper prior

    // Population responsibility vector
    // indexed by population
    MultiVector<float, 1> population_responsibility_vector{};
    
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
    Chain() {};
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
    float get_llik(int sample);
    float get_prior(int sample);
    float get_posterior(int sample);

    void set_llik(float llik);
    void set_temp(float temp);
    float get_temp();
    void set_hot(bool hot) { this->hot = hot; };
    bool is_hot() { return hot; };
};

#endif  // CHAIN_H_
