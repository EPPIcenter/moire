#pragma once

#ifndef CHAIN_H_
#define CHAIN_H_

#include "combination_indices_generator.h"
#include "genotyping_data.h"
#include "parameters.h"
#include "prob_any_missing.h"
#include "sampler.h"

#include <Rcpp.h>

class Chain
{
   private:
    GenotypingData genotyping_data;
    Parameters params;
    Sampler sampler;
    probAnyMissingFunctor probAnyMissing_;
    std::vector<float> prVec_{};

    CombinationIndicesGenerator allele_index_generator_;

    void initialize_latent_genotypes();
    void initialize_p();
    void initialize_m();
    void initialize_eps_neg();
    void initialize_eps_pos();
    void initialize_r();
    void initialize_likelihood();

    float calc_transmission_process(
        std::vector<int> const &allele_index_vec,
        std::vector<float> const &allele_frequencies, int coi,
        float relatedness);

    float calc_observation_process(std::vector<int> const &allele_index_vec,
                                   std::vector<int> const &obs_genotype,
                                   float epsilon_neg, float epsilon_pos);

    float calc_genotype_log_pmf(std::vector<int> const &allele_index_vec,
                                std::vector<int> const &obs_genotype,
                                float epsilon_pos, float epsilon_neg, int coi,
                                float relatedness,
                                std::vector<float> const &allele_frequencies);

    std::vector<float> calc_obs_genotype_lliks(
        std::vector<int> const &obs_genotype,
        std::vector<std::vector<int>> const &true_genotypes, float epsilon_neg,
        float epsilon_pos, int num_genotypes);

    float calculate_llik(int num_samples);
    float calc_old_likelihood();
    float calc_new_likelihood();
    float calc_old_prior();
    float calc_new_prior();
    float calculate_new_posterior();
    float calculate_old_posterior();

    void calculate_genotype_likelihood(int sample_idx, int locux_idx);
    void calculate_eps_neg_likelihood(int sample_idx);
    void calculate_eps_pos_likelihood(int sample_idx);
    void calculate_coi_likelihood(int sample_idx);
    void calculate_relatedness_likelihood(int sample_idx);
    void calculate_mean_coi_likelihood();

    void save_genotype_likelihood(int sample_idx, int locus_idx);
    void save_eps_neg_likelihood(int sample_idx);
    void save_eps_pos_likelihood(int sample_idx);
    void save_coi_likelihood(int sample_idx);
    void save_relatedness_likelihood(int sample_idx);
    void save_mean_coi_likelihood();

    void restore_genotype_likelihood(int sample_idx, int locus_idx);
    void restore_eps_neg_likelihood(int sample_idx);
    void restore_eps_pos_likelihood(int sample_idx);
    void restore_coi_likelihood(int sample_idx);
    void restore_relatedness_likelihood(int sample_idx);
    void restore_mean_coi_likelihood();

   public:
    std::vector<float> genotyping_llik_old{};
    std::vector<float> genotyping_llik_new{};

    std::vector<float> eps_neg_prior_old{};
    std::vector<float> eps_neg_prior_new{};
    std::vector<float> eps_pos_prior_old{};
    std::vector<float> eps_pos_prior_new{};
    std::vector<float> coi_prior_new{};
    std::vector<float> coi_prior_old{};
    std::vector<float> relatedness_prior_new{};
    std::vector<float> relatedness_prior_old{};

    float llik;
    float prior;
    float temp;

    // Latent Genotypes
    std::vector<std::vector<std::vector<int>>> latent_genotypes_old{};
    std::vector<std::vector<std::vector<int>>> latent_genotypes_new{};
    std::vector<std::vector<float>> lg_adj_old{};
    std::vector<std::vector<float>> lg_adj_new{};

    // COI
    std::vector<int> m{};
    std::vector<int> m_accept{};
    float mean_coi;
    float mean_coi_var;
    float mean_coi_accept;
    float mean_coi_hyper_prior_old;
    float mean_coi_hyper_prior_new;

    // Relatedness
    std::vector<float> r{};
    std::vector<int> r_accept{};
    std::vector<float> r_var{};

    std::vector<int> m_r_accept{};
    std::vector<float> m_r_var{};

    // Allele Frequencies
    std::vector<std::vector<float>> p{};
    std::vector<float> prop_p{};
    std::vector<std::vector<float>> p_prop_var{};
    std::vector<std::vector<int>> p_accept{};
    std::vector<std::vector<int>> p_attempt{};

    // Epsilon Positive
    // float eps_pos;
    std::vector<float> eps_pos{};
    std::vector<int> eps_pos_accept{};
    std::vector<float> eps_pos_var{};

    // Epsilon Negative
    // float eps_neg;
    std::vector<float> eps_neg{};
    std::vector<int> eps_neg_accept{};
    std::vector<float> eps_neg_var{};

    std::vector<int> sample_accept{};

    Chain(GenotypingData genotyping_data, Parameters params, float temp = 1.0);
    void update_m(int iteration);
    void update_r(int iteration);
    void update_m_r(int iteration);
    void update_eff_coi(int iteration);
    void update_p(int iteration);
    void update_eps(int iteration);
    void update_eps_pos(int iteration);
    void update_eps_neg(int iteration);
    void update_samples(int iteration);
    void update_mean_coi(int iteration);
    float get_llik();
    float get_prior();
    float get_posterior();

    void set_llik(float llik);
    void set_temp(float temp);
    float get_temp();
};

#endif  // CHAIN_H_
