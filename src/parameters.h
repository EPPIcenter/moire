#ifndef PARAMETERS_H_
#define PARAMETERS_H_

#include <Rcpp.h>
#include <vector>

class Parameters
{
   public:
    // MCMC Parameters
    bool verbose;
    bool simple_verbose;
    bool use_message;
    bool allow_relatedness;
    int chain_number;
    int thin;
    int burnin;
    int samples;
    std::vector<float> pt_chains;
    int pt_num_threads;
    bool adapt_temp;
    int pre_adapt_steps;
    int temp_adapt_steps;
    int max_initialization_tries;
    float max_runtime;

    bool record_latent_genotypes;

    // Model Parameters
    // Complexity of Infection
    std::size_t max_coi;
    float population_coi_p_alpha;
    float population_coi_p_beta;
    float population_coi_r_shape;
    float population_coi_r_rate;

    // False Positive Rate
    float max_eps_pos;    // Max allowed value
    float eps_pos_alpha;  // Alpha parameter prior on beta distribution
    float eps_pos_beta;   // Beta parameter prior on beta distribution

    // False Negative Rate
    float max_eps_neg;    // Max allowed value
    float eps_neg_alpha;  // Alpha parameter prior on beta distribution
    float eps_neg_beta;   // Beta parameter prior on beta distribution

    // Relatedness
    float r_alpha;
    float r_beta;

    // number of populations
    std::size_t num_populations;
    std::vector<float> population_responsibility_vector_alpha;
    
    // Initial allele frequencies from R (optional)
    bool use_initial_allele_frequencies;
    Rcpp::List initial_allele_frequencies;

    // constructors
    Parameters(const Rcpp::List &args);
    Parameters(){};
};

#endif  // PARAMETERS_H_
