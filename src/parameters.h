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
    bool allow_relatedness;
    int chain_number;
    int thin;
    int burnin;
    int samples;
    int pt_chains;
    double pt_grad;
    int pt_num_threads;
    double allele_freq_threshold;

    // Model Parameters
    // Complexity of Infection
    int max_coi;
    double mean_coi_var;
    double mean_coi_prior_shape;
    double mean_coi_prior_scale;

    // False Positive Rate
    double eps_pos_0;      // Initial eps pos
    double max_eps_pos;    // Max allowed value
    double eps_pos_var;    // Variance of sampler
    double eps_pos_alpha;  // Alpha parameter prior on beta distribution
    double eps_pos_beta;   // Beta parameter prior on beta distribution

    // False Negative Rate
    double eps_neg_0;      // Initial eps neg
    double max_eps_neg;    // Max allowed value
    double eps_neg_var;    // Variance of sampler
    double eps_neg_alpha;  // Alpha parameter prior on beta distribution
    double eps_neg_beta;   // Beta parameter prior on beta distribution

    // double allele_freq_var;
    std::vector<double> allele_freq_vars{};
    bool adapt_allele_freq_vars;

    // constructors
    Parameters(const Rcpp::List &args);
};

#endif  // PARAMETERS_H_
