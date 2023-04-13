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
    std::vector<double> pt_chains;
    int pt_num_threads;

    // Model Parameters
    // Complexity of Infection
    int max_coi;
    double mean_coi_shape;
    double mean_coi_scale;

    // False Positive Rate
    double max_eps_pos;    // Max allowed value
    double eps_pos_alpha;  // Alpha parameter prior on beta distribution
    double eps_pos_beta;   // Beta parameter prior on beta distribution

    // False Negative Rate
    double max_eps_neg;    // Max allowed value
    double eps_neg_alpha;  // Alpha parameter prior on beta distribution
    double eps_neg_beta;   // Beta parameter prior on beta distribution

    // Relatedness
    double r_alpha;
    double r_beta;

    // double allele_freq_var;

    // constructors
    Parameters(const Rcpp::List &args);
};

#endif  // PARAMETERS_H_
