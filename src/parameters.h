#pragma once

#ifndef PARAMETERS_H_
#define PARAMETERS_H_

#include <Rcpp.h>

class Parameters {

    public:
    
    // MCMC Parameters
    int thin;
    int burnin;
    int samples;
    int num_chains;
    int importance_sampling_depth;

    // Model Parameters
    // Complexity of Infection
    int max_coi;            // Max allowed value
    int coi_delta;          // Max shift in COI during MCMC
    
    // False Positive Rate
    double eps_pos_0;       // Initial eps pos
    double max_eps_pos;     // Max allowed value
    double eps_pos_var;  // Variance of sampler
    double eps_pos_alpha;   // Alpha parameter prior on beta distribution
    double eps_pos_beta;    // Beta parameter prior on beta distribution
    
    // False Negative Rate
    double eps_neg_0; // Initial eps neg
    double max_eps_neg; // Max allowed value
    double eps_neg_var; // Variance of sampler
    double eps_neg_alpha; // Alpha parameter prior on beta distribution
    double eps_neg_beta; // Beta parameter prior on beta distribution
    

    // Allele Frequencies
    double alpha; // Dirichlet sampling variance (higher values create lower variance)


    // constructors
    Parameters() {};
    Parameters(const Rcpp::List &args);

};

#endif // PARAMETERS_H_
