#pragma once

#ifndef PARAMETERS_H_
#define PARAMETERS_H_

#include <Rcpp.h>

class Parameters {

    public:
    
    // MCMC Parameters
    static int thin;
    static int burnin;
    static int samples;
    static int num_chains;
    static int importance_sampling_depth;

    // Model Parameters
    // Complexity of Infection
    static int max_coi; // Max allowed value
    static int coi_delta; // Max shift in COI during MCMC
    
    // False Positive Rate
    static double eps_pos_0; // Initial eps pos
    static double max_eps_pos; // Max allowed value
    static double eps_pos_var; // Variance of sampler
    
    // False Negative Rate
    static double eps_neg_0; // Initial eps neg
    static double max_eps_neg; // Max allowed value
    static double eps_neg_var; // Variance of sampler

    // Allele Frequencies
    static double alpha; // Dirichlet sampling variance (higher values create lower variance)


    // constructors
    Parameters() {};
    Parameters(const Rcpp::List &args);

};

#endif // PARAMETERS_H_
