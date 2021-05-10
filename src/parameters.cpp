
#include "parameters.h"

#include "mcmc_utils.h"

Parameters::Parameters(const Rcpp::List &args)
{
    // MCMC
    verbose = UtilFunctions::r_to_bool(args["verbose"]);
    thin = UtilFunctions::r_to_int(args["thin"]);
    burnin = UtilFunctions::r_to_int(args["burnin"]);
    samples = UtilFunctions::r_to_int(args["samples"]);
    importance_sampling_depth =
        UtilFunctions::r_to_int(args["importance_sampling_depth"]);
    importance_sampling_scaling_factor =
        UtilFunctions::r_to_double(args["importance_sampling_scaling_factor"]);

    // Model
    mean_coi = UtilFunctions::r_to_int(args["mean_coi"]);
    eps_pos_0 = UtilFunctions::r_to_double(args["eps_pos_0"]);
    max_eps_pos = UtilFunctions::r_to_double(args["max_eps_pos"]);
    eps_pos_alpha = UtilFunctions::r_to_double(args["eps_pos_alpha"]);
    eps_pos_beta = UtilFunctions::r_to_double(args["eps_pos_beta"]);
    eps_pos_var = UtilFunctions::r_to_double(args["eps_pos_var"]);
    eps_neg_0 = UtilFunctions::r_to_double(args["eps_neg_0"]);
    max_eps_neg = UtilFunctions::r_to_double(args["max_eps_neg"]);
    eps_neg_alpha = UtilFunctions::r_to_double(args["eps_neg_alpha"]);
    eps_neg_beta = UtilFunctions::r_to_double(args["eps_neg_beta"]);
    eps_neg_var = UtilFunctions::r_to_double(args["eps_neg_var"]);
    allele_freq_var = UtilFunctions::r_to_double(args["allele_freq_var"]);
};
