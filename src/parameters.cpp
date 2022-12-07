
#include "parameters.h"

#include "mcmc_utils.h"

Parameters::Parameters(const Rcpp::List &args)
{
    // MCMC
    verbose = UtilFunctions::r_to_bool(args["verbose"]);
    simple_verbose = UtilFunctions::r_to_bool(args["simple_verbose"]);
    allow_relatedness = UtilFunctions::r_to_bool(args["allow_relatedness"]);
    chain_number = UtilFunctions::r_to_int(args["chain_number"]);
    thin = UtilFunctions::r_to_int(args["thin"]);
    burnin = UtilFunctions::r_to_int(args["burnin"]);
    samples = UtilFunctions::r_to_int(args["samples"]);
    allele_freq_threshold =
        UtilFunctions::r_to_double(args["allele_freq_threshold"]);
    pt_chains = UtilFunctions::r_to_int(args["pt_chains"]);
    pt_grad = UtilFunctions::r_to_double(args["pt_grad"]);
    pt_num_threads = UtilFunctions::r_to_int(args["pt_num_threads"]);

    // Model
    max_coi = UtilFunctions::r_to_int(args["max_coi"]);
    mean_coi_shape = UtilFunctions::r_to_double(args["mean_coi_shape"]);
    mean_coi_scale = UtilFunctions::r_to_double(args["mean_coi_scale"]);
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

    allele_freq_vars =
        UtilFunctions::r_to_vector_double(args["allele_freq_vars"]);
    adapt_allele_freq_vars =
        UtilFunctions::r_to_bool(args["adapt_allele_freq_vars"]);
};
