
#include "parameters.h"
#include "mcmc_utils.h"

//----------------------------------------
// Declare static member variables

// MCMC
int Parameters::thin;
int Parameters::burnin;
int Parameters::samples;
int Parameters::num_chains;
int Parameters::importance_sampling_depth;

// Model
int Parameters::max_coi;
int Parameters::coi_delta;
double Parameters::eps_pos_0;
double Parameters::max_eps_pos;
double Parameters::eps_pos_var;
double Parameters::eps_neg_0;
double Parameters::max_eps_neg;
double Parameters::eps_neg_var;
double Parameters::alpha;


Parameters::Parameters(const Rcpp::List &args) {
    // MCMC
    thin                      = UtilFunctions::r_to_int(args["thin"]);
    burnin                    = UtilFunctions::r_to_int(args["burnin"]);
    samples                   = UtilFunctions::r_to_int(args["samples"]);
    num_chains                = UtilFunctions::r_to_int(args["num_chains"]);
    importance_sampling_depth = UtilFunctions::r_to_int(args["importance_sampling_depth"]);

    // Model
    max_coi     =   UtilFunctions::r_to_int(args["max_coi"]);
    coi_delta   =   UtilFunctions::r_to_int(args["max_coi_delta"]);
    eps_pos_0   =   UtilFunctions::r_to_double(args["eps_pos_0"]);
    max_eps_pos =   UtilFunctions::r_to_double(args["max_eps_pos"]);
    eps_pos_var =   UtilFunctions::r_to_double(args["eps_pos_var"]);
    eps_neg_0   =   UtilFunctions::r_to_double(args["eps_neg_0"]);
    max_eps_neg =   UtilFunctions::r_to_double(args["max_eps_neg"]);
    eps_neg_var =   UtilFunctions::r_to_double(args["eps_neg_var"]);
    alpha       =   UtilFunctions::r_to_double(args["alpha"]);
};

