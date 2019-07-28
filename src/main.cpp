
#include "main.h"
#include "mcmc_utils.h"
#include "genotyping_data.h"
#include "parameters.h"
#include "lookup.h"

//----------------------------------------------
// [[Rcpp::export(name='run_mcmc_rcpp')]]
Rcpp::List run_mcmc(Rcpp::List args) {
    UtilFunctions::print("Starting Run MCMC");
    UtilFunctions::print("Loading Parameters...");
    Parameters params(args);
    UtilFunctions::print("Loading Genotyping Data...");
    GenotypingData genotyping_data(args);
    UtilFunctions::print("Generating Lookup Tables...");
    Lookup lookup(params.max_coi, genotyping_data.max_alleles);

    MCMC mcmc(genotyping_data, lookup, params);

    UtilFunctions::print("Running Chains:", params.num_chains);

    if(params.burnin > 0) {
        mcmc.burnin();
    }
    
    if(params.samples > 0) {
        mcmc.sample();
    }
    

    Rcpp::List res;
    res.push_back(Rcpp::wrap(mcmc.llik_burnin));
    res.push_back(Rcpp::wrap(mcmc.llik_sample));
    res.push_back(Rcpp::wrap(mcmc.m_store));
    res.push_back(Rcpp::wrap(mcmc.p_store));
    res.push_back(Rcpp::wrap(mcmc.eps_neg_store));
    res.push_back(Rcpp::wrap(mcmc.eps_pos_store));
    res.push_back(Rcpp::wrap(mcmc.mean_coi_store));
    res.push_back(Rcpp::wrap(mcmc.genotyping_data.observed_coi));

    Rcpp::StringVector res_names;
    res_names.push_back("loglike_burnin");
    res_names.push_back("loglike_sample");
    res_names.push_back("m_store");
    res_names.push_back("p_store");
    res_names.push_back("eps_neg_store");
    res_names.push_back("eps_pos_store");
    res_names.push_back("mean_coi_store");
    res_names.push_back("observed_coi");

    res.names() = res_names;
    return res;
}