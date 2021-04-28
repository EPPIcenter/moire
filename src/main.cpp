
#include "main.h"

#include "genotyping_data.h"
#include "lookup.h"
#include "mcmc_utils.h"
#include "parameters.h"

//----------------------------------------------
// [[Rcpp::export(name='run_mcmc_rcpp')]]
Rcpp::List run_mcmc(Rcpp::List args)
{
    Parameters params(args);
    GenotypingData genotyping_data(args);
    Lookup lookup(params.max_coi, genotyping_data.max_alleles);

    MCMC mcmc(genotyping_data, lookup, params);

    if (params.burnin > 0)
    {
        mcmc.burnin();
    }

    if (params.samples > 0)
    {
        mcmc.sample();
    }

    Rcpp::List res;
    res.push_back(Rcpp::wrap(mcmc.llik_burnin));
    res.push_back(Rcpp::wrap(mcmc.llik_sample));
    res.push_back(Rcpp::wrap(mcmc.m_store));
    res.push_back(Rcpp::wrap(mcmc.p_store));
    res.push_back(Rcpp::wrap(mcmc.eps_neg_store));
    res.push_back(Rcpp::wrap(mcmc.eps_pos_store));
    // res.push_back(Rcpp::wrap(mcmc.mean_coi_store));
    res.push_back(Rcpp::wrap(mcmc.genotyping_data.observed_coi));

    Rcpp::StringVector res_names;
    res_names.push_back("llik_burnin");
    res_names.push_back("llik_sample");
    res_names.push_back("coi");
    res_names.push_back("allele_freqs");
    res_names.push_back("eps_neg");
    res_names.push_back("eps_pos");
    // res_names.push_back("mean_coi");
    res_names.push_back("observed_coi");

    res.names() = res_names;
    return res;
}
