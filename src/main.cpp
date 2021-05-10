
#include "main.h"

#include "genotyping_data.h"
#include "lookup.h"
#include "mcmc_progress_bar.h"
#include "mcmc_utils.h"
#include "parameters.h"

#include <progress.hpp>
//----------------------------------------------
// [[Rcpp::export(name='run_mcmc_rcpp')]]
Rcpp::List run_mcmc(Rcpp::List args)
{
    Parameters params(args);
    GenotypingData genotyping_data(args);
    Lookup lookup(params.max_coi, genotyping_data.max_alleles);

    if (params.verbose)
    {
        UtilFunctions::print("Starting MCMC");
        UtilFunctions::print("Total Burnin:", params.burnin);
        UtilFunctions::print("Total Samples:", params.samples);
        UtilFunctions::print("Thinning:", params.thin);
    }

    MCMC mcmc(genotyping_data, lookup, params);
    MCMCProgressBar pb(params.burnin, params.samples);
    Progress p(params.burnin + params.samples, params.verbose, pb);

    int step = 0;
    while (step < params.burnin)
    {
        Rcpp::checkUserInterrupt();
        mcmc.burnin(step);
        ++step;
        pb.set_llik(mcmc.get_llik());
        p.increment();
    }

    step = 0;
    while (step < params.samples)
    {
        Rcpp::checkUserInterrupt();
        mcmc.sample(step);
        ++step;
        pb.set_llik(mcmc.get_llik());
        p.increment();
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
