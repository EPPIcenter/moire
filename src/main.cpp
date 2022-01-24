
#include "main.h"

#include "genotyping_data.h"
#include "lookup.h"
#include "mcmc_progress_bar.h"
#include "mcmc_utils.h"
#include "parameters.h"

#include <Rcpp/utils/tinyformat.h>

#include <progress.hpp>
//----------------------------------------------
// [[Rcpp::export(name='run_mcmc_rcpp')]]
Rcpp::List run_mcmc(Rcpp::List args)
{
    Parameters params(args);
    GenotypingData genotyping_data(args);
    Lookup lookup(genotyping_data.max_alleles);

    if (params.verbose)
    {
        UtilFunctions::print("Starting MCMC");
        UtilFunctions::print("Total Burnin:", params.burnin);
        UtilFunctions::print("Total Samples:", params.samples);
        UtilFunctions::print("Thinning:", params.thin);
    }

    MCMC mcmc(genotyping_data, lookup, params);
    if (std::isnan(mcmc.get_llik()))
    {
        Rcpp::stop("Error: Initial Llik is NaN");
    }

    MCMCProgressBar pb(params.burnin, params.samples);
    Progress p(params.burnin + params.samples, params.verbose, pb);
    pb.set_llik(mcmc.get_llik());

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

    Rcpp::List acceptance_rates;
    acceptance_rates.push_back(Rcpp::wrap(mcmc.chain.p_accept));
    acceptance_rates.push_back(Rcpp::wrap(mcmc.chain.m_accept));
    acceptance_rates.push_back(Rcpp::wrap(mcmc.chain.eps_neg_accept));
    acceptance_rates.push_back(Rcpp::wrap(mcmc.chain.eps_pos_accept));

    Rcpp::StringVector acceptance_rate_names;
    acceptance_rate_names.push_back("allele_freq_accept");
    acceptance_rate_names.push_back("coi_accept");
    acceptance_rate_names.push_back("eps_neg_accept");
    acceptance_rate_names.push_back("eps_pos_accept");
    acceptance_rates.names() = acceptance_rate_names;

    Rcpp::List sampling_variances;
    sampling_variances.push_back(Rcpp::wrap(mcmc.chain.p_prop_var));
    sampling_variances.push_back(Rcpp::wrap(mcmc.chain.eps_neg_var));
    sampling_variances.push_back(Rcpp::wrap(mcmc.chain.eps_pos_var));

    Rcpp::StringVector sampling_variance_names;
    sampling_variance_names.push_back("allele_freq_var");
    sampling_variance_names.push_back("eps_neg_var");
    sampling_variance_names.push_back("eps_pos_var");
    sampling_variances.names() = sampling_variance_names;

    Rcpp::List res;
    res.push_back(Rcpp::wrap(mcmc.llik_burnin));
    res.push_back(Rcpp::wrap(mcmc.llik_sample));
    res.push_back(Rcpp::wrap(mcmc.m_store));
    res.push_back(Rcpp::wrap(mcmc.p_store));
    res.push_back(Rcpp::wrap(mcmc.eps_neg_store));
    res.push_back(Rcpp::wrap(mcmc.eps_pos_store));
    res.push_back(Rcpp::wrap(mcmc.genotyping_data.observed_coi));
    res.push_back(Rcpp::wrap(acceptance_rates));
    res.push_back(Rcpp::wrap(sampling_variances));

    Rcpp::StringVector res_names;
    res_names.push_back("llik_burnin");
    res_names.push_back("llik_sample");
    res_names.push_back("coi");
    res_names.push_back("allele_freqs");
    res_names.push_back("eps_neg");
    res_names.push_back("eps_pos");
    res_names.push_back("observed_coi");
    res_names.push_back("acceptance_rates");
    res_names.push_back("sampling_variances");

    res.names() = res_names;
    return res;
}
