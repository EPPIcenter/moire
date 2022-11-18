
#include "main.h"

#include "genotyping_data.h"
#include "lookup.h"
#include "mcmc_progress_bar.h"
#include "mcmc_utils.h"
#include "parameters.h"

#include <cmath>

#include <Rcpp/utils/tinyformat.h>

#include <progress.hpp>
//----------------------------------------------
// [[Rcpp::export(name='run_mcmc_rcpp')]]
Rcpp::List run_mcmc(Rcpp::List args)
{
    UtilFunctions::print("Running Experimental Relatedness Estimationv4.");
    Parameters params(args);
    GenotypingData genotyping_data(args);
    int chain_number = args["chain_number"];

    UtilFunctions::print("Allow Relatedness:", params.allow_relatedness);

    if (params.verbose && !params.simple_verbose)
    {
        UtilFunctions::print("Starting MCMC -- Chain", params.chain_number);
        UtilFunctions::print("Total Burnin:", params.burnin);
        UtilFunctions::print("Total Samples:", params.samples);
        UtilFunctions::print("Thinning:", params.thin);
    }

    MCMC mcmc(genotyping_data, params);

    // Sometimes when initializing, the likelihood is too extreme and results in
    // a NaN
    bool ill_conditioned = std::isnan(mcmc.get_llik());
    int max_tries = 1000;
    while (ill_conditioned and max_tries != 0)
    {
        mcmc = MCMC(genotyping_data, params);
        max_tries--;
        ill_conditioned = std::isnan(mcmc.get_llik());
    }

    if (ill_conditioned)
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
        if (params.verbose && params.simple_verbose && (step % 100) == 0)
        {
            UtilFunctions::print("Chain", params.chain_number, ":", step, "/",
                                 params.burnin + params.samples, "completed");
        }
        else if (!params.simple_verbose)
        {
            pb.set_llik(mcmc.get_llik());
            p.increment();
        }
    }

    step = 0;
    while (step < params.samples)
    {
        Rcpp::checkUserInterrupt();
        mcmc.sample(step);
        ++step;
        if (params.verbose && params.simple_verbose && (step % 100) == 0)
        {
            UtilFunctions::print("Chain ", params.chain_number, ":", step, "/",
                                 params.burnin + params.samples, "completed");
        }
        else if (!params.simple_verbose)
        {
            pb.set_llik(mcmc.get_llik());
            p.increment();
        }
    }

    Rcpp::List acceptance_rates;
    acceptance_rates.push_back(Rcpp::wrap(mcmc.chain.p_accept));
    acceptance_rates.push_back(Rcpp::wrap(mcmc.chain.m_accept));
    acceptance_rates.push_back(Rcpp::wrap(mcmc.chain.eps_neg_accept));
    acceptance_rates.push_back(Rcpp::wrap(mcmc.chain.eps_pos_accept));
    acceptance_rates.push_back(Rcpp::wrap(mcmc.chain.r_accept));
    acceptance_rates.push_back(Rcpp::wrap(mcmc.chain.sample_accept));

    Rcpp::StringVector acceptance_rate_names;
    acceptance_rate_names.push_back("allele_freq_accept");
    acceptance_rate_names.push_back("coi_accept");
    acceptance_rate_names.push_back("eps_neg_accept");
    acceptance_rate_names.push_back("eps_pos_accept");
    acceptance_rate_names.push_back("r_accept");
    acceptance_rate_names.push_back("full_sample_accept");
    acceptance_rates.names() = acceptance_rate_names;

    Rcpp::List sampling_variances;
    sampling_variances.push_back(Rcpp::wrap(mcmc.chain.p_prop_var));
    sampling_variances.push_back(Rcpp::wrap(mcmc.chain.eps_neg_var));
    sampling_variances.push_back(Rcpp::wrap(mcmc.chain.eps_pos_var));
    sampling_variances.push_back(Rcpp::wrap(mcmc.chain.r_var));

    Rcpp::StringVector sampling_variance_names{"allele_freq_var", "eps_neg_var",
                                               "eps_pos_var", "r_var"};
    sampling_variances.names() = sampling_variance_names;

    Rcpp::List res;
    res.push_back(Rcpp::wrap(mcmc.llik_burnin));
    res.push_back(Rcpp::wrap(mcmc.llik_sample));
    res.push_back(Rcpp::wrap(mcmc.m_store));
    res.push_back(Rcpp::wrap(mcmc.p_store));
    res.push_back(Rcpp::wrap(mcmc.eps_neg_store));
    res.push_back(Rcpp::wrap(mcmc.eps_pos_store));
    res.push_back(Rcpp::wrap(mcmc.r_store));
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
    res_names.push_back("relatedness");
    res_names.push_back("observed_coi");
    res_names.push_back("acceptance_rates");
    res_names.push_back("sampling_variances");

    res.names() = res_names;
    return res;
}
