
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
    Parameters params(args);
    GenotypingData genotyping_data(args);

    if (params.verbose && !params.simple_verbose)
    {
        UtilFunctions::print("-- Starting MCMC --");
        UtilFunctions::print("Total Burnin:", params.burnin);
        UtilFunctions::print("Total Samples:", params.samples);
        UtilFunctions::print("Thinning:", params.thin);
        UtilFunctions::print("Allow Relatedness:",
                             params.allow_relatedness ? "Yes" : "No");
        UtilFunctions::print("Parallel Tempering:",
                             params.pt_chains.size() > 1 ? "Yes" : "No");
        UtilFunctions::print("Adapt Temperature:",
                             params.adapt_temp ? "Yes" : "No");
    }

    MCMC mcmc(genotyping_data, params);
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
            pb.set_llik(mcmc.get_posterior());
            pb.set_hot_chain(mcmc.get_hot_chain());
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
            pb.set_llik(mcmc.get_posterior());
            pb.set_hot_chain(mcmc.get_hot_chain());
            p.increment();
        }
    }

    mcmc.finalize();

    Rcpp::List acceptance_rates;
    Rcpp::List sampling_variances;

    for (const auto &chain : mcmc.chains)
    {
        Rcpp::List chain_acceptance_rates;
        chain_acceptance_rates.push_back(Rcpp::wrap(chain.p_accept));
        chain_acceptance_rates.push_back(Rcpp::wrap(chain.m_accept));
        chain_acceptance_rates.push_back(Rcpp::wrap(chain.eps_neg_accept));
        chain_acceptance_rates.push_back(Rcpp::wrap(chain.eps_pos_accept));
        chain_acceptance_rates.push_back(Rcpp::wrap(chain.r_accept));
        chain_acceptance_rates.push_back(Rcpp::wrap(chain.m_r_accept));
        chain_acceptance_rates.push_back(Rcpp::wrap(chain.sample_accept));
        chain_acceptance_rates.push_back(Rcpp::wrap(chain.mean_coi_accept));

        Rcpp::StringVector acceptance_rate_names;
        acceptance_rate_names.push_back("allele_freq_accept");
        acceptance_rate_names.push_back("coi_accept");
        acceptance_rate_names.push_back("eps_neg_accept");
        acceptance_rate_names.push_back("eps_pos_accept");
        acceptance_rate_names.push_back("r_accept");
        acceptance_rate_names.push_back("m_r_accept");
        acceptance_rate_names.push_back("full_sample_accept");
        acceptance_rate_names.push_back("mean_coi_accept");
        chain_acceptance_rates.names() = acceptance_rate_names;
        acceptance_rates.push_back(chain_acceptance_rates);

        Rcpp::List chain_sampling_variances;
        chain_sampling_variances.push_back(Rcpp::wrap(chain.p_prop_var));
        chain_sampling_variances.push_back(Rcpp::wrap(chain.eps_neg_var));
        chain_sampling_variances.push_back(Rcpp::wrap(chain.eps_pos_var));
        chain_sampling_variances.push_back(Rcpp::wrap(chain.r_var));
        chain_sampling_variances.push_back(Rcpp::wrap(chain.m_r_var));
        chain_sampling_variances.push_back(Rcpp::wrap(chain.mean_coi_var));

        Rcpp::StringVector sampling_variance_names{
            "allele_freq_var", "eps_neg_var", "eps_pos_var",
            "r_var",           "m_r_var",     "mean_coi_var"};
        chain_sampling_variances.names() = sampling_variance_names;
        sampling_variances.push_back(chain_sampling_variances);
    }

    Rcpp::List res;
    res.push_back(Rcpp::wrap(mcmc.llik_burnin));
    res.push_back(Rcpp::wrap(mcmc.llik_sample));
    res.push_back(Rcpp::wrap(mcmc.prior_burnin));
    res.push_back(Rcpp::wrap(mcmc.prior_sample));
    res.push_back(Rcpp::wrap(mcmc.posterior_burnin));
    res.push_back(Rcpp::wrap(mcmc.posterior_sample));
    res.push_back(Rcpp::wrap(mcmc.m_store));
    res.push_back(Rcpp::wrap(mcmc.mean_coi_store));
    res.push_back(Rcpp::wrap(mcmc.p_store));
    res.push_back(Rcpp::wrap(mcmc.eps_neg_store));
    res.push_back(Rcpp::wrap(mcmc.eps_pos_store));
    res.push_back(Rcpp::wrap(mcmc.r_store));
    res.push_back(Rcpp::wrap(mcmc.genotyping_data.observed_coi));
    res.push_back(Rcpp::wrap(mcmc.swap_store));
    res.push_back(Rcpp::wrap(mcmc.swap_acceptances));
    res.push_back(Rcpp::wrap(mcmc.swap_barriers));
    res.push_back(Rcpp::wrap(mcmc.temp_gradient));
    res.push_back(Rcpp::wrap(acceptance_rates));
    res.push_back(Rcpp::wrap(sampling_variances));

    Rcpp::StringVector res_names;
    res_names.push_back("llik_burnin");
    res_names.push_back("llik_sample");
    res_names.push_back("prior_burnin");
    res_names.push_back("prior_sample");
    res_names.push_back("posterior_burnin");
    res_names.push_back("posterior_sample");
    res_names.push_back("coi");
    res_names.push_back("lam_coi");
    res_names.push_back("allele_freqs");
    res_names.push_back("eps_neg");
    res_names.push_back("eps_pos");
    res_names.push_back("relatedness");
    res_names.push_back("observed_coi");
    res_names.push_back("swap_store");
    res_names.push_back("swap_acceptances");
    res_names.push_back("swap_barriers");
    res_names.push_back("temp_gradient");
    res_names.push_back("acceptance_rates");
    res_names.push_back("sampling_variances");

    res.names() = res_names;
    return res;
}
