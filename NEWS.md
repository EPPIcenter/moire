# moire 3.2.0
This is a minor revision that fixes bugs in the adaptive temperature gradient approach and changes the default priors on false positive and false negative rates.

- New default priors
- Fixed bug in adaptive temperature gradient approach
- Numerical stability improvements
- New summarization functions   
    - `calculate_med_allele_freqs()` calculates median allele frequencies
    - `plot_chain_swaps()` creates diagnostic plots for chain swaps when using parallel tempering

# moire 3.1.0
This is a minor revision that introduces new functionality to improve mixing of the MCMC. This also updates the required version of R to 4.0.0. and C++ to C++17.

- Added option to specify prior on within-host relatedness
- Added adaptive approach to tune temperature gradient used during parallel tempering
- Various bugfixes and improvements

# moire 3.0.0
This is a major revision to moire, introducing a simplified API and functionality to infer within host relatedness and effective MOI. This release also introduces a parallel tempering based approach that leverages OpenMP, greatly improving mixing of the MCMC.

# moire 2.2.0
- Added support for running multiple chains simultaneously, then pooling output
- fix bug with missing data
- various other bugfixes and improvements

# moire 2.1.0
- Implemented a new error model that removes sensitivity to total number of alleles at a locus
- Removed option to marginalize out latent genotypes below some complexity threshold
- bug fixes, documentation improvements

# moire 2.0.1
- Minor bugfix when sampling that caused computational slow down

# moire 2.0.0
- Underlying model no longer uses pseudo marginal MH algorithm. Instead, model is augmented with latent genotypes above some user defined complexity level which are then sampled. Loci below the user defined complexity level have the latent state fully marginalized out.
- Removed dependency on GSL
- various bugfixes and speedups

# moire 1.1.1

- fixed overflow bug that would dramatically increase computational costs
- fixed bug in sampling complexity of infection where accepted updates weren't being recorded in some cases

# moire 1.1.0

- Added several new common functions for analyzing data
- Added functions to import data in common formats to the format required
- Made errors an independent parameter across samples
- Changed error model to no longer depend on underlying number of strains contributing alleles
- Removed multiple chain implementations. If multiple chains are desired, use multiprocessing
- Added progress bar for duration of MCMC

# moire 1.0.0

- Initial release of moire.
