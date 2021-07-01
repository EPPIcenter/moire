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
