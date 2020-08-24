
# moire <img src="man/figures/logo.svg" align="right" alt="" height="139" />

moire is a package implementing an MCMC based approach to estimating
complexity of infection and population allele frequencies from
polyallelic genomics data. Details of the method may be found
[here](/refs/thesis.pdf).

## Installation

``` r
# Install development version from Github
remotes::install_github("m-murphy/moire")
```

## Usage

moire takes as input a list of lists of numeric vectors, where each list
element is a collection of observations across samples at a single
genetic locus. A single observation is a vector of 1’s and 0’s
indicating presence or absence of the allele unique to that index.
Consequently, each vector within a locus collection must be the same
length, the order of samples must be the same across locus collections,
and missing data is represented as a vector of 0’s.

``` r
# With data in appropriate format, run MCMC as follows
mcmc_results <- moire::run_mcmc(data)
```
