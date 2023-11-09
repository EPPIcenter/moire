
# moire <img src="man/figures/moire_web_logo.svg" align="right" alt="" width="120"/>

`moire` is a package implementing an MCMC based approach to estimating
complexity of infection (COI), also sometimes referred to as
multiplicity of infection (MOI), and population allele frequencies from
polyallelic genomics data.

## Installation

``` r
# Install development version from Github
remotes::install_github("EPPIcenter/moire")
```

## Usage

moire supports loading data from either a long format `data.frame` using
`load_long_form_data()` or from a wide format `data.frame` using
`load_delimited_data()`.

``` r
df <- read.csv("your_data.csv")
data <- load_long_form_data(df)

# With data in appropriate format, run MCMC as follows
mcmc_results <- moire::run_mcmc(data, is_missing = data$is_missing)
```

## Manuscript
[![DOI](https://zenodo.org/badge/174280517.svg)](https://zenodo.org/doi/10.5281/zenodo.10092402)

The paper describing our method may be found [here](https://doi.org/10.1101/2023.10.03.560769)
