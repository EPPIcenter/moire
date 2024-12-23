
<!-- badges: start -->

[![R-CMD-check](https://github.com/EPPIcenter/moire/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/EPPIcenter/moire/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

# moire <img src="man/figures/logo.svg" align="right"/>

`moire` is a package implementing an MCMC based approach to estimating
complexity of infection (COI), also sometimes referred to as
multiplicity of infection (MOI), population allele frequencies, and
within-host relatedness from polyallelic genomics data.

## Installation

`moire` can be installed either using our r-universe repository
(preferred)

``` r
# Install from r-universe
install.packages("moire", repos = c("https://eppicenter.r-universe.dev", "https://cloud.r-project.org"))
```

or from GitHub using `remotes`

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

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10092403.svg)](https://doi.org/10.5281/zenodo.10092403)

The paper describing our method may be found
[here](https://doi.org/10.1093/bioinformatics/btae619)
