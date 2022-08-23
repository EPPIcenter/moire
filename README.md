
# moire <img src="man/figures/logo.svg" align="right" alt=""/>

`moire` is a package implementing an MCMC based approach to estimating
complexity of infection (COI), also sometimes referred to as
multiplicity of infection (MOI), and population allele frequencies from
polyallelic genomics data.

## Installation

``` r
# Install development version from Github
remotes::install_github("m-murphy/moire")
```

## Usage

moire supports loading data from either a long format `data.frame` using
`load_long_form_data()` or from a wide format `data.frame` using
`load_delimited_data()`.

``` r
df <- read.csv("your_data.csv")
data <- load_long_form_data(df)

# With data in appropriate format, run MCMC as follows
mcmc_results <- moire::run_mcmc(
  data$data, data$sample_ids, data$loci, data$is_missing
)
```
