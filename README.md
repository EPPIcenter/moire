
# moire <img src="man/figures/logo.svg" align="right" alt="" height="139" />

`moire` is a package implementing an MCMC based approach to estimating
complexity of infection and population allele frequencies from
polyallelic genomics data.

## Installation

`moire` requires the GNU Scientific Library (GSL).

``` bash
# Linux
sudo apt-get install libgsl-dev
```

``` bash
# MacOS
brew install gsl
```

For Windows, follow instructions found
[here](https://stackoverflow.com/questions/26939683/linking-gsl-library-to-rcppgsl-on-windows-machine)

After installing GSL, install moire as follows

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
mcmc_results <- moire::run_mcmc(data$data, data$sample_ids, data$loci, data$is_missing)
```
