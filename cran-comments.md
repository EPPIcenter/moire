## R CMD check results

0 errors | 0 warnings | 0 notes

* This is a new submission.

## Possible NOTE: installed package size

The installed size is about 6.5Mb, of which `data` is 4.3Mb. The `data`
directory holds a precomputed `mcmc_results` object (already `xz` compressed)
so that the vignette and the `summarize_*()` and `plot_*()` examples can run
in well under a second each rather than the minutes a real MCMC run takes.

## Compiled code

The package uses OpenMP when available and falls back to serial execution
when it is not (for example on CRAN's macOS builders). All OpenMP flags come
from `$(SHLIB_OPENMP_CXXFLAGS)`. Threading is controlled by the
`pt_num_threads` argument to `run_mcmc()`, which defaults to 1; tests and
examples never use more than 2 cores.

`src/include/spline/spline.h` is a bundled header-only cubic spline library
by Tino Kluge, licensed GPL-2 or later and listed as `cph` in `Authors@R`.
