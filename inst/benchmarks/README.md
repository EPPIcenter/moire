# MCMC baseline benchmarks

Repeatable end-to-end timing for `run_mcmc()`. Baselines are **machine-specific**
(wall-clock depends on core count and thread budget) but profiler hotspots are
useful for spotting regressions in hot paths across machines.

## Quick start

```bash
devtools::install()
just bench-mcmc small          # run once, print timings + profiler
just bench-mcmc-save all       # refresh committed baseline CSV
just bench-mcmc-compare all    # fail if >5% slower than baseline
```

## Parallel modes

| Mode | Flag / env | Chains | Threads | Use for |
|------|------------|--------|---------|---------|
| **single** | default | 1 | auto (cores − 1) | Default regression track; inner C++ parallelism |
| **pt** | `--pt` | `BENCH_PT_CHAINS` (20) | auto | Parallel-tempering throughput (vignette-style) |
| **serial** | `--serial` | 1 | 1 | Algorithmic changes without TBB noise |

Pin threads for reproducibility across machines:

```bash
BENCH_NUM_THREADS=8 Rscript inst/scripts/bench_mcmc_baseline.R vignette --save
```

## Presets

| Preset | Data |
|--------|------|
| `minimal` | 4 samples × 3 loci (long-form) |
| `small` | 20 × 10 simulated |
| `medium` | 40 × 50 simulated |
| `vignette` | 100 × 100 simulated (matches demo scale) |

Non-`single` modes append a suffix in the CSV (`vignette_pt`, `vignette_serial`).

## Environment variables

| Variable | Default | Description |
|----------|---------|-------------|
| `BENCH_PARALLEL_MODE` | `single` | `single`, `pt`, or `serial` |
| `BENCH_NUM_THREADS` | physical cores − 1 | TBB thread budget (`num_threads` in R) |
| `BENCH_PT_CHAINS` | 20 | PT replica count when mode is `pt` |
| `BENCH_BURNIN` | 200 | Burn-in iterations |
| `BENCH_SAMPLES` | 200 | Post-burnin samples per chain |
| `BENCH_REPS` | 3 | Repetitions per preset |
| `BENCH_SEED` | 42 | RNG seed |
| `BENCH_REGRESSION_THRESHOLD` | 5 | Max allowed wall-clock regression (%) |

## Comparing runs

`--compare` checks that parallel mode and thread counts match the saved baseline
(meta rows in the CSV). Mismatched configs print a warning; treat wall-clock
deltas as informational only.

For quick wall-clock without profiler overhead:

```bash
MOIRE_DISABLE_PROFILER_REGISTRY=1 Rscript inst/scripts/bench_wall_clock.R vignette
```

## Interpreting results

* **single** preset times track the common case (one chain, multi-threaded inner work).
* **pt** times measure tempered-chain scheduling; faster per-iteration is not expected
  vs single — PT runs more chain work per MCMC step.
* **serial** is best for isolating algorithmic speedups from parallelism changes.

After changing parallelism behavior, refresh baselines on a representative machine
and record `num_threads` / `physical_cores` in the CSV meta section.
