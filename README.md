# pactR — Posterior Analysis of Coalescent Trees

An R package wrapping the [PACT](https://github.com/trvrb/pact) C++ tool by Trevor Bedford, extended with a suite of R functions for analysing and visualising BEAST DTA (Discrete Trait Analysis) posterior tree files.

## Overview

The package has two layers:

**Core PACT** (C++ via Rcpp) — reads posterior NEWICK tree files from BEAST, Migrate, IM, or LAMARC and computes structured coalescent statistics: TMRCA, coalescent rates, migration rates, nucleotide diversity, FST, Tajima's D, lineage persistence, spatial diffusion coefficients, and time-windowed (skyline) versions of all statistics.

**DTA Extensions** (pure R) — reads BEAST NEXUS multi-tree posterior files and computes: lineage persistence by state (matching PACT's `summary persistence`), proportions of lineages over time, lineage-through-time curves, root-to-tip regression, pairwise TMRCA matrices, and MRCA lookups. Visualization via `ggplot2`; optional tree plotting via `ggtree`.

## Installation

### From source (tarball)

```r
install.packages(
  "/path/to/pact_0.9.4.tar.gz",
  repos = NULL, type = "source"
)
```

### From GitHub

```r
# install.packages("remotes")
remotes::install_github("ldamodaran/pactR")
```

### Dependencies

Required (installed automatically):

```r
install.packages(c("Rcpp", "ape"))
```

Recommended for visualization:

```r
install.packages(c("ggplot2", "scales"))

# Tree plotting (ggtree / treeio) — Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("ggtree", "treeio"))
```

A C++11-capable compiler is required (standard on macOS with Xcode Command Line Tools, and on Linux with `build-essential`).

## Quick start

### Core PACT analysis

```r
library(pactR)

# Use bundled example
ex <- pact_example_files()

# Run analysis — supply a param file, R arguments, or both
results <- run_pact(
  trees_file  = ex$trees,
  param_file  = ex$param,
  output_dir  = tempdir()
)

# Read outputs into data frames
stats    <- read_pact_stats(results$stats)
skylines <- read_pact_skylines(results$skylines)
tips     <- read_pact_tips(results$tips)
```

### DTA posterior analysis

```r
library(ggplot2)

f <- pact_example_dta()   # bundled RSV-A BEAST posterior

# Lineage persistence (dot-and-whisker)
pers <- persistence_dta(f, trait = "region", burnin = 5L)
plot_persistence(pers)

# Proportions over time (ribbon or stacked area)
prop <- dta_proportions(f, trait = "region", step = 0.25, burnin = 5L)
plot_dta_proportions(prop, type = "area")

# Lineage through time (skyline)
ltt <- lineage_through_time(f, step = 0.25, burnin = 5L)
plot_ltt(ltt)
```

If tip names do not contain parseable dates, supply the most recent sampling date:

```r
pers <- persistence_dta(f, mrsd = "2023-06-15", burnin = 5L)
```

## Key functions

| Function | Purpose |
|---|---|
| `run_pact()` | Run full PACT C++ analysis |
| `build_pact_params()` | Build param file content from R arguments |
| `read_pact_stats()` / `read_pact_skylines()` | Parse PACT outputs into data frames |
| `persistence_dta()` | Persistence time per state (PACT-matching algorithm) |
| `plot_persistence()` | Dot-and-whisker plot of persistence |
| `dta_proportions()` | Proportion of lineages in each state over time |
| `plot_dta_proportions()` | Ribbon / stacked area / line plot |
| `lineage_through_time()` | Number of lineages over time (posterior mean + CI) |
| `plot_ltt()` | Lineage-through-time plot |
| `root_to_tip()` | Molecular clock regression |
| `pairwise_tmrca()` | Pairwise TMRCA matrix |
| `find_mrca()` | MRCA of a tip set |
| `plot_tree_beast()` | Annotated tree visualization via ggtree |

See `MANUAL.md` for full parameter documentation and examples.

## Example data

The bundled RSV-A example (`pact_example_dta()`) contains 45 posterior trees with 2103 tips and a `region` discrete trait (12 states, 1956–2023). A ready-to-run analysis script is at:

```
inst/extdata/example_dta_trees/analyse_rsva.R
```

Run it with:

```bash
Rscript analyse_rsva.R
```

## Documentation

- `MANUAL.md` — full function reference and usage guide
- `pact_parameter_reference()` — complete PACT `in.param` command reference
- `pact_manual()` — opens the original C++ PACT PDF user manual
- `vignette("pact-intro")` — introductory vignette

## License

GPL-3. Original C++ PACT by Trevor Bedford.
