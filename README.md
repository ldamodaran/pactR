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

### Core PACT example

The bundled influenza H3N2 example (`pact_example_files()`) is a single structured coalescent tree with 7 geographic states and a time range of 2002–2008. Run the full analysis with:

```r
library(pactR)

ex <- pact_example_files()

results <- run_pact(
  trees_file = ex$trees,
  param_file = ex$param,
  output_dir = tempdir()
)

stats    <- read_pact_stats(results$stats)
skylines <- read_pact_skylines(results$skylines)
tips     <- read_pact_tips(results$tips)
```

A self-contained script that runs the analysis and saves TSV files and figures is at:

```
inst/extdata/example/analyse_example.R
```

Run it with:

```bash
Rscript analyse_example.R
```

### DTA posterior analysis

```r
library(ggplot2)

f <- pact_example_dta()   # bundled global HPAI H5N1 BEAST posterior

# Lineage persistence (dot-and-whisker)
pers <- persistence_dta(f, trait = "geo", burnin = 3L)
plot_persistence(pers)

# Proportions over time (ribbon or stacked area)
prop <- dta_proportions(f, trait = "geo", step = 0.05, burnin = 3L)
plot_dta_proportions(prop, type = "area")

# Lineage through time (skyline)
ltt <- lineage_through_time(f, step = 0.05, burnin = 3L)
plot_ltt(ltt)
```

If tip names do not contain parseable dates, supply the most recent sampling date:

```r
pers <- persistence_dta(f, trait = "geo", mrsd = "2023-06-15", burnin = 3L)
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

### Core PACT (`pact_example_files()`)

Influenza H3N2 structured coalescent tree — 1 posterior tree, 7 geographic states, 2002–2008. Produces TMRCA, coalescent rates, pairwise migration rates, skyline proportions, and time-to-trunk statistics.

### DTA posterior (`pact_example_dta()`)

Global HPAI H5N1 BEAST DTA posterior — 31 trees, 1921 tips, `geo` trait (Asia, Europe, North_America), 2018–2023. Ready-to-run scripts in each folder write both TSV data files and PNG figures.

## Documentation

- `MANUAL.md` — full function reference and usage guide
- `pact_parameter_reference()` — complete PACT `in.param` command reference
- `pact_manual()` — opens the original C++ PACT PDF user manual
- `vignette("pact-intro")` — introductory vignette

## License

GPL-3. Original C++ PACT by Trevor Bedford.
