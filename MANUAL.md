# pactR — User Manual

**Version 0.9.4**

---

## Table of contents

1. [Introduction](#1-introduction)
2. [Installation](#2-installation)
3. [Package architecture](#3-package-architecture)
4. [Core PACT analysis](#4-core-pact-analysis)
   - [run_pact()](#run_pact)
   - [build_pact_params()](#build_pact_params)
   - [Parameter reference](#parameter-reference)
   - [Reading outputs](#reading-outputs)
   - [Helper utilities](#helper-utilities)
5. [DTA posterior extensions](#5-dta-posterior-extensions)
   - [Input format](#input-format)
   - [read_beast_posterior()](#read_beast_posterior)
   - [parse_tip_dates()](#parse_tip_dates)
   - [Persistence analysis](#persistence-analysis)
   - [Proportions over time](#proportions-over-time)
   - [Lineage through time](#lineage-through-time)
   - [Tree utilities](#tree-utilities)
   - [Tree visualization](#tree-visualization)
6. [Bundled example data](#6-bundled-example-data)
7. [Complete workflow examples](#7-complete-workflow-examples)

---

## 1. Introduction

`pactR` is an R package providing two complementary layers for analysing posterior samples of phylogenetic trees from Bayesian coalescent inference.

**Layer 1 — Core PACT** wraps the original C++ [PACT](https://github.com/trvrb/pact) tool (Bedford 2009) via Rcpp. It reads posterior NEWICK tree files from BEAST, Migrate, IM, or LAMARC and computes structured coalescent statistics including TMRCA, coalescent rates, migration rates, nucleotide diversity, FST, Tajima's D, spatial diffusion coefficients, and lineage persistence. Time-windowed (skyline) versions of all statistics are supported.

**Layer 2 — DTA Extensions** provides pure-R functions for BEAST DTA (Discrete Trait Analysis) posterior tree files in NEXUS format. These cover the same analytical ground as the core engine but work directly with annotated BEAST output: persistence time per discrete state (matching PACT's `summary persistence` algorithm exactly), lineage proportions over time, lineage-through-time curves, root-to-tip regression, and pairwise TMRCA matrices.

---

## 2. Installation

### From a source tarball

```r
install.packages("/path/to/pact_0.9.4.tar.gz",
                 repos = NULL, type = "source")
```

### From GitHub

```r
# install.packages("remotes")
remotes::install_github("ldamodaran/pactR")
```

### System requirements

- R >= 3.5.0
- A C++11-capable compiler:
  - **macOS**: Xcode Command Line Tools (`xcode-select --install`)
  - **Linux**: `sudo apt install build-essential` (Debian/Ubuntu)
  - **Windows**: Rtools (https://cran.r-project.org/bin/windows/Rtools/)

### R dependencies

Installed automatically:

```r
install.packages(c("Rcpp", "ape"))
```

Recommended for visualization:

```r
install.packages(c("ggplot2", "scales"))
```

Optional tree visualization (Bioconductor):

```r
BiocManager::install(c("ggtree", "treeio"))
```

---

## 3. Package architecture

```
pactR/
├── src/                   C++ source (original PACT + Rcpp wrapper)
│   ├── rcpp_pact.cpp      Rcpp entry point
│   ├── coaltree.cpp/h     CoalescentTree class (all statistics)
│   ├── io.cpp/h           File I/O, modified to accept path arguments
│   └── param.cpp/h        Parameter parsing
├── R/
│   ├── pact.R             Core wrappers: run_pact(), build_pact_params(), read_*()
│   └── extensions.R       DTA extensions: persistence_dta(), dta_proportions(), etc.
└── inst/extdata/
    ├── example/           PACT structured coalescent example files
    ├── example_dta_trees/ RSV-A BEAST DTA posterior + analysis script
    └── parameters.txt     Full PACT parameter command reference
```

---

## 4. Core PACT analysis

### `run_pact()`

```r
run_pact(
  trees_file,
  param_file    = NULL,
  output_dir    = ".",
  output_prefix = "out",
  ...             # any argument from build_pact_params()
)
```

Runs the PACT C++ analysis engine on a posterior tree file.

**Three usage modes:**

```r
# 1. Param file only
run_pact(trees_file, param_file = "in.param", output_dir = "results/")

# 2. R arguments only
run_pact(trees_file,
         output_dir     = "results/",
         burnin         = 10,
         summary_tmrca  = TRUE,
         skyline_settings   = c(2002, 2008, 0.1),
         skyline_proportions = TRUE)

# 3. Both — param_file as base, R args appended (R args take precedence)
run_pact(trees_file, param_file = "in.param",
         output_dir    = "results/",
         summary_tmrca = TRUE)
```

**Returns** a named list of output file paths:

| Element | File | Contents |
|---|---|---|
| `$stats` | `out.stats` | Summary statistics (posterior mean + CI) |
| `$tips` | `out.tips` | Per-tip statistics |
| `$skylines` | `out.skylines` | Time-windowed statistics |
| `$rules` | `out.rules` | Mathematica-format tree (if `print_tree = TRUE`) |
| `$pairs` | `out.pairs` | Pairwise statistics |

---

### `build_pact_params()`

```r
build_pact_params(burnin = NULL, push_times_back = NULL, ...)
```

Converts named R arguments to PACT `in.param` format lines. Exported for inspection and file saving.

```r
lines <- build_pact_params(
  burnin              = 10,
  push_times_back     = c(2002, 2008),
  summary_tmrca       = TRUE,
  skyline_settings    = c(2002, 2008, 0.1),
  skyline_proportions = TRUE
)
cat(lines, sep = "\n")
# burnin 10
# push times back 2002 2008
# summary tmrca
# skyline settings 2002 2008 0.1
# skyline proportions

# Save to file
writeLines(lines, "my_analysis.param")
```

---

### Parameter reference

All parameters below can be passed either in an `in.param` file or as named arguments to `run_pact()` and `build_pact_params()`.

#### General

| R argument | param file | Description |
|---|---|---|
| `burnin = 10` | `burnin 10` | Discard first *n* trees |

#### Tree manipulation

| R argument | param file | Description |
|---|---|---|
| `push_times_back = c(2002, 2008)` | `push times back 2002 2008` | Rescale times to `[start, stop]` |
| `reduce_tips = 0.5` | `reduce tips 0.5` | Randomly retain fraction of tips |
| `renew_trunk = 1` | `renew trunk 1` | Redefine trunk from samples within *n* years of most-recent tip |
| `prune_to_trunk = TRUE` | `prune to trunk` | Keep trunk lineage only |
| `prune_to_label = "USA"` | `prune to label USA` | Keep only tips with this label |
| `prune_to_tips = c("t1","t2")` | `prune to tips t1 t2` | Keep only named tips |
| `remove_tips = c("t1")` | `remove tips t1` | Remove named tips |
| `prune_to_time = c(2000, 2010)` | `prune to time 2000 2010` | Keep branches in time window |
| `pad_migration_events = TRUE` | `pad migration events` | Add virtual migration events |
| `collapse_labels = TRUE` | `collapse labels` | Treat all tips as one population |
| `trim_ends = c(2000, 2010)` | `trim ends 2000 2010` | Retain only branches in window |
| `section_tree = c(2000, 1, 0.5)` | `section tree 2000 1 0.5` | Slice tree (start, window, step) |
| `time_slice = 2005` | `time slice 2005` | Ancestors alive at this time |
| `rotate = 45` | `rotate 45` | Rotate spatial coordinates (degrees) |
| `accumulate = TRUE` | `accumulate` | Accumulate spatial locations |
| `add_tail = 0.5` | `add tail 0.5` | Add tail branch of this length |
| `ordering = c("A","B")` | `ordering A B` | Set tip output ordering |

#### Tree output

| R argument | param file | Description |
|---|---|---|
| `print_tree = TRUE` | `print rule tree` | Highest-posterior tree → `out.rules` (Mathematica format) |
| `print_circular_tree = TRUE` | `print circular tree` | Circular-layout tree |
| `print_all_trees = TRUE` | `print all trees` | All trees → `trees/` directory |

#### Summary statistics

Each produces lower/mean/upper across the posterior in `out.stats`.

| R argument | param file | Description |
|---|---|---|
| `summary_tmrca = TRUE` | `summary tmrca` | Time to most recent common ancestor |
| `summary_length = TRUE` | `summary length` | Total tree branch length |
| `summary_root_proportions = TRUE` | `summary root proportions` | Label proportions at root |
| `summary_proportions = TRUE` | `summary proportions` | Label proportions on trunk |
| `summary_coal_rates = TRUE` | `summary coal rates` | Per-population coalescent rates |
| `summary_mig_rates = TRUE` | `summary mig rates` | Directional migration rates |
| `summary_sub_rates = TRUE` | `summary sub rates` | Mean substitution rate |
| `summary_diversity = TRUE` | `summary diversity` | Nucleotide diversity (π) |
| `summary_fst = TRUE` | `summary fst` | F_ST between populations |
| `summary_tajima_d = TRUE` | `summary tajima d` | Tajima's D |
| `summary_diffusion_coefficient = TRUE` | `summary diffusion coefficient` | Spatial diffusion coefficient |
| `summary_drift_rate = TRUE` | `summary drift rate` | Spatial drift rate |
| `summary_persistence = TRUE` | `summary persistence` | Mean time from tip to first ancestor with different label |

#### Skyline statistics

`skyline_settings` must be set first. All produce time-series output in `out.skylines`.

| R argument | param file | Description |
|---|---|---|
| `skyline_settings = c(2000, 2010, 0.5)` | `skyline settings 2000 2010 0.5` | Time window: start, end, bin-width |
| `skyline_tmrca = TRUE` | `skyline tmrca` | TMRCA over time |
| `skyline_length = TRUE` | `skyline length` | Branch length over time |
| `skyline_proportions = TRUE` | `skyline proportions` | Lineage proportions over time |
| `skyline_coal_rates = TRUE` | `skyline coal rates` | Coalescent rates over time |
| `skyline_mig_rates = TRUE` | `skyline mig rates` | Migration rates over time |
| `skyline_diversity = TRUE` | `skyline diversity` | Nucleotide diversity over time |
| `skyline_fst = TRUE` | `skyline fst` | F_ST over time |
| `skyline_tajima_d = TRUE` | `skyline tajima d` | Tajima's D over time |
| `skyline_timetofix = TRUE` | `skyline time to fix` | Time to fixation |
| `skyline_xmean = TRUE` | `skyline xmean` | Mean x-coordinate over time |
| `skyline_ymean = TRUE` | `skyline ymean` | Mean y-coordinate over time |

#### Tip statistics

| R argument | param file | Description |
|---|---|---|
| `tips_time_to_trunk = TRUE` | `tips time to trunk` | Each tip's time back to trunk → `out.tips` |
| `x_loc_history = c(2000, 2010, 0.5)` | `tips x loc history 2000 2010 0.5` | x-coordinate history per tip |
| `y_loc_history = c(2000, 2010, 0.5)` | `tips y loc history 2000 2010 0.5` | y-coordinate history per tip |
| `pairs_diversity = 1` | `pairs diversity 1` | Pairwise diversity window → `out.pairs` |

---

### Reading outputs

```r
stats    <- read_pact_stats(results$stats)
# Columns: statistic, lower, mean, upper

skylines <- read_pact_skylines(results$skylines)
# Columns: time, statistic, lower, mean, upper

tips     <- read_pact_tips(results$tips)
# Columns: tip, statistic, lower, mean, upper

pairs    <- read_pact_pairs(results$pairs)
# Columns: tip1, tip2, statistic, lower, mean, upper
```

Convert PACT skyline proportions to the common proportion format for `plot_dta_proportions()`:

```r
prop <- pact_to_proportions(skylines, prefix = "pro_")
# Columns: time, state, lower, mean, upper
plot_dta_proportions(prop)
```

---

### Helper utilities

```r
pact_example_files()        # list(trees=..., param=...) — bundled structured coalescent example
pact_parameter_reference()  # print full in.param command reference
pact_manual()               # open bundled PDF user manual
```

---

## 5. DTA posterior extensions

These functions work with **BEAST NEXUS** `.trees` files containing annotated posterior trees from a Discrete Trait Analysis. No conversion is needed — the BEAST output is read directly.

### Input format

A BEAST `.trees` file has this structure:

```text
#NEXUS
Begin taxa;
  Dimensions ntax=5;
  Taxlist
    1 Sample_A_USA_2014,
    2 Sample_B_UK_2015,
    ...;
End;
Begin trees;
  Translate
    1 Sample_A_USA_2014,
    2 Sample_B_UK_2015,
    ...;
  tree STATE_0 [&lnP=-1234.5] = [&R] (1[&region="USA"]:0.5[&rate=0.001],...)
  tree STATE_1 ...
End;
```

Key features:
- **Translate block**: maps short numeric IDs to full taxon names
- **Tree lines**: each `tree STATE_XXXX` line holds one posterior sample in NEWICK format
- **Tip annotations**: `tipnum[&region="STATE"]` immediately after the tip number
- **Internal node annotations**: `)[&region="STATE"]` after each closing parenthesis, in left-to-right order matching ape's internal node numbering

---

### `read_beast_posterior()`

```r
read_beast_posterior(
  file,
  burnin    = 0L,
  trait     = "region",
  tip_dates = NULL,
  mrsd      = NULL
)
```

Reads the full posterior, returning a list with one edge data frame per tree.

| Argument | Description |
|---|---|
| `file` | Path to BEAST `.trees` NEXUS file |
| `burnin` | Integer — discard first *n* trees |
| `trait` | Annotation key to extract (e.g. `"region"`, `"host"`) |
| `tip_dates` | Named numeric decimal-year dates (names = tip labels). Auto-parsed if `NULL` |
| `mrsd` | Most-recent sampling date for time anchoring. `"YYYY-MM-DD"` string or decimal year |

**Returns** a list:

| Element | Description |
|---|---|
| `$trees` | List of data frames. Columns: `parent`, `child`, `parent_time`, `child_time`, `region`, `branch_length` |
| `$translations` | Named character vector: number → taxon name |
| `$tip_dates` | Named numeric: decimal-year dates per tip |
| `$absolute_time` | `TRUE` if calendar-year times; `FALSE` if root-relative |

**Temporal anchoring** (priority order):

1. `mrsd` — pins the tallest tip to this date; all node times become calendar years
2. Auto-parsed tip dates — regex extracts 4-digit year from tip names; adds 0.5 if no decimal
3. Fallback — root-relative branch lengths (root = 0), matching original PACT's behaviour

```r
# Auto-parsed dates (default)
post <- read_beast_posterior("posterior.trees", trait = "region", burnin = 5L)

# Explicit most-recent date
post <- read_beast_posterior("posterior.trees", mrsd = "2023-06-15", burnin = 5L)

# No date anchoring — root-relative times
post <- read_beast_posterior("posterior.trees", tip_dates = NULL, mrsd = NULL)
```

---

### `parse_tip_dates()`

```r
parse_tip_dates(tips, pattern = "(\\d{4})(?:[._]\\d+)?$", decimal = TRUE)
```

Extracts decimal sampling dates from tip label strings.

| Argument | Description |
|---|---|
| `tips` | Character vector of tip labels |
| `pattern` | Regex with one capture group matching the date portion |
| `decimal` | If `TRUE`, adds 0.5 to whole-year matches (midpoint convention) |

**Returns** a named numeric vector (names = tip labels, values = decimal years, `NA` where no match).

```r
tips <- c("RSV_A_USA_2014", "EPI_ISL_123_2016.75", "Sample_NoDate")
parse_tip_dates(tips)
# RSV_A_USA_2014    EPI_ISL_123_2016.75    Sample_NoDate
#         2014.5               2016.750               NA
```

---

### Persistence analysis

#### `persistence_dta()`

```r
persistence_dta(
  file,
  trait     = "region",
  burnin    = 0L,
  tip_dates = NULL,
  mrsd      = NULL
)
```

Computes lineage persistence by state, matching PACT's `summary persistence` algorithm exactly.

**Algorithm:**
1. For each tip in the tree, walk back through ancestors until the first ancestor with a *different* trait state
2. Persistence = `tip_time − ancestor_time` (how long this lineage has been continuously in its current state)
3. Average across all tips with the same state → one value per state per tree
4. Summarise across the posterior → mean + 95% credible interval

**Returns** a data frame:

| Column | Description |
|---|---|
| `state` | Trait state label |
| `lower` | 2.5th percentile across trees (years) |
| `mean` | Mean across trees (years) |
| `upper` | 97.5th percentile across trees (years) |

```r
pers <- persistence_dta("posterior.trees", trait = "region", burnin = 5L)
#           state lower  mean upper
#  Central_Europe  2.15  2.30  2.43
#   North_America  1.76  1.82  1.87
#          Africa  1.69  1.73  1.78
# ...
```

#### `plot_persistence()`

```r
plot_persistence(
  data,
  palette     = NULL,
  title       = "Lineage persistence by state",
  xlab        = "Mean time in state (years)",
  sort_states = TRUE
)
```

Horizontal dot-and-whisker (pointrange) chart. States sorted by mean (longest at top).

| Argument | Description |
|---|---|
| `data` | Data frame from `persistence_dta()` — columns `state`, `mean`, `lower`, `upper` |
| `palette` | Named colour vector (state → colour). Defaults to Okabe-Ito |
| `sort_states` | Sort by mean persistence if `TRUE` |

**Returns** a `ggplot2` object.

```r
p <- plot_persistence(pers, title = "RSV-A Persistence by Region")
ggsave("persistence.png", p, width = 7, height = 5, dpi = 150)
```

---

### Proportions over time

#### `dta_proportions()`

```r
dta_proportions(
  file,
  trait     = "region",
  times     = NULL,
  step      = 0.1,
  burnin    = 0L,
  tip_dates = NULL,
  mrsd      = NULL
)
```

At each time point, counts what fraction of co-circulating lineages are in each state, averaged across the posterior. This is the DTA equivalent of PACT's `skyline proportions`.

| Argument | Description |
|---|---|
| `times` | Numeric vector of evaluation time points. Auto-generated from data if `NULL` |
| `step` | Time step for the auto-generated grid (years) |

**Returns** a data frame:

| Column | Description |
|---|---|
| `time` | Evaluation time (decimal year or root-relative) |
| `state` | Trait state |
| `lower` | 2.5th percentile proportion across trees |
| `mean` | Mean proportion |
| `upper` | 97.5th percentile |
| `n_lineages` | Mean number of lineages at this time point |

#### `pact_to_proportions()`

```r
pact_to_proportions(skylines, prefix = "pro_")
```

Converts `read_pact_skylines()` output to the same format as `dta_proportions()` so both can be passed to `plot_dta_proportions()`.

#### `plot_dta_proportions()`

```r
plot_dta_proportions(
  data,
  type         = c("ribbon", "area", "line"),
  palette      = NULL,
  title        = "Lineage proportions over time",
  xlab         = "Time",
  ylab         = "Proportion of lineages",
  ribbon_alpha = 0.2,
  smooth       = TRUE,
  smooth_span  = 0.15
)
```

| Argument | Description |
|---|---|
| `type` | `"ribbon"` — lines + 95% CI shading (default); `"area"` — stacked filled area; `"line"` — lines only |
| `smooth` | Apply LOESS smoothing to reduce jaggedness |
| `smooth_span` | LOESS bandwidth (lower = less smoothing) |

**Returns** a `ggplot2` object.

```r
prop <- dta_proportions("posterior.trees", step = 0.25, burnin = 5L)

plot_dta_proportions(prop, type = "ribbon")   # lines + CI
plot_dta_proportions(prop, type = "area")     # stacked
plot_dta_proportions(prop, type = "line")     # no CI
```

---

### Lineage through time

#### `lineage_through_time()`

```r
lineage_through_time(
  file,
  times     = NULL,
  step      = 0.1,
  burnin    = 0L,
  tip_dates = NULL,
  mrsd      = NULL
)
```

Counts the number of phylogenetic lineages alive at each time point, averaged across the posterior.

**Returns** a data frame: `time`, `lower`, `mean`, `upper`.

#### `plot_ltt()`

```r
plot_ltt(
  ltt_data,
  log_y        = FALSE,
  colour       = "steelblue",
  ribbon_alpha = 0.25,
  title        = "Lineage through time"
)
```

| Argument | Description |
|---|---|
| `log_y` | Log10 y-axis |
| `ribbon_alpha` | Transparency of the 95% CI ribbon |

**Returns** a `ggplot2` object.

```r
ltt <- lineage_through_time("posterior.trees", step = 0.25, burnin = 5L)
plot_ltt(ltt, log_y = TRUE)
```

---

### Tree utilities

#### `root_to_tip()`

```r
root_to_tip(
  file,
  tip_dates  = NULL,
  tree_index = -1L,
  plot       = TRUE
)
```

Molecular clock diagnostic. Regresses root-to-tip distances on sampling dates and prints R², slope, and intercept. Uses the last tree by default (`tree_index = -1L`).

**Returns** a data frame: `tip`, `date`, `root_to_tip`. When `plot = TRUE` a scatter with regression line is printed and the data frame is returned invisibly.

#### `pairwise_tmrca()`

```r
pairwise_tmrca(file, tree_index = -1L, tip_dates = NULL)
```

Returns a symmetric matrix of absolute pairwise TMRCAs for all tips. Rows and columns labelled by tip name.

#### `find_mrca()`

```r
find_mrca(file, tips, tree_index = -1L)
```

Finds the MRCA of a specified set of tip labels.

**Returns** a list:

| Element | Description |
|---|---|
| `$node` | ape node index |
| `$n_descendants` | Number of tips under the MRCA |
| `$subtree_tips` | Character vector of those tip labels |

```r
mrca <- find_mrca("posterior.trees",
                   tips = c("Sample_A_USA_2014", "Sample_B_USA_2015"))
cat("MRCA node:", mrca$node, " | Descendants:", mrca$n_descendants, "\n")
```

---

### Tree visualization

#### `plot_tree_beast()`

```r
plot_tree_beast(
  file,
  trait      = NULL,
  tip_labels = FALSE,
  tip_size   = 2,
  palette    = NULL,
  layout     = "rectangular",
  mrsd       = NULL,
  ...
)
```

Plots a BEAST-annotated tree using `ggtree`, optionally colouring branches by a discrete trait. Reads via `treeio::read.beast()`.

**Requires** `ggtree` and `treeio` (Bioconductor):

```r
BiocManager::install(c("ggtree", "treeio"))
```

| Argument | Description |
|---|---|
| `trait` | Annotation key to colour by (e.g. `"region"`). `NULL` for uncoloured |
| `layout` | `"rectangular"`, `"circular"`, `"fan"`, or any other `ggtree` layout |
| `mrsd` | Most-recent sampling date for absolute time axis |
| `tip_labels` | Show tip labels |

**Returns** a `ggtree` / `ggplot2` object.

---

## 6. Bundled example data

| Function | File | Description |
|---|---|---|
| `pact_example_files()` | `inst/extdata/example/` | PACT structured coalescent example: NEWICK trees + `in.param` + expected outputs |
| `pact_example_dta()` | `inst/extdata/example_dta_trees/rsva_resample.trees` | RSV-A BEAST DTA posterior: 45 trees, 2103 tips, `region` trait (12 states), 1956–2023 |
| — | `inst/extdata/example_dta_trees/in.param` | PACT parameter file for the RSV-A trees |
| — | `inst/extdata/example_dta_trees/analyse_rsva.R` | Ready-to-run analysis script |
| `pact_parameter_reference()` | `inst/extdata/parameters.txt` | Full `in.param` command reference |
| `pact_manual()` | `inst/doc/pact_manual.pdf` | Original C++ PACT user manual |

---

## 7. Complete workflow examples

### Structured coalescent analysis (core PACT)

```r
library(pactR)

ex <- pact_example_files()

results <- run_pact(
  trees_file          = ex$trees,
  output_dir          = "output/",
  burnin              = 10,
  push_times_back     = c(2002, 2008),
  renew_trunk         = 1,
  summary_tmrca       = TRUE,
  summary_coal_rates  = TRUE,
  summary_mig_rates   = TRUE,
  summary_diversity   = TRUE,
  summary_persistence = TRUE,
  skyline_settings    = c(2002, 2008, 0.1),
  skyline_proportions = TRUE,
  skyline_diversity   = TRUE,
  tips_time_to_trunk  = TRUE
)

stats    <- read_pact_stats(results$stats)
skylines <- read_pact_skylines(results$skylines)
tips     <- read_pact_tips(results$tips)

# Structured coalescent proportions as a proportion-over-time plot
prop_sc <- pact_to_proportions(skylines)
plot_dta_proportions(prop_sc, type = "ribbon")
```

### DTA posterior analysis (extensions)

```r
library(pactR)
library(ggplot2)

f      <- pact_example_dta()
BURNIN <- 5L

# -- Persistence --
pers <- persistence_dta(f, trait = "region", burnin = BURNIN)
p1   <- plot_persistence(pers, title = "RSV-A Lineage Persistence")
ggsave("fig_persistence.png", p1, width = 7, height = 5, dpi = 150)

# -- Proportions over time --
prop <- dta_proportions(f, trait = "region", step = 0.25, burnin = BURNIN)
p2   <- plot_dta_proportions(prop, type = "ribbon",
                              title = "RSV-A Proportions Over Time")
p3   <- plot_dta_proportions(prop, type = "area")
ggsave("fig_proportions_ribbon.png", p2, width = 11, height = 5.5, dpi = 150)
ggsave("fig_proportions_area.png",   p3, width = 11, height = 5.5, dpi = 150)

# -- Lineage through time --
ltt <- lineage_through_time(f, step = 0.25, burnin = BURNIN)
p4  <- plot_ltt(ltt, title = "RSV-A Lineage Through Time")
ggsave("fig_ltt.png", p4, width = 9, height = 4.5, dpi = 150)

# -- Molecular clock check --
rtt <- root_to_tip(f, plot = TRUE)

# -- No parseable dates: anchor with mrsd --
pers2 <- persistence_dta(f, mrsd = "2023-06-15", burnin = BURNIN)
```

### Using a param file with R argument overrides

```r
# Base settings from file, override burnin via R
results <- run_pact(
  trees_file  = "my.trees",
  param_file  = "base.param",
  output_dir  = "out/",
  burnin      = 20           # overrides burnin in base.param
)
```

---

## Function index

| Function | Layer | Description |
|---|---|---|
| `run_pact()` | Core | Run full PACT C++ analysis |
| `build_pact_params()` | Core | Build param file lines from R arguments |
| `read_pact_stats()` | Core | Parse `out.stats` → data frame |
| `read_pact_tips()` | Core | Parse `out.tips` → data frame |
| `read_pact_skylines()` | Core | Parse `out.skylines` → data frame |
| `read_pact_pairs()` | Core | Parse `out.pairs` → data frame |
| `pact_to_proportions()` | Core | Reshape skylines for `plot_dta_proportions()` |
| `pact_example_files()` | Core | Paths to bundled structured coalescent example |
| `pact_parameter_reference()` | Core | Print full parameter reference |
| `pact_manual()` | Core | Open PDF manual |
| `read_beast_posterior()` | DTA | Parse BEAST NEXUS posterior → edge data frames |
| `parse_tip_dates()` | DTA | Extract decimal dates from tip labels |
| `persistence_dta()` | DTA | Persistence time per state (PACT-matching) |
| `plot_persistence()` | DTA | Dot-and-whisker persistence chart |
| `dta_proportions()` | DTA | Proportion of lineages in each state over time |
| `plot_dta_proportions()` | DTA | Ribbon / area / line proportion chart |
| `lineage_through_time()` | DTA | Lineage count over time (posterior mean + CI) |
| `plot_ltt()` | DTA | Lineage-through-time plot |
| `root_to_tip()` | DTA | Molecular clock regression |
| `pairwise_tmrca()` | DTA | Pairwise TMRCA matrix |
| `find_mrca()` | DTA | MRCA node of a tip set |
| `plot_tree_beast()` | DTA | Annotated tree via ggtree (requires Bioconductor) |
| `pact_example_dta()` | DTA | Path to bundled RSV-A BEAST posterior |
