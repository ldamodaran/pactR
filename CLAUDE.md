# PACT R Package — Developer Guide for Claude

## What this is

`pact` is an R package wrapping the PACT (Posterior Analysis of Coalescent Trees) C++ CLI tool by Trevor Bedford. It provides:

1. **Rcpp wrapper** around the original C++ coalescent statistics engine
2. **Baltic-inspired extension functions** (`R/extensions.R`) for DTA posterior tree analysis and visualization

## Package location

```
/Users/ldamoda/Library/CloudStorage/OneDrive-Emory/Desktop/pact/
```

Built tarball ends up at: `/Users/ldamoda/Library/CloudStorage/OneDrive-Emory/Desktop/pact_0.9.4.tar.gz`

## Build / install cycle

```r
# From R:
Rcpp::compileAttributes("/Users/ldamoda/Library/CloudStorage/OneDrive-Emory/Desktop/pact")
# Only needed if Rcpp-exported C++ functions changed

# From terminal:
cd /Users/ldamoda/Library/CloudStorage/OneDrive-Emory/Desktop
R CMD build pact
R CMD INSTALL pact_0.9.4.tar.gz
```

## Key files

| File | Purpose |
|------|---------|
| `src/rcpp_pact.cpp` | Rcpp entry point — exports `pact_run_cpp(trees, param, prefix)` |
| `src/io.cpp` / `src/io.h` | Modified to accept file path args; `cout` → `Rcpp::Rcout` |
| `src/param.cpp` / `src/param.h` | Constructor takes `const std::string& paramFile` |
| `R/pact.R` | Main wrappers: `run_pact()`, `build_pact_params()`, `read_pact_*()` |
| `R/extensions.R` | DTA posterior analysis + visualization (Baltic-inspired, see below) |
| `NAMESPACE` | All exports declared here — must be updated when new functions added |
| `DESCRIPTION` | `Imports: Rcpp, ape`; `Suggests: ggplot2, scales, ggtree, treeio` |
| `inst/extdata/example/` | PACT example `in.trees`, `in.param`, expected outputs |
| `inst/extdata/example_dta_trees/rsva_resample.trees` | RSV-A DTA posterior (45 trees, 2103 tips, `region` trait) |
| `inst/extdata/parameters.txt` | Full PACT parameter reference |
| `inst/doc/pact_manual.pdf` | Original user manual |

## Exported R functions

### Core PACT (R/pact.R)
- `run_pact(trees_file, param_file, output_dir, output_prefix, ...)` — run analysis
- `build_pact_params(burnin, push_times_back, ...)` — build param file from named R args
- `read_pact_stats(file)`, `read_pact_tips(file)`, `read_pact_skylines(file)`, `read_pact_pairs(file)`
- `pact_example_files()`, `pact_parameter_reference()`, `pact_manual()`

### DTA extensions (R/extensions.R)
- `pact_example_dta()` — path to `rsva_resample.trees`
- `parse_tip_dates(tips, pattern, decimal)` — extract decimal dates from tip names
- `read_beast_posterior(file, burnin, trait, tip_dates, mrsd)` — parse BEAST NEXUS posterior
- `persistence_dta(file, trait, times, step, burnin, tip_dates, mrsd)` — lineage persistence over time
- `pact_to_persistence(skylines, prefix)` — convert PACT skyline output to persistence format
- `plot_persistence(data, type, palette, ...)` — ggplot2 visualization (ribbon/area/line)
- `lineage_through_time(file, times, step, burnin, tip_dates, mrsd)` — LTT computation
- `plot_ltt(ltt_data, log_y, colour, ...)` — LTT visualization
- `plot_tree_beast(file, trait, tip_labels, palette, layout, mrsd, ...)` — ggtree visualization
- `root_to_tip(file, tip_dates, tree_index, plot)` — molecular clock check
- `pairwise_tmrca(file, tree_index, tip_dates)` — pairwise TMRCA matrix
- `find_mrca(file, tips, tree_index)` — MRCA node finder

## Temporal anchoring in extensions.R

The original PACT uses branch lengths directly (BEAST outputs time-calibrated trees in decimal years — no date parsing from tip names needed). Extensions follow the same logic:

1. **`mrsd` parameter** (e.g. `mrsd = "2023-06-15"` or `mrsd = 2023.456`): anchors the tree so the tallest tip = mrsd; all node times become absolute calendar years
2. **Auto-parsed tip dates** (default): regex extracts 4-digit year from tip names; adds 0.5 if no decimal part
3. **Fallback**: if neither mrsd nor parseable dates exist → root-relative branch lengths (root = 0), same as original PACT. A message is printed to inform the user.

`"YYYY-MM-DD"` strings are converted to decimal years via `as.Date()`.

## Internal NEXUS parser (no treeio dependency)

`.parse_nexus_translations(lines)` — reads `Translate` block → named char vector (num → taxon name)

`.parse_beast_tree_annotations(tree_str, translations, trait)` — extracts:
- Tip annotations: `123[&region="X"]` pattern
- Internal node annotations: `)[&region="X"]` in left-to-right closing-paren order (matches ape's n+1, n+2,... internal node numbering)
- Returns clean NEWICK with all `[&...]` blocks stripped

`.compute_node_times(phy, tip_dates, mrsd)` — computes absolute or relative node times using `ape::node.depth.edgelength()`

## plot_persistence() input format

Any data frame with columns `time`, `state`, `mean` (and optionally `lower`, `upper` for CI ribbons). Accepts output from:
- `persistence_dta()` — DTA posterior proportions
- `pact_to_persistence(read_pact_skylines(...))` — PACT structured coalescent proportions
- Custom data for Markov Rewards or other sources

```r
plot_persistence(pers, type = "ribbon")   # lines + CI ribbons
plot_persistence(pers, type = "area")     # stacked filled area
plot_persistence(pers, type = "line")     # lines only, no CI
```

Default palette: Okabe-Ito colour-blind-friendly colours. LOESS smoothing applied by default (`smooth = TRUE`, `smooth_span = 0.15`).

## Typical workflow

```r
library(pactR)

# --- PACT coalescent analysis ---
ex      <- pact_example_files()
results <- run_pact(ex$trees, ex$param, output_dir = tempdir())
sky     <- read_pact_skylines(results$skylines)
pers_sc <- pact_to_persistence(sky)        # structured coalescent proportions
plot_persistence(pers_sc)

# --- DTA posterior analysis ---
f       <- pact_example_dta()              # rsva_resample.trees
pers    <- persistence_dta(f, trait = "region", step = 0.25, burnin = 5L)
plot_persistence(pers, type = "area")

ltt     <- lineage_through_time(f, step = 0.25)
plot_ltt(ltt)

rtt     <- root_to_tip(f)                  # molecular clock check

# --- Without parseable dates: supply mrsd ---
pers2   <- persistence_dta(f, mrsd = "2023-06-30", step = 0.25)
```

## Status (2026-05-26)

- Package builds, installs, runs correctly
- All core PACT functions verified against example data
- All extension functions tested on RSV-A DTA posterior (45 trees, 12 regions)
- `plot_persistence()` ribbon/area/line render cleanly with no warnings
- `mrsd` parameter implemented and tested for `read_beast_posterior`, `persistence_dta`, `lineage_through_time`
- Fallback to root-relative branch lengths when no dates available (matches original PACT behavior)

## Pending / possible future work

- Markov Rewards input support for `plot_persistence()` (currently handles DTA and PACT SC)
- `plot_tree_beast()` requires `ggtree` + `treeio` (Bioconductor) — not yet tested end-to-end
- No `man/` pages (roxygen2 not run); `?function_name` help works from `@export` docstrings in vignette build only
