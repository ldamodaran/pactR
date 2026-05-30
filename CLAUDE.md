# pactR R Package — Developer Guide for Claude

## What this is

`pactR` is an R package wrapping the PACT (Posterior Analysis of Coalescent Trees) C++ CLI tool by Trevor Bedford. It provides:

1. **Rcpp wrapper** around the original C++ coalescent statistics engine
2. **DTA extension functions** (`R/extensions.R`) for BEAST DTA posterior tree analysis and visualization

## Package location

```
/Users/ldamoda/Library/CloudStorage/OneDrive-Emory/Desktop/pact/
```

GitHub: `https://github.com/ldamodaran/pactR`

## Build / install cycle

```bash
# From terminal (fast reinstall):
cd /Users/ldamoda/Library/CloudStorage/OneDrive-Emory/Desktop
R CMD INSTALL pact
```

```r
# Only needed if Rcpp-exported C++ functions changed:
Rcpp::compileAttributes("/Users/ldamoda/Library/CloudStorage/OneDrive-Emory/Desktop/pact")
```

## Key files

| File | Purpose |
|------|---------|
| `src/rcpp_pact.cpp` | Rcpp entry point — exports `pact_run_cpp(trees, param, prefix)` |
| `src/io.cpp` / `src/io.h` | Modified to accept file path args; `cout` → `Rcpp::Rcout` |
| `src/coaltree.cpp` | Core coalescent tree algorithms (persistence, migration, proportions) |
| `R/pact.R` | Core wrappers: `run_pact()`, `build_pact_params()`, `read_pact_*()` |
| `R/extensions.R` | DTA posterior analysis + visualization |
| `NAMESPACE` | All exports declared here — update when adding new functions |
| `DESCRIPTION` | `Package: pactR`; `Imports: Rcpp, ape`; `Suggests: ggplot2, scales, ggtree, treeio` |
| `inst/extdata/example/` | Core PACT example: `in.trees` (influenza H3N2), `in.param`, `analyse_example.R` |
| `inst/extdata/example_dta_trees/` | DTA example: `global_HPAI.trees`, `in.param`, `analyse_global_HPAI.R` |
| `inst/extdata/parameters.txt` | Full PACT parameter reference |
| `inst/doc/pact_manual.pdf` | Original C++ user manual |

## Exported R functions

### Core PACT (R/pact.R)
- `run_pact(trees_file, param_file, output_dir, output_prefix, ...)` — run C++ analysis
- `build_pact_params(burnin, push_times_back, ...)` — build param file from R args
- `read_pact_stats(file)`, `read_pact_tips(file)`, `read_pact_skylines(file)`, `read_pact_pairs(file)`
- `pact_example_files()` — paths to bundled `in.trees` + `in.param`
- `pact_parameter_reference()`, `pact_manual()`

### DTA extensions (R/extensions.R)
- `pact_example_dta()` — path to `global_HPAI.trees`
- `parse_tip_dates(tips, pattern, decimal)` — auto-detects ISO `YYYY-MM-DD` and decimal/year formats
- `read_beast_posterior(file, burnin, trait, tip_dates, mrsd)` — parse BEAST NEXUS posterior
- `persistence_dta(file, trait, burnin, tip_dates, mrsd)` — lineage persistence (matches PACT `summary persistence`)
- `plot_persistence(data, palette, title, sort_states)` — dot-and-whisker persistence plot
- `dta_proportions(file, trait, times, step, burnin, tip_dates, mrsd)` — branch-length-weighted proportions over time (matches PACT `skyline proportions`)
- `pact_to_proportions(skylines, prefix)` — reshape PACT skyline output to proportions format
- `plot_dta_proportions(data, type, palette, ...)` — ribbon / stacked area / line plot
- `lineage_through_time(file, times, step, burnin, tip_dates, mrsd)` — LTT computation
- `plot_ltt(ltt_data, colour, title, ...)` — LTT visualization
- `plot_tree_beast(file, trait, tip_labels, palette, layout, mrsd, ...)` — ggtree visualization
- `root_to_tip(file, tip_dates, tree_index, plot)` — molecular clock regression
- `pairwise_tmrca(file, tree_index, tip_dates)` — pairwise TMRCA matrix
- `find_mrca(file, tips, tree_index)` — MRCA node finder

## Algorithm notes

### persistence_dta()
Matches C++ `getPersistence()` exactly: for each tip, walk ancestors until finding the first node with a *different* label; persistence = `tip_time − ancestor_time`. Per-state mean within each tree; 2.5%/mean/97.5% quantiles across posterior trees.

### dta_proportions()
Matches C++ `skyline proportions` (`trimEnds` + `getLabelPro`): for window `[t, t+step]`, computes `overlap = min(child_time, t+step) − max(parent_time, t)` for each branch, then `sum(overlap for state) / sum(all overlaps)`. Time reported at midpoint `t + step/2`. This is branch-length-weighted, not lineage-count-based.

### lineage_through_time()
Counts branches where `child_time >= t AND parent_time < t` at each time point t (matches C++ `getCoalWeight` logic).

## Temporal anchoring

Priority order in `.compute_node_times()`:
1. `mrsd` (explicit most-recent sampling date, string `"YYYY-MM-DD"` or decimal year)
2. Auto-parsed tip dates via `parse_tip_dates()` — tries ISO `YYYY-MM-DD` first, then trailing year
3. Fallback: root-relative branch lengths (root = 0), matching original PACT behavior

## Internal NEXUS parser

`.parse_nexus_translations(lines)` — reads `Translate` block → named char vector (num → taxon name)

`.parse_beast_tree_annotations(tree_str, translations, trait)` — extracts:
- Tip annotations: `123[&geo="X"]` pattern
- Internal node annotations: `)[&geo="X"]` in left-to-right closing-paren order (matches ape's n+1, n+2,... internal node numbering)
- Returns clean NEWICK with all `[&...]` blocks stripped

`.compute_node_times(phy, tip_dates, mrsd)` — computes absolute or relative node times via `ape::node.depth.edgelength()`

## Bundled examples

### Core PACT (`inst/extdata/example/`)
- `in.trees` — single influenza H3N2 NEWICK tree with migration annotations `[&M state1 state2:rate]`; 7 numeric states; time range 2002–2008
- `in.param` — runs TMRCA, coalescent rates, migration rates, skyline proportions, time-to-trunk
- `analyse_example.R` — runs full analysis, writes `stats.tsv`, `skylines.tsv`, `tips.tsv`, and 5 PNG figures

### DTA posterior (`inst/extdata/example_dta_trees/`)
- `global_HPAI.trees` — global HPAI H5N1 BEAST DTA posterior; 31 trees, 1921 tips, `geo` trait (Asia / Europe / North_America), ISO dates in tip names, time range 2018–2023
- `in.param` — burnin 3, monthly skyline bins 2021–2023
- `analyse_global_HPAI.R` — runs persistence, proportions, LTT; writes 3 TSV + 4 PNG files

## Typical workflow

```r
library(pactR)

# --- Core PACT coalescent analysis ---
ex      <- pact_example_files()
results <- run_pact(ex$trees, ex$param, output_dir = tempdir())
stats   <- read_pact_stats(results$stats)
sky     <- read_pact_skylines(results$skylines)

# --- DTA posterior analysis ---
f    <- pact_example_dta()              # global_HPAI.trees
pers <- persistence_dta(f, trait = "geo", burnin = 3L)
plot_persistence(pers)

prop <- dta_proportions(f, trait = "geo", step = 0.05, burnin = 3L)
plot_dta_proportions(prop, type = "area")

ltt  <- lineage_through_time(f, step = 0.05, burnin = 3L)
plot_ltt(ltt)

# Supply mrsd when tip dates are not in names:
pers2 <- persistence_dta(f, trait = "geo", mrsd = "2023-06-30", burnin = 3L)
```

## Status (2026-05-30)

- Package builds and installs as `pactR`; all functions load correctly
- All three DTA functions (`persistence_dta`, `dta_proportions`, `lineage_through_time`) audited against C++ source and verified correct
- `dta_proportions()` corrected to branch-length-weighted proportions (was count-based); time now reported at window midpoint
- `parse_tip_dates()` auto-detects ISO `YYYY-MM-DD` dates (needed for HPAI tip names)
- Tested on global HPAI H5N1 (31 trees, 1921 tips) and original H3N2 structured coalescent example
- All example scripts run cleanly with no warnings
- **Bug fix (2026-05-30):** `RcppExports.cpp` and `RcppExports.R` had stale `_pact_` symbol prefix and `R_init_pact` from before the `pact` → `pactR` rename. R looks for `R_init_pactR` on package load, so routines were never registered and `run_pact()` failed with "pact_run_cpp not found". Fixed by updating all symbols to `_pactR_` prefix.
- **Note:** `persistence_dta` and `dta_proportions` are pure R and do not use `in.param`. That file is only needed by `run_pact()` (the C++ engine). DTA function parameters (burnin, time window, step) are passed as direct R arguments.
- Pushed to `github.com/ldamodaran/pactR` via SSH
