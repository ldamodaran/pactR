# ── Internal BEAST NEXUS parser ───────────────────────────────────────────────
# Depends only on base R + ape. Does NOT require treeio.

.parse_nexus_translations <- function(lines) {
  trans_start <- grep("(?i)^\\s*translate\\s*$", lines, perl = TRUE)[1]
  if (is.na(trans_start)) return(NULL)
  semi_after <- which(grepl("^\\s*;\\s*$", lines) & seq_along(lines) > trans_start)
  trans_end <- semi_after[1]
  raw <- trimws(lines[(trans_start + 1L):(trans_end - 1L)])
  raw <- sub(",\\s*$", "", raw)
  raw <- raw[nzchar(raw)]
  parts <- strsplit(raw, "\\s+")
  setNames(vapply(parts, `[[`, "", 2L), vapply(parts, `[[`, "", 1L))
}

.parse_beast_tree_annotations <- function(tree_str, translations, trait = "region") {
  # Strip the "tree NAME [...] = [&R] " prefix
  newick_part <- sub("^tree\\s+\\S+[^=]+=\\s*\\[&[A-Za-z]+\\]\\s*", "",
                     tree_str, perl = TRUE)

  # Clean NEWICK (strip all annotations)
  clean <- gsub("\\[&[^\\]]*\\]", "", newick_part)

  # If no trait requested, return early with empty annotation vectors
  if (is.null(trait)) {
    return(list(clean_newick = clean, tip_region_by_num = character(0),
                internal_regions = character(0)))
  }

  # --- Tip annotations: 123[&region="X"] ---
  tip_pat <- sprintf('(\\d+)\\[&[^\\]]*%s="([^"]+)"[^\\]]*\\]', trait)
  tip_hits <- gregexpr(tip_pat, newick_part, perl = TRUE)
  tip_matches <- regmatches(newick_part, tip_hits)[[1L]]

  tip_nums    <- regmatches(tip_matches, regexpr("^\\d+", tip_matches))
  tip_regions <- regmatches(
    tip_matches,
    regexpr(sprintf('(?<=%s=")[^"]+', trait), tip_matches, perl = TRUE)
  )
  tip_region_by_num <- setNames(tip_regions, tip_nums)

  # --- Internal node annotations: )[&region="X"] ---
  int_pat     <- sprintf('\\)\\[&[^\\]]*%s="([^"]+)"[^\\]]*\\]', trait)
  int_hits    <- gregexpr(int_pat, newick_part, perl = TRUE)
  int_matches <- regmatches(newick_part, int_hits)[[1L]]
  internal_regions <- regmatches(
    int_matches,
    regexpr(sprintf('(?<=%s=")[^"]+', trait), int_matches, perl = TRUE)
  )

  # --- Strip all [&...] blocks to get clean NEWICK with numeric tip labels ---
  clean <- gsub("\\[&[^\\]]*\\]", "", newick_part)

  list(
    clean_newick      = clean,
    tip_region_by_num = tip_region_by_num,  # "123" -> "Oceania"
    internal_regions  = internal_regions     # ordered to match ape node n+1, n+2...
  )
}

.compute_node_times <- function(phy, tip_dates, mrsd = NULL) {
  heights     <- ape::node.depth.edgelength(phy)
  n_tips      <- length(phy$tip.label)
  tip_heights <- heights[seq_len(n_tips)]

  if (!is.null(mrsd)) {
    # Anchor: most recent tip height maps to mrsd
    return((mrsd - max(tip_heights)) + heights)
  }

  if (is.null(tip_dates) || all(is.na(tip_dates[phy$tip.label]))) {
    # No date information — return root-relative heights (root = 0)
    return(heights)
  }

  tip_date_vec <- tip_dates[phy$tip.label]
  root_time    <- stats::median(tip_date_vec - tip_heights, na.rm = TRUE)
  root_time + heights
}


# ── Tip date parsing ──────────────────────────────────────────────────────────

#' Parse sampling dates from tip labels
#'
#' Extracts decimal sampling dates from phylogenetic tip names. Automatically
#' detects ISO \code{YYYY-MM-DD} dates (e.g. \code{sample|2022-03-15|host})
#' and converts them to decimal years. Falls back to matching trailing
#' 4-digit years or decimal years (e.g. \code{RSV_USA_2014.75}).
#'
#' @param tips Character vector of tip labels.
#' @param pattern Regex with one capture group matching the date component.
#'   If \code{NULL} (default), the function first tries ISO \code{YYYY-MM-DD}
#'   format, then falls back to trailing year/decimal-year format.
#' @param decimal Logical. If \code{TRUE} (default) and only a whole year is
#'   captured (no decimal point), adds 0.5 to represent the midpoint of the
#'   year. Ignored for ISO dates (which are already precise).
#' @return Named numeric vector of decimal dates. Names are the tip labels.
#'   Tips with no date match receive \code{NA}.
#'
#' @examples
#' tips <- c("RSV_A_Argentina_2014", "RSV_B_USA_2016", "RSV_A_UK_2015")
#' parse_tip_dates(tips)
#'
#' tips2 <- c("epi_isl_123|a/chicken/iowa/2022|2022-03-06|domestic",
#'            "epi_isl_456|a/duck/kansas/2022|2022-04-01|backyard_bird")
#' parse_tip_dates(tips2)
#'
#' @export
parse_tip_dates <- function(tips, pattern = NULL, decimal = TRUE) {
  iso_pat  <- "(\\d{4}-\\d{2}-\\d{2})"
  year_pat <- "(\\d{4})(?:[._]\\d+)?$"

  # Auto-detect: try ISO dates first when no pattern supplied
  if (is.null(pattern)) {
    m_iso <- regexpr(iso_pat, tips, perl = TRUE)
    if (sum(m_iso != -1L) > 0L) {
      # At least some tips have ISO dates — use ISO path
      date_str <- rep(NA_character_, length(tips))
      has_match <- m_iso != -1L
      date_str[has_match] <- regmatches(tips, m_iso)
      # Convert YYYY-MM-DD to decimal year
      dates <- rep(NA_real_, length(tips))
      valid <- !is.na(date_str)
      if (any(valid)) {
        d <- as.Date(date_str[valid])
        yr <- as.integer(format(d, "%Y"))
        jd <- as.integer(format(d, "%j"))
        # Days in year (account for leap years)
        days_in_yr <- as.integer(format(as.Date(paste0(yr, "-12-31")), "%j"))
        dates[valid] <- yr + (jd - 1L) / days_in_yr
      }
      return(setNames(dates, tips))
    }
    pattern <- year_pat
  }

  m         <- regexpr(pattern, tips, perl = TRUE)
  date_str  <- rep(NA_character_, length(tips))
  has_match <- m != -1L
  date_str[has_match] <- regmatches(tips, m)
  dates <- as.numeric(date_str)
  if (decimal) {
    is_year <- !is.na(date_str) & nchar(date_str) == 4L & !grepl("\\.", date_str)
    dates[is_year] <- dates[is_year] + 0.5
  }
  setNames(dates, tips)
}


# ── Read BEAST posterior trees ────────────────────────────────────────────────

#' Read a BEAST posterior tree file
#'
#' Reads a NEXUS file containing one or more BEAST-annotated trees (e.g., a
#' DTA posterior sample). Returns a list of edge-level data frames suitable
#' for further analysis.
#'
#' @param file Path to the BEAST \code{.trees} NEXUS file.
#' @param burnin Integer number of trees to discard from the start.
#' @param trait Character. Discrete trait annotation key to extract (e.g.
#'   \code{"region"}).
#' @param tip_dates Named numeric vector of decimal tip dates, with names
#'   matching tip labels. If \code{NULL} (default), dates are parsed
#'   automatically from tip names using \code{\link{parse_tip_dates}}.
#' @param mrsd Numeric or character. Most recent sampling date used to anchor
#'   absolute time. Accepts a decimal year (e.g. \code{2017.5}) or an ISO
#'   date string (e.g. \code{"2017-06-30"}). When provided, overrides
#'   \code{tip_dates}-based anchoring. If neither \code{tip_dates} nor
#'   \code{mrsd} yield usable dates, node times are root-relative branch
#'   lengths (root = 0), matching PACT's original approach.
#'
#' @return A named list with:
#' \describe{
#'   \item{\code{trees}}{A list of data frames, one per tree. Each has columns
#'     \code{parent}, \code{child}, \code{parent_time}, \code{child_time},
#'     \code{region} (or the trait name), \code{branch_length}.}
#'   \item{\code{translations}}{Named character vector of taxon translations.}
#'   \item{\code{tip_dates}}{Named numeric vector of tip decimal dates.}
#'   \item{\code{absolute_time}}{Logical. \code{TRUE} if times are in calendar
#'     years; \code{FALSE} if root-relative branch lengths.}
#' }
#'
#' @seealso \code{\link{persistence_dta}}, \code{\link{parse_tip_dates}}
#'
#' @export
read_beast_posterior <- function(file, burnin = 0L, trait = "region",
                                  tip_dates = NULL, mrsd = NULL) {
  lines        <- readLines(file, warn = FALSE)
  translations <- .parse_nexus_translations(lines)
  if (is.null(translations))
    stop("No Translate block found in '", file, "'.")

  tree_idx    <- grep("^\\s*tree\\s+", lines, ignore.case = TRUE)
  tree_lines  <- trimws(lines[tree_idx])
  if (burnin > 0L) tree_lines <- tree_lines[seq_len(length(tree_lines))[-seq_len(burnin)]]
  if (length(tree_lines) == 0L) stop("No trees remain after burnin.")

  # Parse mrsd to decimal year if given as a date string
  if (!is.null(mrsd) && is.character(mrsd)) {
    d    <- as.Date(mrsd)
    mrsd <- as.numeric(format(d, "%Y")) +
            (as.numeric(format(d, "%j")) - 1L) / 365.25
  } else if (!is.null(mrsd)) {
    mrsd <- as.numeric(mrsd)
  }

  # Compute tip dates once (same tip names across all trees)
  if (is.null(tip_dates)) {
    tip_dates <- parse_tip_dates(unname(translations))
  }

  # Determine whether we have usable temporal anchors
  has_dates <- !is.null(mrsd) || (!is.null(tip_dates) && !all(is.na(tip_dates)))
  if (!has_dates) {
    message("No tip dates or mrsd found — using root-relative branch lengths ",
            "(matching PACT's original approach). Supply 'mrsd' for absolute time.")
  }

  tree_list <- lapply(tree_lines, function(tl) {
    ann  <- .parse_beast_tree_annotations(tl, translations, trait)
    phy  <- ape::read.tree(text = ann$clean_newick)

    # Map numeric tip labels to taxon names
    phy$tip.label <- unname(translations[phy$tip.label])

    n_tips <- length(phy$tip.label)
    node_times <- tryCatch(
      .compute_node_times(phy, tip_dates, mrsd = mrsd),
      error = function(e) rep(NA_real_, n_tips + phy$Nnode)
    )

    # Build region vector indexed by ape node number
    region_by_node <- rep(NA_character_, n_tips + phy$Nnode)

    # Tip regions: map taxon number -> name -> node index
    for (num in names(ann$tip_region_by_num)) {
      tip_name <- translations[num]
      idx <- match(tip_name, phy$tip.label)
      if (!is.na(idx)) region_by_node[idx] <- ann$tip_region_by_num[num]
    }
    # Internal node regions (ordered n+1, n+2, ...)
    n_internal <- length(ann$internal_regions)
    if (n_internal > 0L) {
      internal_idx <- n_tips + seq_len(min(n_internal, phy$Nnode))
      region_by_node[internal_idx] <- ann$internal_regions[seq_len(length(internal_idx))]
    }

    edges <- phy$edge
    data.frame(
      parent        = edges[, 1L],
      child         = edges[, 2L],
      parent_time   = node_times[edges[, 1L]],
      child_time    = node_times[edges[, 2L]],
      region        = region_by_node[edges[, 2L]],   # child-node annotation
      branch_length = phy$edge.length,
      stringsAsFactors = FALSE
    )
  })

  list(trees = tree_list, translations = translations,
       tip_dates = tip_dates, absolute_time = has_dates)
}


# ── Persistence (PACT-matching: time to first label change) ───────────────────

#' Compute lineage persistence from a BEAST DTA posterior sample
#'
#' Replicates PACT's \code{summarypersistence} statistic for DTA posterior
#' trees. For each tip, walks back through its ancestors to find the first
#' ancestor with a \emph{different} discrete trait state, and records
#' tip_time \eqn{-} ancestor_time. Per-state means are computed across tips
#' within each tree, then summarised across the posterior as mean + 95\%
#' credible interval.
#'
#' This answers: \emph{"On average, how long has a lineage sampled from region X
#' been continuously in region X before the most recent state change?"}
#'
#' @param file Path to the BEAST \code{.trees} NEXUS file.
#' @param trait Character. Discrete trait annotation key (e.g. \code{"region"}).
#' @param burnin Integer. Number of trees to discard from the start.
#' @param tip_dates Named numeric vector of decimal tip dates. Auto-parsed if
#'   \code{NULL}.
#' @param mrsd Numeric or character most-recent sampling date (decimal year or
#'   \code{"YYYY-MM-DD"}). See \code{\link{read_beast_posterior}}.
#'
#' @return A data frame with columns \code{state}, \code{lower} (2.5th
#'   percentile across trees), \code{mean}, and \code{upper} (97.5th
#'   percentile), all in years.
#'
#' @examples
#' \donttest{
#' f    <- system.file("extdata", "example_dta_trees",
#'                     "rsva_resample.trees", package = "pactR")
#' pers <- persistence_dta(f, trait = "region", burnin = 5L)
#' plot_persistence(pers)
#' }
#'
#' @seealso \code{\link{plot_persistence}}, \code{\link{read_beast_posterior}}
#'
#' @export
persistence_dta <- function(file, trait = "region", burnin = 0L,
                             tip_dates = NULL, mrsd = NULL) {
  message("Reading trees from: ", file)
  post  <- read_beast_posterior(file, burnin = burnin, trait = trait,
                                tip_dates = tip_dates, mrsd = mrsd)
  trees <- post$trees

  all_states <- sort(unique(unlist(lapply(trees, `[[`, "region"),
                                   use.names = FALSE)))
  all_states <- all_states[!is.na(all_states)]

  message("Computing persistence for ", length(trees), " trees...")

  # For each tree, compute per-tip persistence then average by state
  per_tree <- lapply(trees, function(df) {
    df <- df[is.finite(df$child_time) & is.finite(df$parent_time), ]
    max_node <- max(c(df$parent, df$child))

    # O(1) lookup tables indexed by node number
    parent_of <- integer(max_node)
    parent_of[df$child] <- df$parent

    region_of <- rep(NA_character_, max_node)
    region_of[df$child] <- df$region    # covers all non-root nodes

    time_of <- numeric(max_node)
    time_of[df$child] <- df$child_time
    root_nd <- setdiff(df$parent, df$child)
    if (length(root_nd) == 1L)
      time_of[root_nd] <- df$parent_time[match(root_nd, df$parent)]

    tips <- setdiff(df$child, df$parent)   # leaf nodes

    # PACT algorithm: walk from tip to first ancestor with different label
    tip_persist <- vapply(tips, function(tip) {
      tip_rgn <- region_of[tip]
      if (is.na(tip_rgn)) return(NA_real_)
      tip_t   <- time_of[tip]
      cur     <- parent_of[tip]
      while (cur != 0L) {
        anc_rgn <- region_of[cur]
        if (!is.na(anc_rgn) && anc_rgn != tip_rgn)
          return(tip_t - time_of[cur])
        cur <- parent_of[cur]
      }
      NA_real_   # root has same label as tip — excluded (matches PACT)
    }, numeric(1L))

    tip_rgns <- region_of[tips]

    # Per-state mean across tips in this tree
    sapply(all_states, function(s) {
      v <- tip_persist[!is.na(tip_rgns) & tip_rgns == s]
      v <- v[!is.na(v)]
      if (length(v) == 0L) NA_real_ else mean(v)
    })
  })

  mat <- do.call(rbind, per_tree)   # n_trees × n_states

  do.call(rbind, lapply(seq_along(all_states), function(j) {
    vals <- mat[, j][!is.na(mat[, j])]
    data.frame(
      state = all_states[j],
      lower = unname(stats::quantile(vals, 0.025)),
      mean  = mean(vals),
      upper = unname(stats::quantile(vals, 0.975)),
      stringsAsFactors = FALSE
    )
  }))
}


#' Plot lineage persistence as a dot-and-whisker chart
#'
#' Draws a horizontal point-range (dot-and-whisker) plot showing the mean and
#' 95\% credible interval of lineage persistence time per state. Matches PACT's
#' \code{summarypersistence} output. Accepts output from
#' \code{\link{persistence_dta}} or any data frame with columns \code{state},
#' \code{mean}, and optionally \code{lower}, \code{upper}.
#'
#' @param data Data frame with columns \code{state}, \code{mean}, and
#'   optionally \code{lower} and \code{upper}.
#' @param palette Named colour vector mapping state names to colours. If
#'   \code{NULL}, uses the Okabe-Ito colour-blind-friendly palette.
#' @param title Character. Plot title.
#' @param xlab X-axis label.
#' @param sort_states Logical. If \code{TRUE} (default), states are sorted by
#'   mean persistence (shortest to longest, so the longest appears at the top).
#'
#' @return A \code{ggplot2} object.
#'
#' @examples
#' \donttest{
#' f    <- system.file("extdata", "example_dta_trees",
#'                     "rsva_resample.trees", package = "pactR")
#' pers <- persistence_dta(f, trait = "region", burnin = 5L)
#' plot_persistence(pers)
#' }
#'
#' @seealso \code{\link{persistence_dta}}
#'
#' @export
plot_persistence <- function(data,
                              palette     = NULL,
                              title       = "Lineage persistence by state",
                              xlab        = "Mean time in state (years)",
                              sort_states = TRUE) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is required. Install with install.packages('ggplot2').")

  required <- all(c("state", "mean") %in% names(data))
  if (!required) stop("'data' must have columns: state, mean (and optionally lower, upper).")

  if (sort_states) data <- data[order(data$mean), ]
  data$state <- factor(data$state, levels = data$state)

  states <- levels(data$state)
  default_pal <- c(
    "#E69F00", "#56B4E9", "#009E73", "#F0E442",
    "#0072B2", "#D55E00", "#CC79A7", "#999999",
    "#44AA99", "#882255", "#117733", "#DDCC77"
  )
  if (is.null(palette))
    palette <- setNames(default_pal[seq_along(states)], states)

  has_ci <- all(c("lower", "upper") %in% names(data))

  p <- ggplot2::ggplot(data, ggplot2::aes(x = mean, y = state, colour = state))

  if (has_ci) {
    p <- p + ggplot2::geom_linerange(
      ggplot2::aes(xmin = lower, xmax = upper),
      linewidth = 0.9, na.rm = TRUE)
  }

  p +
    ggplot2::geom_point(size = 3.5, na.rm = TRUE) +
    ggplot2::scale_colour_manual(values = palette, guide = "none") +
    ggplot2::labs(title = title, x = xlab, y = NULL) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor   = ggplot2::element_blank(),
      axis.text.y        = ggplot2::element_text(hjust = 1)
    )
}


# ── Proportions over time (DTA posterior) ─────────────────────────────────────

#' Compute lineage proportions over time from a BEAST DTA posterior sample
#'
#' At each time slice, counts the proportion of phylogenetic lineages in each
#' discrete trait state, averaged across a posterior sample of trees.
#'
#' This is the DTA equivalent of PACT's \code{skyline proportions} (Structured
#' Coalescent). The result can be passed directly to \code{\link{plot_dta_proportions}}.
#'
#' @param file Path to the BEAST \code{.trees} NEXUS file.
#' @param trait Character. Discrete trait key (e.g. \code{"region"}).
#' @param times Numeric vector of time points at which to evaluate proportions.
#'   If \code{NULL}, a grid from the oldest to the most recent tip date is used
#'   with step size \code{step}.
#' @param step Numeric. Time step for the auto-generated grid when \code{times}
#'   is \code{NULL}. Default \code{0.1}.
#' @param burnin Integer. Number of trees to discard from the start.
#' @param tip_dates Named numeric vector of decimal tip dates. Auto-parsed if
#'   \code{NULL}.
#' @param mrsd Numeric or character. Most recent sampling date for time
#'   anchoring (decimal year or ISO date string). See
#'   \code{\link{read_beast_posterior}}.
#'
#' @return A data frame with columns \code{time}, \code{state}, \code{lower}
#'   (2.5th percentile), \code{mean}, \code{upper} (97.5th percentile), and
#'   \code{n_lineages} (mean number of lineages at that time point).
#'
#' @examples
#' \donttest{
#' f <- system.file("extdata", "example_dta_trees",
#'                  "rsva_resample.trees", package = "pactR")
#' prop <- dta_proportions(f, trait = "region", step = 0.25)
#' plot_dta_proportions(prop)
#' }
#'
#' @seealso \code{\link{plot_dta_proportions}}, \code{\link{read_beast_posterior}}
#'
#' @export
dta_proportions <- function(file, trait = "region", times = NULL, step = 0.1,
                             burnin = 0L, tip_dates = NULL, mrsd = NULL) {
  message("Reading trees from: ", file)
  post <- read_beast_posterior(file, burnin = burnin, trait = trait,
                                tip_dates = tip_dates, mrsd = mrsd)
  trees <- post$trees

  # Determine time grid from available node times
  if (is.null(times)) {
    all_node_times <- unlist(lapply(trees, function(df)
      c(df$parent_time, df$child_time)), use.names = FALSE)
    all_node_times <- all_node_times[is.finite(all_node_times)]
    t_min <- floor(min(all_node_times) * 10) / 10
    t_max <- ceiling(max(all_node_times) * 10) / 10
    times <- seq(t_min, t_max, by = step)
    if (!post$absolute_time)
      message("Time axis: root-relative branch lengths. Supply 'mrsd' for calendar years.")
  }

  message("Computing persistence at ", length(times), " time points across ",
          length(trees), " trees...")

  # For each tree × time: count proportion in each state
  all_states <- unique(unlist(lapply(trees, `[[`, "region")))
  all_states <- sort(all_states[!is.na(all_states)])

  # Results matrix: [n_trees × n_times × n_states]
  # Store proportions per tree per time point
  prop_array <- array(
    NA_real_,
    dim = c(length(trees), length(times), length(all_states)),
    dimnames = list(NULL, as.character(times), all_states)
  )
  n_lineages_mat <- matrix(NA_real_, nrow = length(trees), ncol = length(times))

  for (i in seq_along(trees)) {
    df <- trees[[i]]
    df <- df[!is.na(df$parent_time) & !is.na(df$child_time), ]
    for (j in seq_along(times)) {
      t     <- times[j]
      t_end <- t + step
      # Branch-length overlap with window [t, t+step] — matches C++ trimEnds + getLabelPro
      overlap    <- pmax(0, pmin(df$child_time, t_end) - pmax(df$parent_time, t))
      total_len  <- sum(overlap)
      n_lineages_mat[i, j] <- sum(overlap > 0)
      if (total_len > 0) {
        for (s in seq_along(all_states)) {
          state_len          <- sum(overlap[!is.na(df$region) & df$region == all_states[s]])
          prop_array[i, j, s] <- state_len / total_len
        }
      }
    }
  }

  # Summarise across trees: mean + 95% CI per state × time
  # Report at window midpoint (t + step/2) to match C++ PACT output convention
  rows <- vector("list", length(times) * length(all_states))
  k <- 0L
  for (j in seq_along(times)) {
    for (s in seq_along(all_states)) {
      vals <- prop_array[, j, s]
      vals <- vals[!is.na(vals)]
      k <- k + 1L
      rows[[k]] <- data.frame(
        time       = times[j] + step / 2,
        state      = all_states[s],
        lower      = if (length(vals) > 1L) unname(stats::quantile(vals, 0.025)) else vals[1L],
        mean       = if (length(vals) > 0L) mean(vals)                           else NA_real_,
        upper      = if (length(vals) > 1L) unname(stats::quantile(vals, 0.975)) else vals[1L],
        n_lineages = mean(n_lineages_mat[, j], na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}


#' Convert PACT skyline output to proportion-over-time format
#'
#' Reshapes \code{\link{read_pact_skylines}} output into the data frame format
#' used by \code{\link{plot_dta_proportions}}.
#'
#' @param skylines Data frame from \code{\link{read_pact_skylines}}.
#' @param prefix Statistic name prefix that identifies proportion rows.
#'   Defaults to \code{"pro_"} (PACT's \code{skyline proportions} output).
#'
#' @return A data frame with columns \code{time}, \code{state}, \code{lower},
#'   \code{mean}, \code{upper}.
#'
#' @examples
#' \donttest{
#' ex      <- pact_example_files()
#' results <- run_pact(ex$trees, ex$param, output_dir = tempdir())
#' sky     <- read_pact_skylines(results$skylines)
#' prop    <- pact_to_proportions(sky)
#' plot_dta_proportions(prop)
#' }
#'
#' @export
pact_to_proportions <- function(skylines, prefix = "pro_") {
  rows <- skylines[grepl(paste0("^", prefix), skylines$statistic), ]
  if (nrow(rows) == 0L)
    stop("No rows matching prefix '", prefix, "' found in skylines data frame.")
  data.frame(
    time  = rows$time,
    state = sub(paste0("^", prefix), "", rows$statistic),
    lower = rows$lower,
    mean  = rows$mean,
    upper = rows$upper,
    stringsAsFactors = FALSE
  )
}


# ── Persistence plot ──────────────────────────────────────────────────────────

#' Plot lineage proportions over time
#'
#' Visualises the proportion of lineages in each discrete state (e.g. region)
#' over time, with credible-interval ribbons. Accepts output from
#' \code{\link{dta_proportions}}, \code{\link{pact_to_proportions}}, or any
#' data frame with columns \code{time}, \code{state}, \code{mean}.
#'
#' @param data Data frame with columns \code{time}, \code{state}, \code{mean},
#'   and optionally \code{lower} and \code{upper} for credible intervals.
#' @param type Character. \code{"ribbon"} (default) draws lines + shaded CI;
#'   \code{"area"} draws a filled stacked area; \code{"line"} draws lines only.
#' @param palette Named character vector mapping state names to colours. If
#'   \code{NULL}, uses a colour-blind-friendly palette.
#' @param title Character. Plot title.
#' @param xlab,ylab Axis labels.
#' @param ribbon_alpha Numeric (0–1). Transparency of CI ribbons.
#' @param smooth Logical. If \code{TRUE} (default), applies a mild LOESS
#'   smoother to reduce jaggedness in the mean line.
#' @param smooth_span LOESS span parameter (used when \code{smooth = TRUE}).
#'
#' @return A \code{ggplot2} object.
#'
#' @examples
#' \donttest{
#' f    <- system.file("extdata", "example_dta_trees",
#'                     "rsva_resample.trees", package = "pactR")
#' prop <- dta_proportions(f, trait = "region", step = 0.25)
#' plot_dta_proportions(prop)
#' plot_dta_proportions(prop, type = "area")
#' }
#'
#' @seealso \code{\link{dta_proportions}}, \code{\link{pact_to_proportions}}
#'
#' @export
plot_dta_proportions <- function(data,
                                  type         = c("ribbon", "area", "line"),
                                  palette      = NULL,
                                  title        = "Lineage proportions over time",
                                  xlab         = "Time",
                                  ylab         = "Proportion of lineages",
                                  ribbon_alpha = 0.2,
                                  smooth       = TRUE,
                                  smooth_span  = 0.15) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is required. Install with install.packages('ggplot2').")

  type <- match.arg(type)

  # Validate / rename columns
  required <- "time" %in% names(data) && "state" %in% names(data) &&
              "mean" %in% names(data)
  if (!required)
    stop("'data' must have columns: time, state, mean (and optionally lower, upper).")

  has_ci <- all(c("lower", "upper") %in% names(data))

  states <- sort(unique(data$state))

  # Default palette: colour-blind friendly (Okabe-Ito extended)
  default_pal <- c(
    "#E69F00", "#56B4E9", "#009E73", "#F0E442",
    "#0072B2", "#D55E00", "#CC79A7", "#999999",
    "#44AA99", "#882255", "#117733", "#DDCC77"
  )
  if (is.null(palette)) {
    palette <- setNames(default_pal[seq_along(states)], states)
  }

  # Optional smoothing of mean (and CI) lines
  if (smooth) {
    data <- do.call(rbind, lapply(states, function(s) {
      sub_d <- data[data$state == s, ]
      sub_d <- sub_d[order(sub_d$time), ]
      if (nrow(sub_d) < 5L) return(sub_d)
      sub_d$mean <- pmax(0, pmin(1, stats::predict(
        stats::loess(mean ~ time, data = sub_d, span = smooth_span), sub_d)))
      if (has_ci) {
        sub_d$lower <- pmax(0, stats::predict(
          stats::loess(lower ~ time, data = sub_d, span = smooth_span), sub_d))
        sub_d$upper <- pmin(1, stats::predict(
          stats::loess(upper ~ time, data = sub_d, span = smooth_span), sub_d))
      }
      sub_d
    }))
  }

  p <- ggplot2::ggplot(data, ggplot2::aes(x = time, colour = state, fill = state))

  if (type == "ribbon") {
    if (has_ci) {
      p <- p + ggplot2::geom_ribbon(
        ggplot2::aes(ymin = lower, ymax = upper),
        alpha = ribbon_alpha, colour = NA, na.rm = TRUE)
    }
    p <- p + ggplot2::geom_line(ggplot2::aes(y = mean), linewidth = 0.8, na.rm = TRUE)

  } else if (type == "area") {
    # Stacked area — fill in proportions that sum to 1
    # Drop time points where every state has NA mean (avoids non-finite warnings)
    valid_times <- unique(data$time[!is.na(data$mean)])
    area_data   <- data[data$time %in% valid_times, ]
    area_data$mean[!is.finite(area_data$mean)] <- 0
    p <- ggplot2::ggplot(area_data,
           ggplot2::aes(x = time, y = mean, fill = state)) +
         ggplot2::geom_area(position = "fill", alpha = 0.85, colour = "white",
                            linewidth = 0.3)

  } else {  # "line"
    p <- p + ggplot2::geom_line(ggplot2::aes(y = mean), linewidth = 0.8, na.rm = TRUE)
  }

  out <- p + ggplot2::scale_fill_manual(values = palette, name = "State") +
             ggplot2::scale_y_continuous(limits = c(0, 1),
                                         labels = scales::percent_format()) +
             ggplot2::labs(title = title, x = xlab, y = ylab) +
             ggplot2::theme_minimal(base_size = 12) +
             ggplot2::theme(legend.position  = "right",
                            panel.grid.minor = ggplot2::element_blank())
  # Add colour scale only for plot types that use a colour aesthetic
  if (type != "area") {
    out <- out + ggplot2::scale_colour_manual(values = palette, name = "State")
  }
  out
}


# ── Lineage through time ──────────────────────────────────────────────────────

#' Compute lineage-through-time (LTT) from a BEAST posterior sample
#'
#' Counts the number of phylogenetic lineages alive at each point in time,
#' averaged across a posterior sample of trees. Provides the mean and 95%
#' credible interval.
#'
#' @param file Path to the BEAST \code{.trees} NEXUS file.
#' @param times Numeric vector of time points. Auto-generated if \code{NULL}.
#' @param step Time step for auto-generated grid. Default \code{0.1}.
#' @param burnin Integer. Trees to discard.
#' @param tip_dates Named numeric decimal dates. Auto-parsed if \code{NULL}.
#' @param mrsd Numeric or character most-recent sampling date (decimal year or
#'   \code{"YYYY-MM-DD"}). Used for absolute time anchoring when tip names
#'   contain no parseable dates.
#'
#' @return A data frame with columns \code{time}, \code{lower}, \code{mean},
#'   \code{upper}.
#'
#' @examples
#' \donttest{
#' f   <- system.file("extdata", "example_dta_trees",
#'                    "rsva_resample.trees", package = "pactR")
#' ltt <- lineage_through_time(f, step = 0.25)
#' plot_ltt(ltt)
#' }
#'
#' @seealso \code{\link{plot_ltt}}
#'
#' @export
lineage_through_time <- function(file, times = NULL, step = 0.1,
                                  burnin = 0L, tip_dates = NULL, mrsd = NULL) {
  post  <- read_beast_posterior(file, burnin = burnin, trait = NULL,
                                 tip_dates = tip_dates, mrsd = mrsd)
  trees <- post$trees

  if (is.null(times)) {
    all_node_times <- unlist(lapply(trees, function(df)
      c(df$parent_time, df$child_time)), use.names = FALSE)
    all_node_times <- all_node_times[is.finite(all_node_times)]
    t_min <- floor(min(all_node_times) * 10) / 10
    t_max <- ceiling(max(all_node_times) * 10) / 10
    times <- seq(t_min, t_max, by = step)
  }

  counts <- matrix(NA_real_, nrow = length(trees), ncol = length(times))
  for (i in seq_along(trees)) {
    df <- trees[[i]]
    df <- df[!is.na(df$parent_time) & !is.na(df$child_time), ]
    for (j in seq_along(times)) {
      t <- times[j]
      counts[i, j] <- sum(df$parent_time <= t & df$child_time > t)
    }
  }

  data.frame(
    time  = times,
    lower = unname(apply(counts, 2L, stats::quantile, 0.025, na.rm = TRUE)),
    mean  = unname(colMeans(counts, na.rm = TRUE)),
    upper = unname(apply(counts, 2L, stats::quantile, 0.975, na.rm = TRUE))
  )
}


#' Plot a lineage-through-time (LTT) curve
#'
#' @param ltt_data Data frame from \code{\link{lineage_through_time}} with
#'   columns \code{time}, \code{mean}, and optionally \code{lower}/\code{upper}.
#' @param log_y Logical. If \code{TRUE}, uses a log10 y-axis.
#' @param colour Line colour. Default \code{"steelblue"}.
#' @param ribbon_alpha Ribbon transparency (0–1).
#' @param title Plot title.
#'
#' @return A \code{ggplot2} object.
#'
#' @export
plot_ltt <- function(ltt_data,
                     log_y        = FALSE,
                     colour       = "steelblue",
                     ribbon_alpha = 0.25,
                     title        = "Lineage through time") {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is required.")

  has_ci <- all(c("lower", "upper") %in% names(ltt_data))
  p <- ggplot2::ggplot(ltt_data, ggplot2::aes(x = time, y = mean))
  if (has_ci) {
    p <- p + ggplot2::geom_ribbon(
      ggplot2::aes(ymin = lower, ymax = upper),
      fill = colour, alpha = ribbon_alpha)
  }
  p <- p + ggplot2::geom_line(colour = colour, linewidth = 0.9)
  if (log_y) p <- p + ggplot2::scale_y_log10()
  p + ggplot2::labs(title = title, x = "Time", y = "Number of lineages") +
      ggplot2::theme_minimal(base_size = 12)
}


# ── Tree visualisation (ggtree) ───────────────────────────────────────────────

#' Plot a BEAST-annotated tree with ggtree
#'
#' Reads a BEAST MCC or single tree (NEXUS format) and draws it with
#' \code{ggtree}, optionally colouring branches and tips by a discrete
#' trait annotation.
#'
#' @param file Path to a BEAST NEXUS tree file (e.g., MCC tree from
#'   TreeAnnotator, or a single tree from a posterior).
#' @param trait Character. Branch/node annotation to colour by (e.g.
#'   \code{"region"}). Set to \code{NULL} for an uncoloured tree.
#' @param tip_labels Logical. Show tip labels?
#' @param tip_size Numeric. Tip label size (passed to
#'   \code{ggtree::geom_tiplab}).
#' @param palette Named colour vector for trait states. Auto-generated if
#'   \code{NULL}.
#' @param layout Character passed to \code{ggtree::ggtree}: \code{"rectangular"}
#'   (default), \code{"circular"}, \code{"fan"}, etc.
#' @param mrsd Character or numeric. Most recent sampling date used to display
#'   absolute time on the x-axis. E.g. \code{"2017-12-31"} or \code{2017.997}.
#'   If \code{NULL}, relative time is shown.
#' @param ... Additional arguments passed to \code{ggtree::ggtree()}.
#'
#' @return A \code{ggtree}/\code{ggplot2} object.
#'
#' @examples
#' \donttest{
#' f <- system.file("extdata", "example_dta_trees",
#'                  "rsva_resample.trees", package = "pactR")
#' # Plot the last tree in the file
#' plot_tree_beast(f, trait = "region")
#' }
#'
#' @export
plot_tree_beast <- function(file,
                             trait      = NULL,
                             tip_labels = FALSE,
                             tip_size   = 2,
                             palette    = NULL,
                             layout     = "rectangular",
                             mrsd       = NULL,
                             ...) {
  for (pkg in c("ggtree", "treeio", "ggplot2")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      stop("Package '", pkg, "' is required. ",
           if (pkg %in% c("ggtree","treeio"))
             "Install with: BiocManager::install(c('ggtree','treeio'))"
           else
             paste0("Install with: install.packages('", pkg, "')"))
  }

  td <- treeio::read.beast(file)

  p <- if (!is.null(mrsd)) {
    ggtree::ggtree(td, layout = layout, mrsd = mrsd, ...)
  } else {
    ggtree::ggtree(td, layout = layout, ...)
  }

  if (!is.null(trait)) {
    # Map trait to branch colour
    p <- p + ggtree::aes(colour = .data[[trait]])

    # Build palette if not provided
    if (is.null(palette)) {
      states <- sort(unique(td@data[[trait]]))
      states <- states[!is.na(states)]
      default_pal <- c(
        "#E69F00","#56B4E9","#009E73","#F0E442",
        "#0072B2","#D55E00","#CC79A7","#999999",
        "#44AA99","#882255","#117733","#DDCC77"
      )
      palette <- setNames(default_pal[seq_along(states)], states)
    }
    p <- p + ggplot2::scale_colour_manual(values = palette, na.value = "grey80",
                                           name = trait)
  }

  if (tip_labels) {
    p <- p + ggtree::geom_tiplab(size = tip_size)
  }

  p + ggplot2::theme(legend.position = "right")
}


# ── Root-to-tip regression ────────────────────────────────────────────────────

#' Root-to-tip distance regression (molecular clock check)
#'
#' Computes root-to-tip distances for each tip of a phylogenetic tree and
#' regresses them on sampling dates. Returns a data frame and an optional
#' plot. This is a standard molecular clock diagnostic.
#'
#' @param file Path to a tree file (NEXUS or NEWICK). For BEAST NEXUS files
#'   with multiple trees, the last tree is used unless \code{tree_index} is
#'   specified.
#' @param tip_dates Named numeric vector of decimal sampling dates. If
#'   \code{NULL}, dates are parsed from tip names using
#'   \code{\link{parse_tip_dates}}.
#' @param tree_index Integer. Which tree to use when the file contains
#'   multiple trees. Default: last tree (\code{-1L}).
#' @param plot Logical. If \code{TRUE} (default), prints a ggplot2 scatter with
#'   regression line and returns the data frame invisibly.
#'
#' @return A data frame with columns \code{tip}, \code{date},
#'   \code{root_to_tip}. Regression statistics (R², slope, intercept) are
#'   printed if \code{plot = TRUE}.
#'
#' @examples
#' \donttest{
#' f   <- system.file("extdata", "example_dta_trees",
#'                    "rsva_resample.trees", package = "pactR")
#' rtt <- root_to_tip(f)
#' }
#'
#' @export
root_to_tip <- function(file, tip_dates = NULL, tree_index = -1L, plot = TRUE) {
  if (!requireNamespace("ape", quietly = TRUE))
    stop("Package 'ape' is required.")

  # Read tree(s) — try NEXUS first
  trees <- tryCatch(ape::read.nexus(file), error = function(e) NULL)
  if (is.null(trees)) trees <- ape::read.tree(file)

  # Select single tree
  phy <- if (inherits(trees, "multiPhylo")) {
    idx <- if (tree_index == -1L) length(trees) else tree_index
    trees[[idx]]
  } else {
    trees
  }

  # Root-to-tip distances (one value per tip)
  dist_mat <- ape::dist.nodes(phy)
  root_node <- length(phy$tip.label) + 1L
  rtt_dist  <- dist_mat[seq_along(phy$tip.label), root_node]

  # Tip dates
  if (is.null(tip_dates)) {
    tip_dates <- parse_tip_dates(phy$tip.label)
  }
  dates <- tip_dates[phy$tip.label]

  df <- data.frame(
    tip        = phy$tip.label,
    date       = dates,
    root_to_tip = rtt_dist,
    stringsAsFactors = FALSE
  )
  df <- df[!is.na(df$date), ]

  fit    <- stats::lm(root_to_tip ~ date, data = df)
  r2     <- summary(fit)$r.squared
  slope  <- stats::coef(fit)["date"]
  intercept <- stats::coef(fit)["(Intercept)"]
  message(sprintf("Root-to-tip regression: R² = %.4f | slope = %.2e | intercept = %.4f",
                  r2, slope, intercept))

  if (plot) {
    if (!requireNamespace("ggplot2", quietly = TRUE))
      stop("Package 'ggplot2' is required for plotting.")

    p <- ggplot2::ggplot(df, ggplot2::aes(x = date, y = root_to_tip)) +
      ggplot2::geom_point(alpha = 0.5, colour = "steelblue") +
      ggplot2::geom_smooth(method = "lm", formula = y ~ x,
                            colour = "firebrick", linewidth = 0.8, se = TRUE) +
      ggplot2::annotate("text",
        x = min(df$date, na.rm = TRUE),
        y = max(df$root_to_tip, na.rm = TRUE),
        label = sprintf("R² = %.4f\nslope = %.2e", r2, slope),
        hjust = 0, vjust = 1, size = 3.5) +
      ggplot2::labs(
        title = "Root-to-tip regression",
        x = "Sampling date", y = "Root-to-tip distance") +
      ggplot2::theme_minimal(base_size = 12)
    print(p)
    return(invisible(df))
  }
  df
}


# ── Pairwise TMRCA ────────────────────────────────────────────────────────────

#' Compute pairwise TMRCA matrix
#'
#' Returns the time to most recent common ancestor (TMRCA) for every pair of
#' tips in a tree, expressed in the same units as branch lengths.
#'
#' @param file Path to a tree file (NEXUS or NEWICK). For multi-tree files,
#'   uses the last tree.
#' @param tree_index Integer. Which tree to use (\code{-1L} = last).
#' @param tip_dates Named numeric decimal dates. Used to convert root-to-node
#'   heights to absolute times. Auto-parsed if \code{NULL}.
#'
#' @return A symmetric matrix of pairwise TMRCAs with tip labels as
#'   row/column names.
#'
#' @examples
#' \donttest{
#' f <- system.file("extdata", "example_dta_trees",
#'                  "rsva_resample.trees", package = "pactR")
#' mat <- pairwise_tmrca(f)
#' }
#'
#' @export
pairwise_tmrca <- function(file, tree_index = -1L, tip_dates = NULL) {
  trees <- tryCatch(ape::read.nexus(file), error = function(e) NULL)
  if (is.null(trees)) trees <- ape::read.tree(file)
  phy <- if (inherits(trees, "multiPhylo")) {
    idx <- if (tree_index == -1L) length(trees) else tree_index
    trees[[idx]]
  } else { trees }

  if (is.null(tip_dates)) tip_dates <- parse_tip_dates(phy$tip.label)

  dist_mat  <- ape::dist.nodes(phy)
  n_tips    <- length(phy$tip.label)
  root_node <- n_tips + 1L
  heights   <- dist_mat[seq_len(n_tips), root_node]

  most_recent  <- max(tip_dates[phy$tip.label], na.rm = TRUE)
  max_height   <- max(heights)
  root_abs     <- most_recent - max_height
  tip_abs_time <- root_abs + heights   # absolute time of each tip

  # For two tips, their TMRCA is at:
  # tmrca(i,j) = root_abs + (root_to_mrca distance)
  # dist.nodes gives total distance between nodes; for two tips:
  # dist(tip_i, tip_j) = height(tip_i) + height(tip_j) - 2 * height(MRCA)
  # => height(MRCA) = (height(tip_i) + height(tip_j) - dist(tip_i,tip_j)) / 2
  tip_tip_dist <- dist_mat[seq_len(n_tips), seq_len(n_tips)]

  tmrca_mat <- matrix(NA_real_, n_tips, n_tips,
                      dimnames = list(phy$tip.label, phy$tip.label))
  for (i in seq_len(n_tips)) {
    for (j in seq_len(n_tips)) {
      if (i == j) {
        tmrca_mat[i, j] <- tip_abs_time[i]
      } else {
        mrca_height <- (heights[i] + heights[j] - tip_tip_dist[i, j]) / 2
        tmrca_mat[i, j] <- root_abs + mrca_height
      }
    }
  }
  tmrca_mat
}


# ── Find MRCA ─────────────────────────────────────────────────────────────────

#' Find the MRCA of a set of tips
#'
#' Identifies the most recent common ancestor (MRCA) node of a specified set
#' of tip labels in a tree.
#'
#' @param file Path to a tree file (NEXUS or NEWICK).
#' @param tips Character vector of tip labels whose MRCA to find.
#' @param tree_index Integer. Which tree to use for multi-tree files.
#'
#' @return A named list with:
#' \describe{
#'   \item{\code{node}}{Integer ape node index of the MRCA.}
#'   \item{\code{n_descendants}}{Number of tips descending from the MRCA.}
#'   \item{\code{subtree_tips}}{Tip labels in the subtree rooted at the MRCA.}
#' }
#'
#' @export
find_mrca <- function(file, tips, tree_index = -1L) {
  trees <- tryCatch(ape::read.nexus(file), error = function(e) NULL)
  if (is.null(trees)) trees <- ape::read.tree(file)
  phy <- if (inherits(trees, "multiPhylo")) {
    idx <- if (tree_index == -1L) length(trees) else tree_index
    trees[[idx]]
  } else { trees }

  mrca_node   <- ape::getMRCA(phy, tips)
  subtree     <- ape::extract.clade(phy, mrca_node)

  list(
    node          = mrca_node,
    n_descendants = length(subtree$tip.label),
    subtree_tips  = subtree$tip.label
  )
}


# ── Example DTA file helper ───────────────────────────────────────────────────

#' Get path to the bundled example DTA posterior tree file
#'
#' Returns the path to the RSV-A DTA posterior trees bundled with the package.
#' These trees have a \code{region} discrete trait from a Discrete Trait
#' Analysis (DTA) run in BEAST.
#'
#' @return Character path to \code{rsva_resample.trees}.
#'
#' @examples
#' f <- pact_example_dta()
#' cat("File:", f, "\n")
#' cat("Trees:", length(readLines(f)[grep("^\\s*tree ", readLines(f))]), "\n")
#'
#' @export
pact_example_dta <- function() {
  f <- system.file("extdata", "example_dta_trees",
                   "global_HPAI.trees", package = "pactR")
  if (!nzchar(f)) stop("DTA example file not found. Is the pact package installed?")
  f
}
