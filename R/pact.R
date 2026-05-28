# ── Parameter file builder ────────────────────────────────────────────────────

#' Build a PACT parameter file from R arguments
#'
#' Converts named R arguments into the line-based text format expected by
#' \code{in.param}. Used internally by \code{\link{run_pact}} but also
#' exported so users can inspect or save the generated parameter content.
#'
#' @inheritParams run_pact
#'
#' @return A character vector of parameter lines (one directive per element).
#'
#' @seealso \code{\link{run_pact}}, \code{\link{pact_parameter_reference}}
#'
#' @examples
#' lines <- build_pact_params(
#'   push_times_back        = c(2002, 2008),
#'   renew_trunk            = 1,
#'   summary_tmrca          = TRUE,
#'   summary_coal_rates     = TRUE,
#'   skyline_settings       = c(2002, 2008, 0.01),
#'   skyline_proportions    = TRUE
#' )
#' cat(lines, sep = "\n")
#'
#' @export
build_pact_params <- function(
    # General
    burnin                        = NULL,
    # Tree manipulation
    push_times_back               = NULL,
    reduce_tips                   = NULL,
    renew_trunk                   = NULL,
    prune_to_trunk                = FALSE,
    prune_to_label                = NULL,
    prune_to_tips                 = NULL,
    remove_tips                   = NULL,
    prune_to_time                 = NULL,
    pad_migration_events          = FALSE,
    collapse_labels               = FALSE,
    trim_ends                     = NULL,
    section_tree                  = NULL,
    time_slice                    = NULL,
    rotate                        = NULL,
    accumulate                    = FALSE,
    add_tail                      = NULL,
    ordering                      = NULL,
    # Tree output
    print_tree                    = FALSE,
    print_circular_tree           = FALSE,
    print_all_trees               = FALSE,
    # Summary statistics
    summary_tmrca                 = FALSE,
    summary_length                = FALSE,
    summary_root_proportions      = FALSE,
    summary_proportions           = FALSE,
    summary_coal_rates            = FALSE,
    summary_mig_rates             = FALSE,
    summary_sub_rates             = FALSE,
    summary_diversity             = FALSE,
    summary_fst                   = FALSE,
    summary_tajima_d              = FALSE,
    summary_diffusion_coefficient = FALSE,
    summary_drift_rate            = FALSE,
    summary_persistence           = FALSE,
    # Skyline statistics
    skyline_settings              = NULL,
    skyline_tmrca                 = FALSE,
    skyline_length                = FALSE,
    skyline_proportions           = FALSE,
    skyline_coal_rates            = FALSE,
    skyline_mig_rates             = FALSE,
    skyline_pro_history_from_tips = FALSE,
    skyline_diversity             = FALSE,
    skyline_fst                   = FALSE,
    skyline_tajima_d              = FALSE,
    skyline_timetofix             = FALSE,
    skyline_xmean                 = FALSE,
    skyline_ymean                 = FALSE,
    skyline_xdrift                = FALSE,
    skyline_ratemean              = FALSE,
    skyline_xtrunkdiff            = FALSE,
    skyline_locsample             = FALSE,
    skyline_locgrid               = FALSE,
    skyline_drift_rate_from_tips  = FALSE,
    # Tip statistics
    tips_time_to_trunk            = FALSE,
    x_loc_history                 = NULL,
    y_loc_history                 = NULL,
    # Pair statistics
    pairs_diversity               = NULL
) {
  lines <- character(0)

  # General
  if (!is.null(burnin))
    lines <- c(lines, paste("burnin", burnin[1]))

  # Tree manipulation
  if (!is.null(push_times_back))
    lines <- c(lines, paste("push times back", paste(push_times_back, collapse = " ")))
  if (!is.null(reduce_tips))
    lines <- c(lines, paste("reduce tips", reduce_tips[1]))
  if (!is.null(renew_trunk))
    lines <- c(lines, paste("renew trunk", renew_trunk[1]))
  if (isTRUE(prune_to_trunk))
    lines <- c(lines, "prune to trunk")
  if (!is.null(prune_to_label))
    lines <- c(lines, paste("prune to label", prune_to_label[1]))
  if (!is.null(prune_to_tips))
    lines <- c(lines, paste("prune to tips", paste(prune_to_tips, collapse = " ")))
  if (!is.null(remove_tips))
    lines <- c(lines, paste("remove tips", paste(remove_tips, collapse = " ")))
  if (!is.null(prune_to_time))
    lines <- c(lines, paste("prune to time", paste(prune_to_time[1:2], collapse = " ")))
  if (isTRUE(pad_migration_events))
    lines <- c(lines, "pad migration events")
  if (isTRUE(collapse_labels))
    lines <- c(lines, "collapse labels")
  if (!is.null(trim_ends))
    lines <- c(lines, paste("trim ends", paste(trim_ends[1:2], collapse = " ")))
  if (!is.null(section_tree))
    lines <- c(lines, paste("section tree", paste(section_tree[1:3], collapse = " ")))
  if (!is.null(time_slice))
    lines <- c(lines, paste("time slice", time_slice[1]))
  if (!is.null(rotate))
    lines <- c(lines, paste("rotate", rotate[1]))
  if (isTRUE(accumulate))
    lines <- c(lines, "accumulate")
  if (!is.null(add_tail))
    lines <- c(lines, paste("add tail", add_tail[1]))
  if (!is.null(ordering))
    lines <- c(lines, paste("ordering", paste(ordering, collapse = " ")))

  # Tree output
  if (isTRUE(print_tree))
    lines <- c(lines, "print rule tree")
  if (isTRUE(print_circular_tree))
    lines <- c(lines, "print circular tree")
  if (isTRUE(print_all_trees))
    lines <- c(lines, "print all trees")

  # Summary statistics
  if (isTRUE(summary_tmrca))                 lines <- c(lines, "summary tmrca")
  if (isTRUE(summary_length))                lines <- c(lines, "summary length")
  if (isTRUE(summary_root_proportions))      lines <- c(lines, "summary root proportions")
  if (isTRUE(summary_proportions))           lines <- c(lines, "summary proportions")
  if (isTRUE(summary_coal_rates))            lines <- c(lines, "summary coal rates")
  if (isTRUE(summary_mig_rates))             lines <- c(lines, "summary mig rates")
  if (isTRUE(summary_sub_rates))             lines <- c(lines, "summary sub rates")
  if (isTRUE(summary_diversity))             lines <- c(lines, "summary diversity")
  if (isTRUE(summary_fst))                   lines <- c(lines, "summary fst")
  if (isTRUE(summary_tajima_d))              lines <- c(lines, "summary tajima d")
  if (isTRUE(summary_diffusion_coefficient)) lines <- c(lines, "summary diffusion coefficient")
  if (isTRUE(summary_drift_rate))            lines <- c(lines, "summary drift rate")
  if (isTRUE(summary_persistence))           lines <- c(lines, "summary persistence")

  # Skyline statistics — settings must come before the individual flags
  if (!is.null(skyline_settings))
    lines <- c(lines, paste("skyline settings", paste(skyline_settings[1:3], collapse = " ")))
  if (isTRUE(skyline_tmrca))                 lines <- c(lines, "skyline tmrca")
  if (isTRUE(skyline_length))                lines <- c(lines, "skyline length")
  if (isTRUE(skyline_proportions))           lines <- c(lines, "skyline proportions")
  if (isTRUE(skyline_coal_rates))            lines <- c(lines, "skyline coal rates")
  if (isTRUE(skyline_mig_rates))             lines <- c(lines, "skyline mig rates")
  if (isTRUE(skyline_pro_history_from_tips)) lines <- c(lines, "skyline pro history from tips")
  if (isTRUE(skyline_diversity))             lines <- c(lines, "skyline diversity")
  if (isTRUE(skyline_fst))                   lines <- c(lines, "skyline fst")
  if (isTRUE(skyline_tajima_d))              lines <- c(lines, "skyline tajima d")
  if (isTRUE(skyline_timetofix))             lines <- c(lines, "skyline time to fix")
  if (isTRUE(skyline_xmean))                 lines <- c(lines, "skyline xmean")
  if (isTRUE(skyline_ymean))                 lines <- c(lines, "skyline ymean")
  if (isTRUE(skyline_xdrift))               lines <- c(lines, "skyline xdrift")
  if (isTRUE(skyline_ratemean))              lines <- c(lines, "skyline ratemean")
  if (isTRUE(skyline_xtrunkdiff))            lines <- c(lines, "skyline xtrunkdiff")
  if (isTRUE(skyline_locsample))             lines <- c(lines, "skyline locsample")
  if (isTRUE(skyline_locgrid))               lines <- c(lines, "skyline locgrid")
  if (isTRUE(skyline_drift_rate_from_tips))  lines <- c(lines, "skyline drift rate from tips")

  # Tip statistics
  if (isTRUE(tips_time_to_trunk))
    lines <- c(lines, "tips time to trunk")
  if (!is.null(x_loc_history))
    lines <- c(lines, paste("tips x loc history", paste(x_loc_history[1:3], collapse = " ")))
  if (!is.null(y_loc_history))
    lines <- c(lines, paste("tips y loc history", paste(y_loc_history[1:3], collapse = " ")))

  # Pair statistics
  if (!is.null(pairs_diversity))
    lines <- c(lines, paste("pairs diversity", pairs_diversity[1]))

  lines
}


# ── Main analysis function ────────────────────────────────────────────────────

#' Run PACT analysis
#'
#' Runs the full PACT (Posterior Analysis of Coalescent Trees) analysis pipeline.
#' Parameters can be supplied via a \code{param_file}, as named R arguments, or
#' both. When both are given the param file is used as a base and any named
#' arguments are appended (and therefore take precedence for value parameters
#' such as \code{push_times_back}).
#'
#' @param trees_file Path to the input NEWICK tree file (e.g., produced by BEAST
#'   or Migrate). May contain log-probability annotations.
#' @param param_file Path to a PACT parameter file (\code{in.param} format).
#'   Optional when named arguments are supplied directly. See
#'   \code{\link{pact_parameter_reference}} for the full syntax.
#' @param output_dir Directory for output files. Created if needed. Defaults to
#'   the current working directory.
#' @param output_prefix Filename prefix for outputs. Defaults to \code{"out"},
#'   producing \code{out.stats}, \code{out.tips}, \code{out.skylines},
#'   \code{out.rules}, and \code{out.pairs}.
#'
#' @section Tree manipulation parameters:
#' \describe{
#'   \item{\code{burnin}}{Integer. Discard the first \emph{n} trees.}
#'   \item{\code{push_times_back}}{Numeric vector of length 1 or 2. Rescale node
#'     times so the most recent tip is at \code{stop}, or so tips span
#'     \code{[start, stop]}.}
#'   \item{\code{reduce_tips}}{Numeric in (0,1]. Randomly retain this proportion
#'     of tips.}
#'   \item{\code{renew_trunk}}{Numeric. Redefine the trunk lineage working
#'     backwards from all samples within this many years of the most recent tip.}
#'   \item{\code{prune_to_trunk}}{Logical. Reduce tree to the trunk lineage only.}
#'   \item{\code{prune_to_label}}{Character. Keep only tips with this label.}
#'   \item{\code{prune_to_tips}}{Character vector. Keep only the named tips.}
#'   \item{\code{remove_tips}}{Character vector. Remove the named tips.}
#'   \item{\code{prune_to_time}}{Numeric vector \code{c(start, stop)}. Keep
#'     branches dated within this window.}
#'   \item{\code{pad_migration_events}}{Logical. Add virtual migration events.}
#'   \item{\code{collapse_labels}}{Logical. Treat all tips as one population.}
#'   \item{\code{trim_ends}}{Numeric vector \code{c(start, stop)}. Retain only
#'     branches within this time window.}
#'   \item{\code{section_tree}}{Numeric vector \code{c(start, window, step)}.
#'     Slice tree into sections.}
#'   \item{\code{time_slice}}{Numeric. Reduce tree to ancestors alive at this
#'     time.}
#'   \item{\code{rotate}}{Numeric. Rotate spatial coordinates by this many
#'     degrees.}
#'   \item{\code{accumulate}}{Logical. Accumulate spatial locations.}
#'   \item{\code{add_tail}}{Numeric. Add a tail branch of this length.}
#'   \item{\code{ordering}}{Character vector. Set tip ordering for output.}
#' }
#'
#' @section Tree output parameters:
#' \describe{
#'   \item{\code{print_tree}}{Logical. Write the highest-posterior tree to
#'     \code{out.rules} in Mathematica rule-list format.}
#'   \item{\code{print_circular_tree}}{Logical. Write a circular-layout version.}
#'   \item{\code{print_all_trees}}{Logical. Write all trees to
#'     \code{trees/} directory.}
#' }
#'
#' @section Summary statistic parameters:
#' \describe{
#'   \item{\code{summary_tmrca}}{Logical. Compute TMRCA.}
#'   \item{\code{summary_length}}{Logical. Total branch length.}
#'   \item{\code{summary_root_proportions}}{Logical. Label proportions at root.}
#'   \item{\code{summary_proportions}}{Logical. Label proportions on trunk.}
#'   \item{\code{summary_coal_rates}}{Logical. Per-population coalescent rates.}
#'   \item{\code{summary_mig_rates}}{Logical. Directional migration rates.}
#'   \item{\code{summary_sub_rates}}{Logical. Mean substitution rate.}
#'   \item{\code{summary_diversity}}{Logical. Nucleotide diversity (pi).}
#'   \item{\code{summary_fst}}{Logical. F_ST between populations.}
#'   \item{\code{summary_tajima_d}}{Logical. Tajima's D.}
#'   \item{\code{summary_diffusion_coefficient}}{Logical. Spatial diffusion
#'     coefficient.}
#'   \item{\code{summary_drift_rate}}{Logical. Spatial drift rate.}
#'   \item{\code{summary_persistence}}{Logical. Lineage persistence.}
#' }
#'
#' @section Skyline statistic parameters:
#' \describe{
#'   \item{\code{skyline_settings}}{Numeric vector \code{c(start, stop, step)}.
#'     Required for any skyline output.}
#'   \item{\code{skyline_tmrca}}{Logical.}
#'   \item{\code{skyline_length}}{Logical.}
#'   \item{\code{skyline_proportions}}{Logical. Population proportions per
#'     window.}
#'   \item{\code{skyline_coal_rates}}{Logical.}
#'   \item{\code{skyline_mig_rates}}{Logical.}
#'   \item{\code{skyline_pro_history_from_tips}}{Logical.}
#'   \item{\code{skyline_diversity}}{Logical.}
#'   \item{\code{skyline_fst}}{Logical.}
#'   \item{\code{skyline_tajima_d}}{Logical.}
#'   \item{\code{skyline_timetofix}}{Logical.}
#'   \item{\code{skyline_xmean}}{Logical. Mean X location per window.}
#'   \item{\code{skyline_ymean}}{Logical. Mean Y location per window.}
#'   \item{\code{skyline_xdrift}}{Logical.}
#'   \item{\code{skyline_ratemean}}{Logical.}
#'   \item{\code{skyline_xtrunkdiff}}{Logical.}
#'   \item{\code{skyline_locsample}}{Logical.}
#'   \item{\code{skyline_locgrid}}{Logical.}
#'   \item{\code{skyline_drift_rate_from_tips}}{Logical.}
#' }
#'
#' @section Tip statistic parameters:
#' \describe{
#'   \item{\code{tips_time_to_trunk}}{Logical. Time from each tip to the trunk.}
#'   \item{\code{x_loc_history}}{Numeric vector \code{c(start, stop, step)}.
#'     X-location history per tip.}
#'   \item{\code{y_loc_history}}{Numeric vector \code{c(start, stop, step)}.
#'     Y-location history per tip.}
#' }
#'
#' @section Pair statistic parameters:
#' \describe{
#'   \item{\code{pairs_diversity}}{Numeric. Maximum time difference between tips
#'     for pairwise diversity computation.}
#' }
#'
#' @return A named list of file paths (\code{stats}, \code{tips},
#'   \code{skylines}, \code{rules}, \code{pairs}). Files not requested by the
#'   parameters will not exist on disk.
#'
#' @seealso \code{\link{build_pact_params}}, \code{\link{read_pact_stats}},
#'   \code{\link{read_pact_tips}}, \code{\link{read_pact_skylines}},
#'   \code{\link{pact_example_files}}, \code{\link{pact_parameter_reference}}
#'
#' @examples
#' \donttest{
#' ex <- pact_example_files()
#' outdir <- tempdir()
#'
#' # Using a param file
#' results <- run_pact(ex$trees, param_file = ex$param, output_dir = outdir)
#'
#' # Using named arguments only
#' results <- run_pact(
#'   trees_file          = ex$trees,
#'   output_dir          = outdir,
#'   push_times_back     = c(2002, 2008),
#'   renew_trunk         = 1,
#'   summary_tmrca       = TRUE,
#'   summary_coal_rates  = TRUE,
#'   summary_mig_rates   = TRUE,
#'   skyline_settings    = c(2002, 2008, 0.01),
#'   skyline_proportions = TRUE,
#'   tips_time_to_trunk  = TRUE
#' )
#'
#' # Using a param file as base and overriding one value
#' results <- run_pact(
#'   trees_file      = ex$trees,
#'   param_file      = ex$param,
#'   output_dir      = outdir,
#'   summary_fst     = TRUE   # add FST on top of what param_file requests
#' )
#'
#' stats <- read_pact_stats(results$stats)
#' head(stats)
#' }
#'
#' @export
run_pact <- function(
    trees_file,
    param_file                    = NULL,
    output_dir                    = ".",
    output_prefix                 = "out",
    # General
    burnin                        = NULL,
    # Tree manipulation
    push_times_back               = NULL,
    reduce_tips                   = NULL,
    renew_trunk                   = NULL,
    prune_to_trunk                = FALSE,
    prune_to_label                = NULL,
    prune_to_tips                 = NULL,
    remove_tips                   = NULL,
    prune_to_time                 = NULL,
    pad_migration_events          = FALSE,
    collapse_labels               = FALSE,
    trim_ends                     = NULL,
    section_tree                  = NULL,
    time_slice                    = NULL,
    rotate                        = NULL,
    accumulate                    = FALSE,
    add_tail                      = NULL,
    ordering                      = NULL,
    # Tree output
    print_tree                    = FALSE,
    print_circular_tree           = FALSE,
    print_all_trees               = FALSE,
    # Summary statistics
    summary_tmrca                 = FALSE,
    summary_length                = FALSE,
    summary_root_proportions      = FALSE,
    summary_proportions           = FALSE,
    summary_coal_rates            = FALSE,
    summary_mig_rates             = FALSE,
    summary_sub_rates             = FALSE,
    summary_diversity             = FALSE,
    summary_fst                   = FALSE,
    summary_tajima_d              = FALSE,
    summary_diffusion_coefficient = FALSE,
    summary_drift_rate            = FALSE,
    summary_persistence           = FALSE,
    # Skyline statistics
    skyline_settings              = NULL,
    skyline_tmrca                 = FALSE,
    skyline_length                = FALSE,
    skyline_proportions           = FALSE,
    skyline_coal_rates            = FALSE,
    skyline_mig_rates             = FALSE,
    skyline_pro_history_from_tips = FALSE,
    skyline_diversity             = FALSE,
    skyline_fst                   = FALSE,
    skyline_tajima_d              = FALSE,
    skyline_timetofix             = FALSE,
    skyline_xmean                 = FALSE,
    skyline_ymean                 = FALSE,
    skyline_xdrift                = FALSE,
    skyline_ratemean              = FALSE,
    skyline_xtrunkdiff            = FALSE,
    skyline_locsample             = FALSE,
    skyline_locgrid               = FALSE,
    skyline_drift_rate_from_tips  = FALSE,
    # Tip statistics
    tips_time_to_trunk            = FALSE,
    x_loc_history                 = NULL,
    y_loc_history                 = NULL,
    # Pair statistics
    pairs_diversity               = NULL
) {
  trees_file <- normalizePath(trees_file, mustWork = TRUE)

  # Build keyword-arg lines (empty if all defaults)
  arg_lines <- build_pact_params(
    burnin                        = burnin,
    push_times_back               = push_times_back,
    reduce_tips                   = reduce_tips,
    renew_trunk                   = renew_trunk,
    prune_to_trunk                = prune_to_trunk,
    prune_to_label                = prune_to_label,
    prune_to_tips                 = prune_to_tips,
    remove_tips                   = remove_tips,
    prune_to_time                 = prune_to_time,
    pad_migration_events          = pad_migration_events,
    collapse_labels               = collapse_labels,
    trim_ends                     = trim_ends,
    section_tree                  = section_tree,
    time_slice                    = time_slice,
    rotate                        = rotate,
    accumulate                    = accumulate,
    add_tail                      = add_tail,
    ordering                      = ordering,
    print_tree                    = print_tree,
    print_circular_tree           = print_circular_tree,
    print_all_trees               = print_all_trees,
    summary_tmrca                 = summary_tmrca,
    summary_length                = summary_length,
    summary_root_proportions      = summary_root_proportions,
    summary_proportions           = summary_proportions,
    summary_coal_rates            = summary_coal_rates,
    summary_mig_rates             = summary_mig_rates,
    summary_sub_rates             = summary_sub_rates,
    summary_diversity             = summary_diversity,
    summary_fst                   = summary_fst,
    summary_tajima_d              = summary_tajima_d,
    summary_diffusion_coefficient = summary_diffusion_coefficient,
    summary_drift_rate            = summary_drift_rate,
    summary_persistence           = summary_persistence,
    skyline_settings              = skyline_settings,
    skyline_tmrca                 = skyline_tmrca,
    skyline_length                = skyline_length,
    skyline_proportions           = skyline_proportions,
    skyline_coal_rates            = skyline_coal_rates,
    skyline_mig_rates             = skyline_mig_rates,
    skyline_pro_history_from_tips = skyline_pro_history_from_tips,
    skyline_diversity             = skyline_diversity,
    skyline_fst                   = skyline_fst,
    skyline_tajima_d              = skyline_tajima_d,
    skyline_timetofix             = skyline_timetofix,
    skyline_xmean                 = skyline_xmean,
    skyline_ymean                 = skyline_ymean,
    skyline_xdrift                = skyline_xdrift,
    skyline_ratemean              = skyline_ratemean,
    skyline_xtrunkdiff            = skyline_xtrunkdiff,
    skyline_locsample             = skyline_locsample,
    skyline_locgrid               = skyline_locgrid,
    skyline_drift_rate_from_tips  = skyline_drift_rate_from_tips,
    tips_time_to_trunk            = tips_time_to_trunk,
    x_loc_history                 = x_loc_history,
    y_loc_history                 = y_loc_history,
    pairs_diversity               = pairs_diversity
  )

  # Resolve the effective param file
  if (is.null(param_file) && length(arg_lines) == 0) {
    stop("Supply at least one of: 'param_file', or named parameter arguments.")
  }

  if (length(arg_lines) == 0) {
    # Param file only — use it directly
    effective_param <- normalizePath(param_file, mustWork = TRUE)
  } else {
    # Build a temp param file. If param_file is also given, prepend its content
    # so keyword args (appended after) take precedence for value parameters.
    base_lines <- character(0)
    if (!is.null(param_file)) {
      param_file <- normalizePath(param_file, mustWork = TRUE)
      base_lines <- readLines(param_file)
    }
    effective_param <- tempfile(fileext = ".param")
    writeLines(c(base_lines, arg_lines), effective_param)
    on.exit(unlink(effective_param), add = TRUE)
  }

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  output_dir <- normalizePath(output_dir, mustWork = TRUE)
  prefix <- file.path(output_dir, output_prefix)

  pact_run_cpp(trees_file, effective_param, prefix)

  invisible(list(
    stats    = paste0(prefix, ".stats"),
    tips     = paste0(prefix, ".tips"),
    skylines = paste0(prefix, ".skylines"),
    rules    = paste0(prefix, ".rules"),
    pairs    = paste0(prefix, ".pairs")
  ))
}


# ── Output readers ────────────────────────────────────────────────────────────

#' Read PACT summary statistics output
#'
#' Reads the \code{.stats} file produced by \code{\link{run_pact}} into a
#' data frame with columns \code{statistic}, \code{lower}, \code{mean}, and
#' \code{upper}.
#'
#' @param file Path to the \code{.stats} output file.
#' @return A data frame.
#' @seealso \code{\link{run_pact}}
#' @export
read_pact_stats <- function(file) {
  if (!file.exists(file)) stop("File not found: ", file)
  utils::read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
}

#' Read PACT tip statistics output
#'
#' Reads the \code{.tips} file produced by \code{\link{run_pact}} into a
#' data frame with columns \code{statistic}, \code{name}, \code{label},
#' \code{time}, \code{lower}, \code{mean}, and \code{upper}.
#'
#' @param file Path to the \code{.tips} output file.
#' @return A data frame.
#' @seealso \code{\link{run_pact}}
#' @export
read_pact_tips <- function(file) {
  if (!file.exists(file)) stop("File not found: ", file)
  utils::read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
}

#' Read PACT skyline statistics output
#'
#' Reads the \code{.skylines} file produced by \code{\link{run_pact}} into a
#' data frame with columns \code{statistic}, \code{time}, \code{lower},
#' \code{mean}, and \code{upper}.
#'
#' @param file Path to the \code{.skylines} output file.
#' @return A data frame.
#' @seealso \code{\link{run_pact}}
#' @export
read_pact_skylines <- function(file) {
  if (!file.exists(file)) stop("File not found: ", file)
  utils::read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
}

#' Read PACT pairwise statistics output
#'
#' Reads the \code{.pairs} file produced by \code{\link{run_pact}} into a
#' data frame with columns \code{statistic}, \code{nameA}, \code{nameB},
#' \code{lower}, \code{mean}, and \code{upper}.
#'
#' @param file Path to the \code{.pairs} output file.
#' @return A data frame.
#' @seealso \code{\link{run_pact}}
#' @export
read_pact_pairs <- function(file) {
  if (!file.exists(file)) stop("File not found: ", file)
  utils::read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
}


# ── Helper utilities ──────────────────────────────────────────────────────────

#' Get paths to bundled example input files
#'
#' Returns file paths to the example \code{in.trees} and \code{in.param} files
#' bundled with the package (an influenza phylogenetic analysis).
#'
#' @return A named list with elements \code{trees}, \code{param}, and
#'   \code{dir}.
#'
#' @examples
#' ex <- pact_example_files()
#' cat(readLines(ex$param), sep = "\n")
#'
#' @export
pact_example_files <- function() {
  dir <- system.file("extdata", "example", package = "pactR")
  if (!nzchar(dir)) stop("Example files not found. Is the pactR package installed?")
  list(
    trees = file.path(dir, "in.trees"),
    param = file.path(dir, "in.param"),
    dir   = dir
  )
}

#' Open the PACT parameter reference
#'
#' Prints the bundled \code{parameters.txt} file, which documents all available
#' parameter file options with descriptions and examples.
#'
#' @param print Logical. If \code{TRUE} (default) prints to the console.
#'   If \code{FALSE} returns the lines invisibly.
#'
#' @return The parameter reference as a character vector (invisibly when
#'   \code{print = TRUE}).
#'
#' @examples
#' pact_parameter_reference()
#'
#' @export
pact_parameter_reference <- function(print = TRUE) {
  f <- system.file("extdata", "parameters.txt", package = "pactR")
  if (!nzchar(f)) stop("parameters.txt not found. Is the pactR package installed?")
  lines <- readLines(f)
  if (print) cat(lines, sep = "\n")
  invisible(lines)
}

#' Open the PACT user manual
#'
#' Opens the bundled PDF user manual in the system default PDF viewer.
#'
#' @return The path to the PDF file, invisibly.
#'
#' @export
pact_manual <- function() {
  pdf <- system.file("doc", "pact_manual.pdf", package = "pactR")
  if (!nzchar(pdf)) stop("Manual not found. Is the pactR package installed?")
  if (.Platform$OS.type == "windows") {
    shell.exec(pdf)
  } else {
    system2("open", pdf)
  }
  invisible(pdf)
}
