# ── PACT Core Example Analysis ────────────────────────────────────────────────
# Influenza H3N2 phylogenetic analysis using original PACT C++ engine.
# Run with:  Rscript analyse_example.R
# Or source from RStudio after setwd() to this folder.
# Outputs: TSV data files + PNG figures in the same directory as this script.

library(pactR)
library(ggplot2)

# ── Locate files --------------------------------------------------------------

this_dir <- local({
  cmd <- commandArgs(FALSE)
  f   <- sub("--file=", "", cmd[startsWith(cmd, "--file=")])
  if (length(f) == 1L) dirname(normalizePath(f)) else getwd()
})

trees_file <- file.path(this_dir, "in.trees")
param_file <- file.path(this_dir, "in.param")
stopifnot(file.exists(trees_file), file.exists(param_file))

cat("Trees file :", trees_file, "\n")
cat("Param file :", param_file, "\n")
cat("Output dir :", this_dir, "\n\n")

# ── Run PACT C++ engine -------------------------------------------------------

cat("=== Running PACT ===\n")
results <- run_pact(
  trees_file = trees_file,
  param_file = param_file,
  output_dir = this_dir,
  output_prefix = "out"
)
cat("Done.\n\n")

# ── Read outputs --------------------------------------------------------------

stats    <- read_pact_stats(results$stats)
skylines <- read_pact_skylines(results$skylines)
tips     <- read_pact_tips(results$tips)

cat("Stats rows    :", nrow(stats), "\n")
cat("Skylines rows :", nrow(skylines), "\n")
cat("Tips rows     :", nrow(tips), "\n\n")

# ── Write TSV files -----------------------------------------------------------

write.table(stats,    file.path(this_dir, "stats.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(skylines, file.path(this_dir, "skylines.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(tips,     file.path(this_dir, "tips.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("Saved: stats.tsv, skylines.tsv, tips.tsv\n\n")

# ── Okabe-Ito palette (7 states) ---------------------------------------------

OI7 <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
          "#0072B2", "#D55E00", "#CC79A7")

# ── 1. Summary statistics: migration rates heatmap ---------------------------
# Parse mig_X_Y rows into a from/to matrix

cat("=== 1. Migration rate heatmap ===\n")
mig_rows <- stats[grepl("^mig_[0-9]+_[0-9]+$", stats$statistic), ]
mig_rows$from  <- as.integer(sub("mig_(\\d+)_(\\d+)", "\\1", mig_rows$statistic))
mig_rows$to    <- as.integer(sub("mig_(\\d+)_(\\d+)", "\\2", mig_rows$statistic))

p_mig <- ggplot(mig_rows, aes(x = factor(to), y = factor(from), fill = mean)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = round(mean, 2)), size = 3, colour = "black") +
  scale_fill_gradient(low = "#f7fbff", high = "#084594",
                      name = "Rate") +
  labs(
    title = "Pairwise Migration Rates",
    x     = "To state",
    y     = "From state"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text  = element_text(size = 11)
  )
ggsave(file.path(this_dir, "fig_mig_rates.png"),
       p_mig, width = 6, height = 5.5, dpi = 150)
cat("Saved: fig_mig_rates.png\n\n")

# ── 2. Coalescent rates bar chart ---------------------------------------------

cat("=== 2. Coalescent rates ===\n")
coal_rows <- stats[grepl("^coal_[0-9]+$", stats$statistic), ]
coal_rows$state <- as.integer(sub("coal_(\\d+)", "\\1", coal_rows$statistic))

p_coal <- ggplot(coal_rows, aes(x = factor(state), y = mean, fill = factor(state))) +
  geom_col(width = 0.65, show.legend = FALSE) +
  scale_fill_manual(values = OI7) +
  labs(
    title = "Coalescent Rates by State",
    x     = "State",
    y     = "Coalescent rate"
  ) +
  theme_minimal(base_size = 12)
ggsave(file.path(this_dir, "fig_coal_rates.png"),
       p_coal, width = 5, height = 4, dpi = 150)
cat("Saved: fig_coal_rates.png\n\n")

# ── 3. Skyline: lineage proportions over time ---------------------------------

cat("=== 3. Skyline proportions ===\n")
pro <- skylines[grepl("^pro_[0-9]+$", skylines$statistic), ]
pro$state <- sub("pro_(\\d+)", "\\1", pro$statistic)
cat("Time range:", round(range(pro$time), 2), "| States:", length(unique(pro$state)), "\n")

p_sky_ribbon <- ggplot(pro,
    aes(x = time, y = mean, colour = state, fill = state)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.8) +
  scale_colour_manual(values = OI7, name = "State") +
  scale_fill_manual(values = OI7, name = "State") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  labs(
    title = "Lineage Proportions Over Time",
    x     = "Year",
    y     = "Proportion of lineages"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "right")
ggsave(file.path(this_dir, "fig_skyline_proportions.png"),
       p_sky_ribbon, width = 9, height = 4.5, dpi = 150)
cat("Saved: fig_skyline_proportions.png\n")

# Stacked area version (rescale to sum to 1 per time step)
pro_wide <- reshape(
  pro[, c("time", "state", "mean")],
  idvar     = "time",
  timevar   = "state",
  direction = "wide"
)
state_cols <- grep("^mean\\.", names(pro_wide), value = TRUE)
row_sums   <- rowSums(pro_wide[, state_cols], na.rm = TRUE)
pro_wide[, state_cols] <- pro_wide[, state_cols] / ifelse(row_sums == 0, 1, row_sums)
pro_norm <- reshape(pro_wide, direction = "long",
                    varying = state_cols, v.names = "mean",
                    timevar = "state", times = sub("mean\\.", "", state_cols))
pro_norm <- pro_norm[order(pro_norm$time, pro_norm$state), ]

p_sky_area <- ggplot(pro_norm, aes(x = time, y = mean, fill = state)) +
  geom_area(position = "stack") +
  scale_fill_manual(values = OI7, name = "State") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = 0)) +
  labs(
    title = "Lineage Proportions Over Time",
    x     = "Year",
    y     = "Proportion of lineages"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "right")
ggsave(file.path(this_dir, "fig_skyline_area.png"),
       p_sky_area, width = 9, height = 4.5, dpi = 150)
cat("Saved: fig_skyline_area.png\n\n")

# ── 4. Tips: time-to-trunk scatter -------------------------------------------

cat("=== 4. Time-to-trunk ===\n")
ttt <- tips[tips$statistic == "time_to_trunk", ]
ttt$state <- as.character(ttt$label)
cat("Tips:", nrow(ttt), "| States:", length(unique(ttt$state)), "\n")

p_ttt <- ggplot(ttt, aes(x = time, y = mean, colour = state)) +
  geom_point(size = 1.8, alpha = 0.7) +
  scale_colour_manual(values = OI7, name = "State") +
  labs(
    title = "Time to Trunk",
    x     = "Sampling time (year)",
    y     = "Time to trunk (years)"
  ) +
  theme_minimal(base_size = 12)
ggsave(file.path(this_dir, "fig_time_to_trunk.png"),
       p_ttt, width = 8, height = 4.5, dpi = 150)
cat("Saved: fig_time_to_trunk.png\n\n")

# ── Summary -------------------------------------------------------------------

cat("Done. Files written to:\n")
outputs <- c(
  "stats.tsv", "skylines.tsv", "tips.tsv",
  "fig_mig_rates.png", "fig_coal_rates.png",
  "fig_skyline_proportions.png", "fig_skyline_area.png",
  "fig_time_to_trunk.png"
)
for (f in outputs) cat(" ", file.path(this_dir, f), "\n")
