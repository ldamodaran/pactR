# ── Global HPAI H5N1 DTA Posterior Analysis ───────────────────────────────────
# Run with:  Rscript analyse_global_HPAI.R
# Or source from RStudio after setwd() to this folder.
# All figures and TSV files are saved to the same directory as this script.

library(pactR)
library(ggplot2)

# ── Locate files --------------------------------------------------------------

this_dir <- local({
  cmd <- commandArgs(FALSE)
  f   <- sub("--file=", "", cmd[startsWith(cmd, "--file=")])
  if (length(f) == 1L) dirname(normalizePath(f)) else getwd()
})

trees_file <- file.path(this_dir, "global_HPAI.trees")
stopifnot(file.exists(trees_file))

cat("Trees file :", trees_file, "\n")
cat("Output dir :", this_dir, "\n\n")

BURNIN <- 3L     # discard first 3 trees (~10% of 31)
TRAIT  <- "geo"

# ── 1. Persistence ────────────────────────────────────────────────────────────

cat("=== 1. Persistence ===\n")
pers <- persistence_dta(trees_file, trait = TRAIT, burnin = BURNIN)
print(pers[order(pers$mean, decreasing = TRUE), ], row.names = FALSE)
write.table(pers, file.path(this_dir, "persistence.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("Saved: persistence.tsv\n")

p_pers <- plot_persistence(
  pers,
  title = "Global HPAI H5N1 Lineage Persistence by Region"
)
ggsave(file.path(this_dir, "fig_persistence.png"),
       p_pers, width = 6, height = 3.5, dpi = 150)
cat("Saved: fig_persistence.png\n\n")

# ── 2. Proportions over time ──────────────────────────────────────────────────

cat("=== 2. Proportions over time ===\n")
prop <- dta_proportions(trees_file, trait = TRAIT, step = 0.05, burnin = BURNIN)
cat("Time range:", round(range(prop$time), 2),
    "| Regions:", length(unique(prop$state)), "\n")
write.table(prop, file.path(this_dir, "proportions.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("Saved: proportions.tsv\n")

p_ribbon <- plot_dta_proportions(
  prop,
  type  = "ribbon",
  title = "Global HPAI H5N1 Regional Proportions Over Time"
)
ggsave(file.path(this_dir, "fig_proportions_ribbon.png"),
       p_ribbon, width = 10, height = 5, dpi = 150)
cat("Saved: fig_proportions_ribbon.png\n")

p_area <- plot_dta_proportions(
  prop,
  type  = "area",
  title = "Global HPAI H5N1 Regional Proportions Over Time"
)
ggsave(file.path(this_dir, "fig_proportions_area.png"),
       p_area, width = 10, height = 5, dpi = 150)
cat("Saved: fig_proportions_area.png\n\n")

# ── 3. Skyline (lineage through time) ─────────────────────────────────────────

cat("=== 3. Skyline (lineage through time) ===\n")
ltt <- lineage_through_time(trees_file, step = 0.05, burnin = BURNIN)
cat("Peak lineages:", round(max(ltt$mean), 0),
    "at time", round(ltt$time[which.max(ltt$mean)], 3), "\n")
write.table(ltt, file.path(this_dir, "ltt.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("Saved: ltt.tsv\n")

p_ltt <- plot_ltt(
  ltt,
  title = "Global HPAI H5N1 Lineage Through Time"
)
ggsave(file.path(this_dir, "fig_skyline_ltt.png"),
       p_ltt, width = 9, height = 4.5, dpi = 150)
cat("Saved: fig_skyline_ltt.png\n\n")

# ── Summary ───────────────────────────────────────────────────────────────────

cat("Done. Files written to:\n")
outputs <- c(
  "persistence.tsv", "proportions.tsv", "ltt.tsv",
  "fig_persistence.png", "fig_proportions_ribbon.png",
  "fig_proportions_area.png", "fig_skyline_ltt.png"
)
for (f in outputs) cat(" ", file.path(this_dir, f), "\n")
