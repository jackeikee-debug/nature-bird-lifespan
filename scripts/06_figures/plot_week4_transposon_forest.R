#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
subset_path <- ifelse(length(args) >= 1, args[[1]], "results/tables/week4_full_sequence_transposon_sensitivity.tsv")
leave_one_path <- ifelse(length(args) >= 2, args[[2]], "results/tables/week4_transposon_gene_leave_one_out_pgls.tsv")
out_dir <- ifelse(length(args) >= 3, args[[3]], "results/figures")
report_path <- ifelse(length(args) >= 4, args[[4]], "results/reports/week4_transposon_forest_notes.md")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(report_path), recursive = TRUE, showWarnings = FALSE)

read_tsv <- function(path) {
  read.delim(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
}

subset <- read_tsv(subset_path)
leave_one <- read_tsv(leave_one_path)

subset_key <- subset[
  subset$score_variant == "transposon_sequence_strict" &
    subset$model == "mass_clade_module" &
    (is.na(subset$error) | subset$error == ""),
]
subset_key$estimate <- as.numeric(subset_key$estimate)
subset_key$se <- as.numeric(subset_key$se)
subset_key$p <- as.numeric(subset_key$p)
subset_key$p_bh_by_variant_model <- as.numeric(subset_key$p_bh_by_variant_model)
subset_key$n <- as.numeric(subset_key$n)

subset_labels <- c(
  all_primary = "Main effect",
  exclude_human = "Exclude human",
  exclude_top_abs_residual_5 = "Exclude top 5 outliers",
  leave_out_Mammalia_Chiroptera = "Exclude bats",
  leave_out_Mammalia_nonChiroptera = "Exclude other mammals",
  leave_out_Reptilia = "Exclude reptiles",
  tier1_only = "Tier1 only",
  leave_out_Aves = "Exclude birds"
)
subset_order <- names(subset_labels)
subset_rows <- subset_key[match(subset_order, subset_key$subset), ]
subset_rows <- subset_rows[!is.na(subset_rows$subset), ]
subset_forest <- data.frame(
  group = "Subset sensitivity",
  label = unname(subset_labels[subset_rows$subset]),
  estimate = subset_rows$estimate,
  se = subset_rows$se,
  p = subset_rows$p,
  bh = subset_rows$p_bh_by_variant_model,
  n = subset_rows$n,
  support = ifelse(subset_rows$estimate > 0 & subset_rows$p_bh_by_variant_model < 0.05, "BH significant", "Not supported"),
  stringsAsFactors = FALSE
)

loo_key <- leave_one[
  grepl("_strict$", leave_one$score_variant) &
    leave_one$model == "mass_clade_module" &
    (is.na(leave_one$error) | leave_one$error == ""),
]
loo_key$module_estimate <- as.numeric(loo_key$module_estimate)
loo_key$module_se <- as.numeric(loo_key$module_se)
loo_key$module_p <- as.numeric(loo_key$module_p)
loo_key$module_p_bh_by_variant_model <- as.numeric(loo_key$module_p_bh_by_variant_model)
loo_key$n <- as.numeric(loo_key$n)

loo_order <- c(
  "drop_MOV10_strict",
  "drop_PIWIL1_strict",
  "drop_PIWIL2_strict",
  "drop_SETDB1_strict",
  "drop_TRIM28_strict"
)
loo_labels <- c(
  drop_MOV10_strict = "Drop MOV10",
  drop_PIWIL1_strict = "Drop PIWIL1",
  drop_PIWIL2_strict = "Drop PIWIL2",
  drop_SETDB1_strict = "Drop SETDB1",
  drop_TRIM28_strict = "Drop TRIM28"
)
loo_rows <- loo_key[match(loo_order, loo_key$score_variant), ]
loo_rows <- loo_rows[!is.na(loo_rows$score_variant), ]
loo_forest <- data.frame(
  group = "Gene leave-one-out",
  label = unname(loo_labels[loo_rows$score_variant]),
  estimate = loo_rows$module_estimate,
  se = loo_rows$module_se,
  p = loo_rows$module_p,
  bh = loo_rows$module_p_bh_by_variant_model,
  n = loo_rows$n,
  support = ifelse(loo_rows$module_estimate > 0 & loo_rows$module_p_bh_by_variant_model < 0.05, "BH significant", "Not supported"),
  stringsAsFactors = FALSE
)

forest <- rbind(subset_forest, loo_forest)
forest$ci_low <- forest$estimate - 1.96 * forest$se
forest$ci_high <- forest$estimate + 1.96 * forest$se
forest$label_text <- forest$label

plot_forest <- function() {
  par(mar = c(4.8, 7.5, 3.0, 1.5), las = 1, family = "sans")
  y <- rev(seq_len(nrow(forest)))
  xlim <- range(c(forest$ci_low, forest$ci_high, 0), na.rm = TRUE)
  xpad <- diff(xlim) * 0.08
  xlim <- c(xlim[1] - xpad, xlim[2] + xpad)
  cols <- ifelse(forest$support == "BH significant", "#1b9e77", "#bdbdbd")
  pch <- ifelse(forest$group == "Subset sensitivity", 19, 17)
  plot(
    forest$estimate,
    y,
    xlim = xlim,
    ylim = c(0.5, nrow(forest) + 0.5),
    yaxt = "n",
    ylab = "",
    xlab = "PGLS coefficient for strict transposon score",
    main = "Main effect and sensitivity forest plot",
    pch = pch,
    col = cols,
    cex = 1.2,
    axes = FALSE
  )
  axis(1)
  axis(2, at = y, labels = forest$label_text, tick = FALSE, las = 1, cex.axis = 0.72)
  abline(v = 0, lty = 3, col = "grey55", lwd = 1.1)
  for (i in seq_along(y)) {
    segments(forest$ci_low[i], y[i], forest$ci_high[i], y[i], col = cols[i], lwd = 2)
    points(forest$estimate[i], y[i], pch = pch[i], col = cols[i], bg = cols[i], cex = 1.25)
  }
  group_change <- max(which(forest$group == "Subset sensitivity"))
  abline(h = nrow(forest) - group_change + 0.5, col = "grey85", lwd = 1)
  box(bty = "l", lwd = 1.1)
  legend(
    "bottomright",
    legend = c("Subset sensitivity", "Gene leave-one-out", "BH significant", "Not supported"),
    pch = c(19, 17, 15, 15),
    col = c("grey25", "grey25", "#1b9e77", "#bdbdbd"),
    bty = "n",
    cex = 0.78
  )
}

png_path <- file.path(out_dir, "week4_transposon_main_effect_sensitivity_forest.png")
pdf_path <- file.path(out_dir, "week4_transposon_main_effect_sensitivity_forest.pdf")
png(png_path, width = 9.2, height = 6.8, units = "in", res = 300, type = "cairo")
plot_forest()
dev.off()
pdf(pdf_path, width = 9.2, height = 6.8, useDingbats = FALSE)
plot_forest()
dev.off()

report <- c(
  "# Week 4 Transposon Forest Plot Notes",
  "",
  "Generated a forest plot combining the strict transposon main effect, subset sensitivity checks, and single-gene leave-one-out checks.",
  "",
  "Intervals are approximate estimate +/- 1.96 SE from the corresponding mass+clade PGLS model rows.",
  "",
  "## Files",
  "",
  paste0("- `", png_path, "`"),
  paste0("- `", pdf_path, "`"),
  "",
  "## Rows",
  ""
)
for (i in seq_len(nrow(forest))) {
  report <- c(
    report,
    sprintf(
      "- %s: estimate=%.3f, 95%% CI [%.3f, %.3f], p=%.4g, BH=%.4g, n=%d, %s",
      forest$label[i], forest$estimate[i], forest$ci_low[i], forest$ci_high[i],
      forest$p[i], forest$bh[i], forest$n[i], forest$support[i]
    )
  )
}
writeLines(report, report_path)
cat("Wrote", png_path, pdf_path, "and", report_path, "\n")
