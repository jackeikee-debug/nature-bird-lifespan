#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
pgls_path <- ifelse(length(args) >= 1, args[[1]], "results/tables/phase2_W3_full_background_expanded_score_variant_pgls.tsv")
module_ranking_path <- ifelse(length(args) >= 2, args[[2]], "results/tables/phase2_final_module_ranking.tsv")
annotation_path <- ifelse(length(args) >= 3, args[[3]], "results/tables/phase2_final_annotation_bias_module_ranking.tsv")
sensitivity_path <- ifelse(length(args) >= 4, args[[4]], "results/tables/phase2_final_clade_sensitivity_pgls.tsv")
gates_path <- ifelse(length(args) >= 5, args[[5]], "results/tables/phase2_final_decision_gates.tsv")
out_dir <- ifelse(length(args) >= 6, args[[6]], "results/figures")
manifest_path <- ifelse(length(args) >= 7, args[[7]], "results/figures/phase2_final_figures_manifest.tsv")
report_path <- ifelse(length(args) >= 8, args[[8]], "results/reports/phase2_final_figure_package_report.md")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(manifest_path), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(report_path), recursive = TRUE, showWarnings = FALSE)

read_tsv <- function(path) {
  read.delim(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
}

num <- function(x) suppressWarnings(as.numeric(x))

fmt_p <- function(x) {
  x <- num(x)
  ifelse(is.na(x), "NA", ifelse(x < 0.001, formatC(x, format = "e", digits = 2), sprintf("%.3f", x)))
}

module_labels <- c(
  "transposon_repeat_suppression" = "Transposon/repeat",
  "chromatin_repression_heterochromatin" = "Chromatin repression",
  "DNA_repair_replication_stress" = "DNA repair",
  "proteostasis_autophagy_mitophagy" = "Proteostasis/autophagy",
  "cancer_surveillance_senescence" = "Cancer/senescence",
  "inflammation_innate_immune_restraint" = "Inflammation restraint"
)

status_colors <- c(
  "pass" = "#1b9e77",
  "caution" = "#d95f02",
  "caution-go" = "#7570b3",
  "fail" = "#bdbdbd"
)

module_colors <- c(
  "transposon_repeat_suppression" = "#1b9e77",
  "chromatin_repression_heterochromatin" = "#7570b3",
  "DNA_repair_replication_stress" = "#386cb0",
  "proteostasis_autophagy_mitophagy" = "#4daf4a",
  "cancer_surveillance_senescence" = "#984ea3",
  "inflammation_innate_immune_restraint" = "#e41a1c"
)

write_png_pdf <- function(stem, width = 8.0, height = 5.6, plot_fun) {
  png_path <- file.path(out_dir, paste0(stem, ".png"))
  pdf_path <- file.path(out_dir, paste0(stem, ".pdf"))
  png(png_path, width = width, height = height, units = "in", res = 300, type = "cairo")
  plot_fun()
  dev.off()
  pdf(pdf_path, width = width, height = height, useDingbats = FALSE)
  plot_fun()
  dev.off()
  c(png = png_path, pdf = pdf_path)
}

module_ranking <- read_tsv(module_ranking_path)
annotation <- read_tsv(annotation_path)
sensitivity <- read_tsv(sensitivity_path)
gates <- read_tsv(gates_path)
pgls <- read_tsv(pgls_path)

for (col in c("p", "p_bh_by_model", "rank_by_p", "mean_coverage")) {
  if (col %in% names(module_ranking)) module_ranking[[col]] <- num(module_ranking[[col]])
}
for (col in c("p", "p_bh_by_model", "rank_by_p", "estimate")) {
  if (col %in% names(annotation)) annotation[[col]] <- num(annotation[[col]])
}
for (col in c("estimate", "se", "p", "p_bh_by_subset_model", "rank_by_p", "n")) {
  if (col %in% names(sensitivity)) sensitivity[[col]] <- num(sensitivity[[col]])
}
for (col in c("module_estimate", "module_se", "module_p", "module_p_bh_by_variant_model")) {
  if (col %in% names(pgls)) pgls[[col]] <- num(pgls[[col]])
}

figures <- list()

figures[["phase2_final_module_ranking"]] <- write_png_pdf(
  "phase2_final_module_ranking",
  width = 8.6,
  height = 5.8,
  plot_fun = function() {
    residual <- module_ranking[module_ranking$model == "pgls_clade_residual_module", ]
    mass <- module_ranking[module_ranking$model == "mass_clade_module", ]
    residual <- residual[order(residual$rank_by_p), ]
    modules <- residual$maintenance_module
    y <- rev(seq_along(modules))
    residual_vals <- -log10(residual$p)
    mass_vals <- -log10(mass$p[match(modules, mass$maintenance_module)])
    xlim <- range(c(residual_vals, mass_vals, -log10(0.05)), na.rm = TRUE)
    xlim <- c(0, xlim[2] * 1.18)
    par(mar = c(4.8, 8.2, 3.0, 1.2), las = 1, family = "sans")
    plot(
      residual_vals, y + 0.12,
      xlim = xlim,
      ylim = c(0.5, length(modules) + 0.5),
      yaxt = "n",
      xlab = "-log10(P value)",
      ylab = "",
      main = "Maintenance-module prioritization",
      pch = 19,
      col = "#1b9e77",
      cex = 1.3,
      axes = FALSE
    )
    axis(1)
    axis(2, at = y, labels = module_labels[modules], tick = FALSE, cex.axis = 0.78)
    abline(v = -log10(0.05), col = "grey55", lty = 3, lwd = 1.1)
    points(mass_vals, y - 0.12, pch = 17, col = "#984ea3", cex = 1.25)
    segments(0, y + 0.12, residual_vals, y + 0.12, col = adjustcolor("#1b9e77", alpha.f = 0.55), lwd = 2)
    segments(0, y - 0.12, mass_vals, y - 0.12, col = adjustcolor("#984ea3", alpha.f = 0.55), lwd = 2)
    text(
      pmax(residual_vals, mass_vals) + 0.04,
      y,
      labels = paste0("cov=", sprintf("%.2f", residual$mean_coverage)),
      cex = 0.68,
      col = "grey25",
      adj = 0
    )
    box(bty = "l", lwd = 1.1)
    legend(
      "bottomright",
      legend = c("PGLS clade residual", "Mass + clade"),
      pch = c(19, 17),
      col = c("#1b9e77", "#984ea3"),
      bty = "n",
      cex = 0.82
    )
    text(-log10(0.05) + 0.03, 0.58, labels = "P = 0.05", cex = 0.66, col = "grey35", adj = 0)
  }
)

figures[["phase2_final_annotation_bias_stress_test"]] <- write_png_pdf(
  "phase2_final_annotation_bias_stress_test",
  width = 9.4,
  height = 6.2,
  plot_fun = function() {
    models <- c(
      "residual_base",
      "residual_tier",
      "residual_coverage",
      "residual_tier_coverage",
      "mass_clade_tier",
      "mass_clade_tier_coverage"
    )
    model_labels <- c(
      "Residual",
      "Residual + tier",
      "Residual + coverage",
      "Residual + tier + coverage",
      "Mass/clade + tier",
      "Mass/clade + tier + coverage"
    )
    residual <- annotation[annotation$model == "residual_base", ]
    residual <- residual[order(residual$rank_by_p), ]
    modules <- residual$maintenance_module
    z <- matrix(NA_real_, nrow = length(modules), ncol = length(models), dimnames = list(modules, models))
    ranks <- matrix(NA_integer_, nrow = length(modules), ncol = length(models), dimnames = list(modules, models))
    for (i in seq_along(models)) {
      sub <- annotation[annotation$model == models[i], ]
      z[, i] <- -log10(sub$p[match(modules, sub$maintenance_module)])
      ranks[, i] <- sub$rank_by_p[match(modules, sub$maintenance_module)]
    }
    z_plot <- pmin(z, 3)
    palette <- colorRampPalette(c("#f7f7f7", "#fee08b", "#d95f02", "#7570b3"))(80)
    par(mar = c(8.8, 8.2, 3.0, 5.0), las = 1, family = "sans")
    plot(
      NA,
      xlim = c(0.5, length(models) + 0.5),
      ylim = c(0.5, length(modules) + 0.5),
      xaxt = "n",
      yaxt = "n",
      xlab = "",
      ylab = "",
      main = "Annotation-bias stress test",
      axes = FALSE
    )
    for (i in seq_along(modules)) {
      for (j in seq_along(models)) {
        value <- z_plot[i, j]
        col_idx <- ifelse(is.na(value), 1, max(1, min(80, round(value / 3 * 79) + 1)))
        rect(j - 0.5, length(modules) - i + 0.5, j + 0.5, length(modules) - i + 1.5, col = palette[col_idx], border = "white")
        text(j, length(modules) - i + 1, labels = paste0("r", ranks[i, j]), cex = 0.7, col = "grey15")
      }
    }
    axis(1, at = seq_along(models), labels = model_labels, las = 2, cex.axis = 0.74, tick = FALSE)
    axis(2, at = rev(seq_along(modules)), labels = module_labels[modules], tick = FALSE, cex.axis = 0.76)
    box(lwd = 1.1)
    legend_y <- seq(1, length(palette), length.out = 6)
    legend_vals <- seq(0, 3, length.out = 6)
    x0 <- length(models) + 0.78
    for (k in seq_along(palette)) {
      rect(x0, 0.65 + (k - 1) * 0.045, x0 + 0.18, 0.65 + k * 0.045, col = palette[k], border = NA, xpd = NA)
    }
    text(x0 + 0.34, 0.65 + (legend_y - 1) * 0.045, labels = sprintf("%.1f", legend_vals), cex = 0.66, adj = 0, xpd = NA)
    text(x0, 0.4, labels = "-log10(P)", cex = 0.72, adj = 0, xpd = NA)
    mtext("Cell text shows within-model rank; coverage-adjusted models deliberately stress observability bias.", side = 1, line = 7.9, cex = 0.72, col = "grey35")
  }
)

figures[["phase2_final_transposon_sensitivity_forest"]] <- write_png_pdf(
  "phase2_final_transposon_sensitivity_forest",
  width = 10.2,
  height = 5.9,
  plot_fun = function() {
    sub <- sensitivity[
      sensitivity$maintenance_module == "transposon_repeat_suppression" &
        sensitivity$model == "pgls_clade_residual_module" &
        (is.na(sensitivity$error) | sensitivity$error == ""),
    ]
    subset_order <- c(
      "all_primary",
      "birds_only",
      "no_bats",
      "no_reptiles",
      "no_nonflying_mammals",
      "no_birds",
      "exclude_top_abs_residual_5",
      "exclude_top_abs_residual_10"
    )
    subset_labels <- c(
      "all_primary" = "All primary species",
      "birds_only" = "Birds only",
      "no_bats" = "Remove bats",
      "no_reptiles" = "Remove reptiles",
      "no_nonflying_mammals" = "Remove nonflight mammals",
      "no_birds" = "Remove birds",
      "exclude_top_abs_residual_5" = "Remove top-5 outliers",
      "exclude_top_abs_residual_10" = "Remove top-10 outliers"
    )
    sub <- sub[match(subset_order, sub$subset), ]
    sub <- sub[!is.na(sub$subset), ]
    sub$ci_low <- sub$estimate - 1.96 * sub$se
    sub$ci_high <- sub$estimate + 1.96 * sub$se
    status <- ifelse(sub$p_bh_by_subset_model < 0.05, "BH < 0.05", ifelse(sub$p < 0.05, "P < 0.05 only", "Not significant"))
    cols <- c("BH < 0.05" = "#1b9e77", "P < 0.05 only" = "#d95f02", "Not significant" = "#bdbdbd")[status]
    y <- rev(seq_len(nrow(sub)))
    xlim <- range(c(sub$ci_low, sub$ci_high, 0), na.rm = TRUE)
    xlim <- c(xlim[1] - diff(xlim) * 0.08, xlim[2] + diff(xlim) * 0.25)
    par(mar = c(4.8, 9.6, 3.0, 1.3), las = 1, family = "sans")
    plot(
      sub$estimate,
      y,
      xlim = xlim,
      ylim = c(0.5, nrow(sub) + 0.5),
      yaxt = "n",
      xlab = "PGLS coefficient for transposon/repeat score",
      ylab = "",
      main = "Final transposon clade/outlier sensitivity",
      pch = 19,
      col = cols,
      cex = 1.25,
      axes = FALSE
    )
    axis(1)
    axis(2, at = y, labels = subset_labels[sub$subset], tick = FALSE, cex.axis = 0.72)
    abline(v = 0, col = "grey55", lty = 3, lwd = 1.1)
    for (i in seq_len(nrow(sub))) {
      segments(sub$ci_low[i], y[i], sub$ci_high[i], y[i], col = cols[i], lwd = 2)
      points(sub$estimate[i], y[i], pch = 19, col = cols[i], cex = 1.3)
      text(xlim[2], y[i], labels = paste0("p=", fmt_p(sub$p[i]), ", n=", sub$n[i]), adj = 1, cex = 0.66, col = "grey25")
    }
    box(bty = "l", lwd = 1.1)
    legend(
      "topleft",
      legend = c("BH < 0.05", "P < 0.05 only", "Not significant"),
      pch = 19,
      col = c("#1b9e77", "#d95f02", "#bdbdbd"),
      bty = "n",
      cex = 0.78
    )
  }
)

figures[["phase2_final_decision_gates"]] <- write_png_pdf(
  "phase2_final_decision_gates",
  width = 9.2,
  height = 5.8,
  plot_fun = function() {
    gates$gate_label <- gsub("_", " ", gates$gate)
    gates$gate_label <- sub("^annotation bias", "annotation-bias", gates$gate_label)
    y <- rev(seq_len(nrow(gates)))
    cols <- status_colors[gates$status]
    cols[is.na(cols)] <- "#bdbdbd"
    par(mar = c(4.2, 9.0, 3.0, 1.3), las = 1, family = "sans")
    plot(
      NA,
      xlim = c(0, 1),
      ylim = c(0.5, nrow(gates) + 0.5),
      xaxt = "n",
      yaxt = "n",
      xlab = "",
      ylab = "",
      main = "Phase 2 decision gates",
      axes = FALSE
    )
    for (i in seq_len(nrow(gates))) {
      rect(0.02, y[i] - 0.32, 0.68, y[i] + 0.32, col = adjustcolor(cols[i], alpha.f = 0.18), border = cols[i], lwd = 1.2)
      text(0.04, y[i], labels = gates$gate_label[i], adj = 0, cex = 0.72, col = "grey15")
      text(0.97, y[i], labels = gates$status[i], adj = 1, cex = 0.78, col = cols[i], font = 2)
    }
    axis(2, at = y, labels = rep("", length(y)), tick = FALSE)
    box(bty = "l", lwd = 1.1)
    legend(
      "bottomleft",
      legend = names(status_colors),
      fill = status_colors,
      border = NA,
      horiz = TRUE,
      bty = "n",
      cex = 0.78
    )
    mtext("Interpretation: caution-go, with coverage/observability and outlier sensitivity as explicit limitations.", side = 1, line = 2.6, cex = 0.72, col = "grey35")
  }
)

manifest <- data.frame(
  figure_id = names(figures),
  png_path = vapply(figures, function(x) x[["png"]], character(1)),
  pdf_path = vapply(figures, function(x) x[["pdf"]], character(1)),
  role = c(
    "main_module_ranking",
    "annotation_bias_stress_test",
    "clade_outlier_sensitivity",
    "decision_gate_summary"
  ),
  stringsAsFactors = FALSE
)
write.table(manifest, manifest_path, sep = "\t", quote = FALSE, row.names = FALSE)

trans_ann <- annotation[annotation$maintenance_module == "transposon_repeat_suppression", ]
trans_sens <- sensitivity[
  sensitivity$maintenance_module == "transposon_repeat_suppression" &
    sensitivity$model == "pgls_clade_residual_module",
]

report <- c(
  "# Phase 2 Final Figure Package",
  "",
  "Generated final Phase 2 figure drafts for the manuscript/preprint package.",
  "",
  "## Figures",
  "",
  "- `phase2_final_module_ranking`: final W3 module ranking by PGLS clade-residual and mass+clade P values, with mean module coverage shown next to each module.",
  "- `phase2_final_annotation_bias_stress_test`: all-module annotation-bias stress test showing that tier adjustment preserves the transposon/repeat ranking while coverage adjustment collapses module-score terms.",
  "- `phase2_final_transposon_sensitivity_forest`: transposon/repeat coefficient under clade-removal and residual-outlier sensitivity.",
  "- `phase2_final_decision_gates`: compact status summary for the final decision gates.",
  "",
  "## Key Numbers",
  "",
  paste0("- transposon residual+tier p: ", fmt_p(trans_ann$p[trans_ann$model == "residual_tier"])),
  paste0("- transposon residual+coverage p: ", fmt_p(trans_ann$p[trans_ann$model == "residual_coverage"])),
  paste0("- transposon no-birds sensitivity p: ", fmt_p(trans_sens$p[trans_sens$subset == "no_birds"])),
  paste0("- transposon top-5 outlier exclusion p: ", fmt_p(trans_sens$p[trans_sens$subset == "exclude_top_abs_residual_5"])),
  "",
  "## Interpretation",
  "",
  "The figure package supports a cautious module-ranking story. Transposon/repeat suppression remains prioritized after genome-tier adjustment, but not after direct module-coverage adjustment. The final signal is bird-dependent or avian-enriched and sensitive to top lifespan-residual outliers.",
  "",
  "## Files",
  ""
)
for (i in seq_len(nrow(manifest))) {
  report <- c(report, paste0("- `", manifest$figure_id[i], "`: `", manifest$png_path[i], "`, `", manifest$pdf_path[i], "`"))
}
writeLines(report, report_path)

cat("Wrote", nrow(manifest), "Phase 2 final figures,", manifest_path, "and", report_path, "\n")
