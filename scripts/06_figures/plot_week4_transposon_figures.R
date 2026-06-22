#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
lifespan_path <- ifelse(length(args) >= 1, args[[1]], "data/processed/maintenance_lifespan_week4_full_sequence_validated.tsv")
ortholog_path <- ifelse(length(args) >= 2, args[[2]], "data/processed/ortholog_matrix_primary_week4_full_sequence_validated.tsv")
out_dir <- ifelse(length(args) >= 3, args[[3]], "results/figures")
manifest_path <- ifelse(length(args) >= 4, args[[4]], "results/figures/week4_transposon_figures_manifest.tsv")
report_path <- ifelse(length(args) >= 5, args[[5]], "results/reports/week4_transposon_figure_notes.md")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(manifest_path), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(report_path), recursive = TRUE, showWarnings = FALSE)

read_tsv <- function(path) {
  read.delim(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
}

write_png_pdf <- function(stem, width = 7.0, height = 5.2, plot_fun) {
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

theme_axes <- function() {
  box(bty = "l", lwd = 1.1)
}

clade_colors <- c(
  "Aves" = "#1b9e77",
  "Mammalia_Chiroptera" = "#d95f02",
  "Mammalia_nonChiroptera" = "#7570b3",
  "Reptilia" = "#666666"
)

clade_labels <- c(
  "Aves" = "Birds",
  "Mammalia_Chiroptera" = "Bats",
  "Mammalia_nonChiroptera" = "Other mammals",
  "Reptilia" = "Reptiles"
)

lifespan <- read_tsv(lifespan_path)
strict <- lifespan[lifespan$score_variant == "transposon_sequence_strict", ]
strict$transposon_suppression_score <- as.numeric(strict$transposon_suppression_score)
strict$lifespan_residual_log10 <- as.numeric(strict$lifespan_residual_log10)
strict$log10_body_mass_g <- as.numeric(strict$log10_body_mass_g)
strict <- strict[is.finite(strict$transposon_suppression_score) & is.finite(strict$lifespan_residual_log10), ]
strict$clade <- factor(strict$clade, levels = names(clade_colors))

ortholog <- read_tsv(ortholog_path)
trans <- ortholog[ortholog$maintenance_module == "transposon_suppression", ]

figures <- list()

figures[["week4_lifespan_residual_vs_strict_transposon_score"]] <- write_png_pdf(
  "week4_lifespan_residual_vs_strict_transposon_score",
  width = 7.2,
  height = 5.4,
  plot_fun = function() {
    par(mar = c(4.8, 5.0, 2.7, 1.2), las = 1, family = "sans")
    cols <- clade_colors[as.character(strict$clade)]
    pch <- ifelse(strict$flight_status %in% c("flighted", "powered_flight"), 19, 17)
    plot(
      strict$transposon_suppression_score,
      strict$lifespan_residual_log10,
      pch = pch,
      col = cols,
      bg = cols,
      cex = 1.2,
      xlab = "Strict sequence-supported transposon score",
      ylab = "Lifespan residual (log10 years)",
      main = "Lifespan residual vs strict transposon score",
      axes = FALSE
    )
    axis(1)
    axis(2)
    abline(h = 0, col = "grey65", lty = 3, lwd = 1)
    fit <- lm(lifespan_residual_log10 ~ transposon_suppression_score, data = strict)
    xs <- seq(min(strict$transposon_suppression_score), max(strict$transposon_suppression_score), length.out = 100)
    preds <- predict(fit, newdata = data.frame(transposon_suppression_score = xs), interval = "confidence")
    polygon(c(xs, rev(xs)), c(preds[, "lwr"], rev(preds[, "upr"])), col = adjustcolor("black", alpha.f = 0.12), border = NA)
    lines(xs, preds[, "fit"], lwd = 2.2, col = "black")
    theme_axes()
    legend(
      "topleft",
      legend = clade_labels[names(clade_colors)],
      col = clade_colors,
      pch = 19,
      pt.cex = 1.0,
      bty = "n",
      cex = 0.82
    )
    legend(
      "bottomright",
      legend = c("Powered flight", "Non-flying/flightless"),
      pch = c(19, 17),
      col = "grey25",
      bty = "n",
      cex = 0.82
    )
  }
)

figures[["week4_clade_colored_transposon_score_distribution"]] <- write_png_pdf(
  "week4_clade_colored_transposon_score_distribution",
  width = 7.2,
  height = 5.2,
  plot_fun = function() {
    par(mar = c(5.7, 4.8, 2.7, 1.2), las = 1, family = "sans")
    groups <- split(strict$transposon_suppression_score, strict$clade)
    groups <- groups[names(clade_colors)]
    boxplot(
      groups,
      col = adjustcolor(clade_colors[names(groups)], alpha.f = 0.25),
      border = clade_colors[names(groups)],
      ylab = "Strict sequence-supported transposon score",
      main = "Transposon score distribution by clade",
      xaxt = "n",
      outline = FALSE,
      ylim = c(0, 1.05)
    )
    axis(1, at = seq_along(groups), labels = clade_labels[names(groups)], las = 2, cex.axis = 0.85)
    set.seed(20260611)
    for (i in seq_along(groups)) {
      y <- groups[[i]]
      x <- jitter(rep(i, length(y)), amount = 0.12)
      points(x, y, pch = 19, col = adjustcolor(clade_colors[names(groups)[i]], alpha.f = 0.82), cex = 1.1)
      text(i, 1.03, labels = paste0("n=", length(y)), cex = 0.78, col = "grey25")
    }
    theme_axes()
  }
)

figures[["week4_sequence_validation_waterfall_counts"]] <- write_png_pdf(
  "week4_sequence_validation_waterfall_counts",
  width = 8.0,
  height = 5.4,
  plot_fun = function() {
    par(mar = c(5.8, 4.8, 2.7, 1.2), las = 1, family = "sans")
    direct_supported <- sum(trans$week4_validation_batch_source == "ncbi_direct_batch" & trans$week4_sequence_status == "sequence_supported", na.rm = TRUE)
    rescue_supported <- sum(trans$week4_validation_batch_source == "rescue_or_unresolved_batch" & trans$week4_sequence_status == "sequence_supported", na.rm = TRUE)
    weak_supported <- sum(trans$week4_sequence_status == "sequence_supported_weak", na.rm = TRUE)
    not_supported <- sum(trans$week4_sequence_status == "sequence_not_supported", na.rm = TRUE)
    counts <- c(
      "All transposon\ngene slots" = nrow(trans),
      "NCBI direct\nsequence-supported" = direct_supported,
      "Rescue/unresolved\nsequence-supported" = rescue_supported,
      "Weak\nsequence support" = weak_supported,
      "Not sequence\nsupported" = not_supported
    )
    bar_cols <- c("#4c4c4c", "#1b9e77", "#66a61e", "#e6ab02", "#bdbdbd")
    mids <- barplot(
      counts,
      col = bar_cols,
      border = NA,
      ylim = c(0, max(counts) * 1.16),
      ylab = "Transposon gene-species rows",
      main = "Week 4 sequence validation counts",
      cex.names = 0.78,
      las = 2
    )
    text(mids, counts + max(counts) * 0.035, labels = counts, cex = 0.86)
    axis(2)
    theme_axes()
    legend(
      "topright",
      legend = c("Initial slots", "Supported", "Weak", "Not supported"),
      fill = c("#4c4c4c", "#1b9e77", "#e6ab02", "#bdbdbd"),
      border = NA,
      bty = "n",
      cex = 0.82
    )
  }
)

manifest <- data.frame(
  figure_id = names(figures),
  png_path = vapply(figures, function(x) x[["png"]], character(1)),
  pdf_path = vapply(figures, function(x) x[["pdf"]], character(1)),
  stringsAsFactors = FALSE
)
write.table(manifest, manifest_path, sep = "\t", quote = FALSE, row.names = FALSE)

status_counts <- table(trans$week4_sequence_status)
report_lines <- c(
  "# Week 4 Transposon Figure Notes",
  "",
  "Generated figure drafts for the final Week 4 visualization pass.",
  "",
  "## Figures",
  "",
  "- `week4_lifespan_residual_vs_strict_transposon_score`: lifespan residual against the strict sequence-supported transposon score. Points are colored by clade and shaped by flight status. The fitted line is an ordinary visual trend, not the inferential PGLS model.",
  "- `week4_clade_colored_transposon_score_distribution`: strict transposon score distribution across birds, bats, other mammals, and reptiles.",
  "- `week4_sequence_validation_waterfall_counts`: count summary for all transposon gene-species rows through Week 4 sequence validation.",
  "",
  "## Counts",
  "",
  paste0("- strict-score species: ", nrow(strict)),
  paste0("- transposon gene-species rows: ", nrow(trans)),
  paste0("- sequence_supported: ", ifelse("sequence_supported" %in% names(status_counts), status_counts[["sequence_supported"]], 0)),
  paste0("- sequence_supported_weak: ", ifelse("sequence_supported_weak" %in% names(status_counts), status_counts[["sequence_supported_weak"]], 0)),
  paste0("- sequence_not_supported: ", ifelse("sequence_not_supported" %in% names(status_counts), status_counts[["sequence_not_supported"]], 0)),
  "",
  "## Files",
  ""
)
for (i in seq_len(nrow(manifest))) {
  report_lines <- c(
    report_lines,
    paste0("- `", manifest$figure_id[[i]], "`: `", manifest$png_path[[i]], "`, `", manifest$pdf_path[[i]], "`")
  )
}
writeLines(report_lines, report_path)
cat("Wrote", nrow(manifest), "figures,", manifest_path, "and", report_path, "\n")
