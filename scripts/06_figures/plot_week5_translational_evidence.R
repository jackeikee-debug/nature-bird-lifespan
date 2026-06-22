#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
human_mapping_path <- ifelse(length(args) >= 1, args[[1]], "data/processed/human_mapping.tsv")
translation_path <- ifelse(length(args) >= 2, args[[2]], "data/processed/human_translation_priority.tsv")
out_dir <- ifelse(length(args) >= 3, args[[3]], "results/figures")
plot_table_path <- ifelse(length(args) >= 4, args[[4]], "data/processed/week5_translational_evidence_plot_table.tsv")
report_path <- ifelse(length(args) >= 5, args[[5]], "results/reports/week5_translational_evidence_figure_report.md")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(plot_table_path), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(report_path), recursive = TRUE, showWarnings = FALSE)

read_tsv <- function(path) {
  read.delim(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
}

yes01 <- function(x) {
  ifelse(tolower(as.character(x)) == "yes", 1, 0)
}

keyword01 <- function(x) {
  ifelse(as.character(x) == "keyword_supported", 1, 0)
}

priority_rank <- c(
  "high_translation_priority" = 1,
  "medium_translation_priority" = 2,
  "supporting_translation_context" = 3,
  "low_translation_priority" = 4
)

module_order <- c(
  "transposon_suppression",
  "DNA_repair",
  "cancer_surveillance",
  "inflammation_control",
  "mitochondrial_quality_control",
  "proteostasis",
  "autophagy"
)

module_colors <- c(
  "transposon_suppression" = "#1b9e77",
  "DNA_repair" = "#386cb0",
  "cancer_surveillance" = "#984ea3",
  "inflammation_control" = "#e41a1c",
  "mitochondrial_quality_control" = "#ff7f00",
  "proteostasis" = "#4daf4a",
  "autophagy" = "#a65628"
)

module_legend_labels <- c(
  "transposon_suppression" = "Transposon",
  "DNA_repair" = "DNA repair",
  "cancer_surveillance" = "Cancer surveillance",
  "inflammation_control" = "Inflammation",
  "mitochondrial_quality_control" = "Mito QC",
  "proteostasis" = "Proteostasis",
  "autophagy" = "Autophagy"
)

human <- read_tsv(human_mapping_path)
translation <- read_tsv(translation_path)

merged <- merge(
  human,
  translation[, c(
    "human_gene_symbol",
    "max_open_targets_score",
    "small_molecule_tractability",
    "small_molecule_tractability_count",
    "cancer_top_disease_count",
    "nervous_system_top_disease_count",
    "immune_inflammation_top_disease_count",
    "week5_translation_priority"
  )],
  by = "human_gene_symbol",
  all.x = TRUE
)

merged$module_rank <- match(merged$maintenance_module, module_order)
merged$module_rank[is.na(merged$module_rank)] <- 99
merged$priority_rank <- priority_rank[merged$week5_translation_priority]
merged$priority_rank[is.na(merged$priority_rank)] <- 99
merged$max_open_targets_score <- as.numeric(merged$max_open_targets_score)
merged$max_open_targets_score[is.na(merged$max_open_targets_score)] <- 0
merged$small_molecule_tractability_count <- as.numeric(merged$small_molecule_tractability_count)
merged$small_molecule_tractability_count[is.na(merged$small_molecule_tractability_count)] <- 0

merged$week4_focal_transposon_gene <- ifelse(merged$maintenance_module == "transposon_suppression", 1, 0)
merged$genage <- yes01(merged$genage_human_evidence)
merged$longevitymap <- yes01(merged$longevitymap_evidence)
merged$cellage <- yes01(merged$cellage_evidence)
merged$genome_instability <- keyword01(merged$genome_instability_evidence)
merged$repeat_control <- keyword01(merged$transposon_or_repeat_evidence)
merged$senescence <- keyword01(merged$senescence_evidence)
merged$cancer_disease <- ifelse(as.numeric(merged$cancer_top_disease_count) > 0, 1, 0)
merged$neuro_disease <- ifelse(as.numeric(merged$nervous_system_top_disease_count) > 0, 1, 0)
merged$immune_disease <- ifelse(as.numeric(merged$immune_inflammation_top_disease_count) > 0, 1, 0)
merged$small_molecule <- ifelse(merged$small_molecule_tractability_count > 0, 1, 0)

merged <- merged[order(merged$module_rank, merged$priority_rank, -merged$max_open_targets_score, merged$human_gene_symbol), ]

plot_table <- merged[, c(
  "human_gene_symbol",
  "maintenance_module",
  "week5_translation_priority",
  "week4_focal_transposon_gene",
  "genage",
  "longevitymap",
  "cellage",
  "genome_instability",
  "repeat_control",
  "senescence",
  "cancer_disease",
  "neuro_disease",
  "immune_disease",
  "small_molecule",
  "max_open_targets_score",
  "small_molecule_tractability"
)]
write.table(plot_table, plot_table_path, sep = "\t", quote = FALSE, row.names = FALSE)

binary_cols <- c(
  "week4_focal_transposon_gene",
  "genage",
  "longevitymap",
  "cellage",
  "genome_instability",
  "repeat_control",
  "senescence",
  "cancer_disease",
  "neuro_disease",
  "immune_disease",
  "small_molecule"
)
binary_labels <- c(
  "Week 4\nfocal",
  "GenAge",
  "Longevity\nMap",
  "CellAge",
  "Genome\ninstability",
  "Repeat\ncontrol",
  "Senescence",
  "Cancer\ndisease",
  "Neuro\ndisease",
  "Immune /\ninflammation",
  "Small\nmolecule"
)

score_col_index <- length(binary_cols) + 1
x_labels <- c(binary_labels, "Open Targets\nmax score")

score_color <- function(score) {
  score <- max(0, min(0.9, score))
  breaks <- seq(0, 0.9, length.out = 101)
  palette <- colorRampPalette(c("#f7f7f7", "#fee08b", "#f46d43", "#a50026"))(100)
  palette[findInterval(score, breaks, all.inside = TRUE)]
}

plot_fun <- function() {
  n <- nrow(plot_table)
  x_positions <- seq_along(x_labels)
  y_positions <- rev(seq_len(n))
  par(mar = c(7.2, 5.6, 3.1, 7.6), las = 1, family = "sans")
  plot(
    NA,
    xlim = c(0.5, length(x_labels) + 3.1),
    ylim = c(0.5, n + 0.5),
    xaxt = "n",
    yaxt = "n",
    xlab = "",
    ylab = "",
    main = "Human translational evidence context",
    bty = "n"
  )
  abline(v = seq(1.5, length(x_labels) - 0.5, by = 1), col = "grey92", lwd = 0.8)
  abline(h = seq(1.5, n - 0.5, by = 1), col = "grey95", lwd = 0.6)
  axis(1, at = x_positions, labels = x_labels, las = 2, tick = FALSE, cex.axis = 0.73)
  row_labels <- plot_table$human_gene_symbol
  axis(2, at = y_positions, labels = row_labels, tick = FALSE, cex.axis = 0.56)

  module_col <- module_colors[plot_table$maintenance_module]
  module_col[is.na(module_col)] <- "grey60"
  points(rep(0.65, n), y_positions, pch = 15, col = module_col, cex = 0.68)

  for (i in seq_len(n)) {
    y <- y_positions[i]
    for (j in seq_along(binary_cols)) {
      val <- plot_table[i, binary_cols[j]]
      if (isTRUE(as.numeric(val) > 0)) {
        points(j, y, pch = 21, bg = "#1b9e77", col = "#1b9e77", cex = 1.2)
      } else {
        points(j, y, pch = 21, bg = "white", col = "grey78", cex = 0.8)
      }
    }
    score <- as.numeric(plot_table$max_open_targets_score[i])
    cex_score <- 0.65 + 1.75 * sqrt(max(score, 0))
    points(score_col_index, y, pch = 21, bg = score_color(score), col = "grey35", cex = cex_score)
  }

  legend(
    x = length(x_labels) + 0.55,
    y = n * 0.58,
    legend = unname(module_legend_labels[names(module_colors)]),
    fill = module_colors[names(module_colors)],
    border = NA,
    bty = "n",
    cex = 0.62,
    title = "Module",
    xpd = FALSE
  )
  legend(
    x = length(x_labels) + 0.55,
    y = 3.5,
    legend = c("Evidence present", "No evidence", "Open Targets score"),
    pt.bg = c("#1b9e77", "white", score_color(0.65)),
    col = c("#1b9e77", "grey78", "grey35"),
    pch = 21,
    pt.cex = c(1.1, 0.8, 1.6),
    bty = "n",
    cex = 0.66,
    xpd = FALSE
  )
}

png_path <- file.path(out_dir, "week5_translational_evidence_map.png")
pdf_path <- file.path(out_dir, "week5_translational_evidence_map.pdf")
png(png_path, width = 12.2, height = 12.5, units = "in", res = 300, type = "cairo")
plot_fun()
dev.off()
pdf(pdf_path, width = 12.2, height = 12.5, useDingbats = FALSE)
plot_fun()
dev.off()

trans <- plot_table[plot_table$maintenance_module == "transposon_suppression", ]
report_lines <- c(
  "# Week 5 Translational Evidence Figure Report",
  "",
  "Generated a gene-level evidence map linking the Week 4 focal transposon module to Week 5 human ageing, disease, and tractability evidence.",
  "",
  "## Files",
  "",
  paste0("- `", png_path, "`"),
  paste0("- `", pdf_path, "`"),
  paste0("- `", plot_table_path, "`"),
  "",
  "## Encoded Evidence Columns",
  "",
  "- Week 4 focal: gene belongs to the transposon-suppression module that carried the strict sequence-supported comparative signal.",
  "- GenAge, LongevityMap, CellAge: curated HAGR source hits.",
  "- Genome instability, repeat control, senescence: keyword-supported human gene context.",
  "- Cancer, neuro, immune/inflammation disease: Open Targets top-disease area flags.",
  "- Small molecule: Open Targets small-molecule tractability flag.",
  "- Open Targets max score: maximum disease-target association score among the retrieved top disease rows.",
  "",
  "## Transposon Rows",
  ""
)
for (i in seq_len(nrow(trans))) {
  report_lines <- c(
    report_lines,
    paste0(
      "- ", trans$human_gene_symbol[i],
      ": priority=", trans$week5_translation_priority[i],
      "; OpenTargetsMax=", sprintf("%.3f", trans$max_open_targets_score[i]),
      "; small_molecule=", ifelse(trans$small_molecule[i] == 1, "yes", "no"),
      "; repeat_control=", ifelse(trans$repeat_control[i] == 1, "yes", "no"),
      "; senescence=", ifelse(trans$senescence[i] == 1, "yes", "no")
    )
  )
}
report_lines <- c(
  report_lines,
  "",
  "## Interpretation",
  "",
  "This is a translational triage figure, not a causal ageing model. Its strongest use is to show that the comparative transposon signal has plausible human hooks through repeat control, disease association, and tractability for selected chromatin/repression genes."
)
writeLines(report_lines, report_path)

cat("Wrote", png_path, pdf_path, plot_table_path, "and", report_path, "\n")
