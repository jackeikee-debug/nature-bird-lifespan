#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  project_libs <- c(file.path(getwd(), "env", "R_library"), file.path(getwd(), "env", "R", "library"))
  project_libs <- project_libs[dir.exists(project_libs)]
  if (length(project_libs) > 0) .libPaths(c(project_libs, .libPaths()))
  library(ape)
  library(ggplot2)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)
tree_path <- ifelse(length(args) >= 1, args[[1]], "data/processed/phylogeny_inputs/opentree_induced_subtree.tre")
compact_table_path <- ifelse(length(args) >= 2, args[[2]], "data/processed/figure1_phylo_heatmap_compact_table.tsv")
full_table_path <- ifelse(length(args) >= 3, args[[3]], "data/processed/figure1_phylo_heatmap_table.tsv")
out_prefix <- ifelse(length(args) >= 4, args[[4]], "results/figures/figure1_phylogenetic_heatmap_R")
report_path <- ifelse(length(args) >= 5, args[[5]], "results/reports/figure1_phylogenetic_heatmap_R_report.md")
dated_tree_path <- ifelse(length(args) >= 6, args[[6]], "data/processed/phylogeny_inputs/opentree_datelife_calibrated_primary68.tre")

root_age_ma <- 320

clade_colors <- c(
  Aves = "#0F766E",
  Mammalia_Chiroptera = "#7C3AED",
  Mammalia_nonChiroptera = "#475569",
  Reptilia = "#B7791F"
)
tier_colors <- c(
  tier1_refseq_annotated_chromosome = "#176B87",
  tier1_refseq_annotated = "#2A9D8F",
  tier2_annotated = "#E9A23B",
  tier3_assembly_only = "#9CA3AF"
)

species_label <- function(x) {
  parts <- strsplit(x, " ", fixed = TRUE)
  vapply(parts, function(p) if (length(p) >= 2) paste0(substr(p[[1]], 1, 1), ". ", p[[2]]) else p[[1]], character(1))
}

make_approx_dated_tree <- function(tree, labels, root_age = 320) {
  tr <- drop.tip(tree, setdiff(tree$tip.label, labels))
  tr <- compute.brlen(tr, method = "Grafen", power = 1)
  depths <- node.depth.edgelength(tr)
  max_depth <- max(depths[seq_along(tr$tip.label)], na.rm = TRUE)
  tr$edge.length <- tr$edge.length * (root_age / max_depth)
  tr
}

make_plot_tree <- function(tree, labels, root_age = 320) {
  tr <- drop.tip(tree, setdiff(tree$tip.label, labels))
  if (!is.null(tr$edge.length) && all(is.finite(tr$edge.length)) && all(tr$edge.length > 0) && is.ultrametric(tr, tol = 1e-6)) {
    return(tr)
  }
  make_approx_dated_tree(tree, labels, root_age)
}

tree_segments <- function(tree) {
  tree <- ladderize(tree)
  n_tip <- length(tree$tip.label)
  yy <- numeric(n_tip + tree$Nnode)
  yy[seq_len(n_tip)] <- seq_len(n_tip)
  depths <- node.depth.edgelength(tree)
  for (node in rev(unique(tree$edge[, 1]))) {
    children <- tree$edge[tree$edge[, 1] == node, 2]
    yy[node] <- mean(yy[children])
  }
  edge_df <- as.data.frame(tree$edge)
  names(edge_df) <- c("parent", "child")
  edge_df$x <- depths[edge_df$parent]
  edge_df$xend <- depths[edge_df$child]
  edge_df$y <- yy[edge_df$child]
  edge_df$yend <- yy[edge_df$child]
  vertical_df <- edge_df |>
    group_by(parent) |>
    summarise(x = first(x), ymin = min(y), ymax = max(y), .groups = "drop")
  tips <- tibble(
    opentree_tip_label = tree$tip.label,
    x = depths[seq_len(n_tip)],
    y = yy[seq_len(n_tip)]
  )
  list(edge = edge_df, vertical = vertical_df, tips = tips, tree = tree)
}

plot_one <- function(plot_table, tree, title, subtitle, panel_label, out_prefix) {
  labels <- plot_table$opentree_tip_label
  tr <- make_plot_tree(tree, labels, root_age_ma)
  seg <- tree_segments(tr)
  tips <- seg$tips |>
    left_join(plot_table, by = "opentree_tip_label") |>
    mutate(
      label = species_label(scientific_name),
      clade = factor(clade, levels = names(clade_colors)),
      genome_analysis_tier = factor(genome_analysis_tier, levels = names(tier_colors)),
      heat_x_clade = max(seg$tips$x) + 28,
      heat_x_tier = heat_x_clade + 10,
      heat_x_residual = heat_x_tier + 13,
      heat_x_score = heat_x_residual + 18,
      heat_x_coverage = heat_x_score + 18,
      heat_x_strict = heat_x_coverage + 18,
      heat_x_ambig = heat_x_strict + 18
    )

  heat_long <- bind_rows(
    tips |> transmute(y, x = heat_x_residual, variable = "Residual\nlifespan", value = lifespan_residual_log10, scale = "residual"),
    tips |> transmute(y, x = heat_x_score, variable = "Repeat\nscore", value = transposon_repeat_suppression_score, scale = "score"),
    tips |> transmute(y, x = heat_x_coverage, variable = "Repeat\ncoverage", value = transposon_repeat_suppression_coverage, scale = "coverage"),
    tips |> transmute(y, x = heat_x_strict, variable = "Strict seq.\nrows", value = strict_sequence_rows / 10, raw_label = ifelse(targeted_sequence_validation, as.character(as.integer(strict_sequence_rows)), ""), scale = "strict"),
    tips |> transmute(y, x = heat_x_ambig, variable = "Ambig./\nunresolved", value = unresolved_or_ambiguous_fraction, scale = "ambig")
  )
  x_breaks <- c(
    unique(tips$heat_x_residual), unique(tips$heat_x_score), unique(tips$heat_x_coverage),
    unique(tips$heat_x_strict), unique(tips$heat_x_ambig)
  )
  x_labels <- c("Residual\nlifespan", "Repeat\nscore", "Repeat\ncoverage", "Strict seq.\nrows", "Ambig./\nunresolved")

  p <- ggplot() +
    geom_segment(data = seg$edge, aes(x = x, xend = xend, y = y, yend = yend), color = "#475569", linewidth = 0.22) +
    geom_segment(data = seg$vertical, aes(x = x, xend = x, y = ymin, yend = ymax), color = "#475569", linewidth = 0.22) +
    geom_point(data = tips, aes(x = x + 1.5, y = y, color = clade), size = 1.1) +
    geom_text(data = tips, aes(x = max(seg$tips$x) + 7, y = y, label = label), hjust = 0, size = 1.9, color = "#111827") +
    geom_tile(data = tips, aes(x = heat_x_clade, y = y, fill = clade), width = 5.2, height = 0.72, show.legend = FALSE) +
    scale_fill_manual(values = clade_colors, na.value = "#E5E7EB") +
    ggnewscale::new_scale_fill() +
    geom_tile(data = tips, aes(x = heat_x_tier, y = y, fill = genome_analysis_tier), width = 5.2, height = 0.72, show.legend = FALSE) +
    scale_fill_manual(values = tier_colors, na.value = "#E5E7EB") +
    ggnewscale::new_scale_fill() +
    geom_tile(data = heat_long |> filter(scale %in% c("score", "coverage", "strict")), aes(x = x, y = y, fill = value), width = 13.5, height = 0.72) +
    scale_fill_gradient(low = "#F0FDF4", high = "#006D2C", limits = c(0, 1), na.value = "#F3F4F6", name = "Score / coverage") +
    ggnewscale::new_scale_fill() +
    geom_tile(data = heat_long |> filter(scale == "residual"), aes(x = x, y = y, fill = value), width = 13.5, height = 0.72) +
    scale_fill_gradient2(low = "#A6CEE3", mid = "#FFF7ED", high = "#9F1239", midpoint = 0, limits = c(-0.45, 0.45), na.value = "#F3F4F6", name = "Residual") +
    ggnewscale::new_scale_fill() +
    geom_tile(data = heat_long |> filter(scale == "ambig"), aes(x = x, y = y, fill = value), width = 13.5, height = 0.72) +
    scale_fill_gradient(low = "#F8FAFC", high = "#DC2626", limits = c(0, 1), na.value = "#F3F4F6", name = "Ambig./unresolved") +
    geom_text(data = heat_long |> filter(scale == "strict", raw_label != ""), aes(x = x, y = y, label = raw_label), size = 1.65, color = "#111827") +
    annotate("text", x = 0, y = 0.0, label = panel_label, fontface = "bold", size = 5.2, hjust = 0) +
    annotate("text", x = 12, y = 0.0, label = title, fontface = "bold", size = 4.0, hjust = 0) +
    annotate("text", x = 12, y = 1.0, label = subtitle, size = 2.45, hjust = 0, color = "#475569") +
    annotate("text", x = unique(tips$heat_x_clade), y = 1.0, label = "Clade", angle = 50, hjust = 0, size = 2.15) +
    annotate("text", x = unique(tips$heat_x_tier), y = 1.0, label = "Genome\ntier", angle = 50, hjust = 0, size = 2.15) +
    annotate("text", x = x_breaks, y = 1.0, label = x_labels, angle = 50, hjust = 0, size = 2.15) +
    scale_color_manual(values = clade_colors, name = "Clade", labels = c("Birds", "Bats", "Other mammals", "Reptiles")) +
    scale_y_reverse(limits = c(max(tips$y) + 1, -1.2), expand = c(0, 0)) +
    scale_x_continuous(expand = c(0.005, 0.005)) +
    labs(caption = "Literature-calibrated branch lengths: fixed OpenTree topology with median node ages from peer-reviewed DateLife chronograms; four uncovered internal nodes use an audited Grafen fallback. PGLS sensitivity is reported separately.") +
    theme_void(base_size = 8) +
    theme(
      legend.position = "right",
      plot.margin = margin(6, 6, 6, 6),
      plot.caption = element_text(size = 6.5, color = "#64748B", hjust = 0),
      legend.title = element_text(size = 7, face = "bold"),
      legend.text = element_text(size = 6.5)
    )

  dir.create(dirname(out_prefix), recursive = TRUE, showWarnings = FALSE)
  ggsave(paste0(out_prefix, ".png"), p, width = 13.8, height = max(8.8, nrow(tips) * 0.18 + 1.2), dpi = 320, bg = "white")
  ggsave(paste0(out_prefix, ".pdf"), p, width = 13.8, height = max(8.8, nrow(tips) * 0.18 + 1.2), bg = "white")
  ggsave(paste0(out_prefix, ".svg"), p, width = 13.8, height = max(8.8, nrow(tips) * 0.18 + 1.2), bg = "white")
  tr
}

tree <- if (file.exists(dated_tree_path)) read.tree(dated_tree_path) else read.tree(tree_path)
compact_table <- read_tsv(compact_table_path, show_col_types = FALSE)
full_table <- read_tsv(full_table_path, show_col_types = FALSE)

if (!requireNamespace("ggnewscale", quietly = TRUE)) {
  stop("R package ggnewscale is required. Install it with install.packages('ggnewscale').")
}

compact_tree <- plot_one(
  compact_table,
  tree,
  "R-rendered phylogenetic heatmap of lifespan residuals and repeat-control evidence",
  "Compact manuscript view: all targeted validation birds plus anchor and clade-representative species.",
  "A",
  out_prefix
)
full_prefix <- sub("figure1_phylogenetic_heatmap_R$", "figureS1_phylogenetic_heatmap_full_R", out_prefix)
full_tree <- plot_one(
  full_table,
  tree,
  "R-rendered full 68-species genome-panel phylogeny",
  "Supplementary full-panel view using the same approximate dated branch-length transform.",
  "S1",
  full_prefix
)

if (!file.exists(dated_tree_path)) {
  dir.create(dirname(dated_tree_path), recursive = TRUE, showWarnings = FALSE)
  write.tree(make_approx_dated_tree(tree, full_table$opentree_tip_label, root_age_ma), file = dated_tree_path)
}

lines <- c(
  "# R Figure 1 Phylogenetic Heatmap Report",
  "",
  "Generated an R-rendered Figure 1 candidate and full supplementary tree heatmap.",
  "",
  "## Outputs",
  "",
  paste0("- compact PNG: `", out_prefix, ".png`"),
  paste0("- compact PDF: `", out_prefix, ".pdf`"),
  paste0("- compact SVG: `", out_prefix, ".svg`"),
  paste0("- full supplementary PNG: `", full_prefix, ".png`"),
  paste0("- literature-calibrated Newick: `", dated_tree_path, "`"),
  "",
  "## Branch-Length Caveat",
  "",
  "The plotted tree retains the OpenTree topology and uses median node ages from peer-reviewed chronograms in the DateLife cache. Four uncovered internal nodes use an explicitly audited Grafen fallback. This is a secondary literature-calibrated chronogram, not a newly inferred sequence-based species tree."
)
dir.create(dirname(report_path), recursive = TRUE, showWarnings = FALSE)
writeLines(lines, report_path)

cat("Wrote R Figure 1 outputs and report\n")
