#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ape)
})

root <- normalizePath(getwd())

traits_path <- file.path(root, "data", "processed", "jme_revision2_analysis_table.tsv")
pc1_path <- file.path(root, "results", "tables", "jme_revision2_shared_axis_species.tsv")
tree_path <- file.path(root, "data", "processed", "phylogeny_inputs", "opentree_datelife_calibrated_primary68.tre")
figure_dir <- file.path(root, "results", "figures")
table_dir <- file.path(root, "results", "tables")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

traits <- read.delim(traits_path, check.names = FALSE, stringsAsFactors = FALSE)
pc1 <- read.delim(pc1_path, check.names = FALSE, stringsAsFactors = FALSE)
tree <- read.tree(tree_path)

pc1_column <- if ("shared_module_pc1" %in% names(pc1)) "shared_module_pc1" else "pc1"
traits$shared_module_pc1 <- pc1[[pc1_column]][match(traits$scientific_name, pc1$scientific_name)]
coverage_columns <- grep("_coverage$", names(traits), value = TRUE)
traits$mean_module_coverage <- rowMeans(traits[, coverage_columns], na.rm = TRUE)
traits$log10_max_lifespan <- log10(traits$max_lifespan_years)
traits$assembly_tier <- ifelse(grepl("^tier1", traits$genome_analysis_tier), "Tier 1", "Tier 2")

missing_tips <- setdiff(tree$tip.label, traits$opentree_tip_label)
extra_rows <- setdiff(traits$opentree_tip_label, tree$tip.label)
if (length(missing_tips) || length(extra_rows)) {
  stop(sprintf("Tree/metadata mismatch: %d missing tips and %d extra rows", length(missing_tips), length(extra_rows)))
}

meta <- traits[match(tree$tip.label, traits$opentree_tip_label), ]
stopifnot(identical(meta$opentree_tip_label, tree$tip.label))

clade_colors <- c(
  Aves = "#2F6B9A",
  Mammalia_Chiroptera = "#D97A35",
  Mammalia_nonChiroptera = "#3D8C7C",
  Reptilia = "#B59A34"
)
tier_colors <- c("Tier 1" = "#284B63", "Tier 2" = "#B8C4CC")
missing_color <- "#D9DDE1"

hex_to_rgb <- function(hex) {
  col2rgb(hex) / 255
}

mix_colors <- function(low, high, n = 256) {
  rgb(colorRamp(c(low, high))(seq(0, 1, length.out = n)), maxColorValue = 255)
}

diverging_colors <- colorRampPalette(c("#2C6EAB", "#F7F7F4", "#B9473E"))(256)
sequential_colors <- colorRampPalette(c("#F1F4F1", "#9CC4B8", "#236A73"))(256)
lifespan_colors <- colorRampPalette(c("#F4EFE3", "#E2B85B", "#9D4F2D"))(256)

numeric_colors <- function(values, palette, symmetric = FALSE) {
  out <- rep(missing_color, length(values))
  valid <- is.finite(values)
  if (!any(valid)) return(out)
  if (symmetric) {
    limit <- max(abs(values[valid]))
    range <- c(-limit, limit)
  } else {
    range <- range(values[valid])
  }
  if (diff(range) == 0) {
    index <- rep(ceiling(length(palette) / 2), sum(valid))
  } else {
    scaled <- (values[valid] - range[1]) / diff(range)
    index <- pmax(1, pmin(length(palette), floor(scaled * (length(palette) - 1)) + 1))
  }
  out[valid] <- palette[index]
  out
}

track_spec <- list(
  list(key = "clade", label = "Clade", colors = unname(clade_colors[meta$clade])),
  list(key = "lifespan_residual_log10", label = "Lifespan residual", colors = numeric_colors(meta$lifespan_residual_log10, diverging_colors, TRUE)),
  list(key = "log10_max_lifespan", label = "Maximum lifespan", colors = numeric_colors(meta$log10_max_lifespan, lifespan_colors)),
  list(key = "busco_complete", label = "BUSCO completeness", colors = numeric_colors(meta$busco_complete, sequential_colors)),
  list(key = "assembly_tier", label = "Assembly tier", colors = unname(tier_colors[meta$assembly_tier])),
  list(key = "mean_module_coverage", label = "Mean module coverage", colors = numeric_colors(meta$mean_module_coverage, sequential_colors)),
  list(key = "shared_module_pc1", label = "Shared module PC1", colors = numeric_colors(meta$shared_module_pc1, diverging_colors, TRUE)),
  list(key = "transposon_repeat_suppression_score", label = "Repeat-control score", colors = numeric_colors(meta$transposon_repeat_suppression_score, sequential_colors))
)

draw_wedge <- function(r_inner, r_outer, theta, half_width, fill) {
  angle <- seq(theta - half_width, theta + half_width, length.out = 16)
  polygon(
    c(r_inner * cos(angle), rev(r_outer * cos(angle))),
    c(r_inner * sin(angle), rev(r_outer * sin(angle))),
    col = fill,
    border = "white",
    lwd = 0.18
  )
}

draw_gradient <- function(x0, y0, x1, y1, palette, border = "#B7BEC4") {
  n <- length(palette)
  width <- (x1 - x0) / n
  for (i in seq_len(n)) {
    rect(x0 + (i - 1) * width, y0, x0 + i * width, y1, col = palette[i], border = NA)
  }
  rect(x0, y0, x1, y1, border = border, lwd = 0.6)
}

plot_circular <- function(device) {
  if (device == "png") {
    png(file.path(figure_dir, "jme_revision2_figure1_circular_phylogeny.png"), width = 7200, height = 5500, res = 500, bg = "white")
  } else {
    pdf(file.path(figure_dir, "jme_revision2_figure1_circular_phylogeny.pdf"), width = 14.4, height = 11, useDingbats = FALSE)
  }
  on.exit(dev.off(), add = TRUE)

  layout(matrix(c(1, 2), nrow = 1), widths = c(4.9, 1.55))
  par(oma = c(0.2, 0.3, 1.0, 0.3), family = "sans")
  tree_radius_target <- max(node.depth.edgelength(tree))
  plotting_limit <- tree_radius_target * 1.70
  par(mar = c(0.2, 0.2, 0.2, 0.2), xpd = FALSE)
  plot.phylo(
    tree,
    type = "fan",
    show.tip.label = FALSE,
    edge.color = "#64737C",
    edge.width = 0.72,
    no.margin = TRUE,
    rotate.tree = 90,
    x.lim = c(-plotting_limit, plotting_limit),
    y.lim = c(-plotting_limit, plotting_limit)
  )
  state <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  tip_x <- state$xx[seq_len(Ntip(tree))]
  tip_y <- state$yy[seq_len(Ntip(tree))]
  theta <- atan2(tip_y, tip_x)
  radius <- max(sqrt(tip_x^2 + tip_y^2))
  half_width <- pi / Ntip(tree) * 0.93
  ring_width <- radius * 0.055
  ring_gap <- radius * 0.008
  ring_start <- radius * 1.025

  for (track_index in seq_along(track_spec)) {
    inner <- ring_start + (track_index - 1) * (ring_width + ring_gap)
    outer <- inner + ring_width
    for (tip_index in seq_len(Ntip(tree))) {
      draw_wedge(inner, outer, theta[tip_index], half_width, track_spec[[track_index]]$colors[tip_index])
    }
  }

  outer_radius <- ring_start + length(track_spec) * (ring_width + ring_gap)
  label_radius <- outer_radius + radius * 0.055
  species_labels <- gsub("_ott[0-9]+$", "", tree$tip.label)
  species_labels <- gsub("_", " ", species_labels)
  theta_degrees <- theta * 180 / pi
  for (i in seq_len(Ntip(tree))) {
    angle <- theta_degrees[i]
    flip <- angle > 90 || angle < -90
    rotation <- if (flip) angle + 180 else angle
    alignment <- if (flip) 1 else 0
    text(
      label_radius * cos(theta[i]),
      label_radius * sin(theta[i]),
      labels = species_labels[i],
      srt = rotation,
      adj = c(alignment, 0.5),
      cex = 0.42,
      font = 3,
      col = "#27343A"
    )
  }

  par(mar = c(0.8, 0.4, 0.8, 0.5), xpd = NA)
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))
  text(0, 0.97, "Ring order", adj = c(0, 1), cex = 1.05, font = 2, col = "#1F2D33")
  text(0, 0.925, "Inner to outer", adj = c(0, 1), cex = 0.72, col = "#68767D")

  y <- 0.875
  for (i in seq_along(track_spec)) {
    rect(0, y - 0.017, 0.075, y + 0.017, col = "#EEF1F2", border = "#AEB7BC", lwd = 0.6)
    text(0.0375, y, as.character(i), cex = 0.62, font = 2, col = "#314149")
    text(0.095, y, track_spec[[i]]$label, adj = c(0, 0.5), cex = 0.72, col = "#27343A")
    y <- y - 0.055
  }

  y <- 0.405
  text(0, y, "Clade", adj = c(0, 1), cex = 0.78, font = 2, col = "#1F2D33")
  y <- y - 0.042
  clade_labels <- c(Aves = "Birds", Mammalia_Chiroptera = "Bats", Mammalia_nonChiroptera = "Other mammals", Reptilia = "Reptiles")
  for (key in names(clade_labels)) {
    rect(0, y - 0.012, 0.055, y + 0.012, col = clade_colors[key], border = NA)
    text(0.073, y, clade_labels[key], adj = c(0, 0.5), cex = 0.66, col = "#27343A")
    y <- y - 0.038
  }

  tier_y <- 0.405
  text(0.55, tier_y, "Assembly tier", adj = c(0, 1), cex = 0.78, font = 2, col = "#1F2D33")
  tier_y <- tier_y - 0.045
  for (key in names(tier_colors)) {
    rect(0.55, tier_y - 0.012, 0.605, tier_y + 0.012, col = tier_colors[key], border = NA)
    text(0.623, tier_y, key, adj = c(0, 0.5), cex = 0.62, col = "#27343A")
    tier_y <- tier_y - 0.038
  }

  y <- 0.215
  text(0, y, "Continuous tracks", adj = c(0, 1), cex = 0.78, font = 2, col = "#1F2D33")
  y <- y - 0.05
  draw_gradient(0, y - 0.012, 0.34, y + 0.012, diverging_colors)
  text(0, y - 0.028, "low", adj = c(0, 1), cex = 0.57, col = "#68767D")
  text(0.34, y - 0.028, "high", adj = c(1, 1), cex = 0.57, col = "#68767D")
  y <- y - 0.075
  draw_gradient(0, y - 0.012, 0.34, y + 0.012, sequential_colors)
  text(0, y - 0.028, "low", adj = c(0, 1), cex = 0.57, col = "#68767D")
  text(0.34, y - 0.028, "high", adj = c(1, 1), cex = 0.57, col = "#68767D")
  rect(0.48, y - 0.012, 0.535, y + 0.012, col = missing_color, border = "#AEB7BC", lwd = 0.5)
  text(0.55, y, "missing", adj = c(0, 0.5), cex = 0.57, col = "#68767D")

  text(0, 0.015, "Tree: OpenTree topology with DateLife-calibrated branch lengths.\nFour internal nodes used documented Grafen fallback ages.", adj = c(0, 0), cex = 0.57, col = "#68767D")

  mtext("Phylogeny, lifespan, and gene observability across the 68-species panel", outer = TRUE, side = 3, line = -0.05, adj = 0.02, cex = 1.18, font = 2, col = "#1F2D33")
}

plot_circular("png")
plot_circular("pdf")

source_data <- data.frame(
  tree_tip_order = seq_len(Ntip(tree)),
  tree_tip_label = tree$tip.label,
  scientific_name = meta$scientific_name,
  clade = meta$clade,
  lifespan_residual_log10 = meta$lifespan_residual_log10,
  max_lifespan_years = meta$max_lifespan_years,
  log10_max_lifespan = meta$log10_max_lifespan,
  busco_complete = meta$busco_complete,
  assembly_tier = meta$assembly_tier,
  mean_module_coverage = meta$mean_module_coverage,
  shared_module_pc1 = meta$shared_module_pc1,
  transposon_repeat_suppression_score = meta$transposon_repeat_suppression_score,
  stringsAsFactors = FALSE
)
write.table(
  source_data,
  file.path(table_dir, "jme_revision2_figure1_circular_phylogeny.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE,
  na = "NA"
)

cat("Wrote circular phylogenetic heatmap and source data.\n")
