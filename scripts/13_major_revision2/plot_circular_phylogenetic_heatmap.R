#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ape)
})

root <- normalizePath(getwd())

traits_path <- file.path(root, "data", "processed", "jme_revision2_analysis_table.tsv")
pc1_path <- file.path(root, "results", "tables", "jme_revision2_shared_axis_species.tsv")
samhd1_rows_path <- file.path(root, "results", "tables", "targeted_domain_conservation_rows.tsv")
samhd1_position_path <- file.path(root, "results", "tables", "samhd1_alignment_position_qc.tsv")
samhd1_alignment_path <- file.path(root, "data", "interim", "protein_conservation", "SAMHD1.aligned.faa")
tree_path <- file.path(root, "data", "processed", "phylogeny_inputs", "opentree_datelife_calibrated_primary68.tre")
figure_dir <- file.path(root, "results", "figures")
table_dir <- file.path(root, "results", "tables")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

traits <- read.delim(traits_path, check.names = FALSE, stringsAsFactors = FALSE)
pc1 <- read.delim(pc1_path, check.names = FALSE, stringsAsFactors = FALSE)
domain_rows <- read.delim(samhd1_rows_path, check.names = FALSE, stringsAsFactors = FALSE)
samhd1_positions <- read.delim(samhd1_position_path, check.names = FALSE, stringsAsFactors = FALSE)
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
no_sequence_color <- "#F1F3F3"
residue_colors <- c(identical = "#2B837C", substitution = "#E6B09D", gap = "#E6EAEC")
domain_colors <- c(SAM = "#6AA89E", HD = "#D46A4C", non_domain = "#BFC7CA")
fingerprint_palette <- colorRampPalette(c("#F4F2EE", "#D2A3B2", "#67466F"))(256)
fingerprint_genes <- c("DNMT1", "DNMT3A", "DNMT3B", "HELLS", "MBD2", "MBD3", "MORC3", "SETDB2", "UHRF1")

hex_to_rgb <- function(hex) {
  col2rgb(hex) / 255
}

mix_colors <- function(low, high, n = 256) {
  rgb(colorRamp(c(low, high))(seq(0, 1, length.out = n)), maxColorValue = 255)
}

diverging_colors <- colorRampPalette(c("#2C6EAB", "#F7F7F4", "#B9473E"))(256)
sequential_colors <- colorRampPalette(c("#F1F4F1", "#9CC4B8", "#236A73"))(256)
lifespan_colors <- colorRampPalette(c("#F4EFE3", "#E2B85B", "#9D4F2D"))(256)
track_background <- "#F1F3F2"
bar_lifespan_color <- "#C58A2A"
bar_coverage_color <- "#4F9388"
dot_busco_color <- "#335F7A"
dot_repeat_color <- "#7C4F78"

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

scale01 <- function(values) {
  out <- rep(NA_real_, length(values))
  valid <- is.finite(values)
  if (!any(valid)) return(out)
  value_range <- range(values[valid])
  if (diff(value_range) == 0) {
    out[valid] <- 0.5
  } else {
    out[valid] <- (values[valid] - value_range[1]) / diff(value_range)
  }
  out
}

scale_abs01 <- function(values) {
  out <- rep(NA_real_, length(values))
  valid <- is.finite(values)
  if (!any(valid)) return(out)
  limit <- max(abs(values[valid]))
  out[valid] <- if (limit == 0) 0 else abs(values[valid]) / limit
  out
}

read_fasta_strings <- function(path) {
  lines <- readLines(path, warn = FALSE)
  headers <- which(startsWith(lines, ">"))
  output <- character(length(headers))
  names(output) <- sub("^>([^ ]+).*$", "\\1", lines[headers])
  for (i in seq_along(headers)) {
    start <- headers[i] + 1
    end <- if (i < length(headers)) headers[i + 1] - 1 else length(lines)
    output[i] <- paste(lines[start:end], collapse = "")
  }
  output
}

mix_residue_colors <- function(identical, substitution, gap) {
  weights <- c(identical, substitution, gap)
  palette <- col2rgb(residue_colors) / 255
  mixed <- as.numeric(palette %*% weights)
  rgb(mixed[1], mixed[2], mixed[3])
}

samhd1_alignment <- read_fasta_strings(samhd1_alignment_path)
samhd1_rows <- domain_rows[domain_rows$human_gene_symbol == "SAMHD1", ]
samhd1_rows <- samhd1_rows[match(meta$scientific_name, samhd1_rows$scientific_name), ]
stopifnot(identical(meta$scientific_name, samhd1_rows$scientific_name))
samhd1_rows$sequence_available_flag <- tolower(as.character(samhd1_rows$sequence_available)) %in% c("true", "1", "yes")

samhd1_window_count <- 24L
reference_length <- max(samhd1_positions$human_reference_residue)
samhd1_positions$window_index <- pmin(
  samhd1_window_count,
  ceiling(samhd1_positions$human_reference_residue / reference_length * samhd1_window_count)
)

window_rows <- list()
row_counter <- 0L
for (tip_index in seq_len(nrow(meta))) {
  alignment_id <- samhd1_rows$alignment_record_id[tip_index]
  sequence_available <- samhd1_rows$sequence_available_flag[tip_index] &&
    !is.na(alignment_id) && alignment_id %in% names(samhd1_alignment)
  target <- if (sequence_available) strsplit(samhd1_alignment[[alignment_id]], "", fixed = TRUE)[[1]] else character()
  for (window_index in seq_len(samhd1_window_count)) {
    positions <- samhd1_positions[samhd1_positions$window_index == window_index, ]
    identical_fraction <- substitution_fraction <- gap_fraction <- NA_real_
    plot_color <- no_sequence_color
    if (sequence_available) {
      target_residues <- toupper(target[positions$alignment_column_1based])
      reference_residues <- toupper(positions$reference_amino_acid)
      gaps <- target_residues %in% c("-", ".", "?")
      identical <- !gaps & target_residues == reference_residues
      substitution <- !gaps & !identical
      identical_fraction <- mean(identical)
      substitution_fraction <- mean(substitution)
      gap_fraction <- mean(gaps)
      plot_color <- mix_residue_colors(identical_fraction, substitution_fraction, gap_fraction)
    }
    domain_counts <- table(positions$domain_class)
    dominant_domain <- names(domain_counts)[which.max(domain_counts)]
    row_counter <- row_counter + 1L
    window_rows[[row_counter]] <- data.frame(
      tree_tip_order = tip_index,
      scientific_name = meta$scientific_name[tip_index],
      alignment_record_id = if (sequence_available) alignment_id else NA_character_,
      sequence_available = sequence_available,
      domain_analysis_qualified = sequence_available && is.finite(samhd1_rows$domain_reference_coverage[tip_index]) && samhd1_rows$domain_reference_coverage[tip_index] >= 0.5,
      window_index = window_index,
      human_residue_start = min(positions$human_reference_residue),
      human_residue_end = max(positions$human_reference_residue),
      dominant_domain = dominant_domain,
      identical_fraction = identical_fraction,
      substitution_fraction = substitution_fraction,
      gap_fraction = gap_fraction,
      plot_color = plot_color,
      stringsAsFactors = FALSE
    )
  }
}
samhd1_window_data <- do.call(rbind, window_rows)

fingerprint_rows <- list()
row_counter <- 0L
for (gene_index in seq_along(fingerprint_genes)) {
  gene <- fingerprint_genes[gene_index]
  gene_rows <- domain_rows[domain_rows$human_gene_symbol == gene, ]
  gene_rows <- gene_rows[match(meta$scientific_name, gene_rows$scientific_name), ]
  stopifnot(identical(meta$scientific_name, gene_rows$scientific_name))
  available <- tolower(as.character(gene_rows$sequence_available)) %in% c("true", "1", "yes")
  score <- gene_rows$domain_identity_coverage_product
  score[!available] <- NA_real_
  score <- pmax(0, pmin(1, score))
  color_index <- pmax(1, pmin(length(fingerprint_palette), floor(score * (length(fingerprint_palette) - 1)) + 1))
  plot_color <- rep(no_sequence_color, nrow(gene_rows))
  plot_color[is.finite(score)] <- fingerprint_palette[color_index[is.finite(score)]]
  for (tip_index in seq_len(nrow(meta))) {
    row_counter <- row_counter + 1L
    fingerprint_rows[[row_counter]] <- data.frame(
      tree_tip_order = tip_index,
      scientific_name = meta$scientific_name[tip_index],
      gene_ring_order = gene_index,
      human_gene_symbol = gene,
      sequence_available = available[tip_index],
      domain_reference_coverage = gene_rows$domain_reference_coverage[tip_index],
      domain_aligned_identity = gene_rows$domain_aligned_identity[tip_index],
      domain_identity_coverage_product = score[tip_index],
      plot_color = plot_color[tip_index],
      stringsAsFactors = FALSE
    )
  }
}
fingerprint_data <- do.call(rbind, fingerprint_rows)

track_spec <- list(
  list(key = "clade", label = "Clade", glyph = "band", colors = unname(clade_colors[meta$clade])),
  list(key = "lifespan_residual_log10", label = "Lifespan residual", glyph = "bubble", values = meta$lifespan_residual_log10, scaled = scale_abs01(meta$lifespan_residual_log10), colors = numeric_colors(meta$lifespan_residual_log10, diverging_colors, TRUE)),
  list(key = "log10_max_lifespan", label = "Maximum lifespan", glyph = "bar", values = meta$log10_max_lifespan, scaled = scale01(meta$log10_max_lifespan), color = bar_lifespan_color),
  list(key = "busco_complete", label = "BUSCO completeness", glyph = "dot", values = meta$busco_complete, scaled = scale01(meta$busco_complete), color = dot_busco_color),
  list(key = "assembly_tier", label = "Assembly tier", glyph = "band", colors = unname(tier_colors[meta$assembly_tier])),
  list(key = "mean_module_coverage", label = "Mean module coverage", glyph = "bar", values = meta$mean_module_coverage, scaled = scale01(meta$mean_module_coverage), color = bar_coverage_color),
  list(key = "shared_module_pc1", label = "Shared module PC1", glyph = "bubble", values = meta$shared_module_pc1, scaled = scale_abs01(meta$shared_module_pc1), colors = numeric_colors(meta$shared_module_pc1, diverging_colors, TRUE)),
  list(key = "transposon_repeat_suppression_score", label = "Repeat-control score", glyph = "dot", values = meta$transposon_repeat_suppression_score, scaled = scale01(meta$transposon_repeat_suppression_score), color = dot_repeat_color)
)

draw_wedge <- function(r_inner, r_outer, theta, half_width, fill, border = "white", lwd = 0.18) {
  angle <- seq(theta - half_width, theta + half_width, length.out = 16)
  polygon(
    c(r_inner * cos(angle), rev(r_outer * cos(angle))),
    c(r_inner * sin(angle), rev(r_outer * sin(angle))),
    col = fill,
    border = border,
    lwd = lwd
  )
}

draw_hatched_wedge <- function(r_inner, r_outer, theta, half_width) {
  angle <- seq(theta - half_width, theta + half_width, length.out = 24)
  polygon(
    c(r_inner * cos(angle), rev(r_outer * cos(angle))),
    c(r_inner * sin(angle), rev(r_outer * sin(angle))),
    col = "#AAB4B9",
    border = "#C4CBCF",
    density = 34,
    angle = 45,
    lwd = 0.25
  )
}

draw_residue_composition <- function(r_inner, r_outer, theta, half_width, identical, substitution, gap) {
  fractions <- c(identical = identical, substitution = substitution, gap = gap)
  cursor <- r_inner
  for (state in names(fractions)) {
    segment_outer <- cursor + (r_outer - r_inner) * fractions[state]
    if (segment_outer > cursor) {
      draw_wedge(cursor, segment_outer, theta, half_width, residue_colors[state], border = NA, lwd = 0.05)
    }
    cursor <- segment_outer
  }
  draw_wedge(r_inner, r_outer, theta, half_width, NA, border = "white", lwd = 0.10)
}

draw_radial_circle <- function(radius, color, width = 1) {
  angle <- seq(0, 2 * pi, length.out = 720)
  lines(radius * cos(angle), radius * sin(angle), col = color, lwd = width)
}

draw_circle <- function(radius, theta, circle_radius, fill, border = "white", lwd = 0.35) {
  angle <- seq(0, 2 * pi, length.out = 32)
  polygon(
    radius * cos(theta) + circle_radius * cos(angle),
    radius * sin(theta) + circle_radius * sin(angle),
    col = fill,
    border = border,
    lwd = lwd
  )
}

draw_legend_circle <- function(x, y, radius, fill, border = "white") {
  angle <- seq(0, 2 * pi, length.out = 32)
  polygon(x + radius * cos(angle), y + radius * sin(angle), col = fill, border = border, lwd = 0.5)
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

  layout(matrix(c(1, 2), nrow = 1), widths = c(4.75, 1.75))
  par(oma = c(0.2, 0.3, 1.0, 0.3), family = "sans")
  tree_radius_reference <- max(node.depth.edgelength(tree))
  tree_plot <- tree
  tree_plot$edge.length <- tree$edge.length * 0.72
  plotting_limit <- tree_radius_reference * 2.12
  par(mar = c(0.2, 0.2, 0.2, 0.2), xpd = FALSE)
  plot.phylo(
    tree_plot,
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
  ring_width <- tree_radius_reference * 0.075
  ring_gap <- tree_radius_reference * 0.007
  ring_start <- radius + tree_radius_reference * 0.012

  for (track_index in seq_along(track_spec)) {
    inner <- ring_start + (track_index - 1) * (ring_width + ring_gap)
    outer <- inner + ring_width
    track <- track_spec[[track_index]]
    for (tip_index in seq_len(Ntip(tree))) {
      if (track$glyph == "band") {
        draw_wedge(inner, outer, theta[tip_index], half_width, track$colors[tip_index])
      } else {
        draw_wedge(inner, outer, theta[tip_index], half_width, track_background)
        scaled_value <- track$scaled[tip_index]
        if (!is.finite(scaled_value)) next
        if (track$glyph == "bar") {
          bar_outer <- inner + ring_width * (0.10 + 0.86 * scaled_value)
          draw_wedge(inner + ring_width * 0.06, bar_outer, theta[tip_index], half_width * 0.74, track$color)
        } else if (track$glyph == "bubble") {
          bubble_radius <- ring_width * (0.13 + 0.30 * sqrt(scaled_value))
          draw_circle((inner + outer) / 2, theta[tip_index], bubble_radius, track$colors[tip_index], "white", 0.45)
        } else if (track$glyph == "dot") {
          dot_radius <- inner + ring_width * (0.12 + 0.76 * scaled_value)
          draw_circle(dot_radius, theta[tip_index], ring_width * 0.145, track$color, "white", 0.45)
        }
      }
    }
  }

  fingerprint_gap <- tree_radius_reference * 0.013
  fingerprint_ring_width <- tree_radius_reference * 0.011
  fingerprint_ring_gap <- tree_radius_reference * 0.0015
  fingerprint_start <- ring_start + length(track_spec) * (ring_width + ring_gap) + fingerprint_gap
  for (gene_index in seq_along(fingerprint_genes)) {
    inner <- fingerprint_start + (gene_index - 1) * (fingerprint_ring_width + fingerprint_ring_gap)
    outer <- inner + fingerprint_ring_width
    for (tip_index in seq_len(Ntip(tree))) {
      row_index <- (gene_index - 1) * Ntip(tree) + tip_index
      draw_wedge(inner, outer, theta[tip_index], half_width * 0.88, fingerprint_data$plot_color[row_index], border = "white", lwd = 0.08)
    }
  }
  fingerprint_outer <- fingerprint_start + length(fingerprint_genes) * (fingerprint_ring_width + fingerprint_ring_gap) - fingerprint_ring_gap
  draw_radial_circle(fingerprint_start, "#8C6A8E", 0.8)
  draw_radial_circle(fingerprint_outer, "#8C6A8E", 0.8)

  sequence_gap <- tree_radius_reference * 0.014
  sequence_ring_width <- tree_radius_reference * 0.015
  sequence_ring_gap <- tree_radius_reference * 0.0015
  sequence_start <- fingerprint_outer + sequence_gap
  sequence_outer <- sequence_start + samhd1_window_count * (sequence_ring_width + sequence_ring_gap) - sequence_ring_gap

  for (tip_index in seq_len(Ntip(tree))) {
    if (!samhd1_rows$sequence_available_flag[tip_index]) {
      draw_wedge(sequence_start, sequence_outer, theta[tip_index], half_width * 0.90, no_sequence_color, border = "white", lwd = 0.20)
      draw_hatched_wedge(sequence_start, sequence_outer, theta[tip_index], half_width * 0.82)
    }
  }

  for (window_index in seq_len(samhd1_window_count)) {
    inner <- sequence_start + (window_index - 1) * (sequence_ring_width + sequence_ring_gap)
    outer <- inner + sequence_ring_width
    for (tip_index in seq_len(Ntip(tree))) {
      row_index <- (tip_index - 1) * samhd1_window_count + window_index
      if (samhd1_window_data$sequence_available[row_index]) {
        draw_residue_composition(
          inner,
          outer,
          theta[tip_index],
          half_width * 0.90,
          samhd1_window_data$identical_fraction[row_index],
          samhd1_window_data$substitution_fraction[row_index],
          samhd1_window_data$gap_fraction[row_index]
        )
      }
    }
  }

  sam_inner <- sequence_start + (2 - 1) * (sequence_ring_width + sequence_ring_gap)
  sam_outer <- sequence_start + 5 * (sequence_ring_width + sequence_ring_gap) - sequence_ring_gap
  hd_inner <- sequence_start + (7 - 1) * (sequence_ring_width + sequence_ring_gap)
  hd_outer <- sequence_start + 9 * (sequence_ring_width + sequence_ring_gap) - sequence_ring_gap
  draw_radial_circle(sam_inner, domain_colors["SAM"], 1.15)
  draw_radial_circle(sam_outer, domain_colors["SAM"], 1.15)
  draw_radial_circle(hd_inner, domain_colors["HD"], 1.15)
  draw_radial_circle(hd_outer, domain_colors["HD"], 1.15)

  samhd1_non_gap_coverage <- vapply(seq_len(Ntip(tree)), function(tip_index) {
    rows <- samhd1_window_data$tree_tip_order == tip_index & samhd1_window_data$sequence_available
    if (!any(rows)) return(NA_real_)
    mean(1 - samhd1_window_data$gap_fraction[rows])
  }, numeric(1))
  audit_marker_inner <- sequence_outer + tree_radius_reference * 0.006
  audit_track_width <- tree_radius_reference * 0.030
  for (tip_index in seq_len(Ntip(tree))) {
    if (samhd1_rows$sequence_available_flag[tip_index]) {
      coverage_outer <- audit_marker_inner + audit_track_width * samhd1_non_gap_coverage[tip_index]
      draw_wedge(audit_marker_inner, coverage_outer, theta[tip_index], half_width * 0.66, "#477F87", border = NA, lwd = 0.05)
      if (is.finite(samhd1_rows$domain_reference_coverage[tip_index]) && samhd1_rows$domain_reference_coverage[tip_index] >= 0.5) {
        draw_wedge(coverage_outer, coverage_outer + tree_radius_reference * 0.006, theta[tip_index], half_width * 0.70, "#25343B", border = NA, lwd = 0.05)
      }
    }
  }

  center_busco_n <- sum(is.finite(meta$busco_complete))
  text(0, radius * 0.105, sprintf("%d", Ntip(tree)), cex = 1.35, font = 2, col = "#1F2D33")
  text(0, radius * 0.025, "species", cex = 0.63, col = "#5D6B72")
  text(0, -radius * 0.070, "200 genes | 6 modules", cex = 0.60, font = 2, col = "#33464F")
  text(0, -radius * 0.145, sprintf("BUSCO observed: %d", center_busco_n), cex = 0.52, col = "#68767D")

  outer_radius <- audit_marker_inner + audit_track_width + tree_radius_reference * 0.008
  label_radius <- outer_radius + tree_radius_reference * 0.038
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
      cex = 0.36,
      font = 3,
      col = "#27343A"
    )
  }

  par(mar = c(0.8, 0.4, 0.8, 0.5), xpd = NA)
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))
  text(0, 0.97, "Tracks and sequence audit", adj = c(0, 1), cex = 0.98, font = 2, col = "#1F2D33")
  text(0, 0.932, "Inner to outer", adj = c(0, 1), cex = 0.66, col = "#68767D")

  y <- 0.890
  for (i in seq_along(track_spec)) {
    glyph <- track_spec[[i]]$glyph
    rect(0, y - 0.015, 0.066, y + 0.015, col = "#EEF1F2", border = "#AEB7BC", lwd = 0.6)
    text(0.033, y, as.character(i), cex = 0.57, font = 2, col = "#314149")
    if (glyph == "band") {
      rect(0.083, y - 0.010, 0.125, y + 0.010, col = "#627F8A", border = NA)
    } else if (glyph == "bar") {
      rect(0.083, y - 0.010, 0.118, y + 0.010, col = track_spec[[i]]$color, border = NA)
      rect(0.118, y - 0.010, 0.132, y + 0.010, col = track_background, border = NA)
    } else if (glyph == "bubble") {
      draw_legend_circle(0.104, y, 0.013, "#B9473E")
    } else {
      segments(0.083, y, 0.132, y, col = "#C4CBCE", lwd = 1)
      draw_legend_circle(0.112, y, 0.008, track_spec[[i]]$color)
    }
    text(0.148, y, track_spec[[i]]$label, adj = c(0, 0.5), cex = 0.64, col = "#27343A")
    y <- y - 0.044
  }

  y <- 0.505
  text(0, y, "Panel summary", adj = c(0, 1), cex = 0.76, font = 2, col = "#1F2D33")
  y <- y - 0.039
  clade_labels <- c(Aves = "Birds", Mammalia_Chiroptera = "Bats", Mammalia_nonChiroptera = "Other mammals", Reptilia = "Reptiles")
  for (key in names(clade_labels)) {
    keep <- meta$clade == key
    n_clade <- sum(keep)
    median_lifespan <- median(meta$max_lifespan_years[keep], na.rm = TRUE)
    rect(0, y - 0.009, 0.040, y + 0.009, col = clade_colors[key], border = NA)
    text(0.055, y, sprintf("%-14s n=%d  median %.1f y", clade_labels[key], n_clade, median_lifespan), adj = c(0, 0.5), cex = 0.57, col = "#27343A")
    y <- y - 0.032
  }

  y <- 0.305
  text(0, y, "Visual encoding", adj = c(0, 1), cex = 0.76, font = 2, col = "#1F2D33")
  y <- y - 0.043
  draw_legend_circle(0.018, y, 0.009, "#2C6EAB")
  draw_legend_circle(0.050, y, 0.016, "#B9473E")
  text(0.082, y, "Bubble area = |value|; color = sign", adj = c(0, 0.5), cex = 0.56, col = "#3E4C53")
  y <- y - 0.043
  rect(0, y - 0.009, 0.052, y + 0.009, col = bar_coverage_color, border = NA)
  rect(0.052, y - 0.009, 0.075, y + 0.009, col = track_background, border = NA)
  text(0.092, y, "Bar length = scaled value", adj = c(0, 0.5), cex = 0.56, col = "#3E4C53")
  y <- y - 0.043
  segments(0, y, 0.075, y, col = "#BFC7CA", lwd = 1.2)
  draw_legend_circle(0.052, y, 0.008, dot_repeat_color)
  text(0.092, y, "Dot position = scaled value", adj = c(0, 0.5), cex = 0.56, col = "#3E4C53")
  y <- y - 0.045
  draw_gradient(0, y - 0.009, 0.28, y + 0.009, diverging_colors)
  text(0, y - 0.020, "negative", adj = c(0, 1), cex = 0.51, col = "#68767D")
  text(0.28, y - 0.020, "positive", adj = c(1, 1), cex = 0.51, col = "#68767D")
  rect(0.42, y - 0.009, 0.46, y + 0.009, col = tier_colors["Tier 1"], border = NA)
  rect(0.47, y - 0.009, 0.51, y + 0.009, col = tier_colors["Tier 2"], border = NA)
  text(0.53, y, "Tier 1 / Tier 2", adj = c(0, 0.5), cex = 0.51, col = "#68767D")
  y <- y - 0.058
  rect(0, y - 0.009, 0.040, y + 0.009, col = missing_color, border = "#AEB7BC", lwd = 0.5)
  text(0.055, y, "missing", adj = c(0, 0.5), cex = 0.53, col = "#68767D")

  sequence_x <- 0.56
  sequence_y <- 0.505
  text(sequence_x, sequence_y, "SAMHD1 residue heatmap", adj = c(0, 1), cex = 0.76, font = 2, col = "#1F2D33")
  sequence_y <- sequence_y - 0.038
  text(sequence_x, sequence_y, "626 aa in 24 windows | 56 sequences", adj = c(0, 1), cex = 0.52, col = "#68767D")
  sequence_y <- sequence_y - 0.036
  text(sequence_x, sequence_y, "Inner: aa 1-27  ->  outer: aa 601-626", adj = c(0, 1), cex = 0.44, col = "#68767D")
  sequence_y <- sequence_y - 0.032
  for (state in names(residue_colors)) {
    rect(sequence_x, sequence_y - 0.008, sequence_x + 0.032, sequence_y + 0.008, col = residue_colors[state], border = NA)
    text(sequence_x + 0.043, sequence_y, tools::toTitleCase(state), adj = c(0, 0.5), cex = 0.48, col = "#3E4C53")
    sequence_y <- sequence_y - 0.029
  }
  rect(sequence_x, sequence_y - 0.008, sequence_x + 0.032, sequence_y + 0.008, col = "#AAB4B9", border = "#C4CBCF", density = 22, angle = 45, lwd = 0.4)
  text(sequence_x + 0.043, sequence_y, "No sequence", adj = c(0, 0.5), cex = 0.48, col = "#3E4C53")
  sequence_y <- sequence_y - 0.038
  text(sequence_x, sequence_y, "Human-reference domains", adj = c(0, 0.5), cex = 0.49, font = 2, col = "#33434A")
  sequence_y <- sequence_y - 0.029
  rect(sequence_x, sequence_y - 0.008, sequence_x + 0.032, sequence_y + 0.008, col = domain_colors["SAM"], border = NA)
  text(sequence_x + 0.043, sequence_y, "SAM: 42-107 aa (rings 2-5)", adj = c(0, 0.5), cex = 0.44, col = "#3E4C53")
  sequence_y <- sequence_y - 0.027
  rect(sequence_x, sequence_y - 0.008, sequence_x + 0.032, sequence_y + 0.008, col = domain_colors["HD"], border = NA)
  text(sequence_x + 0.043, sequence_y, "HD: 164-227 aa (rings 7-9)", adj = c(0, 0.5), cex = 0.44, col = "#3E4C53")
  sequence_y <- sequence_y - 0.035
  rect(sequence_x, sequence_y - 0.006, sequence_x + 0.023, sequence_y + 0.006, col = "#477F87", border = NA)
  rect(sequence_x + 0.023, sequence_y - 0.006, sequence_x + 0.032, sequence_y + 0.006, col = "#25343B", border = NA)
  text(sequence_x + 0.043, sequence_y, "Non-gap coverage; black cap = qualified", adj = c(0, 0.5), cex = 0.41, col = "#68767D")

  sequence_y <- sequence_y - 0.042
  text(sequence_x, sequence_y, "Nine-gene domain fingerprint", adj = c(0, 0.5), cex = 0.43, font = 2, col = "#33434A")
  sequence_y <- sequence_y - 0.026
  draw_gradient(sequence_x, sequence_y - 0.006, sequence_x + 0.18, sequence_y + 0.006, fingerprint_palette)
  text(sequence_x + 0.195, sequence_y, "coverage x identity", adj = c(0, 0.5), cex = 0.37, col = "#68767D")
  sequence_y <- sequence_y - 0.024
  text(sequence_x, sequence_y, "Inner -> outer: DNMT1, DNMT3A, DNMT3B", adj = c(0, 0.5), cex = 0.37, col = "#68767D")
  sequence_y <- sequence_y - 0.020
  text(sequence_x, sequence_y, "HELLS, MBD2, MBD3, MORC3, SETDB2, UHRF1", adj = c(0, 0.5), cex = 0.37, col = "#68767D")

  text(0, 0.012, "Continuous values are descriptively scaled.\nOuter rings summarize amino-acid identity, substitution, and gaps; they are not a codon-selection analysis.\nTree: OpenTree + DateLife; four Grafen fallback nodes.", adj = c(0, 0), cex = 0.44, col = "#68767D")

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
  lifespan_residual_bubble_magnitude = scale_abs01(meta$lifespan_residual_log10),
  maximum_lifespan_bar_scaled = scale01(meta$log10_max_lifespan),
  busco_dot_position_scaled = scale01(meta$busco_complete),
  module_coverage_bar_scaled = scale01(meta$mean_module_coverage),
  shared_pc1_bubble_magnitude = scale_abs01(meta$shared_module_pc1),
  repeat_score_dot_position_scaled = scale01(meta$transposon_repeat_suppression_score),
  samhd1_non_gap_coverage = vapply(seq_len(Ntip(tree)), function(tip_index) {
    rows <- samhd1_window_data$tree_tip_order == tip_index & samhd1_window_data$sequence_available
    if (!any(rows)) return(NA_real_)
    mean(1 - samhd1_window_data$gap_fraction[rows])
  }, numeric(1)),
  samhd1_domain_analysis_qualified = vapply(seq_len(Ntip(tree)), function(tip_index) {
    rows <- samhd1_window_data$tree_tip_order == tip_index
    any(samhd1_window_data$domain_analysis_qualified[rows])
  }, logical(1)),
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

write.table(
  samhd1_window_data,
  file.path(table_dir, "jme_revision2_figure1_samhd1_residue_windows.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE,
  na = "NA"
)

write.table(
  fingerprint_data,
  file.path(table_dir, "jme_revision2_figure1_multigene_sequence_fingerprints.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE,
  na = "NA"
)

cat("Wrote circular phylogenetic heatmap, multigene fingerprints, SAMHD1 residue windows, and source data.\n")
