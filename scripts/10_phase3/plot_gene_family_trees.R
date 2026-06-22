#!/usr/bin/env Rscript

# Plot targeted gene-family trees from IQ-TREE outputs.
# Node labels show SH-aLRT/UFBoot support where support is reasonably high.

suppressPackageStartupMessages({
  project_libs <- c(file.path(getwd(), "env", "R_library"), file.path(getwd(), "env", "R", "library"))
  project_libs <- project_libs[dir.exists(project_libs)]
  if (length(project_libs) > 0) .libPaths(c(project_libs, .libPaths()))
  library(ape)
  library(readr)
  library(dplyr)
})

manifest_path <- "data/interim/phase3/gene_family_trees/gene_family_tree_sequence_manifest.tsv"
tree_dir <- "results/trees"
figure_dir <- "results/figures"
report_path <- "results/reports/gene_family_tree_validation_report.md"
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(report_path), recursive = TRUE, showWarnings = FALSE)

manifest <- read_tsv(manifest_path, show_col_types = FALSE)

gene_colors <- c(
  DNMT1 = "#0F766E",
  DNMT3A = "#4F46E5",
  DNMT3B = "#BE123C",
  MBD2 = "#B45309",
  MBD3 = "#15803D",
  HELLS = "#0E7490",
  UHRF1 = "#7C2D12",
  SETDB2 = "#6D28D9",
  SAMHD1 = "#0369A1"
)

abbr_species <- function(species) {
  species <- gsub("_", " ", species)
  parts <- strsplit(species, " ", fixed = TRUE)
  vapply(parts, function(p) {
    if (length(p) >= 2) paste0(substr(p[[1]], 1, 1), ". ", p[[2]]) else p[[1]]
  }, character(1))
}

format_support <- function(labels) {
  vapply(labels, function(x) {
    if (is.na(x) || x == "") return("")
    parts <- strsplit(x, "/", fixed = TRUE)[[1]]
    vals <- suppressWarnings(as.numeric(parts))
    if (length(vals) >= 2 && any(is.finite(vals))) {
      alrt <- vals[[1]]
      ufboot <- vals[[2]]
      if ((is.finite(alrt) && alrt >= 80) || (is.finite(ufboot) && ufboot >= 80)) {
        return(paste0(round(alrt), "/", round(ufboot)))
      }
      return("")
    }
    val <- suppressWarnings(as.numeric(x))
    if (is.finite(val) && val >= 80) return(as.character(round(val)))
    ""
  }, character(1))
}

extract_iqtree_line <- function(path, pattern) {
  if (!file.exists(path)) return(NA_character_)
  lines <- readLines(path, warn = FALSE)
  hit <- grep(pattern, lines, value = TRUE)
  if (length(hit) == 0) return(NA_character_)
  trimws(hit[[1]])
}

extract_convergence_note <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  lines <- readLines(path, warn = FALSE)
  final_warning <- grep("bootstrap analysis did not converge", lines, value = TRUE, ignore.case = TRUE)
  if (length(final_warning) > 0) return(trimws(final_warning[[length(final_warning)]]))
  correlations <- grep("Bootstrap correlation coefficient", lines, value = TRUE, ignore.case = TRUE)
  if (length(correlations) > 0) return(trimws(correlations[[length(correlations)]]))
  NA_character_
}

family_title <- function(family) {
  if (family == "DNMT_family") return("Supplementary DNMT family tree supporting DNMT paralog assignment")
  if (family == "MBD_family") return("Supplementary MBD family tree for MBD2/MBD3 ambiguity review")
  gene <- sub("_support_tree$", "", family)
  paste0("Supplementary ", gene, " support tree for targeted sequence validation")
}

plot_family <- function(family, title = family_title(family), width = 13.2, height = 8.8, prefix_suffix = "") {
  tree_path <- file.path(tree_dir, paste0(family, "_iqtree", prefix_suffix, ".treefile"))
  if (!file.exists(tree_path)) stop(paste("Missing tree:", tree_path))
  tr <- read.tree(tree_path)
  meta <- manifest |> filter(family == !!family)
  meta$pretty_label <- paste0(meta$gene, " | ", abbr_species(meta$species), " | ", meta$source_class)
  label_map <- setNames(meta$pretty_label, meta$output_id)
  gene_map <- setNames(meta$gene, meta$output_id)

  tip_genes <- gene_map[tr$tip.label]
  tip_colors <- gene_colors[tip_genes]
  tip_colors[is.na(tip_colors)] <- "#111827"
  tr$tip.label <- ifelse(tr$tip.label %in% names(label_map), label_map[tr$tip.label], tr$tip.label)
  node_support <- format_support(tr$node.label)

  fig_suffix <- ifelse(prefix_suffix == "", "", prefix_suffix)
  png_path <- file.path(figure_dir, paste0("figureS_gene_tree_", family, fig_suffix, ".png"))
  pdf_path <- file.path(figure_dir, paste0("figureS_gene_tree_", family, fig_suffix, ".pdf"))
  svg_path <- file.path(figure_dir, paste0("figureS_gene_tree_", family, fig_suffix, ".svg"))

  draw <- function() {
    par(mar = c(2.8, 1.2, 4.4, 1.2), xpd = TRUE)
    plot.phylo(
      tr,
      type = "phylogram",
      cex = ifelse(length(tr$tip.label) > 28, 0.55, 0.68),
      tip.color = tip_colors,
      label.offset = 0.01,
      no.margin = FALSE
    )
    nodelabels(
      text = node_support,
      frame = "none",
      cex = 0.48,
      col = "#475569",
      bg = "transparent",
      adj = c(1.1, -0.25)
    )
    title(main = title, adj = 0, cex.main = 1.0, font.main = 2, line = 1.4)
    mtext("Tip labels are colored by annotated gene. Node labels: SH-aLRT/UFBoot, shown when either support metric is >=80.", side = 1, adj = 0, line = 1.1, cex = 0.62, col = "#64748B")
  }

  png(png_path, width = width, height = height, units = "in", res = 320, bg = "white")
  draw()
  dev.off()
  pdf(pdf_path, width = width, height = height, bg = "white")
  draw()
  dev.off()
  svg(svg_path, width = width, height = height, bg = "white")
  draw()
  dev.off()

  iqtree_path <- file.path(tree_dir, paste0(family, "_iqtree", prefix_suffix, ".iqtree"))
  log_path <- file.path(tree_dir, paste0(family, "_iqtree", prefix_suffix, ".log"))
  data.frame(
    family = family,
    tree_variant = ifelse(prefix_suffix == "", "full_alignment", sub("^_", "", prefix_suffix)),
    sequences = length(tr$tip.label),
    genes = paste(sort(unique(meta$gene)), collapse = ", "),
    species = length(unique(meta$species)),
    best_model = {
      model_line <- extract_iqtree_line(iqtree_path, "Best-fit model:")
      if (is.na(model_line)) model_line <- extract_iqtree_line(log_path, "Best-fit model:")
      model_line
    },
    convergence_note = extract_convergence_note(log_path),
    png = png_path,
    pdf = pdf_path,
    svg = svg_path,
    treefile = tree_path
  )
}

families <- c("DNMT_family", "MBD_family", "HELLS_support_tree", "UHRF1_support_tree", "SETDB2_support_tree", "SAMHD1_support_tree")
families <- families[families %in% unique(manifest$family)]
summary_rows <- list()
for (family in families) {
  width <- ifelse(family %in% c("DNMT_family", "MBD_family"), 13.6, 10.8)
  height <- ifelse(family == "DNMT_family", 9.4, ifelse(family == "MBD_family", 8.2, 6.6))
  summary_rows[[length(summary_rows) + 1]] <- plot_family(family, width = width, height = height)
  trimmed_tree <- file.path(tree_dir, paste0(family, "_iqtree_trimmed_gappy70.treefile"))
  if (family %in% c("DNMT_family", "MBD_family") && file.exists(trimmed_tree)) {
    summary_rows[[length(summary_rows) + 1]] <- plot_family(
      family,
      title = paste0(family_title(family), " (trimmed alignment sensitivity)"),
      width = width,
      height = height,
      prefix_suffix = "_trimmed_gappy70"
    )
  }
}
summary <- bind_rows(summary_rows)

write_tsv(summary, "results/tables/gene_family_tree_validation_summary.tsv")

lines <- c(
  "# Gene-Family Tree Validation Report",
  "",
  "Targeted protein gene trees were generated for paralog-prone and high-priority transposon/chromatin module genes. These trees are intended as supplementary sequence-validation evidence, not as species-tree inference.",
  "",
  "## Methods",
  "",
  "- Sequences: phase3 assembly/GFF rescue candidates, phase3 UniProt rescue candidates, CDS-translation rescue candidates, rank-1 human references, and a manual UniProt human DNMT3B anchor.",
  "- Alignment: MAFFT `--auto` protein alignment; DNMT and MBD also include trimmed-alignment sensitivity trees where available.",
  "- Tree inference: IQ-TREE ModelFinder plus SH-aLRT 1000 and UFBoot 1000.",
  "- Visualization: R/ape phylogram; node labels are SH-aLRT/UFBoot and are shown when either metric is at least 80.",
  "",
  "## Outputs",
  "",
  apply(summary, 1, function(row) {
    paste0(
      "- ", row[["family"]], " / ", row[["tree_variant"]], ": ", row[["sequences"]], " sequences, ", row[["species"]],
      " species, genes=", row[["genes"]], "; tree=`", row[["treefile"]],
      "`, figure=`", row[["png"]], "`"
    )
  }),
  "",
  "## Caveats",
  "",
  "- DNMT includes many partial assembly/GFF fragments, so high gap fractions are expected; interpretation should focus on whether candidates fall in the expected DNMT1, DNMT3A, or DNMT3B neighborhood.",
  "- MBD proteins are short and domain-similar; a tree helps flag ambiguity but does not by itself convert short MBD-domain fragments into strict orthology evidence.",
  "- Bootstrap support is for these gene-family alignments only. It should not be reported as support for the OpenTree species topology.",
  "- Single-gene support trees for HELLS, UHRF1, SETDB2, and SAMHD1 are validation aids; they mainly check sequence clustering among target-labeled candidates rather than testing broad gene-family evolution."
)
writeLines(lines, report_path)

cat("Wrote gene-family tree figures, summary table, and report\n")
