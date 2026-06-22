#!/usr/bin/env Rscript

# Build a literature-calibrated chronogram while retaining the OpenTree topology.
# Node ages are medians across peer-reviewed chronograms in the DateLife cache.
# Nodes without usable cross-child coverage retain a scaled Grafen fallback and
# are explicitly identified in the audit table.

suppressPackageStartupMessages({
  project_libs <- c(file.path(getwd(), "env", "R_library"), file.path(getwd(), "env", "R", "library"))
  project_libs <- project_libs[dir.exists(project_libs)]
  if (length(project_libs) > 0) .libPaths(c(project_libs, .libPaths()))
  library(ape)
})

args <- commandArgs(trailingOnly = TRUE)
target_tree_path <- ifelse(length(args) >= 1, args[[1]], "data/processed/phylogeny_inputs/opentree_grafen_scaled_320Ma_primary68.tre")
cache_path <- ifelse(length(args) >= 2, args[[2]], "data/raw/phylogeny/datelife_opentree_chronograms.rda")
output_tree_path <- ifelse(length(args) >= 3, args[[3]], "data/processed/phylogeny_inputs/opentree_datelife_calibrated_primary68.tre")
node_audit_path <- ifelse(length(args) >= 4, args[[4]], "results/tables/datelife_node_calibration_audit.tsv")
species_coverage_path <- ifelse(length(args) >= 5, args[[5]], "results/tables/datelife_species_coverage.tsv")
report_path <- ifelse(length(args) >= 6, args[[6]], "results/reports/datelife_calibrated_tree_report.md")

cache_url <- "https://raw.githubusercontent.com/phylotastic/datelife/master/data/opentree_chronograms.rda"

normalize_tip <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("'", "", x, fixed = TRUE)
  x <- gsub("[[:space:]]+", "_", x)
  x
}

clean_target_tip <- function(x) gsub("_ott[0-9]+$", "", normalize_tip(x))

if (!file.exists(cache_path)) {
  dir.create(dirname(cache_path), recursive = TRUE, showWarnings = FALSE)
  download.file(cache_url, cache_path, mode = "wb", quiet = FALSE)
}

loaded_names <- load(cache_path)
if (!"opentree_chronograms" %in% loaded_names) {
  stop("DateLife cache does not contain opentree_chronograms")
}

target <- read.tree(target_tree_path)
target <- ladderize(target)
n_tip <- length(target$tip.label)
target_clean <- clean_target_tip(target$tip.label)
if (anyDuplicated(target_clean)) stop("Target tree has duplicate normalized tip labels")

all_source_tips <- unique(unlist(lapply(opentree_chronograms$trees, function(tr) normalize_tip(tr$tip.label))))
species_coverage <- data.frame(
  opentree_tip_label = target$tip.label,
  scientific_name = gsub("_", " ", target_clean),
  datelife_exact_tip_match = target_clean %in% all_source_tips,
  stringsAsFactors = FALSE
)

source_records <- list()
source_status <- data.frame(
  source_tree_id = seq_along(opentree_chronograms$trees),
  status = "not_evaluated",
  matched_target_tips = 0L,
  stringsAsFactors = FALSE
)

for (i in seq_along(opentree_chronograms$trees)) {
  tr <- opentree_chronograms$trees[[i]]
  tr$tip.label <- normalize_tip(tr$tip.label)
  duplicated_tips <- which(duplicated(tr$tip.label))
  if (length(duplicated_tips) > 0) tr <- drop.tip(tr, duplicated_tips)
  shared <- intersect(target_clean, tr$tip.label)
  source_status$matched_target_tips[[i]] <- length(shared)
  if (length(shared) < 2) {
    source_status$status[[i]] <- "fewer_than_two_target_tips"
    next
  }
  if (is.null(tr$edge.length) || any(!is.finite(tr$edge.length)) || any(tr$edge.length < 0)) {
    source_status$status[[i]] <- "missing_or_invalid_branch_lengths"
    next
  }
  tip_depths <- node.depth.edgelength(tr)[seq_along(tr$tip.label)]
  root_age <- max(tip_depths)
  if (!is.finite(root_age) || root_age <= 0) {
    source_status$status[[i]] <- "invalid_root_age"
    next
  }
  if ((max(tip_depths) - min(tip_depths)) > max(1e-4, 0.01 * root_age)) {
    source_status$status[[i]] <- "non_ultrametric_source"
    next
  }
  dist_matrix <- cophenetic.phylo(tr)[shared, shared, drop = FALSE] / 2
  source_records[[length(source_records) + 1]] <- list(
    source_tree_id = i,
    distances = dist_matrix,
    study = if (length(opentree_chronograms$studies) >= i) as.character(opentree_chronograms$studies[[i]]) else ""
  )
  source_status$status[[i]] <- "used"
}

if (length(source_records) == 0) stop("No DateLife chronograms passed the source-tree filters")

children_of <- function(node) target$edge[target$edge[, 1] == node, 2]

descendant_tips <- local({
  cache <- new.env(parent = emptyenv())
  function(node) {
    key <- as.character(node)
    if (exists(key, envir = cache, inherits = FALSE)) return(get(key, envir = cache, inherits = FALSE))
    if (node <= n_tip) {
      out <- node
    } else {
      out <- unlist(lapply(children_of(node), descendant_tips), use.names = FALSE)
    }
    assign(key, out, envir = cache)
    out
  }
})

internal_nodes <- (n_tip + 1):(n_tip + target$Nnode)
node_rows <- vector("list", length(internal_nodes))

for (j in seq_along(internal_nodes)) {
  node <- internal_nodes[[j]]
  child_nodes <- children_of(node)
  child_tip_names <- lapply(child_nodes, function(child) target_clean[descendant_tips(child)])
  estimates <- numeric(0)
  source_ids <- integer(0)
  for (record in source_records) {
    available_by_child <- lapply(child_tip_names, function(x) intersect(x, rownames(record$distances)))
    present_children <- which(lengths(available_by_child) > 0)
    if (length(present_children) < 2) next
    child_pairs <- combn(present_children, 2, simplify = FALSE)
    cross_values <- unlist(lapply(child_pairs, function(pair) {
      as.numeric(record$distances[available_by_child[[pair[[1]]]], available_by_child[[pair[[2]]]], drop = FALSE])
    }), use.names = FALSE)
    cross_values <- cross_values[is.finite(cross_values) & cross_values > 0]
    if (length(cross_values) == 0) next
    estimates <- c(estimates, median(cross_values))
    source_ids <- c(source_ids, record$source_tree_id)
  }
  node_rows[[j]] <- data.frame(
    node = node,
    descendant_tip_count = length(descendant_tips(node)),
    child_count = length(child_nodes),
    raw_median_age_ma = if (length(estimates) > 0) median(estimates) else NA_real_,
    raw_q25_age_ma = if (length(estimates) > 0) unname(quantile(estimates, 0.25)) else NA_real_,
    raw_q75_age_ma = if (length(estimates) > 0) unname(quantile(estimates, 0.75)) else NA_real_,
    source_tree_count = length(estimates),
    source_tree_ids = paste(source_ids, collapse = ";"),
    stringsAsFactors = FALSE
  )
}

node_audit <- do.call(rbind, node_rows)
root_node <- n_tip + 1
root_raw <- node_audit$raw_median_age_ma[node_audit$node == root_node]
root_age <- ifelse(length(root_raw) == 1 && is.finite(root_raw), root_raw, 320)

grafen <- compute.brlen(target, method = "Grafen", power = 1)
grafen_depth <- node.depth.edgelength(grafen)
grafen_root_depth <- max(grafen_depth[seq_len(n_tip)])
grafen$edge.length <- grafen$edge.length * (root_age / grafen_root_depth)
grafen_depth <- node.depth.edgelength(grafen)
grafen_node_age <- root_age - grafen_depth

node_age <- rep(NA_real_, n_tip + target$Nnode)
node_age[seq_len(n_tip)] <- 0
raw_by_node <- setNames(node_audit$raw_median_age_ma, node_audit$node)
for (node in internal_nodes) {
  raw <- raw_by_node[[as.character(node)]]
  node_age[[node]] <- if (is.finite(raw)) raw else grafen_node_age[[node]]
}

# Enforce time consistency from younger to older nodes. The adjustment is
# deliberately minimal and fully recorded rather than silently smoothing ages.
desc_counts <- vapply(internal_nodes, function(node) length(descendant_tips(node)), integer(1))
for (node in internal_nodes[order(desc_counts)]) {
  child_age <- node_age[children_of(node)]
  required_age <- max(child_age) + 1e-4
  if (!is.finite(node_age[[node]]) || node_age[[node]] < required_age) node_age[[node]] <- required_age
}

output_tree <- target
output_tree$edge.length <- node_age[output_tree$edge[, 1]] - node_age[output_tree$edge[, 2]]
if (any(!is.finite(output_tree$edge.length)) || any(output_tree$edge.length <= 0)) {
  stop("Calibrated tree contains non-positive or invalid branch lengths")
}

node_audit$grafen_fallback_age_ma <- grafen_node_age[node_audit$node]
node_audit$final_age_ma <- node_age[node_audit$node]
node_audit$age_method <- ifelse(
  is.finite(node_audit$raw_median_age_ma),
  ifelse(abs(node_audit$final_age_ma - node_audit$raw_median_age_ma) > 1e-8, "datelife_median_monotonic_adjusted", "datelife_median"),
  ifelse(abs(node_audit$final_age_ma - node_audit$grafen_fallback_age_ma) > 1e-8, "grafen_fallback_monotonic_adjusted", "grafen_fallback")
)
node_audit$absolute_adjustment_ma <- ifelse(
  is.finite(node_audit$raw_median_age_ma),
  abs(node_audit$final_age_ma - node_audit$raw_median_age_ma),
  abs(node_audit$final_age_ma - node_audit$grafen_fallback_age_ma)
)

dir.create(dirname(output_tree_path), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(node_audit_path), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(species_coverage_path), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(report_path), recursive = TRUE, showWarnings = FALSE)
write.tree(output_tree, file = output_tree_path)
write.table(node_audit, node_audit_path, sep = "\t", row.names = FALSE, quote = FALSE)
write.table(species_coverage, species_coverage_path, sep = "\t", row.names = FALSE, quote = FALSE)

calibrated_nodes <- sum(grepl("^datelife", node_audit$age_method))
fallback_nodes <- sum(grepl("^grafen", node_audit$age_method))
adjusted_nodes <- sum(node_audit$absolute_adjustment_ma > 1e-8)
cache_update <- format(opentree_chronograms$update, "%Y-%m-%d %H:%M:%S %Z")
cache_version <- as.character(opentree_chronograms$version)

report_lines <- c(
  "# DateLife-Calibrated OpenTree Chronogram",
  "",
  "The target topology is the 68-species OpenTree subtree. Internal-node ages are medians across usable peer-reviewed chronograms in the DateLife OpenTree cache. Nodes without source coverage use a Grafen fallback scaled to the DateLife-derived root age.",
  "",
  "## Coverage",
  "",
  paste0("- Target species: ", n_tip),
  paste0("- Exact species matches in the DateLife cache: ", sum(species_coverage$datelife_exact_tip_match), "/", n_tip),
  paste0("- Cache chronograms: ", length(opentree_chronograms$trees)),
  paste0("- Chronograms passing overlap, branch-length, and ultrametric filters: ", length(source_records)),
  paste0("- Internal nodes with DateLife median ages: ", calibrated_nodes, "/", target$Nnode),
  paste0("- Internal nodes using explicit Grafen fallback: ", fallback_nodes, "/", target$Nnode),
  paste0("- Nodes minimally adjusted to enforce temporal ordering: ", adjusted_nodes),
  paste0("- Final root age: ", signif(node_age[[root_node]], 6), " Ma"),
  paste0("- Output tree ultrametric: ", is.ultrametric(output_tree, tol = 1e-6)),
  "",
  "## Provenance",
  "",
  paste0("- DateLife cache URL: `", cache_url, "`"),
  paste0("- Cache update timestamp: ", cache_update),
  paste0("- DateLife package cache version: ", cache_version),
  "- Node-level source-tree identifiers, age quartiles, fallback ages, and adjustments are reported in the audit table.",
  "",
  "## Interpretation Guardrail",
  "",
  "This is a literature-calibrated secondary chronogram on a fixed OpenTree topology. It is stronger than a globally scaled Grafen visualization, but it is not a newly inferred sequence-based species tree and it inherits taxonomic and source-tree uncertainty."
)
writeLines(report_lines, report_path)

cat("Wrote ", output_tree_path, ", ", node_audit_path, ", ", species_coverage_path, " and ", report_path, "\n", sep = "")
