# DateLife-tree PGLS tests for the targeted protein-conservation evidence layer.

project_lib <- file.path(getwd(), "env", "R_library")
if (dir.exists(project_lib)) {
  .libPaths(c(project_lib, .libPaths()))
}

library(ape)
library(caper)

args <- commandArgs(trailingOnly = TRUE)
sequence_path <- ifelse(length(args) >= 1, args[[1]], "results/tables/targeted_protein_conservation_rows.tsv")
trait_path <- ifelse(length(args) >= 2, args[[2]], "data/processed/maintenance_lifespan_phase2_W3_full_background_expanded.tsv")
tree_path <- ifelse(length(args) >= 3, args[[3]], "data/processed/phylogeny_inputs/opentree_datelife_calibrated_primary68.tre")
species_output <- ifelse(length(args) >= 4, args[[4]], "results/tables/targeted_protein_conservation_species.tsv")
model_output <- ifelse(length(args) >= 5, args[[5]], "results/tables/targeted_protein_conservation_pgls.tsv")
report_output <- ifelse(length(args) >= 6, args[[6]], "results/reports/targeted_protein_conservation_pgls_report.md")

sequence_rows <- read.delim(sequence_path, stringsAsFactors = FALSE, check.names = FALSE)
traits <- read.delim(trait_path, stringsAsFactors = FALSE, check.names = FALSE)
traits <- traits[traits$score_variant == "phase2_W3_full_background_sensitivity", ]
traits <- traits[!duplicated(traits$scientific_name), ]
tree_full <- read.tree(tree_path)

numeric_columns <- c(
  "human_reference_coverage", "aligned_identity", "identity_coverage_product",
  "length_completeness", "protein_length_ratio"
)
for (column in numeric_columns) {
  sequence_rows[[column]] <- as.numeric(sequence_rows[[column]])
}
sequence_rows$sequence_available <- as.logical(sequence_rows$sequence_available)
sequence_rows$qualified_sequence <- sequence_rows$sequence_available & sequence_rows$human_reference_coverage >= 0.5

zscore <- function(x) {
  if (sum(is.finite(x)) < 2 || sd(x, na.rm = TRUE) == 0) {
    return(rep(NA_real_, length(x)))
  }
  as.numeric(scale(x))
}

sequence_rows$identity_product_gene_z <- ave(
  sequence_rows$identity_coverage_product,
  sequence_rows$human_gene_symbol,
  FUN = zscore
)
sequence_rows$identity_product_gene_clade_z <- ave(
  sequence_rows$identity_coverage_product,
  interaction(sequence_rows$human_gene_symbol, sequence_rows$clade, drop = TRUE),
  FUN = zscore
)
sequence_rows$identity_product_gene_z[!sequence_rows$qualified_sequence] <- NA_real_
sequence_rows$identity_product_gene_clade_z[!sequence_rows$qualified_sequence] <- NA_real_

species_names <- unique(sequence_rows$scientific_name)
species_rows <- lapply(species_names, function(species) {
  x <- sequence_rows[sequence_rows$scientific_name == species, ]
  qualified <- x[x$qualified_sequence, ]
  data.frame(
    scientific_name = species,
    clade = x$clade[1],
    flight_status = x$flight_status[1],
    focal_genes_total = nrow(x),
    focal_sequences_n = sum(x$sequence_available, na.rm = TRUE),
    focal_sequence_fraction = mean(x$sequence_available, na.rm = TRUE),
    qualified_sequences_n = nrow(qualified),
    aggregate_identity_gene_z = mean(qualified$identity_product_gene_z, na.rm = TRUE),
    aggregate_identity_gene_clade_z = mean(qualified$identity_product_gene_clade_z, na.rm = TRUE),
    mean_identity_coverage_product = mean(qualified$identity_coverage_product, na.rm = TRUE),
    mean_length_completeness = mean(qualified$length_completeness, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
})
species_data <- do.call(rbind, species_rows)
species_data[!is.finite(species_data$aggregate_identity_gene_z), "aggregate_identity_gene_z"] <- NA_real_
species_data[!is.finite(species_data$aggregate_identity_gene_clade_z), "aggregate_identity_gene_clade_z"] <- NA_real_
species_data[!is.finite(species_data$mean_identity_coverage_product), "mean_identity_coverage_product"] <- NA_real_
species_data[!is.finite(species_data$mean_length_completeness), "mean_length_completeness"] <- NA_real_

trait_columns <- c(
  "scientific_name", "opentree_tip_label", "body_mass_g", "log10_body_mass_g",
  "max_lifespan_years", "pgls_model_c_mass_clade_residual", "genome_analysis_tier"
)
species_data <- merge(species_data, traits[, trait_columns], by = "scientific_name", all.x = TRUE)
if ("opentree_tip_label" %in% names(sequence_rows)) {
  sequence_rows$opentree_tip_label <- NULL
}
sequence_rows <- merge(
  sequence_rows,
  traits[, c("scientific_name", "opentree_tip_label", "pgls_model_c_mass_clade_residual")],
  by = "scientific_name",
  all.x = TRUE
)

prepare_comp <- function(data, predictor) {
  needed <- c("opentree_tip_label", "pgls_model_c_mass_clade_residual", predictor)
  data <- data[complete.cases(data[, needed]), ]
  data <- data[, needed]
  data <- data[!duplicated(data$opentree_tip_label), ]
  if (nrow(data) < 10 || sd(data[[predictor]]) == 0) {
    stop("insufficient usable observations or predictor variation")
  }
  data$predictor_z <- zscore(data[[predictor]])
  rownames(data) <- data$opentree_tip_label
  shared <- intersect(tree_full$tip.label, data$opentree_tip_label)
  pruned <- drop.tip(tree_full, setdiff(tree_full$tip.label, shared))
  data <- data[pruned$tip.label, ]
  comparative.data(
    phy = pruned,
    data = data,
    names.col = opentree_tip_label,
    vcv = TRUE,
    warn.dropped = FALSE
  )
}

fit_model <- function(data, predictor, test_family, scope, gene = "module_aggregate") {
  result <- tryCatch({
    comp <- prepare_comp(data, predictor)
    model <- pgls(pgls_model_c_mass_clade_residual ~ predictor_z, data = comp, lambda = "ML")
    term <- summary(model)$coefficients["predictor_z", ]
    data.frame(
      test_family = test_family,
      scope = scope,
      human_gene_symbol = gene,
      predictor = predictor,
      n = length(model$residuals),
      lambda = as.numeric(model$param["lambda"]),
      estimate_per_sd = unname(term["Estimate"]),
      se = unname(term["Std. Error"]),
      conf_low = unname(term["Estimate"] - 1.96 * term["Std. Error"]),
      conf_high = unname(term["Estimate"] + 1.96 * term["Std. Error"]),
      t = unname(term["t value"]),
      p = unname(term["Pr(>|t|)"]),
      error = "",
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    data.frame(
      test_family = test_family, scope = scope, human_gene_symbol = gene,
      predictor = predictor, n = NA, lambda = NA, estimate_per_sd = NA, se = NA,
      conf_low = NA, conf_high = NA, t = NA, p = NA,
      error = conditionMessage(e), stringsAsFactors = FALSE
    )
  })
  result
}

model_rows <- list()

# Predeclared gene-level tests use identity x human-reference coverage and require >=50% coverage.
for (gene in unique(sequence_rows$human_gene_symbol)) {
  gene_data <- sequence_rows[sequence_rows$human_gene_symbol == gene & sequence_rows$qualified_sequence, ]
  model_rows[[length(model_rows) + 1]] <- fit_model(
    gene_data, "identity_coverage_product", "gene_level_primary", "all_species", gene
  )
  model_rows[[length(model_rows) + 1]] <- fit_model(
    gene_data[gene_data$clade == "Aves", ], "identity_coverage_product",
    "gene_level_primary", "aves_only", gene
  )
}

aggregate_predictors <- c(
  "aggregate_identity_gene_z",
  "aggregate_identity_gene_clade_z",
  "mean_length_completeness"
)
aggregate_scopes <- list(
  all_species_min5 = species_data[species_data$qualified_sequences_n >= 5, ],
  all_species_min8 = species_data[species_data$qualified_sequences_n >= 8, ],
  aves_min5 = species_data[species_data$clade == "Aves" & species_data$qualified_sequences_n >= 5, ],
  aves_min8 = species_data[species_data$clade == "Aves" & species_data$qualified_sequences_n >= 8, ]
)
for (scope in names(aggregate_scopes)) {
  for (predictor in aggregate_predictors) {
    model_rows[[length(model_rows) + 1]] <- fit_model(
      aggregate_scopes[[scope]], predictor, "module_aggregate", scope
    )
  }
}

# Annotation-bias diagnostics: sequence count and fraction are tested on all 68 species.
for (predictor in c("focal_sequences_n", "focal_sequence_fraction")) {
  model_rows[[length(model_rows) + 1]] <- fit_model(
    species_data, predictor, "sequence_availability_diagnostic", "all_species"
  )
}

models <- do.call(rbind, model_rows)
models$q <- NA_real_
for (family_scope in unique(paste(models$test_family, models$scope, sep = "::"))) {
  index <- paste(models$test_family, models$scope, sep = "::") == family_scope & is.finite(models$p)
  models$q[index] <- p.adjust(models$p[index], method = "BH")
}
models <- models[order(models$test_family, models$scope, models$p), ]

dir.create(dirname(species_output), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(model_output), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(report_output), recursive = TRUE, showWarnings = FALSE)
write.table(species_data, species_output, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
write.table(models, model_output, sep = "\t", row.names = FALSE, quote = FALSE, na = "")

valid_models <- models[models$error == "", ]
gene_all <- valid_models[valid_models$test_family == "gene_level_primary" & valid_models$scope == "all_species", ]
gene_aves <- valid_models[valid_models$test_family == "gene_level_primary" & valid_models$scope == "aves_only", ]
aggregate <- valid_models[valid_models$test_family == "module_aggregate", ]
availability <- valid_models[valid_models$test_family == "sequence_availability_diagnostic", ]

format_best <- function(x) {
  if (nrow(x) == 0) return("No estimable model.")
  x <- x[order(x$p), ]
  paste0(
    x$human_gene_symbol[1], " / ", x$predictor[1], " / ", x$scope[1],
    ": beta per SD = ", signif(x$estimate_per_sd[1], 3),
    ", P = ", signif(x$p[1], 3), ", q = ", signif(x$q[1], 3),
    ", n = ", x$n[1], "."
  )
}

report <- c(
  "# Targeted Protein-Conservation PGLS",
  "",
  "Protein conservation was evaluated against the body-mass-and-clade-adjusted lifespan residual on the DateLife-calibrated 68-species tree. The primary per-gene metric was aligned amino-acid identity multiplied by human-reference coverage, restricted to sequences covering at least 50% of the human reference. Effects are reported per one standard deviation of the predictor.",
  "",
  paste0("- Species with at least five qualified focal proteins: ", sum(species_data$qualified_sequences_n >= 5), "/68."),
  paste0("- Species with at least eight qualified focal proteins: ", sum(species_data$qualified_sequences_n >= 8), "/68."),
  paste0("- Estimable all-species gene models: ", nrow(gene_all), "/10; BH q < 0.05: ", sum(gene_all$q < 0.05), "."),
  paste0("- Estimable bird-only gene models: ", nrow(gene_aves), "/10; BH q < 0.05: ", sum(gene_aves$q < 0.05), "."),
  paste0("- Best all-species gene result: ", format_best(gene_all)),
  paste0("- Best bird-only gene result: ", format_best(gene_aves)),
  paste0("- Best module-level result: ", format_best(aggregate)),
  "",
  "## Annotation-bias diagnostic",
  ""
)
if (nrow(availability) > 0) {
  for (i in seq_len(nrow(availability))) {
    report <- c(
      report,
      paste0(
        "- ", availability$predictor[i], ": beta per SD = ", signif(availability$estimate_per_sd[i], 3),
        ", P = ", signif(availability$p[i], 3), ", lambda = ", signif(availability$lambda[i], 3), "."
      )
    )
  }
}
report <- c(
  report,
  "",
  "These analyses test protein-sequence conservation, not biochemical activity or genomic repeat burden. A null result sets a useful claim boundary: the annotation-derived maintenance score should not be described as evidence of accelerated protein conservation."
)
writeLines(report, report_output)
cat("Wrote ", species_output, ", ", model_output, ", and ", report_output, "\n", sep = "")
