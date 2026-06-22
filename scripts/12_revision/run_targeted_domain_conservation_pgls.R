# DateLife-tree PGLS for Pfam-domain-aware conservation metrics.

project_lib <- file.path(getwd(), "env", "R_library")
if (dir.exists(project_lib)) .libPaths(c(project_lib, .libPaths()))

library(ape)
library(caper)

args <- commandArgs(trailingOnly = TRUE)
domain_path <- ifelse(length(args) >= 1, args[[1]], "results/tables/targeted_domain_conservation_rows.tsv")
trait_path <- ifelse(length(args) >= 2, args[[2]], "data/processed/maintenance_lifespan_phase3_gff_cds_plus_uniprot_sequence_rescued.tsv")
tree_path <- ifelse(length(args) >= 3, args[[3]], "data/processed/phylogeny_inputs/opentree_datelife_calibrated_primary68.tre")
species_output <- ifelse(length(args) >= 4, args[[4]], "results/tables/targeted_domain_conservation_species.tsv")
model_output <- ifelse(length(args) >= 5, args[[5]], "results/tables/targeted_domain_conservation_pgls.tsv")
report_output <- ifelse(length(args) >= 6, args[[6]], "results/reports/targeted_domain_conservation_pgls.md")

rows <- read.delim(domain_path, stringsAsFactors = FALSE, check.names = FALSE)
traits <- read.delim(trait_path, stringsAsFactors = FALSE, check.names = FALSE)
traits <- traits[traits$score_variant == "phase2_W3_full_background_sensitivity", ]
traits <- traits[!duplicated(traits$scientific_name), ]
tree_full <- read.tree(tree_path)

metric_columns <- c(
  "domain_reference_coverage", "domain_aligned_identity", "domain_identity_coverage_product",
  "nondomain_reference_coverage", "nondomain_aligned_identity", "nondomain_identity_coverage_product",
  "domain_minus_nondomain_identity", "domain_minus_nondomain_product"
)
for (column in metric_columns) rows[[column]] <- as.numeric(rows[[column]])
rows$sequence_available <- tolower(as.character(rows$sequence_available)) %in% c("true", "1", "yes")
rows$qualified_domain_sequence <- rows$sequence_available & rows$domain_reference_coverage >= 0.5

zscore <- function(x) {
  if (sum(is.finite(x)) < 2 || sd(x, na.rm = TRUE) == 0) return(rep(NA_real_, length(x)))
  as.numeric(scale(x))
}
rows$domain_product_gene_z <- ave(rows$domain_identity_coverage_product, rows$human_gene_symbol, FUN = zscore)
rows$domain_product_gene_clade_z <- ave(
  rows$domain_identity_coverage_product,
  interaction(rows$human_gene_symbol, rows$clade, drop = TRUE),
  FUN = zscore
)
rows$domain_contrast_gene_clade_z <- ave(
  rows$domain_minus_nondomain_product,
  interaction(rows$human_gene_symbol, rows$clade, drop = TRUE),
  FUN = zscore
)
rows$domain_product_gene_z[!rows$qualified_domain_sequence] <- NA_real_
rows$domain_product_gene_clade_z[!rows$qualified_domain_sequence] <- NA_real_
rows$domain_contrast_gene_clade_z[!rows$qualified_domain_sequence] <- NA_real_

species_rows <- lapply(unique(rows$scientific_name), function(species) {
  x <- rows[rows$scientific_name == species, ]
  q <- x[x$qualified_domain_sequence, ]
  data.frame(
    scientific_name = species,
    clade = x$clade[1],
    qualified_domain_genes_n = nrow(q),
    aggregate_domain_gene_z = mean(q$domain_product_gene_z, na.rm = TRUE),
    aggregate_domain_gene_clade_z = mean(q$domain_product_gene_clade_z, na.rm = TRUE),
    aggregate_domain_contrast_clade_z = mean(q$domain_contrast_gene_clade_z, na.rm = TRUE),
    mean_domain_identity = mean(q$domain_aligned_identity, na.rm = TRUE),
    mean_domain_coverage = mean(q$domain_reference_coverage, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
})
species <- do.call(rbind, species_rows)
numeric_species <- setdiff(names(species), c("scientific_name", "clade"))
for (column in numeric_species) species[[column]][!is.finite(species[[column]])] <- NA_real_
trait_columns <- c("scientific_name", "opentree_tip_label", "pgls_model_c_mass_clade_residual", "log10_body_mass_g", "max_lifespan_years")
species <- merge(species, traits[, trait_columns], by = "scientific_name", all.x = TRUE)
if ("opentree_tip_label" %in% names(rows)) rows$opentree_tip_label <- NULL
rows <- merge(rows, traits[, c("scientific_name", "opentree_tip_label", "pgls_model_c_mass_clade_residual")], by = "scientific_name", all.x = TRUE)

prepare_comp <- function(input, predictor) {
  needed <- c("opentree_tip_label", "pgls_model_c_mass_clade_residual", predictor)
  x <- input[complete.cases(input[, needed]), needed]
  x <- x[!duplicated(x$opentree_tip_label), ]
  if (nrow(x) < 8 || sd(x[[predictor]]) == 0) stop("insufficient observations or predictor variation")
  x$predictor_z <- zscore(x[[predictor]])
  rownames(x) <- x$opentree_tip_label
  shared <- intersect(tree_full$tip.label, x$opentree_tip_label)
  phy <- drop.tip(tree_full, setdiff(tree_full$tip.label, shared))
  x <- x[phy$tip.label, ]
  comparative.data(phy = phy, data = x, names.col = opentree_tip_label, vcv = TRUE, warn.dropped = FALSE)
}

fit_one <- function(input, predictor, family, scope, gene = "domain_aggregate") {
  tryCatch({
    comp <- prepare_comp(input, predictor)
    formula <- pgls_model_c_mass_clade_residual ~ predictor_z
    ml_attempt <- tryCatch(pgls(formula, data = comp, lambda = "ML"), error = function(e) e)
    if (inherits(ml_attempt, "error")) {
      model <- pgls(formula, data = comp, lambda = 1e-6)
      lambda_method <- "fixed_near_zero_optimizer_fallback"
      fit_note <- conditionMessage(ml_attempt)
    } else {
      model <- ml_attempt
      lambda_method <- "maximum_likelihood"
      fit_note <- ""
    }
    term <- summary(model)$coefficients["predictor_z", ]
    data.frame(
      test_family = family, scope = scope, human_gene_symbol = gene, predictor = predictor,
      n = nrow(comp$data), lambda = as.numeric(model$param["lambda"]), lambda_method = lambda_method,
      estimate_per_sd = unname(term["Estimate"]), se = unname(term["Std. Error"]),
      conf_low = unname(term["Estimate"] - 1.96 * term["Std. Error"]),
      conf_high = unname(term["Estimate"] + 1.96 * term["Std. Error"]),
      t = unname(term["t value"]), p = unname(term["Pr(>|t|)"]), fit_note = fit_note, error = "",
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    data.frame(
      test_family = family, scope = scope, human_gene_symbol = gene, predictor = predictor,
      n = nrow(input), lambda = NA, lambda_method = "failed", estimate_per_sd = NA, se = NA, conf_low = NA,
      conf_high = NA, t = NA, p = NA, fit_note = "", error = conditionMessage(e), stringsAsFactors = FALSE
    )
  })
}

model_rows <- list()
for (gene in unique(rows$human_gene_symbol)) {
  gene_rows <- rows[rows$human_gene_symbol == gene & rows$qualified_domain_sequence, ]
  for (predictor in c("domain_identity_coverage_product", "domain_minus_nondomain_product")) {
    model_rows[[length(model_rows) + 1]] <- fit_one(gene_rows, predictor, "gene_level_domain", "all_species", gene)
    model_rows[[length(model_rows) + 1]] <- fit_one(gene_rows[gene_rows$clade == "Aves", ], predictor, "gene_level_domain", "aves_only", gene)
  }
}

aggregate_scopes <- list(
  all_species_min5 = species[species$qualified_domain_genes_n >= 5, ],
  all_species_min8 = species[species$qualified_domain_genes_n >= 8, ],
  aves_min5 = species[species$clade == "Aves" & species$qualified_domain_genes_n >= 5, ],
  aves_min8 = species[species$clade == "Aves" & species$qualified_domain_genes_n >= 8, ]
)
for (scope in names(aggregate_scopes)) {
  for (predictor in c("aggregate_domain_gene_z", "aggregate_domain_gene_clade_z", "aggregate_domain_contrast_clade_z")) {
    model_rows[[length(model_rows) + 1]] <- fit_one(aggregate_scopes[[scope]], predictor, "domain_aggregate", scope)
  }
}
models <- do.call(rbind, model_rows)
models$q <- NA_real_
for (family_scope_predictor in unique(paste(models$test_family, models$scope, models$predictor, sep = "::"))) {
  index <- paste(models$test_family, models$scope, models$predictor, sep = "::") == family_scope_predictor & is.finite(models$p)
  models$q[index] <- p.adjust(models$p[index], method = "BH")
}
models$q_gene_scope_all_metrics <- NA_real_
for (scope in unique(models$scope[models$test_family == "gene_level_domain"])) {
  index <- models$test_family == "gene_level_domain" & models$scope == scope & is.finite(models$p)
  models$q_gene_scope_all_metrics[index] <- p.adjust(models$p[index], method = "BH")
}
models <- models[order(models$test_family, models$scope, models$predictor, models$p), ]

dir.create(dirname(species_output), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(model_output), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(report_output), recursive = TRUE, showWarnings = FALSE)
write.table(species, species_output, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
write.table(models, model_output, sep = "\t", row.names = FALSE, quote = FALSE, na = "")

valid <- models[models$error == "", ]
best <- valid[order(valid$p), ][1, ]
report <- c(
  "# Pfam-Domain Conservation PGLS",
  "",
  "Human Pfam domain boundaries were projected onto per-gene MAFFT alignments. The primary domain metric combines amino-acid identity with coverage of the human domain residues. Domain-minus-nondomain contrasts test whether any lifespan relationship is concentrated within annotated functional regions.",
  "",
  paste0("- Species with at least five qualified domain-scored genes: ", sum(species$qualified_domain_genes_n >= 5), "/68."),
  paste0("- Species with at least eight qualified domain-scored genes: ", sum(species$qualified_domain_genes_n >= 8), "/68."),
  paste0("- Estimable models: ", nrow(valid), "/", nrow(models), "; predictor-family BH q < 0.05: ", sum(valid$q < 0.05), "."),
  paste0("- Smallest nominal P: ", best$test_family, " / ", best$scope, " / ", best$human_gene_symbol, " / ", best$predictor, ": beta per SD = ", signif(best$estimate_per_sd, 3), ", P = ", signif(best$p, 3), ", predictor-family q = ", signif(best$q, 3), ", joint gene-metric q = ", signif(best$q_gene_scope_all_metrics, 3), ", n = ", best$n, "."),
  "",
  "Gene-level FDR support is treated as targeted exploratory evidence and requires gene-specific influence and clade sensitivity checks. Aggregate null results mean the annotation-derived module score is not broadly corroborated by stronger conservation across all 10 Pfam-defined protein architectures. The analysis does not test domain biochemical activity or regulatory evolution."
)
writeLines(report, report_output)
cat("Wrote ", species_output, ", ", model_output, ", and ", report_output, "\n", sep = "")
