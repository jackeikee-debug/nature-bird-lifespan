# Gene-specific robustness checks for the exploratory SAMHD1 Pfam-domain signal.

project_lib <- file.path(getwd(), "env", "R_library")
if (dir.exists(project_lib)) .libPaths(c(project_lib, .libPaths()))

library(ape)
library(caper)

args <- commandArgs(trailingOnly = TRUE)
domain_path <- ifelse(length(args) >= 1, args[[1]], "results/tables/targeted_domain_conservation_rows.tsv")
protein_path <- ifelse(length(args) >= 2, args[[2]], "results/tables/targeted_protein_conservation_rows.tsv")
trait_path <- ifelse(length(args) >= 3, args[[3]], "data/processed/maintenance_lifespan_phase3_gff_cds_plus_uniprot_sequence_rescued.tsv")
tree_path <- ifelse(length(args) >= 4, args[[4]], "data/processed/phylogeny_inputs/opentree_datelife_calibrated_primary68.tre")
model_output <- ifelse(length(args) >= 5, args[[5]], "results/tables/samhd1_domain_robustness_pgls.tsv")
loo_output <- ifelse(length(args) >= 6, args[[6]], "results/tables/samhd1_domain_leave_one_species_out.tsv")
report_output <- ifelse(length(args) >= 7, args[[7]], "results/reports/samhd1_domain_robustness.md")

domain <- read.delim(domain_path, stringsAsFactors = FALSE, check.names = FALSE)
protein <- read.delim(protein_path, stringsAsFactors = FALSE, check.names = FALSE)
traits <- read.delim(trait_path, stringsAsFactors = FALSE, check.names = FALSE)
traits <- traits[traits$score_variant == "phase2_W3_full_background_sensitivity", ]
traits <- traits[!duplicated(traits$scientific_name), ]
tree_full <- read.tree(tree_path)

domain <- domain[domain$human_gene_symbol == "SAMHD1", ]
protein <- protein[protein$human_gene_symbol == "SAMHD1", c("scientific_name", "identity_coverage_product")]
names(protein)[2] <- "whole_protein_identity_coverage_product"
data <- merge(domain, protein, by = "scientific_name", all.x = TRUE)
data <- merge(
  data,
  traits[, c("scientific_name", "opentree_tip_label", "pgls_model_c_mass_clade_residual")],
  by = "scientific_name",
  all.x = TRUE,
  suffixes = c("", "_trait")
)
if ("opentree_tip_label_trait" %in% names(data)) data$opentree_tip_label <- data$opentree_tip_label_trait

numeric_columns <- c(
  "domain_reference_coverage", "domain_aligned_identity", "domain_identity_coverage_product",
  "nondomain_identity_coverage_product", "domain_minus_nondomain_product",
  "whole_protein_identity_coverage_product", "pgls_model_c_mass_clade_residual"
)
for (column in numeric_columns) data[[column]] <- as.numeric(data[[column]])
data$sequence_available <- tolower(as.character(data$sequence_available)) %in% c("true", "1", "yes")
data <- data[data$sequence_available & data$domain_reference_coverage >= 0.5, ]
data$clade <- factor(data$clade)

scope_definitions <- list(
  all_species = function(x) rep(TRUE, nrow(x)),
  domain_coverage_ge_0.8 = function(x) x$domain_reference_coverage >= 0.8,
  domain_coverage_ge_0.9 = function(x) x$domain_reference_coverage >= 0.9,
  birds_only = function(x) x$clade == "Aves",
  nonbirds = function(x) x$clade != "Aves",
  no_birds = function(x) x$clade != "Aves",
  no_bats = function(x) x$clade != "Mammalia_Chiroptera",
  no_nonflying_mammals = function(x) x$clade != "Mammalia_nonChiroptera",
  no_reptiles = function(x) x$clade != "Reptilia"
)

prepare_comp <- function(input, predictor, fixed_center = NULL, fixed_scale = NULL) {
  needed <- c("opentree_tip_label", "scientific_name", "clade", "pgls_model_c_mass_clade_residual", predictor)
  x <- input[complete.cases(input[, needed]), needed]
  x <- x[!duplicated(x$opentree_tip_label), ]
  center <- ifelse(is.null(fixed_center), mean(x[[predictor]]), fixed_center)
  scale_value <- ifelse(is.null(fixed_scale), sd(x[[predictor]]), fixed_scale)
  if (nrow(x) < 8 || !is.finite(scale_value) || scale_value == 0) stop("insufficient observations or variation")
  x$predictor_z <- (x[[predictor]] - center) / scale_value
  rownames(x) <- x$opentree_tip_label
  shared <- intersect(tree_full$tip.label, x$opentree_tip_label)
  phy <- drop.tip(tree_full, setdiff(tree_full$tip.label, shared))
  x <- x[phy$tip.label, ]
  list(
    comp = comparative.data(phy = phy, data = x, names.col = opentree_tip_label, vcv = TRUE, warn.dropped = FALSE),
    center = center,
    scale = scale_value
  )
}

fit_one <- function(input, scope, predictor, center = NULL, scale_value = NULL) {
  tryCatch({
    prepared <- prepare_comp(input, predictor, center, scale_value)
    formula <- pgls_model_c_mass_clade_residual ~ predictor_z
    ml_attempt <- tryCatch(pgls(formula, data = prepared$comp, lambda = "ML"), error = function(e) e)
    if (inherits(ml_attempt, "error")) {
      model <- pgls(formula, data = prepared$comp, lambda = 1e-6)
      lambda_method <- "fixed_near_zero_optimizer_fallback"
      fit_note <- conditionMessage(ml_attempt)
    } else {
      model <- ml_attempt
      lambda_method <- "maximum_likelihood"
      fit_note <- ""
    }
    term <- summary(model)$coefficients["predictor_z", ]
    data.frame(
      scope = scope, predictor = predictor, n = nrow(prepared$comp$data),
      predictor_center = prepared$center, predictor_scale = prepared$scale,
      lambda = as.numeric(model$param["lambda"]), lambda_method = lambda_method,
      estimate_per_sd = unname(term["Estimate"]), se = unname(term["Std. Error"]),
      conf_low = unname(term["Estimate"] - 1.96 * term["Std. Error"]),
      conf_high = unname(term["Estimate"] + 1.96 * term["Std. Error"]),
      t = unname(term["t value"]), p = unname(term["Pr(>|t|)"]), fit_note = fit_note, error = "",
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    data.frame(
      scope = scope, predictor = predictor, n = nrow(input), predictor_center = NA,
      predictor_scale = NA, lambda = NA, lambda_method = "failed", estimate_per_sd = NA, se = NA,
      conf_low = NA, conf_high = NA, t = NA, p = NA, fit_note = "",
      error = conditionMessage(e), stringsAsFactors = FALSE
    )
  })
}

predictors <- c(
  "domain_identity_coverage_product",
  "domain_minus_nondomain_product",
  "nondomain_identity_coverage_product",
  "whole_protein_identity_coverage_product"
)
model_rows <- list()
for (scope in names(scope_definitions)) {
  subset <- data[scope_definitions[[scope]](data), ]
  for (predictor in predictors) {
    model_rows[[length(model_rows) + 1]] <- fit_one(subset, scope, predictor)
  }
}
models <- do.call(rbind, model_rows)
models$q_within_scope <- NA_real_
for (scope in unique(models$scope)) {
  index <- models$scope == scope & is.finite(models$p)
  models$q_within_scope[index] <- p.adjust(models$p[index], method = "BH")
}
models <- models[order(models$scope, models$p), ]

loo_rows <- list()
for (predictor in c("domain_identity_coverage_product", "domain_minus_nondomain_product")) {
  baseline <- fit_one(data, "all_species", predictor)
  for (species in data$scientific_name) {
    loo <- fit_one(
      data[data$scientific_name != species, ], "leave_one_species_out", predictor,
      baseline$predictor_center, baseline$predictor_scale
    )
    loo$omitted_species <- species
    loo$baseline_n <- baseline$n
    loo$baseline_estimate <- baseline$estimate_per_sd
    loo$baseline_p <- baseline$p
    loo$delta_estimate <- loo$estimate_per_sd - baseline$estimate_per_sd
    loo$direction_preserved <- sign(loo$estimate_per_sd) == sign(baseline$estimate_per_sd)
    loo$nominal_support_preserved <- loo$p < 0.05
    loo_rows[[length(loo_rows) + 1]] <- loo
  }
}
loo <- do.call(rbind, loo_rows)
loo <- loo[order(loo$predictor, -abs(loo$delta_estimate)), ]

dir.create(dirname(model_output), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(loo_output), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(report_output), recursive = TRUE, showWarnings = FALSE)
write.table(models, model_output, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
write.table(loo, loo_output, sep = "\t", row.names = FALSE, quote = FALSE, na = "")

fmt <- function(x) ifelse(is.finite(x), format(signif(x, 4), scientific = FALSE), "NA")
report <- c(
  "# SAMHD1 Pfam-Domain Robustness",
  "",
  "The exploratory SAMHD1 signal was compared with non-domain and whole-protein metrics, high-domain-coverage subsets, taxonomic exclusions, and leave-one-species refits on the DateLife-calibrated tree.",
  "",
  "## Scope models",
  ""
)
for (scope in names(scope_definitions)) {
  subset <- models[models$scope == scope, ]
  for (i in seq_len(nrow(subset))) {
    row <- subset[i, ]
    report <- c(report, paste0(
      "- ", scope, " / ", row$predictor, ": beta per SD = ", fmt(row$estimate_per_sd),
      " (95% CI ", fmt(row$conf_low), " to ", fmt(row$conf_high),
      "), P = ", fmt(row$p), ", q-within-scope = ", fmt(row$q_within_scope), ", n = ", row$n, "."
    ))
  }
}
report <- c(report, "", "## Leave-one-species stability", "")
for (predictor in unique(loo$predictor)) {
  x <- loo[loo$predictor == predictor & loo$error == "", ]
  top <- x[which.max(abs(x$delta_estimate)), ]
  report <- c(report, paste0(
    "- ", predictor, ": estimate range ", fmt(min(x$estimate_per_sd)), " to ", fmt(max(x$estimate_per_sd)),
    "; positive direction retained ", sum(x$direction_preserved), "/", nrow(x),
    "; nominal P < 0.05 retained ", sum(x$nominal_support_preserved), "/", nrow(x),
    "; largest shift after omitting ", top$omitted_species, "."
  ))
}
report <- c(
  report,
  "",
  "A stable domain-specific association would require positive effects in high-coverage and taxonomic-exclusion models, resistance to single-species omission, and stronger support than the corresponding whole-protein or non-domain metric. These checks do not establish SAMHD1 biochemical activity or causality."
)
writeLines(report, report_output)
cat("Wrote ", model_output, ", ", loo_output, ", and ", report_output, "\n", sep = "")
