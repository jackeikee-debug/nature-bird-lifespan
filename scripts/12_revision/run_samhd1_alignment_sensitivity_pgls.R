# PGLS sensitivity of the SAMHD1 signal to removal of gap-heavy alignment columns.

project_lib <- file.path(getwd(), "env", "R_library")
if (dir.exists(project_lib)) .libPaths(c(project_lib, .libPaths()))

library(ape)
library(caper)

args <- commandArgs(trailingOnly = TRUE)
data_path <- ifelse(length(args) >= 1, args[[1]], "results/tables/samhd1_alignment_species_qc.tsv")
tree_path <- ifelse(length(args) >= 2, args[[2]], "data/processed/phylogeny_inputs/opentree_datelife_calibrated_primary68.tre")
output_path <- ifelse(length(args) >= 3, args[[3]], "results/tables/samhd1_alignment_sensitivity_pgls.tsv")
report_path <- ifelse(length(args) >= 4, args[[4]], "results/reports/samhd1_alignment_sensitivity_pgls.md")

data_all <- read.delim(data_path, stringsAsFactors = FALSE, check.names = FALSE)
tree_full <- read.tree(tree_path)
numeric_columns <- c(
  "original_domain_reference_coverage", "domain_identity_coverage_product",
  "nondomain_identity_coverage_product", "domain_minus_nondomain_product",
  "pgls_model_c_mass_clade_residual"
)
for (column in numeric_columns) data_all[[column]] <- as.numeric(data_all[[column]])
data_all$qualified_original_domain_coverage_ge_0.5 <- tolower(as.character(data_all$qualified_original_domain_coverage_ge_0.5)) %in% c("true", "1", "yes")
data_all <- data_all[data_all$qualified_original_domain_coverage_ge_0.5, ]

fit_one <- function(input, variant, scope, predictor) {
  tryCatch({
    x <- input[complete.cases(input[, c("opentree_tip_label", "pgls_model_c_mass_clade_residual", predictor)]), ]
    x <- x[!duplicated(x$opentree_tip_label), ]
    x$predictor_z <- as.numeric(scale(x[[predictor]]))
    rownames(x) <- x$opentree_tip_label
    shared <- intersect(tree_full$tip.label, x$opentree_tip_label)
    phy <- drop.tip(tree_full, setdiff(tree_full$tip.label, shared))
    x <- x[phy$tip.label, ]
    comp <- comparative.data(phy = phy, data = x, names.col = opentree_tip_label, vcv = TRUE, warn.dropped = FALSE)
    formula <- pgls_model_c_mass_clade_residual ~ predictor_z
    ml <- tryCatch(pgls(formula, data = comp, lambda = "ML"), error = function(e) e)
    if (inherits(ml, "error")) {
      model <- pgls(formula, data = comp, lambda = 1e-6)
      lambda_method <- "fixed_near_zero_optimizer_fallback"
      fit_note <- conditionMessage(ml)
    } else {
      model <- ml
      lambda_method <- "maximum_likelihood"
      fit_note <- ""
    }
    term <- summary(model)$coefficients["predictor_z", ]
    data.frame(
      alignment_variant = variant, scope = scope, predictor = predictor,
      n = nrow(comp$data), lambda = as.numeric(model$param["lambda"]), lambda_method = lambda_method,
      estimate_per_sd = unname(term["Estimate"]), se = unname(term["Std. Error"]),
      conf_low = unname(term["Estimate"] - 1.96 * term["Std. Error"]),
      conf_high = unname(term["Estimate"] + 1.96 * term["Std. Error"]),
      t = unname(term["t value"]), p = unname(term["Pr(>|t|)"]), fit_note = fit_note, error = "",
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    data.frame(
      alignment_variant = variant, scope = scope, predictor = predictor,
      n = nrow(input), lambda = NA, lambda_method = "failed", estimate_per_sd = NA,
      se = NA, conf_low = NA, conf_high = NA, t = NA, p = NA, fit_note = "",
      error = conditionMessage(e), stringsAsFactors = FALSE
    )
  })
}

rows <- list()
predictors <- c("domain_identity_coverage_product", "domain_minus_nondomain_product")
for (variant in unique(data_all$alignment_variant)) {
  x <- data_all[data_all$alignment_variant == variant, ]
  scopes <- list(
    all_species = x,
    original_domain_coverage_ge_0.8 = x[x$original_domain_reference_coverage >= 0.8, ]
  )
  for (scope in names(scopes)) {
    for (predictor in predictors) {
      rows[[length(rows) + 1]] <- fit_one(scopes[[scope]], variant, scope, predictor)
    }
  }
}
models <- do.call(rbind, rows)
models$q_within_variant_scope <- NA_real_
families <- unique(paste(models$alignment_variant, models$scope, sep = "::"))
for (family in families) {
  index <- paste(models$alignment_variant, models$scope, sep = "::") == family & is.finite(models$p)
  models$q_within_variant_scope[index] <- p.adjust(models$p[index], method = "BH")
}
models <- models[order(models$scope, models$predictor, models$alignment_variant), ]

dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(report_path), recursive = TRUE, showWarnings = FALSE)
write.table(models, output_path, sep = "\t", row.names = FALSE, quote = FALSE, na = "")

fmt <- function(x) ifelse(is.finite(x), format(signif(x, 4), scientific = FALSE), "NA")
report <- c(
  "# SAMHD1 Alignment-Column Sensitivity PGLS",
  "",
  "The original SAMHD1 metrics were compared with metrics recalculated after removal of human-reference columns having target-sequence gap fractions greater than 0.70 or 0.50. The same species qualifying under the original domain-coverage threshold were retained for direct comparison.",
  ""
)
for (i in seq_len(nrow(models))) {
  row <- models[i, ]
  report <- c(report, paste0(
    "- ", row$scope, " / ", row$predictor, " / ", row$alignment_variant,
    ": beta per SD = ", fmt(row$estimate_per_sd),
    " (95% CI ", fmt(row$conf_low), " to ", fmt(row$conf_high),
    "), P = ", fmt(row$p), ", q = ", fmt(row$q_within_variant_scope), ", n = ", row$n, "."
  ))
}
writeLines(report, report_path)
cat("Wrote ", output_path, " and ", report_path, "\n", sep = "")
