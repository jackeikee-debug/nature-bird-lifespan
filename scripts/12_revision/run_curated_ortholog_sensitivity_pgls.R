# Diagnostic DateLife PGLS for stricter full-matrix ortholog rules.

project_lib <- file.path(getwd(), "env", "R_library")
if (dir.exists(project_lib)) .libPaths(c(project_lib, .libPaths()))

library(ape)
library(caper)

args <- commandArgs(trailingOnly = TRUE)
score_path <- ifelse(length(args) >= 1, args[[1]], "results/tables/curated_ortholog_module_scores.tsv")
tree_path <- ifelse(length(args) >= 2, args[[2]], "data/processed/phylogeny_inputs/opentree_datelife_calibrated_primary68.tre")
output_path <- ifelse(length(args) >= 3, args[[3]], "results/tables/curated_ortholog_sensitivity_pgls.tsv")
report_path <- ifelse(length(args) >= 4, args[[4]], "results/reports/curated_ortholog_sensitivity_pgls.md")

data_all <- read.delim(score_path, stringsAsFactors = FALSE, check.names = FALSE)
tree_full <- read.tree(tree_path)
for (column in c("curated_fraction", "log10_body_mass_g", "max_lifespan_years", "pgls_model_c_mass_clade_residual")) data_all[[column]] <- as.numeric(data_all[[column]])
data_all$log10_max_lifespan_years <- log10(data_all$max_lifespan_years)
data_all$clade <- factor(data_all$clade)

prepare_comp <- function(input) {
  needed <- c("opentree_tip_label", "clade", "curated_fraction", "log10_body_mass_g", "log10_max_lifespan_years", "pgls_model_c_mass_clade_residual")
  x <- input[complete.cases(input[, needed]), needed]
  x <- x[!duplicated(x$opentree_tip_label), ]
  if (nrow(x) < 8 || sd(x$curated_fraction) == 0) stop("insufficient observations or ortholog-score variation")
  x$curated_z <- as.numeric(scale(x$curated_fraction))
  x$clade <- factor(x$clade)
  rownames(x) <- x$opentree_tip_label
  shared <- intersect(tree_full$tip.label, x$opentree_tip_label)
  phy <- drop.tip(tree_full, setdiff(tree_full$tip.label, shared))
  x <- x[phy$tip.label, ]
  comparative.data(phy = phy, data = x, names.col = opentree_tip_label, vcv = TRUE, warn.dropped = FALSE)
}

fit_one <- function(input, rule, module, scope, model_name, formula_text) {
  tryCatch({
    comp <- prepare_comp(input)
    model <- pgls(as.formula(formula_text), data = comp, lambda = "ML")
    term <- summary(model)$coefficients["curated_z", ]
    data.frame(
      ortholog_rule = rule, maintenance_module = module, scope = scope, model = model_name,
      formula = formula_text, n = nrow(comp$data), mean_curated_fraction = mean(comp$data$curated_fraction),
      lambda = as.numeric(model$param["lambda"]), estimate_per_sd = unname(term["Estimate"]),
      se = unname(term["Std. Error"]), conf_low = unname(term["Estimate"] - 1.96 * term["Std. Error"]),
      conf_high = unname(term["Estimate"] + 1.96 * term["Std. Error"]),
      t = unname(term["t value"]), p = unname(term["Pr(>|t|)"]), error = "",
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    data.frame(
      ortholog_rule = rule, maintenance_module = module, scope = scope, model = model_name,
      formula = formula_text, n = nrow(input), mean_curated_fraction = mean(input$curated_fraction),
      lambda = NA, estimate_per_sd = NA, se = NA, conf_low = NA, conf_high = NA,
      t = NA, p = NA, error = conditionMessage(e), stringsAsFactors = FALSE
    )
  })
}

rows <- list()
for (rule in unique(data_all$ortholog_rule)) {
  for (module in unique(data_all$maintenance_module)) {
    x <- data_all[data_all$ortholog_rule == rule & data_all$maintenance_module == module, ]
    scopes <- list(
      all_species = x,
      birds_only = x[x$clade == "Aves", ],
      curated_fraction_ge_0.5 = x[x$curated_fraction >= 0.5, ]
    )
    for (scope in names(scopes)) {
      subset <- scopes[[scope]]
      rows[[length(rows) + 1]] <- fit_one(subset, rule, module, scope, "residual_curated_fraction", "pgls_model_c_mass_clade_residual ~ curated_z")
      if (scope == "all_species") {
        rows[[length(rows) + 1]] <- fit_one(subset, rule, module, scope, "mass_clade_curated_fraction", "log10_max_lifespan_years ~ log10_body_mass_g + clade + curated_z")
      }
    }
  }
}
models <- do.call(rbind, rows)
models$q <- NA_real_
for (family in unique(paste(models$ortholog_rule, models$scope, models$model, sep = "::"))) {
  index <- paste(models$ortholog_rule, models$scope, models$model, sep = "::") == family & is.finite(models$p)
  models$q[index] <- p.adjust(models$p[index], method = "BH")
}
models <- models[order(models$ortholog_rule, models$scope, models$model, models$p), ]

dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(report_path), recursive = TRUE, showWarnings = FALSE)
write.table(models, output_path, sep = "\t", row.names = FALSE, quote = FALSE, na = "")

valid <- models[models$error == "", ]
report <- c(
  "# Curated Ortholog-Rule Diagnostic PGLS",
  "",
  "These models are diagnostic because no strict full-matrix rule passed the prespecified combination of clade balance, minimum clade coverage, and non-degenerate species variation.",
  "",
  paste0("- Estimable models: ", nrow(valid), "/", nrow(models), "."),
  paste0("- Models with BH q < 0.05: ", sum(valid$q < 0.05), "."),
  "",
  "A significant strict-rule model cannot by itself rescue the primary claim when the rule is clade-imbalanced or extremely sparse. Interpretation must follow the ortholog-rule eligibility audit."
)
writeLines(report, report_path)
cat("Wrote ", output_path, " and ", report_path, "\n", sep = "")
