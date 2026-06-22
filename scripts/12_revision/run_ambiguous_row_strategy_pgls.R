# DateLife-tree PGLS sensitivity for ambiguous-row score encodings.

project_lib <- file.path(getwd(), "env", "R_library")
if (dir.exists(project_lib)) .libPaths(c(project_lib, .libPaths()))

library(ape)
library(caper)

args <- commandArgs(trailingOnly = TRUE)
data_path <- ifelse(length(args) >= 1, args[[1]], "results/tables/ambiguous_row_scoring_strategies_species.tsv")
tree_path <- ifelse(length(args) >= 2, args[[2]], "data/processed/phylogeny_inputs/opentree_datelife_calibrated_primary68.tre")
output_path <- ifelse(length(args) >= 3, args[[3]], "results/tables/ambiguous_row_scoring_strategy_pgls.tsv")
report_path <- ifelse(length(args) >= 4, args[[4]], "results/reports/ambiguous_row_scoring_strategy_pgls.md")

data_all <- read.delim(data_path, stringsAsFactors = FALSE, check.names = FALSE)
tree_full <- read.tree(tree_path)
for (column in c("strategy_score", "strategy_coverage", "log10_body_mass_g", "max_lifespan_years", "pgls_model_c_mass_clade_residual")) {
  data_all[[column]] <- as.numeric(data_all[[column]])
}
data_all$log10_max_lifespan_years <- log10(data_all$max_lifespan_years)
data_all$clade <- factor(data_all$clade)

prepare_comp <- function(input, needed) {
  columns <- unique(c("opentree_tip_label", needed))
  x <- input[, columns]
  x <- x[complete.cases(x), ]
  x <- x[!duplicated(x$opentree_tip_label), ]
  x$score_z <- as.numeric(scale(x$strategy_score))
  x$coverage_z <- as.numeric(scale(x$strategy_coverage))
  x$clade <- factor(x$clade)
  rownames(x) <- x$opentree_tip_label
  shared <- intersect(tree_full$tip.label, x$opentree_tip_label)
  phy <- drop.tip(tree_full, setdiff(tree_full$tip.label, shared))
  x <- x[phy$tip.label, ]
  comparative.data(phy = phy, data = x, names.col = opentree_tip_label, vcv = TRUE, warn.dropped = FALSE)
}

fit_one <- function(input, strategy, scope, model_name, formula_text, needed) {
  tryCatch({
    comp <- prepare_comp(input, needed)
    model <- pgls(as.formula(formula_text), data = comp, lambda = "ML")
    term <- summary(model)$coefficients["score_z", ]
    data.frame(
      ambiguity_strategy = strategy, scope = scope, model = model_name,
      formula = formula_text, n = nrow(comp$data), lambda = as.numeric(model$param["lambda"]),
      estimate_per_score_sd = unname(term["Estimate"]), se = unname(term["Std. Error"]),
      conf_low = unname(term["Estimate"] - 1.96 * term["Std. Error"]),
      conf_high = unname(term["Estimate"] + 1.96 * term["Std. Error"]),
      t = unname(term["t value"]), p = unname(term["Pr(>|t|)"]), error = "",
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    data.frame(
      ambiguity_strategy = strategy, scope = scope, model = model_name,
      formula = formula_text, n = nrow(input), lambda = NA, estimate_per_score_sd = NA,
      se = NA, conf_low = NA, conf_high = NA, t = NA, p = NA,
      error = conditionMessage(e), stringsAsFactors = FALSE
    )
  })
}

rows <- list()
for (strategy in unique(data_all$ambiguity_strategy)) {
  all_species <- data_all[data_all$ambiguity_strategy == strategy, ]
  birds <- all_species[all_species$clade == "Aves", ]
  specs <- list(
    list(all_species, "all_species", "residual_score", "pgls_model_c_mass_clade_residual ~ score_z", c("pgls_model_c_mass_clade_residual", "strategy_score", "strategy_coverage", "clade")),
    list(all_species, "all_species", "residual_score_coverage", "pgls_model_c_mass_clade_residual ~ score_z + coverage_z", c("pgls_model_c_mass_clade_residual", "strategy_score", "strategy_coverage", "clade")),
    list(all_species, "all_species", "mass_clade_score", "log10_max_lifespan_years ~ log10_body_mass_g + clade + score_z", c("log10_max_lifespan_years", "log10_body_mass_g", "strategy_score", "strategy_coverage", "clade")),
    list(birds, "birds_only", "residual_score", "pgls_model_c_mass_clade_residual ~ score_z", c("pgls_model_c_mass_clade_residual", "strategy_score", "strategy_coverage", "clade")),
    list(birds, "birds_only", "residual_score_coverage", "pgls_model_c_mass_clade_residual ~ score_z + coverage_z", c("pgls_model_c_mass_clade_residual", "strategy_score", "strategy_coverage", "clade"))
  )
  for (spec in specs) {
    rows[[length(rows) + 1]] <- fit_one(spec[[1]], strategy, spec[[2]], spec[[3]], spec[[4]], spec[[5]])
  }
}
models <- do.call(rbind, rows)
models$q <- NA_real_
for (scope_model in unique(paste(models$scope, models$model, sep = "::"))) {
  index <- paste(models$scope, models$model, sep = "::") == scope_model & is.finite(models$p)
  models$q[index] <- p.adjust(models$p[index], method = "BH")
}
models <- models[order(models$scope, models$model, models$ambiguity_strategy), ]

dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(report_path), recursive = TRUE, showWarnings = FALSE)
write.table(models, output_path, sep = "\t", row.names = FALSE, quote = FALSE, na = "")

fmt <- function(x) ifelse(is.finite(x), format(signif(x, 4), scientific = FALSE), "NA")
report <- c(
  "# Ambiguous-Row Strategy PGLS",
  "",
  "Partial/family-ambiguous rows were encoded as fixed-denominator missing, conservative 0.5-weight present-like, or species-specific denominator excluded. All models use the DateLife-calibrated OpenTree chronogram.",
  ""
)
for (i in seq_len(nrow(models))) {
  row <- models[i, ]
  report <- c(report, paste0(
    "- ", row$scope, " / ", row$model, " / ", row$ambiguity_strategy,
    ": beta per score SD = ", fmt(row$estimate_per_score_sd),
    " (95% CI ", fmt(row$conf_low), " to ", fmt(row$conf_high),
    "), P = ", fmt(row$p), ", q = ", fmt(row$q), ", n = ", row$n,
    ifelse(row$error == "", ".", paste0(", error: ", row$error))
  ))
}
writeLines(report, report_path)
cat("Wrote ", output_path, " and ", report_path, "\n", sep = "")
