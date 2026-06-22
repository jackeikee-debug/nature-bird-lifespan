# PGLS sensitivity to evidence-confidence weights in the final 200-gene matrix.

project_lib <- file.path(getwd(), "env", "R_library")
if (dir.exists(project_lib)) .libPaths(c(project_lib, .libPaths()))

library(ape)
library(caper)

args <- commandArgs(trailingOnly = TRUE)
data_path <- ifelse(length(args) >= 1, args[[1]], "results/tables/module_weight_sensitivity_species.tsv")
tree_path <- ifelse(length(args) >= 2, args[[2]], "data/processed/phylogeny_inputs/opentree_datelife_calibrated_primary68.tre")
output_path <- ifelse(length(args) >= 3, args[[3]], "results/tables/module_weight_sensitivity_pgls.tsv")
report_path <- ifelse(length(args) >= 4, args[[4]], "results/reports/module_weight_sensitivity_pgls.md")

data_all <- read.delim(data_path, stringsAsFactors = FALSE, check.names = FALSE)
tree_full <- read.tree(tree_path)
for (column in c("module_score", "coverage_fraction", "log10_body_mass_g", "log10_max_lifespan_years", "pgls_model_c_mass_clade_residual")) {
  data_all[[column]] <- as.numeric(data_all[[column]])
}
data_all$clade <- factor(data_all$clade)

fit_one <- function(input, scheme, module, endpoint, formula_text, needed) {
  tryCatch({
    x <- input[complete.cases(input[, needed]), unique(c("opentree_tip_label", needed))]
    x <- x[!duplicated(x$opentree_tip_label), ]
    x$score_z <- as.numeric(scale(x$module_score))
    x$clade <- factor(x$clade)
    rownames(x) <- x$opentree_tip_label
    shared <- intersect(tree_full$tip.label, x$opentree_tip_label)
    phy <- drop.tip(tree_full, setdiff(tree_full$tip.label, shared))
    x <- x[phy$tip.label, ]
    comp <- comparative.data(phy = phy, data = x, names.col = opentree_tip_label, vcv = TRUE, warn.dropped = FALSE)
    model <- pgls(as.formula(formula_text), data = comp, lambda = "ML")
    term <- summary(model)$coefficients["score_z", ]
    data.frame(
      weight_scheme = scheme, maintenance_module = module, endpoint = endpoint,
      formula = formula_text, n = nrow(comp$data), lambda = as.numeric(model$param["lambda"]),
      estimate_per_score_sd = unname(term["Estimate"]), se = unname(term["Std. Error"]),
      conf_low = unname(term["Estimate"] - 1.96 * term["Std. Error"]),
      conf_high = unname(term["Estimate"] + 1.96 * term["Std. Error"]),
      t = unname(term["t value"]), p = unname(term["Pr(>|t|)"]), error = "",
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    data.frame(
      weight_scheme = scheme, maintenance_module = module, endpoint = endpoint,
      formula = formula_text, n = nrow(input), lambda = NA, estimate_per_score_sd = NA,
      se = NA, conf_low = NA, conf_high = NA, t = NA, p = NA,
      error = conditionMessage(e), stringsAsFactors = FALSE
    )
  })
}

rows <- list()
for (scheme in unique(data_all$weight_scheme)) {
  for (module in unique(data_all$maintenance_module)) {
    x <- data_all[data_all$weight_scheme == scheme & data_all$maintenance_module == module, ]
    rows[[length(rows) + 1]] <- fit_one(
      x, scheme, module, "primary_mass_clade",
      "log10_max_lifespan_years ~ log10_body_mass_g + clade + score_z",
      c("log10_max_lifespan_years", "log10_body_mass_g", "clade", "module_score")
    )
    rows[[length(rows) + 1]] <- fit_one(
      x, scheme, module, "parallel_residual",
      "pgls_model_c_mass_clade_residual ~ score_z",
      c("pgls_model_c_mass_clade_residual", "clade", "module_score")
    )
  }
}
models <- do.call(rbind, rows)
models$q_within_scheme_endpoint <- NA_real_
families <- unique(paste(models$weight_scheme, models$endpoint, sep = "::"))
for (family in families) {
  index <- paste(models$weight_scheme, models$endpoint, sep = "::") == family & is.finite(models$p)
  models$q_within_scheme_endpoint[index] <- p.adjust(models$p[index], method = "BH")
}
models <- models[order(models$endpoint, models$maintenance_module, models$weight_scheme), ]

dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(report_path), recursive = TRUE, showWarnings = FALSE)
write.table(models, output_path, sep = "\t", row.names = FALSE, quote = FALSE, na = "")

fmt <- function(x) ifelse(is.finite(x), format(signif(x, 4), scientific = FALSE), "NA")
report <- c(
  "# Module Score-Weight Sensitivity PGLS",
  "",
  "The primary confidence weights were compared with unweighted local presence and a more conservative alternative. External sensitivity rows remained excluded from primary scoring in every scheme.",
  ""
)
for (i in seq_len(nrow(models))) {
  row <- models[i, ]
  report <- c(report, paste0(
    "- ", row$endpoint, " / ", row$maintenance_module, " / ", row$weight_scheme,
    ": beta per score SD = ", fmt(row$estimate_per_score_sd),
    " (95% CI ", fmt(row$conf_low), " to ", fmt(row$conf_high),
    "), P = ", fmt(row$p), ", q = ", fmt(row$q_within_scheme_endpoint), ", n = ", row$n, "."
  ))
}
writeLines(report, report_path)
cat("Wrote ", output_path, " and ", report_path, "\n", sep = "")
