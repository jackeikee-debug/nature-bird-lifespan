# Species-level influence diagnostics for the focal transposon/repeat PGLS models.

project_lib <- file.path(getwd(), "env", "R_library")
if (dir.exists(project_lib)) {
  .libPaths(c(project_lib, .libPaths()))
}

library(ape)
library(caper)

args <- commandArgs(trailingOnly = TRUE)
data_path <- ifelse(length(args) >= 1, args[[1]], "data/processed/maintenance_lifespan_phase2_W3_full_background_expanded.tsv")
tree_path <- ifelse(length(args) >= 2, args[[2]], "data/processed/phylogeny_inputs/opentree_induced_subtree.tre")
influence_path <- ifelse(length(args) >= 3, args[[3]], "results/tables/communications_biology_species_influence_pgls.tsv")
diagnostics_path <- ifelse(length(args) >= 4, args[[4]], "results/tables/communications_biology_pgls_residual_diagnostics.tsv")
report_path <- ifelse(length(args) >= 5, args[[5]], "results/reports/communications_biology_species_influence_report.md")

variant <- "phase2_W3_full_background_sensitivity"
score_col <- "transposon_repeat_suppression_score"
tree_full <- read.tree(tree_path)
data_all <- read.delim(data_path, stringsAsFactors = FALSE)
data_all <- data_all[data_all$score_variant == variant, ]
data_all$log10_max_lifespan_years <- log10(as.numeric(data_all$max_lifespan_years))
data_all$log10_body_mass_g <- as.numeric(data_all$log10_body_mass_g)
data_all$pgls_model_c_mass_clade_residual <- as.numeric(data_all$pgls_model_c_mass_clade_residual)
data_all[[score_col]] <- as.numeric(data_all[[score_col]])
data_all$clade <- factor(data_all$clade)

prepare_comp <- function(data, needed_cols) {
  metadata_cols <- c("scientific_name", "clade", "genome_analysis_tier", "transposon_repeat_suppression_coverage")
  cols <- unique(c("opentree_tip_label", metadata_cols, needed_cols))
  data <- data[, cols]
  data <- data[complete.cases(data[, unique(c("opentree_tip_label", needed_cols))]), ]
  rownames(data) <- data$opentree_tip_label
  shared <- intersect(tree_full$tip.label, data$opentree_tip_label)
  pruned <- drop.tip(tree_full, setdiff(tree_full$tip.label, shared))
  if (is.null(pruned$edge.length) || any(!is.finite(pruned$edge.length)) || any(pruned$edge.length <= 0)) {
    pruned <- compute.brlen(pruned, method = "Grafen")
  }
  data <- data[pruned$tip.label, ]
  data$clade <- factor(data$clade)
  comparative.data(
    phy = pruned,
    data = data,
    names.col = opentree_tip_label,
    vcv = TRUE,
    warn.dropped = FALSE
  )
}

fit_spec <- function(data, formula_text, needed_cols) {
  comp <- prepare_comp(data, needed_cols)
  model <- pgls(as.formula(formula_text), data = comp, lambda = "ML")
  coefs <- summary(model)$coefficients
  term <- coefs[score_col, ]
  list(
    model = model,
    comp = comp,
    estimate = unname(term["Estimate"]),
    se = unname(term["Std. Error"]),
    t = unname(term["t value"]),
    p = unname(term["Pr(>|t|)"]),
    lambda = as.numeric(model$param["lambda"]),
    n = length(model$residuals)
  )
}

specs <- list(
  mass_clade_module = list(
    formula = paste("log10_max_lifespan_years ~ log10_body_mass_g + clade +", score_col),
    cols = c("log10_max_lifespan_years", "log10_body_mass_g", "clade", score_col)
  ),
  pgls_clade_residual_module = list(
    formula = paste("pgls_model_c_mass_clade_residual ~", score_col),
    cols = c("pgls_model_c_mass_clade_residual", score_col)
  )
)

influence_rows <- list()
baseline_fits <- list()
for (model_name in names(specs)) {
  spec <- specs[[model_name]]
  baseline <- fit_spec(data_all, spec$formula, spec$cols)
  baseline_fits[[model_name]] <- baseline
  ordered_data <- baseline$comp$data
  threshold <- 2 / sqrt(baseline$n)

  for (tip in rownames(ordered_data)) {
    fit <- tryCatch(
      fit_spec(data_all[data_all$opentree_tip_label != tip, ], spec$formula, spec$cols),
      error = function(e) NULL
    )
    meta <- ordered_data[tip, ]
    if (is.null(fit)) {
      influence_rows[[length(influence_rows) + 1]] <- data.frame(
        model = model_name, omitted_tip_label = tip, scientific_name = meta$scientific_name,
        clade = as.character(meta$clade), baseline_n = baseline$n, baseline_estimate = baseline$estimate,
        baseline_se = baseline$se, baseline_p = baseline$p, baseline_lambda = baseline$lambda,
        loo_n = NA, loo_estimate = NA, loo_se = NA, loo_p = NA, loo_lambda = NA,
        loo_conf_low = NA, loo_conf_high = NA, delta_estimate = NA, standardized_estimate_shift = NA,
        influence_threshold = threshold, high_influence = NA, direction_preserved = NA,
        nominal_support_preserved = NA, error = "model_fit_failed", stringsAsFactors = FALSE
      )
    } else {
      delta <- fit$estimate - baseline$estimate
      standardized <- delta / baseline$se
      influence_rows[[length(influence_rows) + 1]] <- data.frame(
        model = model_name, omitted_tip_label = tip, scientific_name = meta$scientific_name,
        clade = as.character(meta$clade), baseline_n = baseline$n, baseline_estimate = baseline$estimate,
        baseline_se = baseline$se, baseline_p = baseline$p, baseline_lambda = baseline$lambda,
        loo_n = fit$n, loo_estimate = fit$estimate, loo_se = fit$se, loo_p = fit$p, loo_lambda = fit$lambda,
        loo_conf_low = fit$estimate - 1.96 * fit$se, loo_conf_high = fit$estimate + 1.96 * fit$se,
        delta_estimate = delta, standardized_estimate_shift = standardized,
        influence_threshold = threshold, high_influence = abs(standardized) > threshold,
        direction_preserved = sign(fit$estimate) == sign(baseline$estimate),
        nominal_support_preserved = fit$p < 0.05, error = "", stringsAsFactors = FALSE
      )
    }
  }
}

influence <- do.call(rbind, influence_rows)
influence <- as.data.frame(influence, stringsAsFactors = FALSE)
numeric_columns <- c(
  "baseline_n", "baseline_estimate", "baseline_se", "baseline_p", "baseline_lambda",
  "loo_n", "loo_estimate", "loo_se", "loo_p", "loo_lambda", "loo_conf_low", "loo_conf_high",
  "delta_estimate", "standardized_estimate_shift", "influence_threshold"
)
for (column in numeric_columns) {
  influence[[column]] <- as.numeric(influence[[column]])
}
influence[["model"]] <- as.character(unlist(influence[["model"]]))
influence <- influence[order(influence[["model"]], -abs(influence[["standardized_estimate_shift"]])), ]

# Generalized-leverage and Cook-distance approximation for the mass-plus-clade PGLS model.
baseline <- baseline_fits[["mass_clade_module"]]
model <- baseline$model
X <- model$x
Vinv <- solve(model$Vt)
XtVinvX_inv <- solve(t(X) %*% Vinv %*% X)
H <- X %*% XtVinvX_inv %*% t(X) %*% Vinv
leverage <- pmax(0, Re(diag(H)))
residual <- as.numeric(model$residuals)
fitted_value <- as.numeric(model$fitted)
n <- length(residual)
p <- ncol(X)
mse <- as.numeric(t(residual) %*% Vinv %*% residual) / (n - p)
observation_variance <- pmax(.Machine$double.eps, Re(diag(model$Vt)))
standardized_residual <- residual / sqrt(pmax(.Machine$double.eps, mse * observation_variance * (1 - leverage)))
cook_gls <- (standardized_residual^2 * leverage) / (p * pmax(.Machine$double.eps, 1 - leverage))
tip_labels <- rownames(X)
ordered_data <- baseline$comp$data[tip_labels, ]
diagnostics <- data.frame(
  opentree_tip_label = tip_labels,
  scientific_name = ordered_data$scientific_name,
  clade = as.character(ordered_data$clade),
  genome_analysis_tier = ordered_data$genome_analysis_tier,
  transposon_repeat_suppression_score = ordered_data[[score_col]],
  transposon_repeat_suppression_coverage = ordered_data$transposon_repeat_suppression_coverage,
  fitted_log10_lifespan = fitted_value,
  raw_residual = residual,
  standardized_residual = standardized_residual,
  generalized_leverage = leverage,
  cook_distance_gls_approx = cook_gls,
  leverage_threshold = 2 * p / n,
  cook_threshold = 4 / n,
  high_leverage = leverage > 2 * p / n,
  high_cook = cook_gls > 4 / n,
  large_standardized_residual = abs(standardized_residual) > 2,
  stringsAsFactors = FALSE
)
diagnostics <- diagnostics[order(-diagnostics$cook_distance_gls_approx), ]

dir.create(dirname(influence_path), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(diagnostics_path), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(report_path), recursive = TRUE, showWarnings = FALSE)
write.table(influence, influence_path, sep = "\t", row.names = FALSE, quote = FALSE)
write.table(diagnostics, diagnostics_path, sep = "\t", row.names = FALSE, quote = FALSE)

report <- c(
  "# Species-Level PGLS Influence Diagnostics",
  "",
  "The focal transposon/repeat models were refitted after omitting each of the 68 species in turn. A standardized estimate shift larger than 2/sqrt(n) is used as a screening heuristic, not as a formal hypothesis test. Generalized leverage and Cook distance are GLS approximations based on the fitted PGLS covariance matrix.",
  ""
)
for (model_name in names(specs)) {
  x <- influence[influence$model == model_name & influence$error == "", ]
  top <- x[which.max(abs(x$standardized_estimate_shift)), ]
  report <- c(
    report,
    paste0("## ", model_name),
    "",
    paste0("- Baseline estimate: ", signif(top$baseline_estimate[1], 4), "; P = ", signif(top$baseline_p[1], 4), "; lambda = ", signif(top$baseline_lambda[1], 4), "."),
    paste0("- Leave-one-species estimate range: ", signif(min(x$loo_estimate), 4), " to ", signif(max(x$loo_estimate), 4), "."),
    paste0("- Positive direction retained: ", sum(x$direction_preserved), "/", nrow(x), "; nominal P < 0.05 retained: ", sum(x$nominal_support_preserved), "/", nrow(x), "."),
    paste0("- Heuristic high-influence omissions: ", sum(x$high_influence), "/", nrow(x), "."),
    paste0("- Largest standardized shift: ", top$scientific_name, " (", signif(top$standardized_estimate_shift, 3), "; omitted estimate = ", signif(top$loo_estimate, 4), ")."),
    ""
  )
}
report <- c(
  report,
  "## Residual and leverage screen",
  "",
  paste0("- High generalized leverage: ", sum(diagnostics$high_leverage), "/", nrow(diagnostics), "."),
  paste0("- Cook-distance approximation above 4/n: ", sum(diagnostics$high_cook), "/", nrow(diagnostics), "."),
  paste0("- Absolute standardized residual above 2: ", sum(diagnostics$large_standardized_residual), "/", nrow(diagnostics), "."),
  paste0("- Largest Cook-distance approximation: ", diagnostics$scientific_name[1], " (", signif(diagnostics$cook_distance_gls_approx[1], 3), ").")
)
writeLines(report, report_path)
cat("Wrote ", influence_path, ", ", diagnostics_path, ", and ", report_path, "\n", sep = "")
