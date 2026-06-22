# Annotation-bias adjusted models for Phase 2 W2-expanded scores.

project_lib <- file.path(getwd(), "env", "R", "library")
if (dir.exists(project_lib)) {
  .libPaths(c(project_lib, .libPaths()))
}

library(ape)
library(caper)

args <- commandArgs(trailingOnly = TRUE)
data_path <- ifelse(length(args) >= 1, args[[1]], "data/processed/maintenance_lifespan_phase2_W2_expanded.tsv")
bias_path <- ifelse(length(args) >= 2, args[[2]], "data/processed/phase2_annotation_bias_covariates.tsv")
table_path <- ifelse(length(args) >= 3, args[[3]], "results/tables/phase2_W2_annotation_bias_models.tsv")
report_path <- ifelse(length(args) >= 4, args[[4]], "results/reports/phase2_W2_annotation_bias_models_report.md")
tree_path <- "data/processed/phylogeny_inputs/opentree_induced_subtree.tre"

variant_name <- "phase2_W2_crossdb_sensitivity"
target_modules <- c("transposon_repeat_suppression", "chromatin_repression_heterochromatin")

tree_full <- read.tree(tree_path)
data_all <- read.delim(data_path, stringsAsFactors = FALSE)
bias <- read.delim(bias_path, stringsAsFactors = FALSE)
data <- merge(
  data_all[data_all$score_variant == variant_name, ],
  bias[, c("scientific_name", "tier_numeric", "genome_quality_risk", "has_annotation_report", "busco_complete", "busco_missing")],
  by = "scientific_name",
  all.x = TRUE
)

for (col in c(
  "max_lifespan_years", "log10_body_mass_g", "lifespan_residual_log10",
  "pgls_model_c_mass_clade_residual", "tier_numeric", "busco_complete", "busco_missing"
)) {
  data[[col]] <- as.numeric(data[[col]])
}
data$log10_max_lifespan_years <- log10(data$max_lifespan_years)
data$clade <- factor(data$clade)
data$genome_quality_risk <- factor(data$genome_quality_risk)

prepare_comp <- function(data, needed_cols) {
  cols <- unique(c("opentree_tip_label", needed_cols))
  model_data <- data[, cols]
  model_data <- model_data[complete.cases(model_data), ]
  rownames(model_data) <- model_data$opentree_tip_label
  shared <- intersect(tree_full$tip.label, model_data$opentree_tip_label)
  pruned <- drop.tip(tree_full, setdiff(tree_full$tip.label, shared))
  pruned <- compute.brlen(pruned, method = "Grafen")
  model_data <- model_data[shared, ]
  comparative.data(
    phy = pruned,
    data = model_data,
    names.col = opentree_tip_label,
    vcv = TRUE,
    warn.dropped = FALSE
  )
}

extract_term <- function(model, predictor) {
  coef_table <- summary(model)$coefficients
  if (!predictor %in% rownames(coef_table)) {
    return(c(estimate = NA, se = NA, t = NA, p = NA))
  }
  c(
    estimate = coef_table[predictor, "Estimate"],
    se = coef_table[predictor, "Std. Error"],
    t = coef_table[predictor, "t value"],
    p = coef_table[predictor, "Pr(>|t|)"]
  )
}

fit_one <- function(module, model_name, formula_text, needed_cols, score_col) {
  result <- tryCatch(
    {
      comp <- prepare_comp(data, needed_cols)
      model <- pgls(as.formula(formula_text), data = comp, lambda = "ML")
      term <- extract_term(model, score_col)
      data.frame(
        score_variant = variant_name,
        maintenance_module = module,
        score_col = score_col,
        model = model_name,
        formula = formula_text,
        n = length(model$residuals),
        lambda = as.numeric(model$param["lambda"]),
        logLik = as.numeric(logLik(model)),
        AIC = AIC(model),
        estimate = term["estimate"],
        se = term["se"],
        t = term["t"],
        p = term["p"],
        error = ""
      )
    },
    error = function(e) {
      data.frame(
        score_variant = variant_name,
        maintenance_module = module,
        score_col = score_col,
        model = model_name,
        formula = formula_text,
        n = NA,
        lambda = NA,
        logLik = NA,
        AIC = NA,
        estimate = NA,
        se = NA,
        t = NA,
        p = NA,
        error = conditionMessage(e)
      )
    }
  )
  result
}

rows <- list()
for (module in target_modules) {
  score_col <- paste0(module, "_score")
  coverage_col <- paste0(module, "_coverage")
  missing_col <- paste0(module, "_missing_matrix_rows")
  for (col in c(score_col, coverage_col, missing_col)) {
    data[[col]] <- as.numeric(data[[col]])
  }

  specs <- list(
    residual_base = list(
      formula = paste("pgls_model_c_mass_clade_residual ~", score_col),
      cols = c("pgls_model_c_mass_clade_residual", score_col)
    ),
    residual_tier = list(
      formula = paste("pgls_model_c_mass_clade_residual ~", score_col, "+ tier_numeric"),
      cols = c("pgls_model_c_mass_clade_residual", score_col, "tier_numeric")
    ),
    residual_missing = list(
      formula = paste("pgls_model_c_mass_clade_residual ~", score_col, "+", missing_col),
      cols = c("pgls_model_c_mass_clade_residual", score_col, missing_col)
    ),
    residual_coverage = list(
      formula = paste("pgls_model_c_mass_clade_residual ~", score_col, "+", coverage_col),
      cols = c("pgls_model_c_mass_clade_residual", score_col, coverage_col)
    ),
    residual_tier_missing = list(
      formula = paste("pgls_model_c_mass_clade_residual ~", score_col, "+ tier_numeric +", missing_col),
      cols = c("pgls_model_c_mass_clade_residual", score_col, "tier_numeric", missing_col)
    ),
    residual_tier_coverage = list(
      formula = paste("pgls_model_c_mass_clade_residual ~", score_col, "+ tier_numeric +", coverage_col),
      cols = c("pgls_model_c_mass_clade_residual", score_col, "tier_numeric", coverage_col)
    ),
    lifespan_mass_clade_tier_missing = list(
      formula = paste("log10_max_lifespan_years ~ log10_body_mass_g + clade +", score_col, "+ tier_numeric +", missing_col),
      cols = c("log10_max_lifespan_years", "log10_body_mass_g", "clade", score_col, "tier_numeric", missing_col)
    )
  )

  for (model_name in names(specs)) {
    rows[[length(rows) + 1]] <- fit_one(
      module,
      model_name,
      specs[[model_name]]$formula,
      specs[[model_name]]$cols,
      score_col
    )
  }
}

results <- do.call(rbind, rows)
idx <- which(is.finite(results$p))
results$p_bh <- NA
results$p_bh[idx] <- p.adjust(results$p[idx], method = "BH")
results <- results[order(results$maintenance_module, results$model), ]

dir.create(dirname(table_path), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(report_path), recursive = TRUE, showWarnings = FALSE)
write.table(results, table_path, sep = "\t", row.names = FALSE, quote = FALSE)

fmt_row <- function(module, model_name) {
  row <- results[results$maintenance_module == module & results$model == model_name, ][1, ]
  paste0(
    "- ", module, " / ", model_name,
    ": estimate=", signif(row$estimate, 4),
    ", p=", signif(row$p, 4),
    ", BH=", signif(row$p_bh, 4),
    ", n=", row$n,
    ", lambda=", signif(row$lambda, 4)
  )
}

report_lines <- c(
  "# Phase 2 W2 Annotation-Bias Models",
  "",
  "## Key Results",
  "",
  fmt_row("transposon_repeat_suppression", "residual_base"),
  fmt_row("transposon_repeat_suppression", "residual_tier"),
  fmt_row("transposon_repeat_suppression", "residual_missing"),
  fmt_row("transposon_repeat_suppression", "residual_coverage"),
  fmt_row("transposon_repeat_suppression", "residual_tier_missing"),
  fmt_row("chromatin_repression_heterochromatin", "residual_base"),
  fmt_row("chromatin_repression_heterochromatin", "residual_tier"),
  fmt_row("chromatin_repression_heterochromatin", "residual_missing"),
  fmt_row("chromatin_repression_heterochromatin", "residual_coverage"),
  fmt_row("chromatin_repression_heterochromatin", "residual_tier_missing"),
  "",
  "## Interpretation",
  "",
  "Tier and missingness-adjusted models ask whether the repeat/chromatin signal survives genome annotation quality controls. Coverage-adjusted models are intentionally stringent because confidence-weighted scores are mathematically tied to module coverage; attenuation there should be interpreted as a specificity stress test rather than a simple failure."
)
writeLines(report_lines, report_path)

cat("Wrote ", table_path, " and ", report_path, "\n", sep = "")
