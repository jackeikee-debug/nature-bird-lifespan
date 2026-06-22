# Final W3 all-module annotation-bias adjusted ranking.

project_lib <- file.path(getwd(), "env", "R", "library")
if (dir.exists(project_lib)) {
  .libPaths(c(project_lib, .libPaths()))
}

library(ape)
library(caper)

args <- commandArgs(trailingOnly = TRUE)
data_path <- ifelse(length(args) >= 1, args[[1]], "data/processed/maintenance_lifespan_phase2_W3_full_background_expanded.tsv")
bias_path <- ifelse(length(args) >= 2, args[[2]], "data/processed/phase2_annotation_bias_covariates.tsv")
table_path <- ifelse(length(args) >= 3, args[[3]], "results/tables/phase2_final_annotation_bias_module_ranking.tsv")
report_path <- ifelse(length(args) >= 4, args[[4]], "results/reports/phase2_final_annotation_bias_module_ranking_report.md")
tree_path <- "data/processed/phylogeny_inputs/opentree_induced_subtree.tre"
variant_name <- "phase2_W3_full_background_sensitivity"

tree_full <- read.tree(tree_path)
data_all <- read.delim(data_path, stringsAsFactors = FALSE)
bias <- read.delim(bias_path, stringsAsFactors = FALSE)

data <- merge(
  data_all[data_all$score_variant == variant_name, ],
  bias[, c("scientific_name", "tier_numeric", "genome_quality_risk", "has_annotation_report", "busco_complete", "busco_missing")],
  by = "scientific_name",
  all.x = TRUE
)

score_cols <- grep("_score$", names(data), value = TRUE)
module_names <- sub("_score$", "", score_cols)

for (col in c(
  "max_lifespan_years", "log10_body_mass_g", "lifespan_residual_log10",
  "pgls_model_c_mass_clade_residual", "tier_numeric", "busco_complete", "busco_missing",
  score_cols
)) {
  data[[col]] <- as.numeric(data[[col]])
}
for (module in module_names) {
  coverage_col <- paste0(module, "_coverage")
  missing_col <- paste0(module, "_missing_matrix_rows")
  if (coverage_col %in% names(data)) data[[coverage_col]] <- as.numeric(data[[coverage_col]])
  if (missing_col %in% names(data)) data[[missing_col]] <- as.numeric(data[[missing_col]])
}
data$log10_max_lifespan_years <- log10(data$max_lifespan_years)
data$clade <- factor(data$clade)

prepare_comp <- function(data, needed_cols) {
  cols <- unique(c("opentree_tip_label", needed_cols))
  model_data <- data[, cols]
  model_data <- model_data[complete.cases(model_data), ]
  if (nrow(model_data) < 12) stop("too_few_complete_rows")
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

fit_one <- function(module, score_col, model_name, formula_text, needed_cols) {
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
for (idx in seq_along(score_cols)) {
  score_col <- score_cols[[idx]]
  module <- module_names[[idx]]
  coverage_col <- paste0(module, "_coverage")
  missing_col <- paste0(module, "_missing_matrix_rows")

  specs <- list(
    residual_base = list(
      formula = paste("pgls_model_c_mass_clade_residual ~", score_col),
      cols = c("pgls_model_c_mass_clade_residual", score_col)
    ),
    residual_tier = list(
      formula = paste("pgls_model_c_mass_clade_residual ~", score_col, "+ tier_numeric"),
      cols = c("pgls_model_c_mass_clade_residual", score_col, "tier_numeric")
    ),
    residual_coverage = list(
      formula = paste("pgls_model_c_mass_clade_residual ~", score_col, "+", coverage_col),
      cols = c("pgls_model_c_mass_clade_residual", score_col, coverage_col)
    ),
    residual_tier_coverage = list(
      formula = paste("pgls_model_c_mass_clade_residual ~", score_col, "+ tier_numeric +", coverage_col),
      cols = c("pgls_model_c_mass_clade_residual", score_col, "tier_numeric", coverage_col)
    ),
    mass_clade_tier = list(
      formula = paste("log10_max_lifespan_years ~ log10_body_mass_g + clade +", score_col, "+ tier_numeric"),
      cols = c("log10_max_lifespan_years", "log10_body_mass_g", "clade", score_col, "tier_numeric")
    ),
    mass_clade_tier_coverage = list(
      formula = paste("log10_max_lifespan_years ~ log10_body_mass_g + clade +", score_col, "+ tier_numeric +", coverage_col),
      cols = c("log10_max_lifespan_years", "log10_body_mass_g", "clade", score_col, "tier_numeric", coverage_col)
    )
  )

  for (model_name in names(specs)) {
    rows[[length(rows) + 1]] <- fit_one(
      module,
      score_col,
      model_name,
      specs[[model_name]]$formula,
      specs[[model_name]]$cols
    )
  }
}

results <- do.call(rbind, rows)
results$p_bh_by_model <- NA
results$rank_by_p <- NA
for (model_name in unique(results$model)) {
  idx <- which(results$model == model_name & is.finite(results$p))
  results$p_bh_by_model[idx] <- p.adjust(results$p[idx], method = "BH")
  ordered <- idx[order(results$p[idx])]
  results$rank_by_p[ordered] <- seq_along(ordered)
}
results <- results[order(results$model, results$p), ]

dir.create(dirname(table_path), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(report_path), recursive = TRUE, showWarnings = FALSE)
write.table(results, table_path, sep = "\t", row.names = FALSE, quote = FALSE)

fmt <- function(x) {
  if (is.na(as.numeric(x))) return("NA")
  signif(as.numeric(x), 4)
}

lines <- c(
  "# Phase 2 Final Annotation-Bias Module Ranking",
  "",
  "## Top Modules by Annotation-Adjusted Model",
  ""
)
for (model_name in c("residual_base", "residual_tier", "residual_coverage", "residual_tier_coverage", "mass_clade_tier", "mass_clade_tier_coverage")) {
  lines <- c(lines, paste0("### ", model_name))
  sub <- results[results$model == model_name & is.finite(results$p), ]
  sub <- head(sub[order(sub$p), ], 6)
  for (idx in seq_len(nrow(sub))) {
    row <- sub[idx, ]
    lines <- c(
      lines,
      paste0(
        "- rank ", row$rank_by_p, ": ", row$maintenance_module,
        ", estimate=", fmt(row$estimate),
        ", p=", fmt(row$p),
        ", BH=", fmt(row$p_bh_by_model),
        ", n=", row$n,
        ", lambda=", fmt(row$lambda)
      )
    )
  }
  lines <- c(lines, "")
}
lines <- c(
  lines,
  "## Interpretation",
  "",
  "Genome-tier adjusted models test whether module ranking survives broad annotation-quality differences. Coverage-adjusted models are deliberately stringent because module scores are confidence-weighted presence fractions; strong attenuation after coverage adjustment indicates that gene observability remains a major vulnerability."
)
writeLines(lines, report_path)

cat("Wrote ", table_path, " and ", report_path, "\n", sep = "")
