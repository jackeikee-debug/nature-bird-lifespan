# Joint transposon/chromatin models for Phase 2 W2-expanded scores.
#
# This tests whether the W2 transposon signal survives after explicitly
# adjusting for the newly available chromatin-repression module.

project_lib <- file.path(getwd(), "env", "R", "library")
if (dir.exists(project_lib)) {
  .libPaths(c(project_lib, .libPaths()))
}

library(ape)
library(caper)

args <- commandArgs(trailingOnly = TRUE)
data_path <- ifelse(length(args) >= 1, args[[1]], "data/processed/maintenance_lifespan_phase2_W2_expanded.tsv")
table_path <- ifelse(length(args) >= 2, args[[2]], "results/tables/phase2_W2_transposon_chromatin_joint_models.tsv")
diag_path <- ifelse(length(args) >= 3, args[[3]], "results/tables/phase2_W2_transposon_chromatin_collinearity.tsv")
report_path <- ifelse(length(args) >= 4, args[[4]], "results/reports/phase2_W2_transposon_chromatin_joint_models_report.md")
tree_path <- "data/processed/phylogeny_inputs/opentree_induced_subtree.tre"

variant_name <- "phase2_W2_crossdb_sensitivity"
transposon_col <- "transposon_repeat_suppression_score"
chromatin_col <- "chromatin_repression_heterochromatin_score"

tree_full <- read.tree(tree_path)
data_all <- read.delim(data_path, stringsAsFactors = FALSE)
data <- data_all[data_all$score_variant == variant_name, ]

numeric_cols <- c(
  "body_mass_g",
  "max_lifespan_years",
  "log10_body_mass_g",
  "lifespan_residual_log10",
  "pgls_model_c_mass_clade_residual",
  transposon_col,
  chromatin_col,
  "transposon_repeat_suppression_coverage",
  "chromatin_repression_heterochromatin_coverage"
)
for (col in numeric_cols) {
  data[[col]] <- as.numeric(data[[col]])
}
data$log10_max_lifespan_years <- log10(data$max_lifespan_years)
data$clade <- factor(data$clade)
data$flight_status <- factor(data$flight_status)

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

extract_term <- function(model, term) {
  coef_table <- summary(model)$coefficients
  if (!term %in% rownames(coef_table)) {
    return(c(estimate = NA, se = NA, t = NA, p = NA))
  }
  c(
    estimate = coef_table[term, "Estimate"],
    se = coef_table[term, "Std. Error"],
    t = coef_table[term, "t value"],
    p = coef_table[term, "Pr(>|t|)"]
  )
}

fit_pgls_terms <- function(model_name, formula_text, needed_cols, terms) {
  result <- tryCatch(
    {
      comp <- prepare_comp(data, needed_cols)
      model <- pgls(as.formula(formula_text), data = comp, lambda = "ML")
      rows <- list()
      for (term in terms) {
        term_stats <- extract_term(model, term)
        rows[[length(rows) + 1]] <- data.frame(
          score_variant = variant_name,
          model = model_name,
          formula = formula_text,
          term = term,
          n = length(model$residuals),
          lambda = as.numeric(model$param["lambda"]),
          logLik = as.numeric(logLik(model)),
          AIC = AIC(model),
          estimate = term_stats["estimate"],
          se = term_stats["se"],
          t = term_stats["t"],
          p = term_stats["p"],
          error = ""
        )
      }
      do.call(rbind, rows)
    },
    error = function(e) {
      data.frame(
        score_variant = variant_name,
        model = model_name,
        formula = formula_text,
        term = terms,
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

vif_for_predictor <- function(df, y, predictors, predictor) {
  other_predictors <- setdiff(predictors, predictor)
  formula_text <- paste(predictor, "~", paste(other_predictors, collapse = " + "))
  fit <- lm(as.formula(formula_text), data = df[, c(predictor, other_predictors), drop = FALSE])
  r2 <- summary(fit)$r.squared
  1 / (1 - r2)
}

complete_diag <- data[complete.cases(data[, c(transposon_col, chromatin_col, "log10_body_mass_g", "clade")]), ]
predictors <- c(transposon_col, chromatin_col, "log10_body_mass_g", "clade")
diag_rows <- data.frame(
  score_variant = variant_name,
  n = nrow(complete_diag),
  pearson_r_scores = cor(complete_diag[[transposon_col]], complete_diag[[chromatin_col]], method = "pearson"),
  spearman_r_scores = cor(complete_diag[[transposon_col]], complete_diag[[chromatin_col]], method = "spearman"),
  pearson_r_coverage = cor(
    complete_diag[["transposon_repeat_suppression_coverage"]],
    complete_diag[["chromatin_repression_heterochromatin_coverage"]],
    method = "pearson",
    use = "complete.obs"
  ),
  vif_transposon_score = vif_for_predictor(complete_diag, "unused", predictors, transposon_col),
  vif_chromatin_score = vif_for_predictor(complete_diag, "unused", predictors, chromatin_col)
)

model_specs <- list(
  mass_clade_joint = list(
    formula = paste(
      "log10_max_lifespan_years ~ log10_body_mass_g + clade +",
      transposon_col, "+", chromatin_col
    ),
    cols = c("log10_max_lifespan_years", "log10_body_mass_g", "clade", transposon_col, chromatin_col)
  ),
  residual_joint = list(
    formula = paste("lifespan_residual_log10 ~", transposon_col, "+", chromatin_col),
    cols = c("lifespan_residual_log10", transposon_col, chromatin_col)
  ),
  pgls_clade_residual_joint = list(
    formula = paste("pgls_model_c_mass_clade_residual ~", transposon_col, "+", chromatin_col),
    cols = c("pgls_model_c_mass_clade_residual", transposon_col, chromatin_col)
  )
)

rows <- list()
for (model_name in names(model_specs)) {
  rows[[length(rows) + 1]] <- fit_pgls_terms(
    model_name,
    model_specs[[model_name]]$formula,
    model_specs[[model_name]]$cols,
    c(transposon_col, chromatin_col)
  )
}
results <- do.call(rbind, rows)
idx <- which(is.finite(results$p))
results$p_bh <- NA
results$p_bh[idx] <- p.adjust(results$p[idx], method = "BH")
results <- results[order(results$model, results$term), ]

dir.create(dirname(table_path), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(diag_path), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(report_path), recursive = TRUE, showWarnings = FALSE)
write.table(results, table_path, sep = "\t", row.names = FALSE, quote = FALSE)
write.table(diag_rows, diag_path, sep = "\t", row.names = FALSE, quote = FALSE)

term_line <- function(model_name, term) {
  row <- results[results$model == model_name & results$term == term, ][1, ]
  paste0(
    "- ", model_name, " / ", term,
    ": estimate=", signif(row$estimate, 4),
    ", p=", signif(row$p, 4),
    ", BH=", signif(row$p_bh, 4),
    ", n=", row$n,
    ", lambda=", signif(row$lambda, 4)
  )
}

report_lines <- c(
  "# Phase 2 W2 Transposon-Chromatin Joint Models",
  "",
  "## Collinearity Diagnostics",
  "",
  paste0("- Pearson r between scores: ", signif(diag_rows$pearson_r_scores, 4)),
  paste0("- Spearman r between scores: ", signif(diag_rows$spearman_r_scores, 4)),
  paste0("- Pearson r between module coverage fractions: ", signif(diag_rows$pearson_r_coverage, 4)),
  paste0("- VIF transposon score adjusted for chromatin, body mass, and clade: ", signif(diag_rows$vif_transposon_score, 4)),
  paste0("- VIF chromatin score adjusted for transposon, body mass, and clade: ", signif(diag_rows$vif_chromatin_score, 4)),
  "",
  "## Joint Model Terms",
  "",
  term_line("mass_clade_joint", transposon_col),
  term_line("mass_clade_joint", chromatin_col),
  term_line("residual_joint", transposon_col),
  term_line("residual_joint", chromatin_col),
  term_line("pgls_clade_residual_joint", transposon_col),
  term_line("pgls_clade_residual_joint", chromatin_col),
  "",
  "## Interpretation",
  "",
  "These models explicitly ask whether transposon/repeat suppression and chromatin repression retain independent associations after entering the same model. If both terms attenuate, the current evidence supports a broader repeat/chromatin maintenance axis rather than a transposon-only claim."
)
writeLines(report_lines, report_path)

cat("Wrote ", table_path, ", ", diag_path, " and ", report_path, "\n", sep = "")
