# PGLS sensitivity across Week 4 maintenance-score variants.
#
# This tests whether module signals depend on rescued or lower-confidence
# ortholog candidates. It keeps the same feasibility-stage tree caveat as Week 3:
# OpenTree synthetic topology with Grafen branch lengths.

project_lib <- file.path(getwd(), "env", "R", "library")
if (dir.exists(project_lib)) {
  .libPaths(c(project_lib, .libPaths()))
}

library(ape)
library(caper)

args <- commandArgs(trailingOnly = TRUE)
data_path <- ifelse(length(args) >= 1, args[[1]], "data/processed/maintenance_lifespan_week4_variants.tsv")
table_path <- ifelse(length(args) >= 2, args[[2]], "results/tables/week4_score_variant_pgls.tsv")
report_path <- ifelse(length(args) >= 3, args[[3]], "results/reports/week4_score_variant_pgls_report.md")
tree_path <- ifelse(length(args) >= 4, args[[4]], "data/processed/phylogeny_inputs/opentree_induced_subtree.tre")

tree_full <- read.tree(tree_path)
data_all <- read.delim(data_path, stringsAsFactors = FALSE)

score_cols <- grep("_score$", names(data_all), value = TRUE)
module_names <- sub("_score$", "", score_cols)

data_all$log10_body_mass_g <- as.numeric(data_all$log10_body_mass_g)
data_all$log10_max_lifespan_years <- log10(as.numeric(data_all$max_lifespan_years))
data_all$lifespan_residual_log10 <- as.numeric(data_all$lifespan_residual_log10)
data_all$pgls_model_c_mass_clade_residual <- as.numeric(data_all$pgls_model_c_mass_clade_residual)
data_all$flight_status <- factor(data_all$flight_status)
data_all$clade <- factor(data_all$clade)

prepare_comp <- function(data, needed_cols) {
  cols <- unique(c("opentree_tip_label", needed_cols))
  data <- data[, cols]
  data <- data[complete.cases(data), ]
  rownames(data) <- data$opentree_tip_label
  shared <- intersect(tree_full$tip.label, data$opentree_tip_label)
  pruned <- drop.tip(tree_full, setdiff(tree_full$tip.label, shared))
  if (is.null(pruned$edge.length) || any(!is.finite(pruned$edge.length)) || any(pruned$edge.length <= 0)) {
    pruned <- compute.brlen(pruned, method = "Grafen")
  }
  data <- data[shared, ]
  comparative.data(
    phy = pruned,
    data = data,
    names.col = opentree_tip_label,
    vcv = TRUE,
    warn.dropped = FALSE
  )
}

extract_predictor <- function(model, predictor) {
  coef_table <- summary(model)$coefficients
  if (!predictor %in% rownames(coef_table)) {
    return(c(estimate = NA, p = NA, se = NA, t = NA))
  }
  c(
    estimate = coef_table[predictor, "Estimate"],
    p = coef_table[predictor, "Pr(>|t|)"],
    se = coef_table[predictor, "Std. Error"],
    t = coef_table[predictor, "t value"]
  )
}

fit_model <- function(data, variant, module, score_col, model_name, formula_text, needed_cols) {
  result <- tryCatch(
    {
      comp <- prepare_comp(data, needed_cols)
      model <- pgls(as.formula(formula_text), data = comp, lambda = "ML")
      term <- extract_predictor(model, score_col)
      data.frame(
        score_variant = variant,
        maintenance_module = module,
        score_col = score_col,
        model = model_name,
        formula = formula_text,
        n = length(model$residuals),
        lambda = as.numeric(model$param["lambda"]),
        logLik = as.numeric(logLik(model)),
        AIC = AIC(model),
        module_estimate = term["estimate"],
        module_se = term["se"],
        module_t = term["t"],
        module_p = term["p"],
        error = ""
      )
    },
    error = function(e) {
      data.frame(
        score_variant = variant,
        maintenance_module = module,
        score_col = score_col,
        model = model_name,
        formula = formula_text,
        n = NA,
        lambda = NA,
        logLik = NA,
        AIC = NA,
        module_estimate = NA,
        module_se = NA,
        module_t = NA,
        module_p = NA,
        error = conditionMessage(e)
      )
    }
  )
  result
}

rows <- list()

for (variant in unique(data_all$score_variant)) {
  data_variant <- data_all[data_all$score_variant == variant, ]
  data_variant$flight_status <- factor(data_variant$flight_status)
  data_variant$clade <- factor(data_variant$clade)

  for (idx in seq_along(score_cols)) {
    score_col <- score_cols[[idx]]
    module <- module_names[[idx]]
    data_variant[[score_col]] <- as.numeric(data_variant[[score_col]])

    specs <- list(
      mass_module = list(
        formula = paste("log10_max_lifespan_years ~ log10_body_mass_g +", score_col),
        cols = c("log10_max_lifespan_years", "log10_body_mass_g", score_col)
      ),
      mass_clade_module = list(
        formula = paste("log10_max_lifespan_years ~ log10_body_mass_g + clade +", score_col),
        cols = c("log10_max_lifespan_years", "log10_body_mass_g", "clade", score_col)
      ),
      residual_module = list(
        formula = paste("lifespan_residual_log10 ~", score_col),
        cols = c("lifespan_residual_log10", score_col)
      ),
      pgls_clade_residual_module = list(
        formula = paste("pgls_model_c_mass_clade_residual ~", score_col),
        cols = c("pgls_model_c_mass_clade_residual", score_col)
      )
    )

    for (model_name in names(specs)) {
      rows[[length(rows) + 1]] <- fit_model(
        data_variant,
        variant,
        module,
        score_col,
        model_name,
        specs[[model_name]]$formula,
        specs[[model_name]]$cols
      )
    }
  }
}

results <- do.call(rbind, rows)
results$module_p_bh_by_variant_model <- NA
for (variant in unique(results$score_variant)) {
  for (model_name in unique(results$model)) {
    idx <- which(
      results$score_variant == variant &
        results$model == model_name &
        is.finite(results$module_p)
    )
    results$module_p_bh_by_variant_model[idx] <- p.adjust(results$module_p[idx], method = "BH")
  }
}
results$module_p_bh_all <- NA
idx_all <- which(is.finite(results$module_p))
results$module_p_bh_all[idx_all] <- p.adjust(results$module_p[idx_all], method = "BH")
results <- results[order(results$module_p), ]

dir.create(dirname(table_path), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(report_path), recursive = TRUE, showWarnings = FALSE)
write.table(results, file = table_path, sep = "\t", row.names = FALSE, quote = FALSE)

top <- head(results[is.na(results$error) | results$error == "", ], 16)
variant_counts <- table(results$score_variant)
summary_lines <- c(
  "# Week 4 Score-Variant PGLS Report",
  "",
  paste0("PGLS sensitivity scan across score variants using tree `", tree_path, "`. Existing positive branch lengths are retained; Grafen lengths are used only as a fallback."),
  "",
  paste0("Rows fitted: ", nrow(results)),
  paste0("Score variants: ", length(unique(results$score_variant))),
  paste0("Modules tested: ", length(module_names)),
  "",
  "## Rows by Variant",
  paste0("- ", names(variant_counts), ": ", as.integer(variant_counts)),
  "",
  "## Top Module Terms",
  apply(top, 1, function(row) {
    paste0(
      "- ", row[["score_variant"]], " / ", row[["maintenance_module"]], " / ", row[["model"]],
      ": estimate=", signif(as.numeric(row[["module_estimate"]]), 4),
      ", p=", signif(as.numeric(row[["module_p"]]), 4),
      ", BH(variant-model)=", signif(as.numeric(row[["module_p_bh_by_variant_model"]]), 4),
      ", n=", row[["n"]]
    )
  }),
  "",
  "## Interpretation",
  "Signals that survive `ncbi_only` and `high_confidence_only` are less likely to be artifacts of annotation rescue. Signals that only appear in `all_validated` should be treated as rescue-sensitive feasibility leads."
)
writeLines(summary_lines, report_path)

cat("Wrote ", table_path, " and ", report_path, "\n", sep = "")
