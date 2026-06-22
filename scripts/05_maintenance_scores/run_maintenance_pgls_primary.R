# PGLS scan for primary maintenance module scores.
#
# Feasibility-stage caveat: this uses the OpenTree synthetic topology with
# Grafen branch lengths. Replace with dated clade-specific trees before final
# inference.

project_lib <- file.path(getwd(), "env", "R", "library")
if (dir.exists(project_lib)) {
  .libPaths(c(project_lib, .libPaths()))
}

library(ape)
library(caper)

tree_path <- "data/processed/phylogeny_inputs/opentree_induced_subtree.tre"
data_path <- "data/processed/maintenance_lifespan_primary.tsv"

tree_full <- read.tree(tree_path)
data_full <- read.delim(data_path, stringsAsFactors = FALSE)

score_cols <- grep("_score$", names(data_full), value = TRUE)
module_names <- sub("_score$", "", score_cols)

data_full$log10_body_mass_g <- as.numeric(data_full$log10_body_mass_g)
data_full$log10_max_lifespan_years <- log10(as.numeric(data_full$max_lifespan_years))
data_full$lifespan_residual_log10 <- as.numeric(data_full$lifespan_residual_log10)
data_full$pgls_model_b_mass_flight_residual <- as.numeric(data_full$pgls_model_b_mass_flight_residual)
data_full$pgls_model_c_mass_clade_residual <- as.numeric(data_full$pgls_model_c_mass_clade_residual)
data_full$flight_status <- factor(data_full$flight_status)
data_full$clade <- factor(data_full$clade)

prepare_comp <- function(data, needed_cols) {
  cols <- unique(c("opentree_tip_label", needed_cols))
  data <- data[, cols]
  data <- data[complete.cases(data), ]
  rownames(data) <- data$opentree_tip_label
  shared <- intersect(tree_full$tip.label, data$opentree_tip_label)
  pruned <- drop.tip(tree_full, setdiff(tree_full$tip.label, shared))
  pruned <- compute.brlen(pruned, method = "Grafen")
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

fit_model <- function(module, score_col, model_name, formula_text, needed_cols) {
  comp <- prepare_comp(data_full, needed_cols)
  result <- tryCatch(
    {
      model <- pgls(as.formula(formula_text), data = comp, lambda = "ML")
      term <- extract_predictor(model, score_col)
      data.frame(
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

for (idx in seq_along(score_cols)) {
  score_col <- score_cols[[idx]]
  module <- module_names[[idx]]
  data_full[[score_col]] <- as.numeric(data_full[[score_col]])

  specs <- list(
    mass_module = list(
      formula = paste("log10_max_lifespan_years ~ log10_body_mass_g +", score_col),
      cols = c("log10_max_lifespan_years", "log10_body_mass_g", score_col)
    ),
    mass_flight_module = list(
      formula = paste("log10_max_lifespan_years ~ log10_body_mass_g + flight_status +", score_col),
      cols = c("log10_max_lifespan_years", "log10_body_mass_g", "flight_status", score_col)
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
      module,
      score_col,
      model_name,
      specs[[model_name]]$formula,
      specs[[model_name]]$cols
    )
  }
}

results <- do.call(rbind, rows)
results$module_p_bh_by_model <- NA
for (model_name in unique(results$model)) {
  idx <- which(results$model == model_name & is.finite(results$module_p))
  results$module_p_bh_by_model[idx] <- p.adjust(results$module_p[idx], method = "BH")
}
results$module_p_bh_all <- NA
idx_all <- which(is.finite(results$module_p))
results$module_p_bh_all[idx_all] <- p.adjust(results$module_p[idx_all], method = "BH")
results <- results[order(results$module_p), ]

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/reports", recursive = TRUE, showWarnings = FALSE)
write.table(
  results,
  file = "results/tables/maintenance_pgls_primary.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

top <- head(results[is.na(results$error) | results$error == "", ], 10)
summary_lines <- c(
  "# Maintenance PGLS Primary Report",
  "",
  "PGLS scan for primary maintenance module scores using OpenTree synthetic topology and Grafen branch lengths.",
  "",
  paste0("Models fitted: ", nrow(results)),
  paste0("Modules tested: ", length(module_names)),
  "",
  "## Top Module Terms",
  apply(top, 1, function(row) {
    paste0(
      "- ", row[["maintenance_module"]], " / ", row[["model"]],
      ": estimate=", signif(as.numeric(row[["module_estimate"]]), 4),
      ", p=", signif(as.numeric(row[["module_p"]]), 4),
      ", BH(model)=", signif(as.numeric(row[["module_p_bh_by_model"]]), 4),
      ", n=", row[["n"]],
      ", lambda=", signif(as.numeric(row[["lambda"]]), 4)
    )
  }),
  "",
  "## Interpretation",
  "This is a feasibility-stage phylogenetic screen. Results are not final because the tree uses Grafen fallback branch lengths and the module scores are candidate-ortholog coverage scores, not direct molecular rate or functional measurements."
)
writeLines(summary_lines, "results/reports/maintenance_pgls_primary_report.md")

cat("Wrote results/tables/maintenance_pgls_primary.tsv and results/reports/maintenance_pgls_primary_report.md\n")
