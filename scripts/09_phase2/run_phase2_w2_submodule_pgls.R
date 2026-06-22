# PGLS for Phase 2 W2 repeat/chromatin submodule scores.

project_lib <- file.path(getwd(), "env", "R", "library")
if (dir.exists(project_lib)) {
  .libPaths(c(project_lib, .libPaths()))
}

library(ape)
library(caper)

args <- commandArgs(trailingOnly = TRUE)
data_path <- ifelse(length(args) >= 1, args[[1]], "data/processed/maintenance_lifespan_phase2_W2_submodules.tsv")
table_path <- ifelse(length(args) >= 2, args[[2]], "results/tables/phase2_W2_submodule_pgls.tsv")
report_path <- ifelse(length(args) >= 3, args[[3]], "results/reports/phase2_W2_submodule_pgls_report.md")
tree_path <- "data/processed/phylogeny_inputs/opentree_induced_subtree.tre"

tree_full <- read.tree(tree_path)
data <- read.delim(data_path, stringsAsFactors = FALSE)
score_cols <- grep("_score$", names(data), value = TRUE)

for (col in c("max_lifespan_years", "log10_body_mass_g", "pgls_model_c_mass_clade_residual", "lifespan_residual_log10", score_cols)) {
  data[[col]] <- as.numeric(data[[col]])
}
data$log10_max_lifespan_years <- log10(data$max_lifespan_years)
data$clade <- factor(data$clade)

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

fit_one <- function(submodule, score_col, model_name, formula_text, needed_cols) {
  result <- tryCatch(
    {
      comp <- prepare_comp(data, needed_cols)
      model <- pgls(as.formula(formula_text), data = comp, lambda = "ML")
      term <- extract_term(model, score_col)
      data.frame(
        submodule_v2 = submodule,
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
        submodule_v2 = submodule,
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
for (score_col in score_cols) {
  submodule <- sub("_score$", "", score_col)
  specs <- list(
    residual = list(
      formula = paste("pgls_model_c_mass_clade_residual ~", score_col),
      cols = c("pgls_model_c_mass_clade_residual", score_col)
    ),
    mass_clade = list(
      formula = paste("log10_max_lifespan_years ~ log10_body_mass_g + clade +", score_col),
      cols = c("log10_max_lifespan_years", "log10_body_mass_g", "clade", score_col)
    )
  )
  for (model_name in names(specs)) {
    rows[[length(rows) + 1]] <- fit_one(
      submodule,
      score_col,
      model_name,
      specs[[model_name]]$formula,
      specs[[model_name]]$cols
    )
  }
}

results <- do.call(rbind, rows)
idx <- which(is.finite(results$p))
results$p_bh <- NA
results$p_bh[idx] <- p.adjust(results$p[idx], method = "BH")
results <- results[order(results$p), ]

dir.create(dirname(table_path), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(report_path), recursive = TRUE, showWarnings = FALSE)
write.table(results, table_path, sep = "\t", row.names = FALSE, quote = FALSE)

top <- results[is.na(results$error) | results$error == "", ]
top <- head(top, 12)
lines <- c(
  "# Phase 2 W2 Submodule PGLS",
  "",
  "## Top Submodule Terms",
  ""
)
for (idx in seq_len(nrow(top))) {
  row <- top[idx, ]
  lines <- c(
    lines,
    paste0(
      "- ", row$submodule_v2, " / ", row$model,
      ": estimate=", signif(row$estimate, 4),
      ", p=", signif(row$p, 4),
      ", BH=", signif(row$p_bh, 4),
      ", n=", row$n,
      ", lambda=", signif(row$lambda, 4)
    )
  )
}
lines <- c(
  lines,
  "",
  "## Interpretation",
  "",
  "Submodule models are exploratory and help identify whether the W2 repeat/chromatin axis is driven by one biological subset or by shared annotation/coverage structure across several related subsets."
)
writeLines(lines, report_path)
cat("Wrote ", table_path, " and ", report_path, "\n", sep = "")
