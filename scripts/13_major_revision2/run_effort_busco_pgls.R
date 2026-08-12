# Study effort, genome completeness, and shared module-axis PGLS analyses.

project_lib <- file.path(getwd(), "env", "R_library")
if (dir.exists(project_lib)) .libPaths(c(project_lib, .libPaths()))

suppressPackageStartupMessages({
  library(ape)
  library(caper)
})

args <- commandArgs(trailingOnly = TRUE)
data_path <- ifelse(length(args) >= 1, args[[1]], "data/processed/jme_revision2_analysis_table.tsv")
tree_path <- ifelse(length(args) >= 2, args[[2]], "data/processed/phylogeny_inputs/opentree_datelife_calibrated_primary68.tre")
model_output <- ifelse(length(args) >= 3, args[[3]], "results/tables/jme_revision2_effort_busco_pgls.tsv")
cor_output <- ifelse(length(args) >= 4, args[[4]], "results/tables/jme_revision2_module_correlations.tsv")
pc_output <- ifelse(length(args) >= 5, args[[5]], "results/tables/jme_revision2_shared_axis_species.tsv")
report_output <- ifelse(length(args) >= 6, args[[6]], "results/reports/jme_revision2_effort_busco_pgls.md")

data <- read.delim(data_path, stringsAsFactors = FALSE, check.names = FALSE)
tree <- read.tree(tree_path)
score_cols <- grep("_score$", names(data), value = TRUE)
module_names <- sub("_score$", "", score_cols)

numeric_cols <- unique(c(
  "max_lifespan_years", "log10_body_mass_g", "pgls_model_c_mass_clade_residual",
  "busco_complete", "tier_numeric", "contig_n50", "scaffold_n50", "coverage",
  "log1p_anage_reference_count", "anage_sample_size_rank", "anage_data_quality_rank",
  "log1p_pubmed_exact_binomial_count", score_cols
))
for (column in intersect(numeric_cols, names(data))) data[[column]] <- as.numeric(data[[column]])
data$log10_max_lifespan_years <- log10(data$max_lifespan_years)
data$log10_contig_n50 <- log10(data$contig_n50)
data$clade <- factor(data$clade)

# The six module scores are strongly correlated. PC1 is reported as a shared
# gene-observability axis, with sign fixed so larger values indicate larger mean scores.
scaled_scores <- scale(data[, score_cols])
pca <- prcomp(scaled_scores, center = FALSE, scale. = FALSE)
pc1 <- pca$x[, 1]
flip_pc1 <- cor(pc1, rowMeans(scaled_scores), use = "complete.obs") < 0
if (flip_pc1) pc1 <- -pc1
data$shared_module_pc1 <- pc1
loadings <- pca$rotation[, 1]
if (flip_pc1) loadings <- -loadings

cor_matrix <- cor(data[, score_cols], method = "spearman", use = "pairwise.complete.obs")
cor_rows <- do.call(rbind, lapply(seq_along(score_cols), function(i) {
  do.call(rbind, lapply(seq_along(score_cols), function(j) {
    data.frame(
      module_1 = module_names[[i]], module_2 = module_names[[j]],
      spearman_rho = cor_matrix[i, j], stringsAsFactors = FALSE
    )
  }))
}))
cor_rows$pc1_loading <- NA_real_
for (i in seq_along(score_cols)) {
  idx <- cor_rows$module_1 == module_names[[i]] & cor_rows$module_2 == module_names[[i]]
  cor_rows$pc1_loading[idx] <- loadings[[i]]
}
pc_species <- data[, c("scientific_name", "clade", score_cols)]
pc_species$shared_module_pc1 <- data$shared_module_pc1

prepare_comp <- function(input, needed) {
  columns <- unique(c("opentree_tip_label", needed))
  x <- input[, columns, drop = FALSE]
  x <- x[complete.cases(x), , drop = FALSE]
  x <- x[!duplicated(x$opentree_tip_label), , drop = FALSE]
  shared <- intersect(tree$tip.label, x$opentree_tip_label)
  if (length(shared) < 20) stop("fewer than 20 complete tree-matched species")
  phy <- drop.tip(tree, setdiff(tree$tip.label, shared))
  rownames(x) <- x$opentree_tip_label
  x <- x[phy$tip.label, , drop = FALSE]
  x$clade <- droplevels(factor(x$clade))
  list(phy = phy, data = x)
}

fit_one <- function(predictor, model_name, formula_text, needed) {
  tryCatch({
    prepared <- prepare_comp(data, needed)
    x <- prepared$data
    standardize <- intersect(
      c(predictor, "busco_complete", "tier_numeric", "log10_contig_n50",
        "log1p_anage_reference_count", "log1p_pubmed_exact_binomial_count",
        "anage_sample_size_rank", "anage_data_quality_rank"),
      names(x)
    )
    for (column in standardize) {
      if (sd(x[[column]], na.rm = TRUE) > 0) x[[paste0(column, "_z")]] <- as.numeric(scale(x[[column]]))
    }
    comp <- comparative.data(phy = prepared$phy, data = x, names.col = opentree_tip_label, vcv = TRUE, warn.dropped = FALSE)
    model <- pgls(as.formula(formula_text), data = comp, lambda = "ML")
    term_name <- paste0(predictor, "_z")
    term <- summary(model)$coefficients[term_name, ]
    data.frame(
      predictor = predictor, model = model_name, formula = formula_text,
      n = length(model$residuals), lambda = as.numeric(model$param["lambda"]),
      estimate_per_predictor_sd = unname(term["Estimate"]),
      se = unname(term["Std. Error"]),
      conf_low = unname(term["Estimate"] - 1.96 * term["Std. Error"]),
      conf_high = unname(term["Estimate"] + 1.96 * term["Std. Error"]),
      t = unname(term["t value"]), p = unname(term["Pr(>|t|)"]),
      AIC = AIC(model), error = "", stringsAsFactors = FALSE
    )
  }, error = function(e) {
    data.frame(
      predictor = predictor, model = model_name, formula = formula_text,
      n = NA, lambda = NA, estimate_per_predictor_sd = NA, se = NA,
      conf_low = NA, conf_high = NA, t = NA, p = NA, AIC = NA,
      error = conditionMessage(e), stringsAsFactors = FALSE
    )
  })
}

specs <- list(
  base = list(
    suffix = "",
    needed = c("log10_max_lifespan_years", "log10_body_mass_g", "clade")
  ),
  anage_effort = list(
    suffix = "+ log1p_anage_reference_count_z",
    needed = c("log10_max_lifespan_years", "log10_body_mass_g", "clade", "log1p_anage_reference_count")
  ),
  pubmed_effort = list(
    suffix = "+ log1p_pubmed_exact_binomial_count_z",
    needed = c("log10_max_lifespan_years", "log10_body_mass_g", "clade", "log1p_pubmed_exact_binomial_count")
  ),
  dual_effort = list(
    suffix = "+ log1p_anage_reference_count_z + log1p_pubmed_exact_binomial_count_z",
    needed = c("log10_max_lifespan_years", "log10_body_mass_g", "clade", "log1p_anage_reference_count", "log1p_pubmed_exact_binomial_count")
  ),
  busco_complete_case = list(
    suffix = "+ busco_complete_z",
    needed = c("log10_max_lifespan_years", "log10_body_mass_g", "clade", "busco_complete")
  ),
  busco_and_dual_effort = list(
    suffix = "+ busco_complete_z + log1p_anage_reference_count_z + log1p_pubmed_exact_binomial_count_z",
    needed = c("log10_max_lifespan_years", "log10_body_mass_g", "clade", "busco_complete", "log1p_anage_reference_count", "log1p_pubmed_exact_binomial_count")
  ),
  assembly_and_dual_effort = list(
    suffix = "+ tier_numeric_z + log10_contig_n50_z + log1p_anage_reference_count_z + log1p_pubmed_exact_binomial_count_z",
    needed = c("log10_max_lifespan_years", "log10_body_mass_g", "clade", "tier_numeric", "log10_contig_n50", "log1p_anage_reference_count", "log1p_pubmed_exact_binomial_count")
  )
)

predictors <- c(score_cols, "shared_module_pc1")
rows <- list()
for (predictor in predictors) {
  for (model_name in names(specs)) {
    specification <- specs[[model_name]]
    formula_text <- paste0(
      "log10_max_lifespan_years ~ log10_body_mass_g + clade + ",
      predictor, "_z ", specification$suffix
    )
    needed <- unique(c(specification$needed, predictor))
    rows[[length(rows) + 1]] <- fit_one(predictor, model_name, formula_text, needed)
  }
}
models <- do.call(rbind, rows)
models$q_within_model <- NA_real_
for (model_name in unique(models$model)) {
  idx <- models$model == model_name & is.finite(models$p) & models$predictor %in% score_cols
  models$q_within_model[idx] <- p.adjust(models$p[idx], method = "BH")
}

dir.create(dirname(model_output), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(report_output), recursive = TRUE, showWarnings = FALSE)
write.table(models, model_output, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
write.table(cor_rows, cor_output, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
write.table(pc_species, pc_output, sep = "\t", row.names = FALSE, quote = FALSE, na = "")

off_diag <- cor_matrix[upper.tri(cor_matrix)]
fmt <- function(x) ifelse(is.finite(x), format(signif(x, 4), scientific = FALSE), "NA")
report <- c(
  "# Study Effort, BUSCO, and Shared Module Axis",
  "",
  paste0("The six module scores had pairwise Spearman correlations from ", fmt(min(off_diag)), " to ", fmt(max(off_diag)), "."),
  paste0("PC1 explained ", fmt(100 * summary(pca)$importance[2, 1]), "% of their standardized variance."),
  paste0("BUSCO completeness was available for ", sum(is.finite(data$busco_complete)), " of ", nrow(data), " species."),
  "",
  "## Transposon/repeat score",
  ""
)
focal <- models[models$predictor == "transposon_repeat_suppression_score", ]
for (i in seq_len(nrow(focal))) {
  row <- focal[i, ]
  report <- c(report, paste0(
    "- ", row$model, ": beta per score SD = ", fmt(row$estimate_per_predictor_sd),
    " (95% CI ", fmt(row$conf_low), " to ", fmt(row$conf_high), "), P = ",
    fmt(row$p), ", n = ", row$n, ", lambda = ", fmt(row$lambda), "."
  ))
}
report <- c(report, "", "## Shared PC1", "")
shared <- models[models$predictor == "shared_module_pc1", ]
for (i in seq_len(nrow(shared))) {
  row <- shared[i, ]
  report <- c(report, paste0(
    "- ", row$model, ": beta per PC1 SD = ", fmt(row$estimate_per_predictor_sd),
    " (95% CI ", fmt(row$conf_low), " to ", fmt(row$conf_high), "), P = ",
    fmt(row$p), ", n = ", row$n, "."
  ))
}
writeLines(report, report_output)
cat("Wrote study-effort, BUSCO, and shared-axis outputs.\n")
