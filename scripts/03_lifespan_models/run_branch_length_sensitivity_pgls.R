#!/usr/bin/env Rscript

# Branch-length sensitivity scan for the transposon/repeat suppression PGLS signal.
# The topology is the OpenTree synthetic subtree; branch lengths are deliberately
# varied to test whether the effect depends on one arbitrary branch-length choice.

suppressPackageStartupMessages({
  project_libs <- c(file.path(getwd(), "env", "R_library"), file.path(getwd(), "env", "R", "library"))
  project_libs <- project_libs[dir.exists(project_libs)]
  if (length(project_libs) > 0) .libPaths(c(project_libs, .libPaths()))
  library(ape)
  library(caper)
  library(ggplot2)
  library(readr)
  library(dplyr)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)
data_path <- ifelse(length(args) >= 1, args[[1]], "data/processed/maintenance_lifespan_phase3_gff_cds_plus_uniprot_sequence_rescued.tsv")
tree_path <- ifelse(length(args) >= 2, args[[2]], "data/processed/phylogeny_inputs/opentree_induced_subtree.tre")
table_path <- ifelse(length(args) >= 3, args[[3]], "results/tables/branch_length_sensitivity_pgls.tsv")
figure_prefix <- ifelse(length(args) >= 4, args[[4]], "results/figures/branch_length_sensitivity_forest")
report_path <- ifelse(length(args) >= 5, args[[5]], "results/reports/branch_length_sensitivity_pgls_report.md")
dated_tree_path <- ifelse(length(args) >= 6, args[[6]], "data/processed/phylogeny_inputs/opentree_datelife_calibrated_primary68.tre")

score_col <- "transposon_repeat_suppression_score"
root_age_ma <- 320

tree_full <- read.tree(tree_path)
dated_tree_full <- if (file.exists(dated_tree_path)) read.tree(dated_tree_path) else NULL
data_all <- read.delim(data_path, stringsAsFactors = FALSE)

data_all$log10_body_mass_g <- as.numeric(data_all$log10_body_mass_g)
data_all$log10_max_lifespan_years <- log10(as.numeric(data_all$max_lifespan_years))
data_all$lifespan_residual_log10 <- as.numeric(data_all$lifespan_residual_log10)
data_all$pgls_model_c_mass_clade_residual <- as.numeric(data_all$pgls_model_c_mass_clade_residual)
data_all[[score_col]] <- as.numeric(data_all[[score_col]])
data_all$clade <- factor(data_all$clade)
data_all$flight_status <- factor(data_all$flight_status)

branch_variants <- c(
  "opentree_original_or_unit",
  "datelife_calibrated",
  "equal_branch_lengths",
  "grafen_power_0.5",
  "grafen_power_1",
  "grafen_power_2",
  "grafen_scaled_320Ma"
)

primary_score_variants <- c(
  "phase2_W3_full_background_sensitivity",
  "phase2_strict48_sequence_updated"
)

make_branch_tree <- function(base_tree, labels, branch_variant) {
  if (branch_variant == "datelife_calibrated") {
    if (is.null(dated_tree_full)) stop("dated_tree_missing")
    if (!all(labels %in% dated_tree_full$tip.label)) stop("dated_tree_tip_mismatch")
    tr <- drop.tip(dated_tree_full, setdiff(dated_tree_full$tip.label, labels))
    if (is.null(tr$edge.length) || any(!is.finite(tr$edge.length)) || any(tr$edge.length <= 0)) {
      stop("dated_tree_invalid_branch_lengths")
    }
    return(ladderize(tr))
  }
  tr <- drop.tip(base_tree, setdiff(base_tree$tip.label, labels))
  tr <- ladderize(tr)

  if (branch_variant == "opentree_original_or_unit") {
    if (is.null(tr$edge.length) || any(!is.finite(tr$edge.length)) || sum(tr$edge.length, na.rm = TRUE) <= 0) {
      tr$edge.length <- rep(1, nrow(tr$edge))
    }
    return(tr)
  }
  if (branch_variant == "equal_branch_lengths") {
    tr$edge.length <- rep(1, nrow(tr$edge))
    return(tr)
  }
  if (branch_variant == "grafen_power_0.5") {
    return(compute.brlen(tr, method = "Grafen", power = 0.5))
  }
  if (branch_variant == "grafen_power_1") {
    return(compute.brlen(tr, method = "Grafen", power = 1))
  }
  if (branch_variant == "grafen_power_2") {
    return(compute.brlen(tr, method = "Grafen", power = 2))
  }
  if (branch_variant == "grafen_scaled_320Ma") {
    tr <- compute.brlen(tr, method = "Grafen", power = 1)
    depths <- node.depth.edgelength(tr)
    max_tip_depth <- max(depths[seq_along(tr$tip.label)], na.rm = TRUE)
    tr$edge.length <- tr$edge.length * (root_age_ma / max_tip_depth)
    return(tr)
  }
  stop(paste("Unknown branch variant:", branch_variant))
}

prepare_comp <- function(data, branch_variant, needed_cols) {
  cols <- unique(c("opentree_tip_label", needed_cols))
  data <- data[, cols]
  data <- data[complete.cases(data), ]
  if (nrow(data) < 12) stop("too_few_complete_rows")
  if (length(unique(data[[score_col]])) < 2) stop("zero_score_variance")
  shared <- intersect(tree_full$tip.label, data$opentree_tip_label)
  if (length(shared) < 12) stop("too_few_tree_matched_rows")
  tr <- make_branch_tree(tree_full, shared, branch_variant)
  data <- data[match(tr$tip.label, data$opentree_tip_label), ]
  rownames(data) <- data$opentree_tip_label
  comparative.data(
    phy = tr,
    data = data,
    names.col = opentree_tip_label,
    vcv = TRUE,
    warn.dropped = FALSE
  )
}

extract_predictor <- function(model, predictor) {
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

fit_one <- function(data, score_variant, branch_variant, model_name, formula_text, needed_cols) {
  tryCatch(
    {
      comp <- prepare_comp(data, branch_variant, needed_cols)
      model <- pgls(as.formula(formula_text), data = comp, lambda = "ML")
      term <- extract_predictor(model, score_col)
      data.frame(
        score_variant = score_variant,
        branch_variant = branch_variant,
        model = model_name,
        formula = formula_text,
        n = length(model$residuals),
        clade_count = length(unique(comp$data$clade)),
        lambda = as.numeric(model$param["lambda"]),
        logLik = as.numeric(logLik(model)),
        AIC = AIC(model),
        estimate = term["estimate"],
        se = term["se"],
        t = term["t"],
        p = term["p"],
        conf_low = term["estimate"] - 1.96 * term["se"],
        conf_high = term["estimate"] + 1.96 * term["se"],
        error = ""
      )
    },
    error = function(e) {
      data.frame(
        score_variant = score_variant,
        branch_variant = branch_variant,
        model = model_name,
        formula = formula_text,
        n = nrow(data),
        clade_count = length(unique(data$clade)),
        lambda = NA,
        logLik = NA,
        AIC = NA,
        estimate = NA,
        se = NA,
        t = NA,
        p = NA,
        conf_low = NA,
        conf_high = NA,
        error = conditionMessage(e)
      )
    }
  )
}

model_specs <- list(
  mass_score = list(
    formula = paste("log10_max_lifespan_years ~ log10_body_mass_g +", score_col),
    cols = c("log10_max_lifespan_years", "log10_body_mass_g", score_col)
  ),
  mass_clade_score = list(
    formula = paste("log10_max_lifespan_years ~ log10_body_mass_g + clade +", score_col),
    cols = c("log10_max_lifespan_years", "log10_body_mass_g", "clade", score_col)
  ),
  clade_residual_score = list(
    formula = paste("pgls_model_c_mass_clade_residual ~", score_col),
    cols = c("pgls_model_c_mass_clade_residual", score_col)
  )
)

rows <- list()
for (score_variant in unique(data_all$score_variant)) {
  variant_data <- data_all[data_all$score_variant == score_variant, ]
  variant_data$clade <- factor(variant_data$clade)
  for (branch_variant in branch_variants) {
    for (model_name in names(model_specs)) {
      rows[[length(rows) + 1]] <- fit_one(
        variant_data,
        score_variant,
        branch_variant,
        model_name,
        model_specs[[model_name]]$formula,
        model_specs[[model_name]]$cols
      )
    }
  }
}

results <- do.call(rbind, rows)
results$p_bh_by_score_model <- NA
for (score_variant in unique(results$score_variant)) {
  for (model_name in unique(results$model)) {
    idx <- which(results$score_variant == score_variant & results$model == model_name & is.finite(results$p))
    results$p_bh_by_score_model[idx] <- p.adjust(results$p[idx], method = "BH")
  }
}
results <- results[order(results$score_variant, results$model, results$branch_variant), ]

dir.create(dirname(table_path), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(figure_prefix), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(report_path), recursive = TRUE, showWarnings = FALSE)
write.table(results, file = table_path, sep = "\t", row.names = FALSE, quote = FALSE)

plot_df <- results |>
  filter(
    score_variant %in% primary_score_variants,
    model %in% c("mass_clade_score", "clade_residual_score"),
    is.finite(estimate),
    is.finite(se)
  ) |>
  mutate(
    branch_variant = factor(branch_variant, levels = rev(branch_variants)),
    model_label = recode(
      model,
      mass_clade_score = "Lifespan ~ mass + clade + score",
      clade_residual_score = "Clade-adjusted residual ~ score"
    ),
    score_label = recode(
      score_variant,
      phase2_W3_full_background_sensitivity = "Full background score",
      phase2_strict48_sequence_updated = "Strict48 sequence-updated score"
    )
  )

forest <- ggplot(plot_df, aes(x = estimate, y = branch_variant, color = score_label)) +
  geom_vline(xintercept = 0, linewidth = 0.35, color = "#94A3B8") +
  geom_errorbarh(aes(xmin = conf_low, xmax = conf_high), height = 0.22, linewidth = 0.55, alpha = 0.85, position = position_dodge(width = 0.45)) +
  geom_point(size = 2.0, position = position_dodge(width = 0.45)) +
  facet_wrap(~model_label, scales = "free_x") +
  scale_color_manual(values = c("Full background score" = "#0F766E", "Strict48 sequence-updated score" = "#9F1239")) +
  scale_y_discrete(labels = c(
    opentree_original_or_unit = "OpenTree original / unit",
    datelife_calibrated = "DateLife calibrated",
    equal_branch_lengths = "Equal branch lengths",
    grafen_power_0.5 = "Grafen power 0.5",
    grafen_power_1 = "Grafen power 1",
    grafen_power_2 = "Grafen power 2",
    grafen_scaled_320Ma = "Grafen scaled to 320 Ma"
  )) +
  labs(
    x = "PGLS estimate for transposon/repeat suppression score",
    y = "Branch-length assumption",
    color = "Score variant",
    title = "Branch-length sensitivity of transposon/repeat suppression PGLS signal",
    subtitle = "OpenTree topology with DateLife-calibrated, unit, equal, and alternative Grafen branch lengths"
  ) +
  theme_minimal(base_size = 9) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 11),
    strip.text = element_text(face = "bold")
  )

ggsave(paste0(figure_prefix, ".png"), forest, width = 10.4, height = 5.6, dpi = 320, bg = "white")
ggsave(paste0(figure_prefix, ".pdf"), forest, width = 10.4, height = 5.6, bg = "white")
ggsave(paste0(figure_prefix, ".svg"), forest, width = 10.4, height = 5.6, bg = "white")

report_core <- results |>
  filter(score_variant %in% primary_score_variants, model == "mass_clade_score") |>
  arrange(score_variant, branch_variant)

report_lines <- c(
  "# Branch-Length Sensitivity PGLS Report",
  "",
  "This scan keeps the OpenTree synthetic topology fixed and varies branch-length assumptions.",
  "It is designed to flag whether the transposon/repeat suppression association depends on an arbitrary branch-length transform.",
  "",
  "## Outputs",
  "",
  paste0("- table: `", table_path, "`"),
  paste0("- forest PNG: `", figure_prefix, ".png`"),
  paste0("- forest PDF: `", figure_prefix, ".pdf`"),
  paste0("- forest SVG: `", figure_prefix, ".svg`"),
  "",
  "## Branch-Length Variants",
  "",
  "- `datelife_calibrated`: fixed OpenTree topology with node ages summarized from peer-reviewed chronograms in the DateLife cache; uncovered nodes use an audited Grafen fallback.",
  "- `opentree_original_or_unit`: OpenTree branch lengths if usable; unit lengths otherwise.",
  "- `equal_branch_lengths`: every edge set to 1.",
  "- `grafen_power_0.5`, `grafen_power_1`, `grafen_power_2`: Grafen branch lengths with different powers.",
  paste0("- `grafen_scaled_320Ma`: Grafen power 1 scaled to a ", root_age_ma, " Ma amniote root; global scaling should often be similar to unscaled Grafen under lambda-ML PGLS."),
  "",
  "## Primary Mass + Clade Results",
  "",
  apply(report_core, 1, function(row) {
    paste0(
      "- ", row[["score_variant"]], " / ", row[["branch_variant"]],
      ": estimate=", signif(as.numeric(row[["estimate"]]), 4),
      ", 95% CI=", signif(as.numeric(row[["conf_low"]]), 4), " to ", signif(as.numeric(row[["conf_high"]]), 4),
      ", p=", signif(as.numeric(row[["p"]]), 4),
      ", lambda=", signif(as.numeric(row[["lambda"]]), 4),
      ", n=", row[["n"]],
      ifelse(row[["error"]] == "", "", paste0(", error=", row[["error"]]))
    )
  }),
  "",
  "## Interpretation Guardrail",
  "",
  "Agreement in effect direction across branch-length variants supports robustness to branch-length specification. The DateLife-calibrated tree is a secondary chronogram on the fixed OpenTree topology, not a newly inferred sequence-based species tree."
)
writeLines(report_lines, report_path)

cat("Wrote ", table_path, ", ", figure_prefix, ".[png/pdf/svg], and ", report_path, "\n", sep = "")
