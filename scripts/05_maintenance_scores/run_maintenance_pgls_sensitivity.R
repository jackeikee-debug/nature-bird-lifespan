# Sensitivity PGLS scans for primary maintenance module scores.

project_lib <- file.path(getwd(), "env", "R", "library")
if (dir.exists(project_lib)) {
  .libPaths(c(project_lib, .libPaths()))
}

library(ape)
library(caper)

tree_path <- "data/processed/phylogeny_inputs/opentree_induced_subtree.tre"
data_path <- "data/processed/maintenance_lifespan_primary.tsv"
coverage_path <- "results/tables/ortholog_coverage_primary_diamond_validated_by_species.tsv"

tree_full <- read.tree(tree_path)
data_full <- read.delim(data_path, stringsAsFactors = FALSE)
coverage <- read.delim(coverage_path, stringsAsFactors = FALSE)

coverage$coverage_fraction <- as.numeric(coverage$coverage_fraction)
names(coverage)[names(coverage) == "coverage_fraction"] <- "overall_ortholog_coverage"
data_full <- merge(
  data_full,
  coverage[, c("scientific_name", "overall_ortholog_coverage")],
  by = "scientific_name",
  all.x = TRUE
)

score_cols <- grep("_score$", names(data_full), value = TRUE)
module_names <- sub("_score$", "", score_cols)

data_full$log10_body_mass_g <- as.numeric(data_full$log10_body_mass_g)
data_full$log10_max_lifespan_years <- log10(as.numeric(data_full$max_lifespan_years))
data_full$lifespan_residual_log10 <- as.numeric(data_full$lifespan_residual_log10)
data_full$pgls_model_c_mass_clade_residual <- as.numeric(data_full$pgls_model_c_mass_clade_residual)
data_full$flight_status <- factor(data_full$flight_status)
data_full$clade <- factor(data_full$clade)

make_subset <- function(name, data) {
  if (name == "all_primary") {
    return(data)
  }
  if (name == "tier1_only") {
    return(subset(data, genome_analysis_tier == "tier1_refseq_annotated_chromosome"))
  }
  if (name == "coverage_ge_0.75") {
    return(subset(data, overall_ortholog_coverage >= 0.75))
  }
  if (name == "coverage_ge_0.60") {
    return(subset(data, overall_ortholog_coverage >= 0.60))
  }
  if (name == "exclude_human") {
    return(subset(data, scientific_name != "Homo sapiens"))
  }
  if (name == "exclude_top_abs_residual_5") {
    data$abs_resid <- abs(data$lifespan_residual_log10)
    keep <- order(data$abs_resid, decreasing = TRUE)[-(1:5)]
    return(data[sort(keep), ])
  }
  if (startsWith(name, "leave_out_")) {
    clade <- sub("^leave_out_", "", name)
    return(subset(data, clade != clade))
  }
  stop(paste("Unknown subset", name))
}

# The base R scoping above is awkward for leave-one-clade, so handle explicitly.
make_subset <- function(name, data) {
  if (name == "all_primary") return(data)
  if (name == "tier1_only") return(subset(data, genome_analysis_tier == "tier1_refseq_annotated_chromosome"))
  if (name == "coverage_ge_0.75") return(subset(data, overall_ortholog_coverage >= 0.75))
  if (name == "coverage_ge_0.60") return(subset(data, overall_ortholog_coverage >= 0.60))
  if (name == "exclude_human") return(subset(data, scientific_name != "Homo sapiens"))
  if (name == "exclude_top_abs_residual_5") {
    data$abs_resid <- abs(data$lifespan_residual_log10)
    drop_idx <- order(data$abs_resid, decreasing = TRUE)[1:5]
    return(data[-drop_idx, ])
  }
  if (name == "leave_out_Aves") return(subset(data, clade != "Aves"))
  if (name == "leave_out_Mammalia_Chiroptera") return(subset(data, clade != "Mammalia_Chiroptera"))
  if (name == "leave_out_Mammalia_nonChiroptera") return(subset(data, clade != "Mammalia_nonChiroptera"))
  if (name == "leave_out_Reptilia") return(subset(data, clade != "Reptilia"))
  stop(paste("Unknown subset", name))
}

prepare_comp <- function(data, needed_cols) {
  cols <- unique(c("opentree_tip_label", needed_cols))
  data <- data[, cols]
  data <- data[complete.cases(data), ]
  if (nrow(data) < 12) {
    stop("too_few_complete_rows")
  }
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

fit_one <- function(subset_name, data, module, score_col, model_name, formula_text, needed_cols) {
  result <- tryCatch(
    {
      comp <- prepare_comp(data, needed_cols)
      model <- pgls(as.formula(formula_text), data = comp, lambda = "ML")
      term <- extract_predictor(model, score_col)
      data.frame(
        subset = subset_name,
        maintenance_module = module,
        score_col = score_col,
        model = model_name,
        formula = formula_text,
        n = length(model$residuals),
        clade_count = length(unique(comp$data$clade)),
        lambda = as.numeric(model$param["lambda"]),
        module_estimate = term["estimate"],
        module_se = term["se"],
        module_t = term["t"],
        module_p = term["p"],
        error = ""
      )
    },
    error = function(e) {
      data.frame(
        subset = subset_name,
        maintenance_module = module,
        score_col = score_col,
        model = model_name,
        formula = formula_text,
        n = nrow(data),
        clade_count = length(unique(data$clade)),
        lambda = NA,
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

subset_names <- c(
  "all_primary",
  "tier1_only",
  "coverage_ge_0.75",
  "coverage_ge_0.60",
  "exclude_human",
  "exclude_top_abs_residual_5",
  "leave_out_Aves",
  "leave_out_Mammalia_Chiroptera",
  "leave_out_Mammalia_nonChiroptera",
  "leave_out_Reptilia"
)

rows <- list()

for (subset_name in subset_names) {
  subset_data <- make_subset(subset_name, data_full)
  subset_data$flight_status <- factor(subset_data$flight_status)
  subset_data$clade <- factor(subset_data$clade)
  for (idx in seq_along(score_cols)) {
    score_col <- score_cols[[idx]]
    module <- module_names[[idx]]
    subset_data[[score_col]] <- as.numeric(subset_data[[score_col]])

    # Always include mass and module. Add clade only where at least two clades remain.
    rows[[length(rows) + 1]] <- fit_one(
      subset_name,
      subset_data,
      module,
      score_col,
      "mass_module",
      paste("log10_max_lifespan_years ~ log10_body_mass_g +", score_col),
      c("log10_max_lifespan_years", "log10_body_mass_g", "clade", score_col)
    )

    if (length(unique(subset_data$clade)) > 1) {
      rows[[length(rows) + 1]] <- fit_one(
        subset_name,
        subset_data,
        module,
        score_col,
        "mass_clade_module",
        paste("log10_max_lifespan_years ~ log10_body_mass_g + clade +", score_col),
        c("log10_max_lifespan_years", "log10_body_mass_g", "clade", score_col)
      )
    }

    rows[[length(rows) + 1]] <- fit_one(
      subset_name,
      subset_data,
      module,
      score_col,
      "pgls_clade_residual_module",
      paste("pgls_model_c_mass_clade_residual ~", score_col),
      c("pgls_model_c_mass_clade_residual", "clade", score_col)
    )
  }
}

results <- do.call(rbind, rows)
results$module_p_bh_by_subset_model <- NA
for (subset_name in unique(results$subset)) {
  for (model_name in unique(results$model)) {
    idx <- which(results$subset == subset_name & results$model == model_name & is.finite(results$module_p))
    results$module_p_bh_by_subset_model[idx] <- p.adjust(results$module_p[idx], method = "BH")
  }
}
results <- results[order(results$subset, results$model, results$module_p), ]

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/reports", recursive = TRUE, showWarnings = FALSE)
write.table(
  results,
  file = "results/tables/maintenance_pgls_sensitivity.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

trans <- subset(results, maintenance_module == "transposon_suppression" & model == "mass_clade_module")
trans <- trans[order(trans$module_p), ]

summary_lines <- c(
  "# Maintenance PGLS Sensitivity Report",
  "",
  "Sensitivity scan for primary maintenance module scores using OpenTree synthetic topology and Grafen branch lengths.",
  "",
  paste0("Rows: ", nrow(results)),
  "",
  "## Subsets",
  paste0("- ", subset_names),
  "",
  "## Transposon Suppression, Mass + Clade Model",
  apply(trans, 1, function(row) {
    paste0(
      "- ", row[["subset"]],
      ": estimate=", signif(as.numeric(row[["module_estimate"]]), 4),
      ", p=", signif(as.numeric(row[["module_p"]]), 4),
      ", BH=", signif(as.numeric(row[["module_p_bh_by_subset_model"]]), 4),
      ", n=", row[["n"]],
      ", lambda=", signif(as.numeric(row[["lambda"]]), 4),
      ifelse(row[["error"]] == "", "", paste0(", error=", row[["error"]]))
    )
  }),
  "",
  "## Caveat",
  "This is a feasibility sensitivity check, not final inference. Leave-one-clade subsets can become underpowered and still use Grafen fallback branch lengths."
)
writeLines(summary_lines, "results/reports/maintenance_pgls_sensitivity_report.md")

cat("Wrote results/tables/maintenance_pgls_sensitivity.tsv and results/reports/maintenance_pgls_sensitivity_report.md\n")
