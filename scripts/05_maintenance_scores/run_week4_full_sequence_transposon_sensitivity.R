# Sensitivity PGLS for full sequence-validated transposon scores.

project_lib <- file.path(getwd(), "env", "R", "library")
if (dir.exists(project_lib)) {
  .libPaths(c(project_lib, .libPaths()))
}

library(ape)
library(caper)

args <- commandArgs(trailingOnly = TRUE)
data_path <- ifelse(length(args) >= 1, args[[1]], "data/processed/maintenance_lifespan_week4_full_sequence_validated.tsv")
table_path <- ifelse(length(args) >= 2, args[[2]], "results/tables/week4_full_sequence_transposon_sensitivity.tsv")
report_path <- ifelse(length(args) >= 3, args[[3]], "results/reports/week4_full_sequence_transposon_sensitivity_report.md")
tree_path <- "data/processed/phylogeny_inputs/opentree_induced_subtree.tre"

tree_full <- read.tree(tree_path)
data_all <- read.delim(data_path, stringsAsFactors = FALSE)

data_all$log10_body_mass_g <- as.numeric(data_all$log10_body_mass_g)
data_all$log10_max_lifespan_years <- log10(as.numeric(data_all$max_lifespan_years))
data_all$lifespan_residual_log10 <- as.numeric(data_all$lifespan_residual_log10)
data_all$pgls_model_c_mass_clade_residual <- as.numeric(data_all$pgls_model_c_mass_clade_residual)
data_all$transposon_suppression_score <- as.numeric(data_all$transposon_suppression_score)
data_all$flight_status <- factor(data_all$flight_status)
data_all$clade <- factor(data_all$clade)

make_subset <- function(name, data) {
  if (name == "all_primary") return(data)
  if (name == "tier1_only") return(subset(data, genome_analysis_tier == "tier1_refseq_annotated_chromosome"))
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
  if (nrow(data) < 12) stop("too_few_complete_rows")
  if (length(unique(data$transposon_suppression_score)) < 2) stop("zero_score_variance")
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

fit_one <- function(data, variant, subset_name, model_name, formula_text, needed_cols) {
  result <- tryCatch(
    {
      comp <- prepare_comp(data, needed_cols)
      model <- pgls(as.formula(formula_text), data = comp, lambda = "ML")
      term <- extract_predictor(model, "transposon_suppression_score")
      data.frame(
        score_variant = variant,
        subset = subset_name,
        model = model_name,
        formula = formula_text,
        n = length(model$residuals),
        clade_count = length(unique(comp$data$clade)),
        lambda = as.numeric(model$param["lambda"]),
        estimate = term["estimate"],
        se = term["se"],
        t = term["t"],
        p = term["p"],
        error = ""
      )
    },
    error = function(e) {
      data.frame(
        score_variant = variant,
        subset = subset_name,
        model = model_name,
        formula = formula_text,
        n = nrow(data),
        clade_count = length(unique(data$clade)),
        lambda = NA,
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

subset_names <- c(
  "all_primary",
  "tier1_only",
  "exclude_human",
  "exclude_top_abs_residual_5",
  "leave_out_Aves",
  "leave_out_Mammalia_Chiroptera",
  "leave_out_Mammalia_nonChiroptera",
  "leave_out_Reptilia"
)
variants <- c("transposon_sequence_strict", "transposon_sequence_weak_inclusive")

rows <- list()
for (variant in variants) {
  variant_data <- subset(data_all, score_variant == variant)
  for (subset_name in subset_names) {
    subset_data <- make_subset(subset_name, variant_data)
    subset_data$flight_status <- factor(subset_data$flight_status)
    subset_data$clade <- factor(subset_data$clade)
    specs <- list(
      mass_module = list(
        formula = "log10_max_lifespan_years ~ log10_body_mass_g + transposon_suppression_score",
        cols = c("log10_max_lifespan_years", "log10_body_mass_g", "clade", "transposon_suppression_score")
      ),
      pgls_clade_residual_module = list(
        formula = "pgls_model_c_mass_clade_residual ~ transposon_suppression_score",
        cols = c("pgls_model_c_mass_clade_residual", "clade", "transposon_suppression_score")
      )
    )
    if (length(unique(subset_data$clade)) > 1) {
      specs$mass_clade_module <- list(
        formula = "log10_max_lifespan_years ~ log10_body_mass_g + clade + transposon_suppression_score",
        cols = c("log10_max_lifespan_years", "log10_body_mass_g", "clade", "transposon_suppression_score")
      )
    }
    for (model_name in names(specs)) {
      rows[[length(rows) + 1]] <- fit_one(
        subset_data,
        variant,
        subset_name,
        model_name,
        specs[[model_name]]$formula,
        specs[[model_name]]$cols
      )
    }
  }
}

results <- do.call(rbind, rows)
results$p_bh_by_variant_model <- NA
for (variant in unique(results$score_variant)) {
  for (model_name in unique(results$model)) {
    idx <- which(results$score_variant == variant & results$model == model_name & is.finite(results$p))
    results$p_bh_by_variant_model[idx] <- p.adjust(results$p[idx], method = "BH")
  }
}
results <- results[order(results$score_variant, results$model, results$p), ]

dir.create(dirname(table_path), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(report_path), recursive = TRUE, showWarnings = FALSE)
write.table(results, file = table_path, sep = "\t", row.names = FALSE, quote = FALSE)

mass_clade <- subset(results, model == "mass_clade_module" & score_variant == "transposon_sequence_strict")
mass_clade <- mass_clade[order(mass_clade$p), ]
summary_lines <- c(
  "# Week 4 Full Sequence Transposon Sensitivity Report",
  "",
  "Sensitivity PGLS for full sequence-validated transposon scores.",
  "",
  paste0("Rows fitted: ", nrow(results)),
  "",
  "## Strict Score, Mass + Clade",
  apply(mass_clade, 1, function(row) {
    paste0(
      "- ", row[["subset"]],
      ": estimate=", signif(as.numeric(row[["estimate"]]), 4),
      ", p=", signif(as.numeric(row[["p"]]), 4),
      ", BH=", signif(as.numeric(row[["p_bh_by_variant_model"]]), 4),
      ", n=", row[["n"]],
      ifelse(row[["error"]] == "", "", paste0(", error=", row[["error"]]))
    )
  }),
  "",
  "## Interpretation",
  "Subsets that keep a positive estimate support robustness; subsets that lose significance identify likely dependence on clade composition or annotation tier."
)
writeLines(summary_lines, report_path)
cat("Wrote ", table_path, " and ", report_path, "\n", sep = "")
