#!/usr/bin/env Rscript

project_lib <- file.path(getwd(), "env", "R", "library")
if (dir.exists(project_lib)) {
  .libPaths(c(project_lib, .libPaths()))
}

library(ape)
library(caper)

args <- commandArgs(trailingOnly = TRUE)
data_path <- ifelse(length(args) >= 1, args[[1]], "data/processed/maintenance_lifespan_phase2_W3_full_background_expanded.tsv")
subset_table_path <- ifelse(length(args) >= 2, args[[2]], "results/tables/phase3_high_coverage_subset_pgls.tsv")
bird_model_path <- ifelse(length(args) >= 3, args[[3]], "results/tables/phase3_bird_enriched_models.tsv")
species_table_path <- ifelse(length(args) >= 4, args[[4]], "data/processed/phase3_high_coverage_subset_species.tsv")
report_path <- ifelse(length(args) >= 5, args[[5]], "results/reports/phase3_high_coverage_bird_models_report.md")
tree_path <- "data/processed/phylogeny_inputs/opentree_induced_subtree.tre"
variant_name <- "phase2_W3_full_background_sensitivity"

tree_full <- read.tree(tree_path)
data_all <- read.delim(data_path, stringsAsFactors = FALSE, check.names = FALSE)
data <- data_all[data_all$score_variant == variant_name, ]

score_col <- "transposon_repeat_suppression_score"
coverage_col <- "transposon_repeat_suppression_coverage"
needed_numeric <- c(
  "max_lifespan_years",
  "log10_body_mass_g",
  "lifespan_residual_log10",
  "pgls_model_c_mass_clade_residual",
  score_col,
  coverage_col
)
for (col in needed_numeric) {
  data[[col]] <- as.numeric(data[[col]])
}
data$log10_max_lifespan_years <- log10(data$max_lifespan_years)
data$clade <- factor(data$clade)
data$flight_status <- factor(data$flight_status)
data$is_bird <- ifelse(data$clade == "Aves", 1, 0)
data$bird_status <- factor(ifelse(data$clade == "Aves", "bird", "non_bird"), levels = c("non_bird", "bird"))

prepare_comp <- function(data, needed_cols) {
  cols <- unique(c("opentree_tip_label", needed_cols))
  model_data <- data[, cols]
  model_data <- model_data[complete.cases(model_data), ]
  if (nrow(model_data) < 12) stop("too_few_complete_rows")
  rownames(model_data) <- model_data$opentree_tip_label
  shared <- intersect(tree_full$tip.label, model_data$opentree_tip_label)
  if (length(shared) < 12) stop("too_few_tree_matched_rows")
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

extract_term <- function(model, term_name) {
  coef_table <- summary(model)$coefficients
  if (!term_name %in% rownames(coef_table)) {
    return(c(estimate = NA, se = NA, t = NA, p = NA))
  }
  c(
    estimate = coef_table[term_name, "Estimate"],
    se = coef_table[term_name, "Std. Error"],
    t = coef_table[term_name, "t value"],
    p = coef_table[term_name, "Pr(>|t|)"]
  )
}

fit_term <- function(label, subset_data, model_name, formula_text, term_name, needed_cols) {
  tryCatch(
    {
      comp <- prepare_comp(subset_data, needed_cols)
      model <- pgls(as.formula(formula_text), data = comp, lambda = "ML")
      term <- extract_term(model, term_name)
      data.frame(
        subset = label,
        model = model_name,
        formula = formula_text,
        term = term_name,
        n = length(model$residuals),
        birds = sum(subset_data$clade == "Aves", na.rm = TRUE),
        bats = sum(subset_data$clade == "Mammalia_Chiroptera", na.rm = TRUE),
        reptiles = sum(subset_data$clade == "Reptilia", na.rm = TRUE),
        nonflying_mammals = sum(subset_data$clade == "Mammalia_nonChiroptera", na.rm = TRUE),
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
        subset = label,
        model = model_name,
        formula = formula_text,
        term = term_name,
        n = nrow(subset_data),
        birds = sum(subset_data$clade == "Aves", na.rm = TRUE),
        bats = sum(subset_data$clade == "Mammalia_Chiroptera", na.rm = TRUE),
        reptiles = sum(subset_data$clade == "Reptilia", na.rm = TRUE),
        nonflying_mammals = sum(subset_data$clade == "Mammalia_nonChiroptera", na.rm = TRUE),
        lambda = NA,
        estimate = NA,
        se = NA,
        t = NA,
        p = NA,
        error = conditionMessage(e)
      )
    }
  )
}

subsets <- list(
  all_primary = data,
  transposon_coverage_ge_0_25 = data[data[[coverage_col]] >= 0.25, ],
  transposon_coverage_ge_0_50 = data[data[[coverage_col]] >= 0.50, ],
  transposon_coverage_ge_0_70 = data[data[[coverage_col]] >= 0.70, ],
  all_module_coverage_ge_0_50 = data[
    data$transposon_repeat_suppression_coverage >= 0.50 &
      data$chromatin_repression_heterochromatin_coverage >= 0.50 &
      data$DNA_repair_replication_stress_coverage >= 0.50 &
      data$proteostasis_autophagy_mitophagy_coverage >= 0.50 &
      data$cancer_surveillance_senescence_coverage >= 0.50 &
      data$inflammation_innate_immune_restraint_coverage >= 0.50,
  ],
  birds_only_transposon_coverage_ge_0_50 = data[data$clade == "Aves" & data[[coverage_col]] >= 0.50, ]
)

subset_rows <- list()
species_rows <- list()
for (subset_name in names(subsets)) {
  d <- subsets[[subset_name]]
  d$clade <- factor(d$clade)
  d$bird_status <- factor(d$bird_status, levels = c("non_bird", "bird"))
  species_rows[[length(species_rows) + 1]] <- data.frame(
    subset = subset_name,
    scientific_name = d$scientific_name,
    clade = d$clade,
    flight_status = d$flight_status,
    transposon_coverage = d[[coverage_col]],
    transposon_score = d[[score_col]],
    lifespan_residual = d$pgls_model_c_mass_clade_residual,
    stringsAsFactors = FALSE
  )
  subset_rows[[length(subset_rows) + 1]] <- fit_term(
    subset_name,
    d,
    "pgls_clade_residual_transposon",
    paste("pgls_model_c_mass_clade_residual ~", score_col),
    score_col,
    c("pgls_model_c_mass_clade_residual", score_col)
  )
  subset_rows[[length(subset_rows) + 1]] <- fit_term(
    subset_name,
    d,
    "mass_clade_transposon",
    paste("log10_max_lifespan_years ~ log10_body_mass_g + clade +", score_col),
    score_col,
    c("log10_max_lifespan_years", "log10_body_mass_g", "clade", score_col)
  )
}
subset_results <- do.call(rbind, subset_rows)
subset_results$p_bh_by_model <- NA
for (model_name in unique(subset_results$model)) {
  idx <- which(subset_results$model == model_name & is.finite(subset_results$p))
  subset_results$p_bh_by_model[idx] <- p.adjust(subset_results$p[idx], method = "BH")
}
species_table <- do.call(rbind, species_rows)

bird_specs <- list(
  residual_bird_interaction = list(
    formula = paste("pgls_model_c_mass_clade_residual ~", score_col, "* bird_status"),
    terms = c(score_col, "bird_statusbird", paste0(score_col, ":bird_statusbird")),
    cols = c("pgls_model_c_mass_clade_residual", score_col, "bird_status")
  ),
  residual_bird_interaction_plus_coverage = list(
    formula = paste("pgls_model_c_mass_clade_residual ~", score_col, "* bird_status +", coverage_col),
    terms = c(score_col, "bird_statusbird", coverage_col, paste0(score_col, ":bird_statusbird")),
    cols = c("pgls_model_c_mass_clade_residual", score_col, "bird_status", coverage_col)
  ),
  mass_bird_interaction = list(
    formula = paste("log10_max_lifespan_years ~ log10_body_mass_g +", score_col, "* bird_status"),
    terms = c("log10_body_mass_g", score_col, "bird_statusbird", paste0(score_col, ":bird_statusbird")),
    cols = c("log10_max_lifespan_years", "log10_body_mass_g", score_col, "bird_status")
  ),
  mass_bird_interaction_plus_coverage = list(
    formula = paste("log10_max_lifespan_years ~ log10_body_mass_g +", score_col, "* bird_status +", coverage_col),
    terms = c("log10_body_mass_g", score_col, "bird_statusbird", coverage_col, paste0(score_col, ":bird_statusbird")),
    cols = c("log10_max_lifespan_years", "log10_body_mass_g", score_col, "bird_status", coverage_col)
  )
)

bird_rows <- list()
for (model_name in names(bird_specs)) {
  spec <- bird_specs[[model_name]]
  for (term_name in spec$terms) {
    bird_rows[[length(bird_rows) + 1]] <- fit_term(
      "all_primary",
      data,
      model_name,
      spec$formula,
      term_name,
      spec$cols
    )
  }
}
bird_results <- do.call(rbind, bird_rows)
bird_results$p_bh_by_model <- NA
for (model_name in unique(bird_results$model)) {
  idx <- which(bird_results$model == model_name & is.finite(bird_results$p))
  bird_results$p_bh_by_model[idx] <- p.adjust(bird_results$p[idx], method = "BH")
}

dir.create(dirname(subset_table_path), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(bird_model_path), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(species_table_path), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(report_path), recursive = TRUE, showWarnings = FALSE)
write.table(subset_results, subset_table_path, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(bird_results, bird_model_path, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(species_table, species_table_path, sep = "\t", quote = FALSE, row.names = FALSE)

fmt <- function(x) {
  if (is.na(as.numeric(x))) return("NA")
  signif(as.numeric(x), 4)
}

subset_key <- subset_results[subset_results$model == "pgls_clade_residual_transposon", ]
bird_key <- bird_results[grepl(":bird_statusbird", bird_results$term), ]

lines <- c(
  "# Phase 3 High-Coverage and Bird-Enriched Model Report",
  "",
  "## High-Coverage Subset Models",
  ""
)
for (idx in seq_len(nrow(subset_key))) {
  row <- subset_key[idx, ]
  lines <- c(
    lines,
    paste0(
      "- ", row$subset,
      ": estimate=", fmt(row$estimate),
      ", p=", fmt(row$p),
      ", BH=", fmt(row$p_bh_by_model),
      ", n=", row$n,
      ", birds=", row$birds,
      ", bats=", row$bats,
      ", reptiles=", row$reptiles
    )
  )
}
lines <- c(lines, "", "## Bird Interaction Terms", "")
for (idx in seq_len(nrow(bird_key))) {
  row <- bird_key[idx, ]
  lines <- c(
    lines,
    paste0(
      "- ", row$model,
      ": term=", row$term,
      ", estimate=", fmt(row$estimate),
      ", p=", fmt(row$p),
      ", BH=", fmt(row$p_bh_by_model),
      ", n=", row$n
    )
  )
}
lines <- c(
  lines,
  "",
  "## Interpretation",
  "",
  "These Phase 3 starter models test whether the final transposon/repeat signal survives removal of low-coverage species and whether an explicit bird interaction is supported. High-coverage subset results should be treated as sensitivity analyses because removing low-coverage species changes clade balance."
)
writeLines(lines, report_path)
cat("Wrote ", subset_table_path, ", ", bird_model_path, " and ", report_path, "\n", sep = "")
