# DateLife-tree phylogenetic signal and clade-specific sensitivity models.

project_lib <- file.path(getwd(), "env", "R_library")
if (dir.exists(project_lib)) {
  .libPaths(c(project_lib, .libPaths()))
}

library(ape)
library(caper)
library(phytools)

args <- commandArgs(trailingOnly = TRUE)
data_path <- ifelse(length(args) >= 1, args[[1]], "data/processed/maintenance_lifespan_phase3_gff_cds_plus_uniprot_sequence_rescued.tsv")
tree_path <- ifelse(length(args) >= 2, args[[2]], "data/processed/phylogeny_inputs/opentree_datelife_calibrated_primary68.tre")
signal_output <- ifelse(length(args) >= 3, args[[3]], "results/tables/communications_biology_phylogenetic_signal.tsv")
model_output <- ifelse(length(args) >= 4, args[[4]], "results/tables/communications_biology_clade_specific_datelife_pgls.tsv")
report_output <- ifelse(length(args) >= 5, args[[5]], "results/reports/communications_biology_phylogenetic_signal_and_clade_models.md")

variant <- "phase2_W3_full_background_sensitivity"
score_col <- "transposon_repeat_suppression_score"
coverage_col <- "transposon_repeat_suppression_coverage"

tree_full <- read.tree(tree_path)
data_all <- read.delim(data_path, stringsAsFactors = FALSE, check.names = FALSE)
data <- data_all[data_all$score_variant == variant, ]
data <- data[!duplicated(data$scientific_name), ]
data$log10_max_lifespan_years <- log10(as.numeric(data$max_lifespan_years))
data$log10_body_mass_g <- as.numeric(data$log10_body_mass_g)
data$pgls_model_c_mass_clade_residual <- as.numeric(data$pgls_model_c_mass_clade_residual)
data[[score_col]] <- as.numeric(data[[score_col]])
data[[coverage_col]] <- as.numeric(data[[coverage_col]])
data$clade <- factor(data$clade)

subset_definitions <- list(
  all_primary = function(x) rep(TRUE, nrow(x)),
  birds_only = function(x) x$clade == "Aves",
  bats_only = function(x) x$clade == "Mammalia_Chiroptera",
  nonflying_mammals_only = function(x) x$clade == "Mammalia_nonChiroptera",
  mammals_all = function(x) x$clade %in% c("Mammalia_Chiroptera", "Mammalia_nonChiroptera"),
  reptiles_only = function(x) x$clade == "Reptilia",
  no_birds = function(x) x$clade != "Aves",
  no_bats = function(x) x$clade != "Mammalia_Chiroptera",
  no_nonflying_mammals = function(x) x$clade != "Mammalia_nonChiroptera",
  no_reptiles = function(x) x$clade != "Reptilia"
)

prepare_tree_data <- function(input, columns) {
  keep <- unique(c("opentree_tip_label", columns))
  x <- input[, keep]
  x <- x[complete.cases(x), ]
  x <- x[!duplicated(x$opentree_tip_label), ]
  shared <- intersect(tree_full$tip.label, x$opentree_tip_label)
  if (length(shared) < 4) stop("too few tree-matched observations")
  phy <- drop.tip(tree_full, setdiff(tree_full$tip.label, shared))
  rownames(x) <- x$opentree_tip_label
  x <- x[phy$tip.label, ]
  list(phy = phy, data = x)
}

signal_traits <- c(
  log10_max_lifespan_years = "log10 maximum lifespan",
  pgls_model_c_mass_clade_residual = "mass-and-clade-adjusted lifespan residual",
  transposon_repeat_suppression_score = "transposon/repeat score",
  transposon_repeat_suppression_coverage = "transposon/repeat coverage"
)

set.seed(20260621)
signal_rows <- list()
for (trait_col in names(signal_traits)) {
  prepared <- prepare_tree_data(data, trait_col)
  values <- prepared$data[[trait_col]]
  names(values) <- prepared$data$opentree_tip_label
  k_result <- tryCatch(
    phytools::phylosig(prepared$phy, values, method = "K", test = TRUE, nsim = 999),
    error = function(e) e
  )
  lambda_result <- tryCatch(
    phytools::phylosig(prepared$phy, values, method = "lambda", test = TRUE),
    error = function(e) e
  )
  signal_rows[[length(signal_rows) + 1]] <- data.frame(
    trait = trait_col,
    trait_label = signal_traits[[trait_col]],
    n = length(values),
    blomberg_k = if (inherits(k_result, "error")) NA else as.numeric(k_result$K),
    blomberg_k_p = if (inherits(k_result, "error")) NA else as.numeric(k_result$P),
    blomberg_k_permutations = 999,
    pagel_lambda = if (inherits(lambda_result, "error")) NA else as.numeric(lambda_result$lambda),
    pagel_lambda_p = if (inherits(lambda_result, "error")) NA else as.numeric(lambda_result$P),
    k_error = if (inherits(k_result, "error")) conditionMessage(k_result) else "",
    lambda_error = if (inherits(lambda_result, "error")) conditionMessage(lambda_result) else "",
    stringsAsFactors = FALSE
  )
}
signal_table <- do.call(rbind, signal_rows)

fit_pgls <- function(input, subset_name, model_name, formula_text, needed_cols) {
  result <- tryCatch({
    prepared <- prepare_tree_data(input, needed_cols)
    x <- prepared$data
    x$score_z <- as.numeric(scale(x[[score_col]]))
    x$clade <- factor(x$clade)
    comp <- comparative.data(
      phy = prepared$phy,
      data = x,
      names.col = opentree_tip_label,
      vcv = TRUE,
      warn.dropped = FALSE
    )
    model <- pgls(as.formula(formula_text), data = comp, lambda = "ML")
    term <- summary(model)$coefficients["score_z", ]
    data.frame(
      score_variant = variant,
      subset = subset_name,
      model = model_name,
      formula = formula_text,
      n = nrow(x),
      clade_count = nlevels(factor(x$clade)),
      score_mean = mean(x[[score_col]]),
      score_sd = sd(x[[score_col]]),
      lambda = as.numeric(model$param["lambda"]),
      estimate_per_score_sd = unname(term["Estimate"]),
      se = unname(term["Std. Error"]),
      conf_low = unname(term["Estimate"] - 1.96 * term["Std. Error"]),
      conf_high = unname(term["Estimate"] + 1.96 * term["Std. Error"]),
      t = unname(term["t value"]),
      p = unname(term["Pr(>|t|)"]),
      error = "",
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    data.frame(
      score_variant = variant, subset = subset_name, model = model_name,
      formula = formula_text, n = nrow(input), clade_count = length(unique(input$clade)),
      score_mean = NA, score_sd = NA, lambda = NA, estimate_per_score_sd = NA,
      se = NA, conf_low = NA, conf_high = NA, t = NA, p = NA,
      error = conditionMessage(e), stringsAsFactors = FALSE
    )
  })
  result
}

model_rows <- list()
for (subset_name in names(subset_definitions)) {
  subset_data <- data[subset_definitions[[subset_name]](data), ]
  subset_data$clade <- droplevels(factor(subset_data$clade))
  model_rows[[length(model_rows) + 1]] <- fit_pgls(
    subset_data,
    subset_name,
    "clade_adjusted_residual",
    "pgls_model_c_mass_clade_residual ~ score_z",
    c("pgls_model_c_mass_clade_residual", score_col, "clade")
  )
  model_rows[[length(model_rows) + 1]] <- fit_pgls(
    subset_data,
    subset_name,
    "mass_adjusted_lifespan",
    "log10_max_lifespan_years ~ log10_body_mass_g + score_z",
    c("log10_max_lifespan_years", "log10_body_mass_g", score_col, "clade")
  )
  if (nlevels(subset_data$clade) > 1) {
    model_rows[[length(model_rows) + 1]] <- fit_pgls(
      subset_data,
      subset_name,
      "mass_and_subset_clade_adjusted_lifespan",
      "log10_max_lifespan_years ~ log10_body_mass_g + clade + score_z",
      c("log10_max_lifespan_years", "log10_body_mass_g", score_col, "clade")
    )
  }
}
models <- do.call(rbind, model_rows)
models$q_by_model <- NA_real_
for (model_name in unique(models$model)) {
  index <- models$model == model_name & is.finite(models$p)
  models$q_by_model[index] <- p.adjust(models$p[index], method = "BH")
}
models <- models[order(models$model, models$p), ]

dir.create(dirname(signal_output), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(model_output), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(report_output), recursive = TRUE, showWarnings = FALSE)
write.table(signal_table, signal_output, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
write.table(models, model_output, sep = "\t", row.names = FALSE, quote = FALSE, na = "")

fmt <- function(x) ifelse(is.finite(x), format(signif(x, 4), scientific = FALSE), "NA")
report <- c(
  "# Phylogenetic Signal and Clade-Specific DateLife PGLS",
  "",
  "Blomberg K used 999 tip-label permutations. Pagel lambda and all PGLS lambda values were estimated on the DateLife-calibrated OpenTree chronogram. Clade-only models are sensitivity analyses and are interpreted using effect sizes and confidence intervals because the bat and reptile samples are small.",
  "",
  "## Phylogenetic signal",
  ""
)
for (i in seq_len(nrow(signal_table))) {
  row <- signal_table[i, ]
  report <- c(
    report,
    paste0(
      "- ", row$trait_label, ": K = ", fmt(row$blomberg_k),
      " (permutation P = ", fmt(row$blomberg_k_p), "); lambda = ",
      fmt(row$pagel_lambda), " (P = ", fmt(row$pagel_lambda_p), "); n = ", row$n, "."
    )
  )
}
report <- c(report, "", "## Clade-adjusted residual models", "")
residual_models <- models[models$model == "clade_adjusted_residual", ]
residual_models <- residual_models[match(names(subset_definitions), residual_models$subset), ]
for (i in seq_len(nrow(residual_models))) {
  row <- residual_models[i, ]
  report <- c(
    report,
    paste0(
      "- ", row$subset, ": beta per score SD = ", fmt(row$estimate_per_score_sd),
      " (95% CI ", fmt(row$conf_low), " to ", fmt(row$conf_high),
      "), P = ", fmt(row$p), ", q = ", fmt(row$q_by_model),
      ", lambda = ", fmt(row$lambda), ", n = ", row$n,
      ifelse(row$error == "", ".", paste0(", error: ", row$error))
    )
  )
}
report <- c(
  report,
  "",
  "Clade-specific non-significance is not treated as evidence of no association when confidence intervals are wide. The primary purpose is to show whether direction and plausible effect magnitude are consistent across taxonomic partitions."
)
writeLines(report, report_output)
cat("Wrote ", signal_output, ", ", model_output, ", and ", report_output, "\n", sep = "")
