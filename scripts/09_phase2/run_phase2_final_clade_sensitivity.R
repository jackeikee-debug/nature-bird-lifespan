# Final Phase 2 clade and outlier sensitivity for W3 full-background scores.

project_lib <- file.path(getwd(), "env", "R", "library")
if (dir.exists(project_lib)) {
  .libPaths(c(project_lib, .libPaths()))
}

library(ape)
library(caper)

args <- commandArgs(trailingOnly = TRUE)
data_path <- ifelse(length(args) >= 1, args[[1]], "data/processed/maintenance_lifespan_phase2_W3_full_background_expanded.tsv")
table_path <- ifelse(length(args) >= 2, args[[2]], "results/tables/phase2_final_clade_sensitivity_pgls.tsv")
drop_path <- ifelse(length(args) >= 3, args[[3]], "results/tables/phase2_final_clade_sensitivity_dropped_species.tsv")
report_path <- ifelse(length(args) >= 4, args[[4]], "results/reports/phase2_final_clade_sensitivity_report.md")
tree_path <- "data/processed/phylogeny_inputs/opentree_induced_subtree.tre"
variant_name <- "phase2_W3_full_background_sensitivity"

tree_full <- read.tree(tree_path)
data_all <- read.delim(data_path, stringsAsFactors = FALSE)
data_full <- data_all[data_all$score_variant == variant_name, ]

score_cols <- grep("_score$", names(data_full), value = TRUE)
module_names <- sub("_score$", "", score_cols)

for (col in c(
  "max_lifespan_years", "log10_body_mass_g", "lifespan_residual_log10",
  "pgls_model_c_mass_clade_residual", score_cols
)) {
  data_full[[col]] <- as.numeric(data_full[[col]])
}
data_full$log10_max_lifespan_years <- log10(data_full$max_lifespan_years)
data_full$clade <- factor(data_full$clade)
data_full$flight_status <- factor(data_full$flight_status)

top_abs_names <- function(data, n) {
  ordered <- data[order(abs(data$pgls_model_c_mass_clade_residual), decreasing = TRUE), ]
  head(ordered$scientific_name, n)
}

dropped_records <- list()
subset_data <- function(subset_name, data) {
  out <- data
  dropped <- character(0)
  if (subset_name == "all_primary") {
    out <- data
  } else if (subset_name == "birds_only") {
    out <- subset(data, clade == "Aves")
    dropped <- setdiff(data$scientific_name, out$scientific_name)
  } else if (subset_name == "no_birds") {
    out <- subset(data, clade != "Aves")
    dropped <- setdiff(data$scientific_name, out$scientific_name)
  } else if (subset_name == "no_bats") {
    out <- subset(data, clade != "Mammalia_Chiroptera")
    dropped <- setdiff(data$scientific_name, out$scientific_name)
  } else if (subset_name == "no_reptiles") {
    out <- subset(data, clade != "Reptilia")
    dropped <- setdiff(data$scientific_name, out$scientific_name)
  } else if (subset_name == "no_nonflying_mammals") {
    out <- subset(data, clade != "Mammalia_nonChiroptera")
    dropped <- setdiff(data$scientific_name, out$scientific_name)
  } else if (subset_name == "exclude_top_abs_residual_5") {
    dropped <- top_abs_names(data, 5)
    out <- data[!(data$scientific_name %in% dropped), ]
  } else if (subset_name == "exclude_top_abs_residual_10") {
    dropped <- top_abs_names(data, 10)
    out <- data[!(data$scientific_name %in% dropped), ]
  } else {
    stop(paste("Unknown subset", subset_name))
  }
  dropped_records[[subset_name]] <<- dropped
  out$clade <- factor(out$clade)
  out$flight_status <- factor(out$flight_status)
  out
}

prepare_comp <- function(data, needed_cols) {
  cols <- unique(c("opentree_tip_label", needed_cols))
  model_data <- data[, cols]
  model_data <- model_data[complete.cases(model_data), ]
  if (nrow(model_data) < 12) {
    stop("too_few_complete_rows")
  }
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

fit_one <- function(subset_name, data, module, score_col, model_name, formula_text, needed_cols) {
  result <- tryCatch(
    {
      comp <- prepare_comp(data, needed_cols)
      model <- pgls(as.formula(formula_text), data = comp, lambda = "ML")
      term <- extract_predictor(model, score_col)
      data.frame(
        score_variant = variant_name,
        subset = subset_name,
        maintenance_module = module,
        score_col = score_col,
        model = model_name,
        formula = formula_text,
        n = length(model$residuals),
        clade_count = length(unique(data$clade)),
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
        score_variant = variant_name,
        subset = subset_name,
        maintenance_module = module,
        score_col = score_col,
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
  "birds_only",
  "no_birds",
  "no_bats",
  "no_reptiles",
  "no_nonflying_mammals",
  "exclude_top_abs_residual_5",
  "exclude_top_abs_residual_10"
)

rows <- list()
for (subset_name in subset_names) {
  d <- subset_data(subset_name, data_full)
  for (idx in seq_along(score_cols)) {
    score_col <- score_cols[[idx]]
    module <- module_names[[idx]]
    specs <- list(
      pgls_clade_residual_module = list(
        formula = paste("pgls_model_c_mass_clade_residual ~", score_col),
        cols = c("pgls_model_c_mass_clade_residual", score_col)
      ),
      mass_module = list(
        formula = paste("log10_max_lifespan_years ~ log10_body_mass_g +", score_col),
        cols = c("log10_max_lifespan_years", "log10_body_mass_g", score_col)
      )
    )
    if (length(unique(d$clade)) > 1) {
      specs$mass_clade_module <- list(
        formula = paste("log10_max_lifespan_years ~ log10_body_mass_g + clade +", score_col),
        cols = c("log10_max_lifespan_years", "log10_body_mass_g", "clade", score_col)
      )
    }
    for (model_name in names(specs)) {
      rows[[length(rows) + 1]] <- fit_one(
        subset_name,
        d,
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
results$p_bh_by_subset_model <- NA
for (subset_name in unique(results$subset)) {
  for (model_name in unique(results$model)) {
    idx <- which(results$subset == subset_name & results$model == model_name & is.finite(results$p))
    results$p_bh_by_subset_model[idx] <- p.adjust(results$p[idx], method = "BH")
  }
}
results$rank_by_p <- NA
for (subset_name in unique(results$subset)) {
  for (model_name in unique(results$model)) {
    idx <- which(results$subset == subset_name & results$model == model_name & is.finite(results$p))
    ordered <- idx[order(results$p[idx])]
    results$rank_by_p[ordered] <- seq_along(ordered)
  }
}
results <- results[order(results$subset, results$model, results$p), ]

drop_rows <- list()
for (subset_name in names(dropped_records)) {
  dropped <- dropped_records[[subset_name]]
  if (length(dropped) == 0) {
    drop_rows[[length(drop_rows) + 1]] <- data.frame(subset = subset_name, scientific_name = "", reason = "none")
  } else {
    for (name in dropped) {
      reason <- ifelse(startsWith(subset_name, "exclude_top_abs"), "top_abs_pgls_residual", "clade_exclusion")
      drop_rows[[length(drop_rows) + 1]] <- data.frame(subset = subset_name, scientific_name = name, reason = reason)
    }
  }
}
dropped_df <- do.call(rbind, drop_rows)

dir.create(dirname(table_path), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(drop_path), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(report_path), recursive = TRUE, showWarnings = FALSE)
write.table(results, table_path, sep = "\t", row.names = FALSE, quote = FALSE)
write.table(dropped_df, drop_path, sep = "\t", row.names = FALSE, quote = FALSE)

trans <- subset(results, maintenance_module == "transposon_repeat_suppression" & model == "pgls_clade_residual_module")
trans <- trans[match(subset_names, trans$subset), ]
key <- subset(results, model == "pgls_clade_residual_module" & rank_by_p <= 3)
key <- key[order(key$subset, key$rank_by_p), ]

fmt <- function(x) {
  if (is.na(as.numeric(x))) return("NA")
  signif(as.numeric(x), 4)
}

lines <- c(
  "# Phase 2 Final Clade and Outlier Sensitivity",
  "",
  "## Transposon/Repeat PGLS Clade-Residual Sensitivity",
  ""
)
for (idx in seq_len(nrow(trans))) {
  row <- trans[idx, ]
  lines <- c(
    lines,
    paste0(
      "- ", row$subset,
      ": rank=", row$rank_by_p,
      ", estimate=", fmt(row$estimate),
      ", p=", fmt(row$p),
      ", BH=", fmt(row$p_bh_by_subset_model),
      ", n=", row$n,
      ", clades=", row$clade_count,
      ", lambda=", fmt(row$lambda),
      ifelse(row$error == "", "", paste0(", error=", row$error))
    )
  )
}
lines <- c(lines, "", "## Top Three Modules by Subset", "")
for (subset_name in subset_names) {
  lines <- c(lines, paste0("### ", subset_name))
  sub <- subset(key, subset == subset_name)
  for (idx in seq_len(nrow(sub))) {
    row <- sub[idx, ]
    lines <- c(
      lines,
      paste0(
        "- rank ", row$rank_by_p, ": ", row$maintenance_module,
        ", estimate=", fmt(row$estimate),
        ", p=", fmt(row$p),
        ", BH=", fmt(row$p_bh_by_subset_model),
        ", n=", row$n
      )
    )
  }
  lines <- c(lines, "")
}
lines <- c(
  lines,
  "## Dropped Outlier Species",
  "",
  paste0("- top 5 absolute residual: ", paste(dropped_records[["exclude_top_abs_residual_5"]], collapse = "; ")),
  paste0("- top 10 absolute residual: ", paste(dropped_records[["exclude_top_abs_residual_10"]], collapse = "; ")),
  "",
  "## Interpretation",
  "",
  "This sensitivity layer tests whether final W3 module ranking is clade- or outlier-dependent. Birds-only and no-birds subsets are the most decisive for whether the result is avian-dependent rather than broadly vertebrate-wide."
)
writeLines(lines, report_path)

cat("Wrote ", table_path, ", ", drop_path, " and ", report_path, "\n", sep = "")
