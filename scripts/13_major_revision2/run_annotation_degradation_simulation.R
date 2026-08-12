# Interventional masking of observed gene rows in high quality genomes.

project_lib <- file.path(getwd(), "env", "R_library")
if (dir.exists(project_lib)) .libPaths(c(project_lib, .libPaths()))

suppressPackageStartupMessages({
  library(ape)
  library(caper)
})

args <- commandArgs(trailingOnly = TRUE)
matrix_path <- ifelse(length(args) >= 1, args[[1]], "data/processed/jme_revision2_degradation_matrix.tsv")
analysis_path <- ifelse(length(args) >= 2, args[[2]], "data/processed/jme_revision2_analysis_table.tsv")
tree_path <- ifelse(length(args) >= 3, args[[3]], "data/processed/phylogeny_inputs/opentree_datelife_calibrated_primary68.tre")
rep_output <- ifelse(length(args) >= 4, args[[4]], "results/tables/jme_revision2_degradation_replicates.tsv")
summary_output <- ifelse(length(args) >= 5, args[[5]], "results/tables/jme_revision2_degradation_summary.tsv")
report_output <- ifelse(length(args) >= 6, args[[6]], "results/reports/jme_revision2_degradation_simulation.md")
n_replicates <- ifelse(length(args) >= 7, as.integer(args[[7]]), 500L)

matrix <- read.delim(matrix_path, stringsAsFactors = FALSE, check.names = FALSE)
traits <- read.delim(analysis_path, stringsAsFactors = FALSE, check.names = FALSE)
tree_full <- read.tree(tree_path)
for (column in c("max_lifespan_years", "log10_body_mass_g", "busco_complete", "log1p_pubmed_exact_binomial_count")) {
  traits[[column]] <- as.numeric(traits[[column]])
}
traits$log10_max_lifespan_years <- log10(traits$max_lifespan_years)
traits$clade <- factor(traits$clade)

high_quality <- traits$genome_analysis_tier == "tier1_refseq_annotated_chromosome" &
  is.finite(traits$busco_complete) & traits$busco_complete >= 0.97
hq_traits <- droplevels(traits[high_quality, ])
hq_species <- hq_traits$scientific_name
hq_matrix <- matrix[matrix$scientific_name %in% hq_species, ]
focal_module <- "transposon_repeat_suppression"
focal <- hq_matrix[hq_matrix$maintenance_module_v2 == focal_module, ]

score_species <- function(weights) {
  aggregate(weights, by = list(scientific_name = focal$scientific_name), FUN = mean)
}

baseline_scores <- score_species(focal$observed_weight)
names(baseline_scores)[2] <- "score"
hq_traits <- merge(hq_traits, baseline_scores, by = "scientific_name", all.x = TRUE)

prepare_comp <- function(input, response, predictor = "score") {
  columns <- unique(c("opentree_tip_label", response, "log10_body_mass_g", "clade", predictor))
  x <- input[, columns, drop = FALSE]
  x <- x[complete.cases(x), , drop = FALSE]
  shared <- intersect(tree_full$tip.label, x$opentree_tip_label)
  phy <- drop.tip(tree_full, setdiff(tree_full$tip.label, shared))
  rownames(x) <- x$opentree_tip_label
  x <- x[phy$tip.label, , drop = FALSE]
  x$clade <- droplevels(factor(x$clade))
  x$score_z <- as.numeric(scale(x[[predictor]]))
  comparative.data(phy = phy, data = x, names.col = opentree_tip_label, vcv = TRUE, warn.dropped = FALSE)
}

fit_pgls <- function(input, response) {
  tryCatch({
    comp <- prepare_comp(input, response)
    model <- pgls(as.formula(paste(response, "~ log10_body_mass_g + clade + score_z")), data = comp, lambda = "ML")
    term <- summary(model)$coefficients["score_z", ]
    c(beta = unname(term["Estimate"]), se = unname(term["Std. Error"]), p = unname(term["Pr(>|t|)"]))
  }, error = function(e) c(beta = NA_real_, se = NA_real_, p = NA_real_))
}

baseline_fit <- fit_pgls(hq_traits, "log10_max_lifespan_years")
loss_rates <- c(0.05, 0.10, 0.20, 0.30)
mechanisms <- c("random", "low_research_attention", "outcome_related")
set.seed(20260812)

phy_hq <- drop.tip(tree_full, setdiff(tree_full$tip.label, hq_traits$opentree_tip_label))
ordered_traits <- hq_traits[match(phy_hq$tip.label, hq_traits$opentree_tip_label), ]
vcv_hq <- vcv(phy_hq, corr = TRUE) + diag(1e-8, length(phy_hq$tip.label))
inv_vcv <- solve(vcv_hq)
L <- chol(vcv_hq)

gls_fit <- function(response, predictor, effort = NULL) {
  predictor_z <- as.numeric(scale(predictor))
  X <- cbind(1, predictor_z)
  if (!is.null(effort)) X <- cbind(X, as.numeric(scale(effort)))
  xtvix <- t(X) %*% inv_vcv %*% X
  beta <- solve(xtvix, t(X) %*% inv_vcv %*% response)
  residual <- response - X %*% beta
  sigma2 <- as.numeric(t(residual) %*% inv_vcv %*% residual) / (length(response) - ncol(X))
  covariance <- sigma2 * solve(xtvix)
  se <- sqrt(diag(covariance))
  t_value <- beta[2] / se[2]
  p_value <- 2 * pt(abs(t_value), df = length(response) - ncol(X), lower.tail = FALSE)
    c(beta = unname(beta[2]), se = unname(se[2]), p = unname(p_value))
}

degraded_scores <- function(probability, loss_rate) {
  probability[focal$observed_weight <= 0] <- 0
  probability <- probability / mean(probability[focal$observed_weight > 0]) * loss_rate
  probability <- pmin(probability, 0.95)
  mask <- runif(nrow(focal)) < probability & focal$observed_weight > 0
  weights <- focal$observed_weight
  weights[mask] <- 0
  scores <- score_species(weights)
  names(scores)[2] <- "score"
  score_by_species <- setNames(scores$score, scores$scientific_name)
  list(
    scores = score_by_species[ordered_traits$scientific_name],
    mask = mask,
    achieved_loss = sum(mask) / sum(focal$observed_weight > 0)
  )
}

rows <- list()
for (mechanism in mechanisms) {
  for (loss_rate in loss_rates) {
    for (iteration in seq_len(n_replicates)) {
      true_score <- ordered_traits$score
      true_score_z <- as.numeric(scale(true_score))
      noise <- as.numeric(t(L) %*% rnorm(nrow(vcv_hq)))
      noise <- as.numeric(scale(noise)) * 0.20
      effort_z <- as.numeric(scale(ordered_traits$log1p_pubmed_exact_binomial_count))
      responses <- list(
        null_independent = noise,
        attention_confounded = 0.10 * effort_z + noise,
        planted = 0.10 * true_score_z + noise
      )

      for (world in names(responses)) {
        response <- responses[[world]]
        true_beta <- ifelse(world == "planted", 0.10, 0)
        effort_by_species <- setNames(ordered_traits$log1p_pubmed_exact_binomial_count, ordered_traits$scientific_name)
        response_by_species <- setNames(response, ordered_traits$scientific_name)
        effort_scaled <- as.numeric(scale(effort_by_species[focal$scientific_name]))
        response_scaled <- as.numeric(scale(response_by_species[focal$scientific_name]))
        if (mechanism == "random") raw_probability <- rep(1, nrow(focal))
        if (mechanism == "low_research_attention") raw_probability <- exp(-effort_scaled)
        if (mechanism == "outcome_related") raw_probability <- exp(response_scaled)
        degraded <- degraded_scores(raw_probability, loss_rate)
        fit <- gls_fit(response, degraded$scores)
        fit_adjusted <- gls_fit(response, degraded$scores, ordered_traits$log1p_pubmed_exact_binomial_count)
        lower <- fit["beta"] - 1.96 * fit["se"]
        upper <- fit["beta"] + 1.96 * fit["se"]
        lower_adjusted <- fit_adjusted["beta"] - 1.96 * fit_adjusted["se"]
        upper_adjusted <- fit_adjusted["beta"] + 1.96 * fit_adjusted["se"]
        rows[[length(rows) + 1]] <- data.frame(
          mechanism = mechanism, world = world, true_beta = true_beta,
          loss_rate = loss_rate, replicate = iteration, n_species = length(response),
          masked_rows = sum(degraded$mask), achieved_loss_fraction = degraded$achieved_loss,
          score_bias = mean(degraded$scores - true_score),
          score_rmse = sqrt(mean((degraded$scores - true_score)^2)),
          estimated_beta = unname(fit["beta"]), estimated_se = unname(fit["se"]), p = unname(fit["p"]),
          beta_bias = unname(fit["beta"] - true_beta),
          ci_covers_truth = lower <= true_beta && upper >= true_beta,
          rejected_p_lt_0.05 = is.finite(fit["p"]) && fit["p"] < 0.05,
          adjusted_estimated_beta = unname(fit_adjusted["beta"]),
          adjusted_beta_bias = unname(fit_adjusted["beta"] - true_beta),
          adjusted_ci_covers_truth = lower_adjusted <= true_beta && upper_adjusted >= true_beta,
          adjusted_rejected_p_lt_0.05 = is.finite(fit_adjusted["p"]) && fit_adjusted["p"] < 0.05,
          stringsAsFactors = FALSE
        )
      }
    }
  }
}
replicates <- do.call(rbind, rows)
if (nrow(replicates) == 0) stop("No degradation replicate rows were generated")
if (all(!complete.cases(replicates[, c("achieved_loss_fraction", "score_bias", "score_rmse", "estimated_beta", "beta_bias")]))) {
  stop("All degradation replicate metrics are missing")
}
summary <- aggregate(
  cbind(achieved_loss_fraction, score_bias, score_rmse, estimated_beta, beta_bias,
        ci_covers_truth, rejected_p_lt_0.05, adjusted_estimated_beta, adjusted_beta_bias,
        adjusted_ci_covers_truth, adjusted_rejected_p_lt_0.05) ~ mechanism + world + true_beta + loss_rate,
  data = replicates, FUN = mean, na.rm = TRUE
)

dir.create(dirname(rep_output), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(report_output), recursive = TRUE, showWarnings = FALSE)
write.table(replicates, rep_output, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
write.table(summary, summary_output, sep = "\t", row.names = FALSE, quote = FALSE, na = "")

fmt <- function(x) format(signif(x, 4), scientific = FALSE)
report <- c(
  "# Annotation Matrix Degradation Experiment",
  "",
  paste0("High quality starting panel: ", nrow(hq_traits), " species (Tier 1 and BUSCO completeness >= 0.97)."),
  paste0("Baseline transposon/repeat PGLS effect: beta per score SD = ", fmt(baseline_fit["beta"]), ", P = ", fmt(baseline_fit["p"]), "."),
  paste0("Replicates per mechanism and loss rate: ", n_replicates, "."),
  "",
  "The experiment masks previously observed gene rows and then recalculates the module score. It does not fragment genome FASTA files or rerun structural annotation, so it is an intervention on the observed annotation matrix rather than a physical assembly degradation benchmark.",
  "",
  "## Results",
  ""
)
for (i in seq_len(nrow(summary))) {
  row <- summary[i, ]
  report <- c(report, paste0(
    "- ", row$mechanism, ", target loss ", 100 * row$loss_rate, "%: achieved loss = ",
    fmt(100 * row$achieved_loss_fraction), "%, score RMSE = ", fmt(row$score_rmse),
    ", world = ", row$world, ", beta bias = ", fmt(row$beta_bias),
    ", 95% CI coverage = ", fmt(row$ci_covers_truth),
    ifelse(row$world == "planted", ", power = ", ", false positive rate = "),
    fmt(row$rejected_p_lt_0.05),
    "; effort-adjusted beta bias = ", fmt(row$adjusted_beta_bias),
    ", effort-adjusted rejection rate = ", fmt(row$adjusted_rejected_p_lt_0.05), "."
  ))
}
writeLines(report, report_output)
cat("Wrote annotation degradation simulation outputs.\n")
