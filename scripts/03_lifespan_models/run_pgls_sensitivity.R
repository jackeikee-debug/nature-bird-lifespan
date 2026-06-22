# PGLS sensitivity analyses for the OpenTree-induced 233-tip subset.

project_lib <- file.path(getwd(), "env", "R", "library")
if (dir.exists(project_lib)) {
  .libPaths(c(project_lib, .libPaths()))
}

library(ape)
library(caper)

tree_path <- "data/processed/phylogeny_inputs/opentree_induced_subtree.tre"
traits_path <- "data/processed/pgls_trait_table.tsv"

tree_full <- read.tree(tree_path)
traits_full <- read.delim(traits_path, stringsAsFactors = FALSE)

if (!"log10_body_mass_g" %in% names(traits_full)) {
  traits_full$log10_body_mass_g <- log10(as.numeric(traits_full$body_mass_g))
}
if (!"log10_max_lifespan_years" %in% names(traits_full)) {
  traits_full$log10_max_lifespan_years <- log10(as.numeric(traits_full$max_lifespan_years))
}
traits_full$body_mass_g <- as.numeric(traits_full$body_mass_g)
traits_full$max_lifespan_years <- as.numeric(traits_full$max_lifespan_years)

make_subset <- function(name, data) {
  if (name == "all_233") {
    return(data)
  }
  if (name == "no_human") {
    return(subset(data, scientific_name != "Homo sapiens"))
  }
  if (name == "no_human_no_huge") {
    return(subset(data, scientific_name != "Homo sapiens" & body_mass_g <= 1000000))
  }
  if (name == "small_body_le_1000g") {
    return(subset(data, body_mass_g <= 1000))
  }
  if (name == "wild_only") {
    return(subset(data, specimen_origin == "wild"))
  }
  if (name == "exclude_questionable_quality") {
    return(subset(data, data_quality != "questionable"))
  }
  stop(paste("Unknown subset", name))
}

prepare_comp <- function(data) {
  data <- data[, c(
    "opentree_tip_label",
    "scientific_name",
    "log10_body_mass_g",
    "log10_max_lifespan_years",
    "flight_status",
    "clade"
  )]
  data <- data[complete.cases(data), ]
  rownames(data) <- data$opentree_tip_label
  shared <- intersect(tree_full$tip.label, data$opentree_tip_label)
  pruned <- drop.tip(tree_full, setdiff(tree_full$tip.label, shared))
  data <- data[shared, ]
  pruned <- compute.brlen(pruned, method = "Grafen")
  data$flight_status <- factor(data$flight_status)
  data$clade <- factor(data$clade)
  comparative.data(
    phy = pruned,
    data = data,
    names.col = opentree_tip_label,
    vcv = TRUE,
    warn.dropped = FALSE
  )
}

extract_term <- function(model, term) {
  coef_table <- summary(model)$coefficients
  if (!term %in% rownames(coef_table)) {
    return(c(estimate = NA, p = NA))
  }
  c(estimate = coef_table[term, "Estimate"], p = coef_table[term, "Pr(>|t|)"])
}

fit_one <- function(comp, formula_text) {
  model <- pgls(as.formula(formula_text), data = comp, lambda = "ML")
  list(
    model = model,
    n = length(model$residuals),
    lambda = as.numeric(model$param["lambda"]),
    logLik = as.numeric(logLik(model)),
    AIC = AIC(model)
  )
}

subset_names <- c(
  "all_233",
  "no_human",
  "no_human_no_huge",
  "small_body_le_1000g",
  "wild_only",
  "exclude_questionable_quality"
)

rows <- list()

for (subset_name in subset_names) {
  subset_data <- make_subset(subset_name, traits_full)
  comp <- prepare_comp(subset_data)
  n <- nrow(comp$data)
  flight_levels <- length(levels(comp$data$flight_status))
  clade_levels <- length(levels(comp$data$clade))

  model_specs <- list(
    model_a_mass = "log10_max_lifespan_years ~ log10_body_mass_g"
  )
  if (flight_levels > 1) {
    model_specs$model_b_mass_flight <- "log10_max_lifespan_years ~ log10_body_mass_g + flight_status"
  }
  if (clade_levels > 1) {
    model_specs$model_c_mass_clade <- "log10_max_lifespan_years ~ log10_body_mass_g + clade"
  }

  for (model_name in names(model_specs)) {
    result <- tryCatch(
      fit_one(comp, model_specs[[model_name]]),
      error = function(e) list(error = conditionMessage(e))
    )
    if (!is.null(result$error)) {
      rows[[length(rows) + 1]] <- data.frame(
        subset = subset_name,
        model = model_name,
        n = n,
        lambda = NA,
        logLik = NA,
        AIC = NA,
        powered_flight_estimate = NA,
        powered_flight_p = NA,
        bat_clade_estimate = NA,
        bat_clade_p = NA,
        reptile_clade_estimate = NA,
        reptile_clade_p = NA,
        error = result$error
      )
      next
    }

    powered <- extract_term(result$model, "flight_statuspowered_flight")
    bat <- extract_term(result$model, "cladeMammalia_Chiroptera")
    reptile <- extract_term(result$model, "cladeReptilia")
    rows[[length(rows) + 1]] <- data.frame(
      subset = subset_name,
      model = model_name,
      n = result$n,
      lambda = result$lambda,
      logLik = result$logLik,
      AIC = result$AIC,
      powered_flight_estimate = powered["estimate"],
      powered_flight_p = powered["p"],
      bat_clade_estimate = bat["estimate"],
      bat_clade_p = bat["p"],
      reptile_clade_estimate = reptile["estimate"],
      reptile_clade_p = reptile["p"],
      error = ""
    )
  }
}

results <- do.call(rbind, rows)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/reports", recursive = TRUE, showWarnings = FALSE)
write.table(
  results,
  file = "results/tables/pgls_sensitivity_models.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

summary_lines <- c(
  "# PGLS Sensitivity Report",
  "",
  "First-pass PGLS sensitivity analyses using OpenTree synthetic topology and Grafen branch lengths.",
  "",
  "## Subsets",
  "- all_233",
  "- no_human",
  "- no_human_no_huge: excludes human and species > 1,000,000 g",
  "- small_body_le_1000g",
  "- wild_only",
  "- exclude_questionable_quality",
  "",
  "## Signal Columns",
  "- `powered_flight_estimate`: coefficient from mass + flight model",
  "- `bat_clade_estimate`: coefficient from mass + clade model",
  "- `reptile_clade_estimate`: coefficient from mass + clade model",
  "",
  "## Caveat",
  "These are feasibility-stage PGLS checks. They use Grafen fallback branch lengths and should be replaced by dated trees before final inference."
)
writeLines(summary_lines, "results/reports/pgls_sensitivity_report.md")

cat("Wrote results/tables/pgls_sensitivity_models.tsv and results/reports/pgls_sensitivity_report.md\n")

