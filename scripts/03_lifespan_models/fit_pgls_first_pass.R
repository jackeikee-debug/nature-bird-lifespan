# First-pass PGLS template.
#
# Requires R packages: ape, caper, nlme.
# The OpenTree synthetic tree is topology-only and should be given branch
# lengths before PGLS. Grafen branch lengths are a pragmatic feasibility-stage
# fallback, not a final manuscript method.

project_lib <- file.path(getwd(), "env", "R", "library")
if (dir.exists(project_lib)) {
  .libPaths(c(project_lib, .libPaths()))
}

library(ape)
library(caper)

tree_path <- "data/processed/phylogeny_inputs/opentree_induced_subtree.tre"
traits_path <- "data/processed/pgls_trait_table.tsv"

tree <- read.tree(tree_path)
traits <- read.delim(traits_path, stringsAsFactors = FALSE)

rownames(traits) <- traits$opentree_tip_label
shared <- intersect(tree$tip.label, traits$opentree_tip_label)
tree <- drop.tip(tree, setdiff(tree$tip.label, shared))
traits <- traits[shared, ]

# OpenTree induced subtrees generally lack branch lengths.
tree <- compute.brlen(tree, method = "Grafen")

if (!"log10_body_mass_g" %in% names(traits)) {
  traits$log10_body_mass_g <- log10(as.numeric(traits$body_mass_g))
}
if (!"log10_max_lifespan_years" %in% names(traits)) {
  traits$log10_max_lifespan_years <- log10(as.numeric(traits$max_lifespan_years))
}
traits$log10_body_mass_g <- as.numeric(traits$log10_body_mass_g)
traits$log10_max_lifespan_years <- as.numeric(traits$log10_max_lifespan_years)
traits$flight_status <- factor(traits$flight_status)
traits$clade <- factor(traits$clade)

pgls_traits <- traits[, c(
  "opentree_tip_label",
  "log10_body_mass_g",
  "log10_max_lifespan_years",
  "flight_status",
  "clade"
)]

comp <- comparative.data(
  phy = tree,
  data = pgls_traits,
  names.col = opentree_tip_label,
  vcv = TRUE,
  warn.dropped = TRUE
)

model_a <- pgls(log10_max_lifespan_years ~ log10_body_mass_g, data = comp, lambda = "ML")
model_b <- pgls(log10_max_lifespan_years ~ log10_body_mass_g + flight_status, data = comp, lambda = "ML")
model_c <- pgls(log10_max_lifespan_years ~ log10_body_mass_g + clade, data = comp, lambda = "ML")

dir.create("results/models", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
sink("results/models/pgls_first_pass_summary.txt")
cat("Model A\n")
print(summary(model_a))
cat("\nModel B\n")
print(summary(model_b))
cat("\nModel C\n")
print(summary(model_c))
sink()

model_table <- data.frame(
  model = c("model_a_mass", "model_b_mass_flight", "model_c_mass_clade"),
  n = c(length(model_a$residuals), length(model_b$residuals), length(model_c$residuals)),
  lambda = c(model_a$param["lambda"], model_b$param["lambda"], model_c$param["lambda"]),
  logLik = c(as.numeric(logLik(model_a)), as.numeric(logLik(model_b)), as.numeric(logLik(model_c))),
  AIC = c(AIC(model_a), AIC(model_b), AIC(model_c))
)
write.table(
  model_table,
  file = "results/tables/pgls_first_pass_models.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

residual_table <- data.frame(
  opentree_tip_label = rownames(comp$data),
  pgls_model_b_mass_flight_residual = residuals(model_b),
  pgls_model_c_mass_clade_residual = residuals(model_c)
)
write.table(
  residual_table,
  file = "data/processed/pgls_first_pass_residuals.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
