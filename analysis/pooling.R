# Compare independent genotype fits with hierarchical pooling on real data.
# Keep the same Gaussian likelihood and beta / genotype prior specification
# for the hierarchy. Independent fits have proper N(0, I) coefficient priors.
# Thin training genotype/environment TRAJECTORIES, not synthetic observations.
library(data.table)
for (file in c("01_data.R", "02_design.R", "03_gibbs.R", "04_posterior.R")) source(file.path("analysis", file))
p <- readRDS("results/prepared.rds"); intervals <- p$intervals
selected <- sprintf("g%03d", c(8, 9, 52, 100, 105, 112, 127, 140, 185, 2, 10, 178))
intervals <- intervals[genotype %in% selected]
training_years <- 1991:2006
scarce_genotypes <- c("g008", "g052", "g185")
training <- intervals$year %in% training_years
dir.create("results", showWarnings = FALSE)
for (density in c("one", "three")) {
  scarce_environments <- if (density == "one") "Yanco_1997" else c("Yanco_1991", "Emerald_1997", "Narrabri_2003")
  keep <- training & (!intervals$genotype %in% scarce_genotypes |
                      intervals$environment %in% scarce_environments)
  design <- prepare_design(p$weather, intervals, keep)
  for (model in c("independent", "hierarchical")) {
    name <- paste0("pooling_", density, "_", model); destination <- paste0("results/", name, ".rds")
    if (file.exists(destination)) next
    fit <- fit_gibbs(design$X[keep, ], intervals$dy[keep], intervals$genotype[keep], mode = model,
                     warmup = 800, samples = 1200)
    evaluation <- evaluate_fit(fit, design, intervals,
      intervals$year %in% 2007:2013 & intervals$genotype %in% scarce_genotypes)
    diagnostics <- diagnose_fit(fit)
    saveRDS(list(fit = fit, design = design, diagnostics = diagnostics,
                 evaluation = evaluation, training_years = training_years,
                 scarce_genotypes = scarce_genotypes, scarce_environments = scarce_environments),
            destination, compress = FALSE)
    print(evaluation$metrics)
  }
}
