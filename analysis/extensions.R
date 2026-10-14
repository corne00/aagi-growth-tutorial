# Optional fitted extensions: known-location responses and SVD compression.
# Both extensions use all 199 genotypes. Evaluate on the same first year block.
library(data.table)
for (file in c("01_data.R", "02_design.R", "03_gibbs.R", "04_posterior.R")) source(file.path("analysis", file))
prepared <- readRDS("results/prepared.rds"); weather <- prepared$weather; intervals <- prepared$intervals
training <- !intervals$year %in% 1983:1990
for (extension in c("location", "svd")) {
  destination <- paste0("results/", extension, "_1983_1990.rds")
  if (file.exists(destination)) next
  design <- prepare_design(weather, intervals, training, pcr = extension == "svd")
  fit <- fit_gibbs(design$X[training, , drop = FALSE], intervals$dy[training], intervals$genotype[training],
    location = if (extension == "location") intervals$location[training] else NULL,
    warmup = as.integer(Sys.getenv("APSIM_WARMUP", "1200")),
    samples = as.integer(Sys.getenv("APSIM_SAMPLES", "1200")))
  diagnostics <- diagnose_fit(fit)
  evaluation <- evaluate_fit(fit, design, intervals, !training)
  print(evaluation$metrics)
  saveRDS(list(fit = fit, design = design, diagnostics = diagnostics,
    evaluation = evaluation, holdout_years = 1983:1990), destination, compress = FALSE)
}
