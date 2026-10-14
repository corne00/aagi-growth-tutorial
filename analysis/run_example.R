# A faster REAL APSIM example. Its results are not the full-data website results.
library(data.table)
for (file in c("01_data.R", "02_design.R", "03_gibbs.R", "04_posterior.R")) source(file.path("analysis", file))
raw <- read_apsim("data/example.tsv")
intervals <- observation_intervals(raw$observations)
training <- intervals$year != 2010
design <- prepare_design(raw$weather, intervals, training)
fit <- fit_gibbs(design$X[training, ], intervals$dy[training], intervals$genotype[training],
                 warmup = 600, samples = 1000, chains = 4)
diagnostics <- diagnose_fit(fit)
evaluation <- evaluate_fit(fit, design, intervals, !training)
dir.create("results", showWarnings = FALSE)
saveRDS(list(fit = fit, design = design, diagnostics = diagnostics,
             evaluation = evaluation), "results/example.rds")
print(evaluation$metrics)
print(diagnostics[, .(max_rhat = max(rhat), min_bulk_ess = min(ess_bulk), min_tail_ess = min(ess_tail))])
