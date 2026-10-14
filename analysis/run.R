# Run from the repository root: Rscript analysis/run.R full
# Dependencies: data.table, jsonlite (MPL-2.0 and MIT; everything else is bundled with R).
library(data.table)
library(jsonlite)
for (file in c("01_data.R", "02_design.R", "03_gibbs.R", "04_posterior.R")) source(file.path("analysis", file))
args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) args[1] else "example"
dir.create("results", showWarnings = FALSE)
path <- if (mode == "example") "data/example.csv" else "data/SimulatedData124Envs_Biomass.txt"
if (!file.exists(path)) stop("Data missing. Run Rscript analysis/download_data.R, or place the attached APSIM biomass file in data/.")
if (mode == "example") {
  stop("The example entry point is analysis/run_example.R; use Rscript analysis/run_example.R.")
}
raw <- read_apsim(path)
intervals <- observation_intervals(raw$observations, 14L)
saveRDS(list(weather = raw$weather, intervals = intervals), "results/prepared.rds", compress = FALSE)
years <- sort(unique(intervals$year))
folds <- split(years, ceiling(seq_along(years) / 8))
jobs <- c(list(full = integer()), setNames(folds, vapply(folds, function(v) paste0(min(v), "_", max(v)), character(1))))
warmup <- as.integer(Sys.getenv("APSIM_WARMUP", "1200")); samples <- as.integer(Sys.getenv("APSIM_SAMPLES", "1200"))
for (job in names(jobs)) {
  holdout <- jobs[[job]]; training <- !intervals$year %in% holdout
  for (model in c("shared", "hierarchical")) {
    name <- paste(model, job, sep = "_"); destination <- paste0("results/", name, ".rds")
    if (file.exists(destination)) { message("Cached: ", name); next }
    message("Fitting ", name, " on ", sum(training), " intervals")
    design <- prepare_design(raw$weather, intervals, training)
    fit <- fit_gibbs(design$X[training, , drop = FALSE], intervals$dy[training],
                     intervals$genotype[training], mode = model, warmup = warmup, samples = samples)
    diagnostics <- diagnose_fit(fit)
    print(diagnostics[, .(max_rhat = max(rhat), min_bulk_ess = min(ess_bulk), min_tail_ess = min(ess_tail))])
    evaluation <- if (length(holdout)) evaluate_fit(fit, design, intervals, !training) else NULL
    if (!is.null(evaluation)) print(evaluation$metrics)
    saveRDS(list(fit = fit, design = design, diagnostics = diagnostics, evaluation = evaluation,
                 holdout_years = holdout), destination, compress = FALSE)
  }
}
message("Main fitting complete. Run Rscript analysis/extensions.R and Rscript analysis/export.R.")
