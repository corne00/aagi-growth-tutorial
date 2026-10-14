# Iterative optimisation audit against the preserved baseline implementation.
# Uses real APSIM data; records user CPU and elapsed time separately.
library(data.table)
source("analysis/02_design.R")
prepared <- readRDS("results/prepared.rds"); train <- prepared$intervals$year > 1990
design <- prepare_design(prepared$weather, prepared$intervals, train)
baseline <- new.env(parent = globalenv()); sys.source("results/gibbs-baseline.R", baseline)
optimised <- new.env(parent = globalenv()); sys.source("analysis/03_gibbs.R", optimised)
X <- design$X[train, ]; y <- prepared$intervals$dy[train]; g <- prepared$intervals$genotype[train]
measure <- function(environment, name, site = FALSE) {
  invisible(environment$fit_gibbs(X, y, g, warmup = 2, samples = 2, chains = 1))
  Rprof(paste0("results/profile-audit-", name, ".out"), interval = .01)
  time <- system.time(fit <- environment$fit_gibbs(X, y, g,
    location = if (site) prepared$intervals$location[train] else NULL,
    warmup = 100, samples = 100, chains = 1))
  Rprof(NULL)
  list(time = time, fit = fit, profile = summaryRprof(paste0("results/profile-audit-", name, ".out")))
}
results <- list()
for (site in c(FALSE, TRUE)) {
  name <- if (site) "location" else "genotype"
  a <- measure(baseline, paste0(name, "-baseline"), site)
  b <- measure(optimised, paste0(name, "-optimised"), site)
  error <- max(abs(a$fit$chains[[1]]$effective - b$fit$chains[[1]]$effective))
  if (!site) stopifnot(error < 1e-6)
  results[[name]] <- list(baseline = as.list(a$time), optimised = as.list(b$time),
    max_same_seed_coefficient_difference = if (!site) error else NULL,
    posterior_check = if (site) "Joint location kernel: dense conditional moment check in check_location.R" else "Same-seed draws agree to numerical precision",
    baseline_profile = head(a$profile$by.total, 12), optimised_profile = head(b$profile$by.total, 12))
  print(results[[name]])
}
jsonlite::write_json(results, "results/performance.json", auto_unbox = TRUE, pretty = TRUE)
message("Two-model timing and posterior-preservation audit passed.")
