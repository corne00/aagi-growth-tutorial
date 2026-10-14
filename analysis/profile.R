# Repeatable profiling on real full-data sufficient statistics.
# Rscript analysis/profile.R before|diagonal|location
library(data.table)
source("analysis/02_design.R"); source("analysis/03_gibbs.R")
prepared <- readRDS("results/prepared.rds")
train <- prepared$intervals$year > 1990
design <- prepare_design(prepared$weather, prepared$intervals, train)
args <- commandArgs(trailingOnly = TRUE); label <- if (length(args)) args[1] else "current"
# Warm up R's byte-code compiler before timed comparisons.
invisible(fit_gibbs(design$X[train, ], prepared$intervals$dy[train],
  prepared$intervals$genotype[train], warmup = 2, samples = 2, chains = 1))
Rprof(paste0("results/profile-", label, ".out"), interval = .005)
time <- system.time(fit <- fit_gibbs(design$X[train, ], prepared$intervals$dy[train],
  prepared$intervals$genotype[train],
  location = if (label == "location") prepared$intervals$location[train] else NULL,
  warmup = 100, samples = 100, chains = 1))
Rprof(NULL)
profile <- summaryRprof(paste0("results/profile-", label, ".out"))
saveRDS(list(time = time, profile = profile), paste0("results/profile-", label, ".rds"))
print(time); print(head(profile$by.total, 15))
