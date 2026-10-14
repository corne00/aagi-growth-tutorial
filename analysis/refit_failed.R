# Extend chains whose rank-normalised R-hat / ESS miss the stated thresholds.
# Continuation starts at saved states with an independent documented RNG seed.
library(data.table)
for (file in c("01_data.R", "02_design.R", "03_gibbs.R", "04_posterior.R")) source(file.path("analysis", file))
prepared <- readRDS("results/prepared.rds"); intervals <- prepared$intervals
for (file in list.files("results", pattern = "^(shared|hierarchical|svd)_(full|[0-9]).*rds$", full.names = TRUE)) {
  result <- readRDS(file)
  for (round in 1:4) {
    d <- result$diagnostics
    if (max(d$rhat) < 1.01 && min(d$ess_bulk) >= 400 && min(d$ess_tail) >= 400) break
    message("Extending ", basename(file), " (max Rhat ", round(max(d$rhat), 4), ")")
    training <- !intervals$year %in% result$holdout_years
    old <- result$fit; new_seed <- old$seed + 100000L + old$samples
    extra <- fit_gibbs(result$design$X[training, ], intervals$dy[training], intervals$genotype[training],
      mode = old$mode, warmup = 0, samples = 600, chains = length(old$chains), seed = new_seed, resume_fit = old)
    combine_array <- function(a, b) {
      dimensions <- dim(a); na <- dimensions[1]; nb <- dim(b)[1]
      if (length(a) == 0) return(array(numeric(), c(na + nb, dimensions[-1])))
      array(rbind(matrix(a, nrow = na), matrix(b, nrow = nb)), c(na + nb, dimensions[-1]))
    }
    for (chain in seq_along(old$chains)) {
      for (name in c("beta", "lambdaG", "lambda0")) old$chains[[chain]][[name]] <-
        rbind(old$chains[[chain]][[name]], extra$chains[[chain]][[name]])
      for (name in c("effective", "site")) old$chains[[chain]][[name]] <-
        combine_array(old$chains[[chain]][[name]], extra$chains[[chain]][[name]])
      old$chains[[chain]]$sigma2 <- c(old$chains[[chain]]$sigma2, extra$chains[[chain]]$sigma2)
    }
    old$samples <- old$samples + extra$samples
    old$refinements <- c(old$refinements, list(list(seed = new_seed, additional_samples = 600)))
    result$fit <- old; result$diagnostics <- diagnose_fit(old)
    if (length(result$holdout_years)) result$evaluation <- evaluate_fit(old, result$design, intervals, !training)
    saveRDS(result, file, compress = FALSE)
  }
  if (max(result$diagnostics$rhat) >= 1.01 || min(result$diagnostics$ess_bulk) < 400 || min(result$diagnostics$ess_tail) < 400)
    stop("Sampling diagnostic threshold still missed: ", file)
}
