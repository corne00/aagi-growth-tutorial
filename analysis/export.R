# Export only summaries and selected joint posterior draws for the static site.
# Heavy posterior RDS files remain local and are reproducible from run.R.
library(data.table)
library(jsonlite)
for (file in c("01_data.R", "02_design.R", "03_gibbs.R", "04_posterior.R")) source(file.path("analysis", file))
dir.create("assets", showWarnings = FALSE)
prepared <- readRDS("results/prepared.rds"); weather <- prepared$weather; intervals <- prepared$intervals
raw <- read_apsim("data/SimulatedData124Envs_Biomass.txt")
selected_genotypes <- sprintf("g%03d", c(8, 9, 52, 100, 105, 112, 127, 140, 185))
selected_years <- c(1986L, 1994L, 2002L, 2010L)
# Fast exercise: 12 actual genotypes, 4 locations, 3 years; native daily rows.
exercise <- fread("data/SimulatedData124Envs_Biomass.txt",
  select = c("Env", "loc", "year", "geno", "das", "biomass", "MinT", "MaxT", "Radn", "Rain", "vp", "day_length"))
exercise <- exercise[geno %in% sprintf("g%03d", c(8, 9, 52, 100, 105, 112, 127, 140, 185, 2, 10, 178)) & year %in% c(1991, 1997, 2010)]
fwrite(exercise, "data/example.tsv", sep = "\t")
rm(exercise); gc()

summarise_draws <- function(values, level = .95) {
  limits <- apply(values, 1, quantile, probs = c((1 - level)/2, .5, (1 + level)/2))
  list(lower = unname(limits[1, ]), median = unname(limits[2, ]), upper = unname(limits[3, ]))
}
diagnostic_summary <- function(d) list(max_rhat = max(d$rhat), min_bulk_ess = min(d$ess_bulk),
  min_tail_ess = min(d$ess_tail), parameters = nrow(d),
  passed = max(d$rhat) < 1.01 && min(d$ess_bulk) >= 400 && min(d$ess_tail) >= 400)

data <- list(genotypes = selected_genotypes, years = selected_years, locations = sort(unique(weather$location)),
  counts = list(records = nrow(raw$observations), trajectories = uniqueN(raw$observations[, .(environment, genotype)]),
                genotypes = uniqueN(intervals$genotype), environments = uniqueN(intervals$environment), intervals = nrow(intervals)),
  source = list(doi = "10.6084/m9.figshare.24843705.v1", md5 = unname(tools::md5sum("data/SimulatedData124Envs_Biomass.txt")),
                license = "CC BY 4.0"), trajectories = list(), metrics = list(), diagnostics = list(), ale = list(), pooling = list())

# Motivation: native trajectories with similar terminal biomass, contrasting timing.
timing <- intervals[, .(endpoint = sum(dy),
  timing = sum((t0 + t1)/2 * pmax(dy, 0)) / sum(pmax(dy, 0))), by = .(environment, genotype)]
near <- timing[abs(endpoint - 8) <= .04 & is.finite(timing)]
setorder(near, timing)
choices <- near[unique(c(1, round(nrow(near)/2), nrow(near)))]
data$introduction <- lapply(seq_len(nrow(choices)), function(i) {
  row <- choices[i]; d <- raw$observations[environment == row$environment & genotype == row$genotype]
  list(id = paste(row$environment, row$genotype), time = d$time, biomass = d$biomass, endpoint = row$endpoint)
})

full <- readRDS("results/hierarchical_full.rds")
data$specification <- list(features = colnames(full$design$X), spacing_days = 14,
  integration = "daily trapezoidal quadrature", warmup = full$fit$warmup,
  samples_per_chain = full$fit$samples, chains = length(full$fit$chains), seed = full$fit$seed,
  response_unit = "t/ha", measurement_noise_added = FALSE,
  genotype_constraint = "Gaussian zero-mean prior; genotype deviations are not empirically constrained to sum to zero",
  refinements = full$fit$refinements,
  priors = full$fit$priors, R = R.version.string,
  packages = list(data.table = as.character(packageVersion("data.table")), jsonlite = as.character(packageVersion("jsonlite"))))
data$diagnostics$hierarchical_full <- diagnostic_summary(full$diagnostics)
fwrite(full$diagnostics, "results/diagnostics_hierarchical_full.csv")
data$traces <- lapply(seq_along(full$fit$chains), function(i) {
  chain <- full$fit$chains[[i]]; list(chain = i, beta = chain$beta[, 1],
     precision = chain$lambdaG[, 1], sigma = sqrt(chain$sigma2))
})

# Response decomposition evaluated along a real weather path, not invented curves.
path <- weather[environment == "Yanco_2010"]
Phi <- rate_features(path, full$design)
mean_coefficients <- posterior_mean(full$fit)
rawPhi <- sweep(evaluate_basis(path, full$design$basis), 2, full$design$scale, "/")
families <- c("intercept", "temperature", "thermal_time", "Radn", "rain14", "day_length", "thermal_radiation", "thermal_rain")
data$basis <- list(time = path$time, features = colnames(rawPhi), values = lapply(seq_len(ncol(rawPhi)), function(j) rawPhi[, j]),
  contributions = lapply(seq_len(ncol(rawPhi)), function(j) rawPhi[, j] * mean_coefficients$beta[j]),
  coefficients = mean_coefficients$beta, rate = drop(Phi %*% mean_coefficients$beta),
  environment = "Yanco_2010")
data$hierarchy <- list(time = path$time, shared = data$basis$rate,
  genotypes = setNames(lapply(selected_genotypes, function(g) {
    draw <- posterior_draws(full$fit, g, "Yanco")
    rates <- Phi %*% t(draw$theta)
    shared_draws <- do.call(rbind, lapply(full$fit$chains, function(c)
      c$beta[unique(round(seq(1, nrow(c$beta), length.out = 40))), ]))
    list(response = summarise_draws(rates), deviation = summarise_draws(rates - Phi %*% t(shared_draws)))
  }), selected_genotypes))
example_rows <- which(intervals$environment == "Yanco_2010" & intervals$genotype == "g008")
data$design <- list(time0 = intervals$t0[example_rows], time1 = intervals$t1[example_rows], dy = intervals$dy[example_rows],
  columns = colnames(full$design$X), X = lapply(example_rows, function(i) unname(full$design$X[i, ])),
  scale = full$design$scale,
  training_rms = sqrt(colMeans(full$design$X^2)))

# Full-fit ALE: environmental days have equal weight, not repeated per genotype.
ale_weather <- weather[unique(round(seq(1, .N, length.out = min(.N, 6000))))]
population_draws <- do.call(rbind, lapply(full$fit$chains, function(c) c$beta[unique(round(seq(1, nrow(c$beta), length.out = 40))), ]))
for (genotype in c("population", selected_genotypes)) {
  theta <- if (genotype == "population") population_draws else posterior_draws(full$fit, genotype, "Yanco")$theta
  for (driver in c("temperature", "thermal_time", "Radn", "rain14")) for (bins in c(6L, 10L, 15L)) {
    ale <- ale_curve(ale_weather, full$design, theta, driver, bins)
    data$ale[[paste(genotype, driver, bins, sep = "|")]] <- c(
      list(x = ale$x, counts = ale$counts, boundaries = ale$boundaries), summarise_draws(ale$draws))
  }
}
data$ale_reference <- list(days = nrow(ale_weather), rule = "6000 evenly indexed environment-days, equal weights; full fit")

oof <- list(); metric_rows <- list()
for (file in list.files("results", pattern = "^(shared|hierarchical)_[0-9].*rds$", full.names = TRUE)) {
  result <- readRDS(file); model <- result$fit$mode
  fold <- paste(range(result$holdout_years), collapse = "–"); key <- sub(".rds$", "", basename(file))
  data$diagnostics[[key]] <- diagnostic_summary(result$diagnostics)
  fwrite(result$diagnostics, paste0("results/diagnostics_", key, ".csv"))
  metric <- copy(result$evaluation$metrics); metric[, `:=`(model = model, fold = fold)]
  metric_rows[[key]] <- metric
  if (model == "hierarchical") oof[[key]] <- result$evaluation$predictions[held_out == TRUE]
  for (year in intersect(selected_years, result$holdout_years)) for (location in data$locations) for (genotype in selected_genotypes) {
    env <- paste(location, year, sep = "_")
    ii <- which(intervals$environment == env & intervals$genotype == genotype)
    selected_genotype <- genotype
    d <- raw$observations[.(env, selected_genotype)]
    w <- weather[environment == env]
    rate_weather <- w[time <= max(d$time)]
    draws <- posterior_draws(result$fit, genotype, location)
    increments <- result$design$X[ii, , drop = FALSE] %*% t(draws$theta)
    biomass <- rbind(0, apply(increments, 2, cumsum))
    set.seed(20261014 + year + match(location, data$locations) * 100 + match(genotype, selected_genotypes))
    # Independent discrepancy increments, accumulating on the SAME joint path.
    repeated <- rep(seq_len(ncol(increments)), each = 10)
    noise <- matrix(rnorm(nrow(increments) * length(repeated)), nrow(increments))
    noise <- sweep(noise, 2, sqrt(draws$sigma2[repeated]), "*")
    predictive <- rbind(0, apply(increments[, repeated, drop = FALSE] + noise, 2, cumsum))
    intervals_by_level <- lapply(c(.5, .8, .95, .99), function(level)
      list(level = level, credible = summarise_draws(biomass, level), predictive = summarise_draws(predictive, level)))
    curve <- list(genotype = genotype, location = location, year = year, fold = fold,
       model = model, time = c(0, intervals$t1[ii]), observed = c(0, cumsum(intervals$dy[ii])), mean_biomass = rowMeans(biomass),
       daily_time = d$time, daily_biomass = d$biomass,
       draws = lapply(seq_len(min(12L, ncol(biomass))), function(j) biomass[, j]),
       intervals = intervals_by_level,
       weather = list(time = w$time, temperature = w$temperature, Radn = w$Radn, Rain = w$Rain,
                      thermal_time = w$thermal_time, rain14 = w$rain14),
       rate = list(time = rate_weather$time, mean = drop(rate_features(rate_weather, result$design) %*% colMeans(draws$theta))))
    data$trajectories[[paste(model, location, year, genotype, sep = "|")]] <- curve
  }
  rm(result); gc()
}
metrics <- rbindlist(metric_rows); data$metrics <- metrics
fwrite(metrics, "results/heldout_metrics.csv")

# Coverage evaluated on all HELD-OUT terminal trajectories, not the displayed subset.
coverage <- list()
for (file in list.files("results", pattern = "^hierarchical_[0-9].*rds$", full.names = TRUE)) {
  result <- readRDS(file); keep <- intervals$year %in% result$holdout_years
  groups <- split(which(keep), paste(intervals$environment[keep], intervals$genotype[keep], sep = "|"))
  nominal <- c(.5, .8, .95, .99); covered_c <- covered_p <- rep(0L, length(nominal))
  set.seed(20261014 + min(result$holdout_years))
  for (ii in groups) {
    draw <- posterior_draws(result$fit, intervals$genotype[ii[1]], intervals$location[ii[1]])
    endpoint <- drop(colSums(result$design$X[ii, , drop = FALSE]) %*% t(draw$theta))
    # Sum of independent Gaussian increments has variance n_intervals * sigma2.
    predictive <- rep(endpoint, each = 10) + rnorm(length(endpoint) * 10) *
                  rep(sqrt(length(ii) * draw$sigma2), each = 10)
    truth <- sum(intervals$dy[ii])
    for (j in seq_along(nominal)) {
      q <- c((1 - nominal[j])/2, (1 + nominal[j])/2)
      c_limits <- quantile(endpoint, q); p_limits <- quantile(predictive, q)
      covered_c[j] <- covered_c[j] + (truth >= c_limits[1] && truth <= c_limits[2])
      covered_p[j] <- covered_p[j] + (truth >= p_limits[1] && truth <= p_limits[2])
    }
  }
  coverage[[basename(file)]] <- list(fold = paste(range(result$holdout_years), collapse = "–"),
    nominal = nominal, credible = covered_c / length(groups), predictive = covered_p / length(groups),
    trajectories = length(groups))
  rm(result); gc()
}
data$coverage <- coverage

# Balanced G/E/GxE decomposition at common DAS checkpoints. Native records only.
all_oof <- rbindlist(oof)
gxe <- list()
for (time in seq(14, 112, 14)) {
  panel <- all_oof[t1 == time]
  stopifnot(nrow(panel) == 199 * 124)
  genotypes <- sort(unique(panel$genotype)); environments <- sort(unique(panel$environment))
  matrices <- lapply(c("observed", "predicted"), function(response) {
    M <- matrix(NA_real_, length(genotypes), length(environments))
    M[cbind(match(panel$genotype, genotypes), match(panel$environment, environments))] <- panel[[response]]
    mean <- mean(M); G <- rowMeans(M) - mean; E <- colMeans(M) - mean
    interaction <- sweep(sweep(M, 1, G + mean, "-"), 2, E, "-")
    fractions <- c(G = mean(G^2), E = mean(E^2), GxE = mean(interaction^2))
    list(G = G, E = E, interaction = interaction, fractions = fractions / sum(fractions))
  })
  gxe[[as.character(time)]] <- list(time = time,
    observed = matrices[[1]]$fractions, predicted = matrices[[2]]$fractions,
    correlation = c(G = cor(matrices[[1]]$G, matrices[[2]]$G), E = cor(matrices[[1]]$E, matrices[[2]]$E),
                    GxE = cor(c(matrices[[1]]$interaction), c(matrices[[2]]$interaction))))
}
data$gxe <- gxe

for (file in list.files("results", pattern = "^pooling_.*rds$", full.names = TRUE)) {
  result <- readRDS(file); key <- sub(".rds$", "", basename(file))
  data$pooling[[key]] <- list(metrics = result$evaluation$metrics, diagnostics = diagnostic_summary(result$diagnostics),
                            training_years = result$training_years)
  fwrite(result$diagnostics, paste0("results/diagnostics_", key, ".csv"))
}
for (extension in c("location", "svd")) {
  result <- readRDS(paste0("results/", extension, "_1983_1990.rds"))
  data[[extension]] <- list(metrics = result$evaluation$metrics, diagnostics = diagnostic_summary(result$diagnostics),
     features = result$fit$p, fold = "1983–1990", energy = result$design$energy)
  fwrite(result$diagnostics, paste0("results/diagnostics_", extension, "_1983_1990.csv"))
  data$diagnostics[[paste0(extension, "_1983_1990")]] <- diagnostic_summary(result$diagnostics)
}
spectral <- svd(full$design$X[intervals$year > 1990, -1], nu = 0, nv = 0)
data$svd$full_energy <- cumsum(spectral$d^2) / sum(spectral$d^2)
  code_functions <- c("read_apsim", "observation_intervals", "make_basis", "evaluate_basis", "integrated_features", "prepare_design",
                    "draw_gaussian", "joint_genotype_update", "fit_gibbs", "posterior_draws", "evaluate_fit", "ale_curve", "diagnose_matrix")
data$code <- setNames(lapply(code_functions, function(name) paste(deparse(get(name), width.cutoff = 95), collapse = "\n")), code_functions)
manifest <- list(created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE), specification = data$specification,
                 source = data$source, counts = data$counts, diagnostics = data$diagnostics,
                 main_model = "hierarchical", interpretation_fit = "full",
                 evaluation = "four held-out year blocks; all 199 genotypes, 124 environments",
                 notes = c("No artificial trajectories or added measurement noise", "No markers or APSIM crop-state predictors",
                   "Only first-block CV for optional SVD/location extensions", "No positivity constraint on predicted growth rate"))
write_json(manifest, "results/manifest.json", pretty = TRUE, auto_unbox = TRUE, digits = 9)
writeLines(paste0("window.WORKSHOP_DATA=", toJSON(data, auto_unbox = TRUE, digits = 7, dataframe = "rows", null = "null"), ";"), "assets/results.js")
message("Exported ", length(data$trajectories), " held-out trajectories; all figures from R / APSIM.")
