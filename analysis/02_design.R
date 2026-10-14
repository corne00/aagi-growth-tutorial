library(splines) # Included with R; natural cubic splines.

make_basis <- function(training_weather) {
  drivers <- c("temperature", "thermal_time", "Radn", "rain14")
  specifications <- lapply(drivers, function(name) {
    x <- training_weather[[name]]
    B <- ns(x, df = 3)
    list(name = name, knots = attr(B, "knots"), boundary = attr(B, "Boundary.knots"))
  })
  scales <- lapply(c("thermal_time", "Radn", "rain14", "day_length"), function(name)
    c(mean = mean(training_weather[[name]]), sd = sd(training_weather[[name]])))
  names(scales) <- c("thermal_time", "Radn", "rain14", "day_length")
  list(specifications = specifications, scales = scales)
}

evaluate_basis <- function(weather, basis) {
  blocks <- lapply(basis$specifications, function(s) {
    B <- ns(weather[[s$name]], knots = s$knots, Boundary.knots = s$boundary)
    colnames(B) <- paste0(s$name, "_s", seq_len(ncol(B)))
    B
  })
  z <- function(name) (weather[[name]] - basis$scales[[name]]["mean"]) /
                       basis$scales[[name]]["sd"]
  cbind(intercept = 1, do.call(cbind, blocks), day_length = z("day_length"),
        thermal_radiation = z("thermal_time") * z("Radn"),
        thermal_rain = z("thermal_time") * z("rain14"))
}

integrated_features <- function(weather, intervals, basis) {
  Phi <- evaluate_basis(weather, basis)
  # Trapezoidal quadrature on the observed daily weather grid.
  cumulative <- lapply(split(seq_len(nrow(weather)), weather$environment), function(ii) {
    P <- Phi[ii, , drop = FALSE]
    daily <- (P[-1, , drop = FALSE] + P[-nrow(P), , drop = FALSE]) / 2
    rbind(0, apply(daily, 2, cumsum))
  })
  X <- matrix(NA_real_, nrow(intervals), ncol(Phi), dimnames = list(NULL, colnames(Phi)))
  for (env in unique(intervals$environment)) {
    ii <- which(intervals$environment == env)
    C <- cumulative[[env]]
    X[ii, ] <- C[intervals$t1[ii] + 1, , drop = FALSE] -
               C[intervals$t0[ii] + 1, , drop = FALSE]
  }
  stopifnot(all(is.finite(X)))
  X
}

prepare_design <- function(weather, intervals, training, pcr = FALSE) {
  basis <- make_basis(weather[year %in% unique(intervals$year[training])])
  X <- integrated_features(weather, intervals, basis)
  scale <- sqrt(colMeans(X[training, , drop = FALSE]^2))
  stopifnot(all(scale > 0))
  X <- sweep(X, 2, scale, "/")
  projection <- diag(ncol(X))
  energy <- NULL
  if (pcr) {
    # Keep the intercept separate; compress the nonconstant feature family.
    spectral <- svd(X[training, -1, drop = FALSE], nu = 0)
    energy <- cumsum(spectral$d^2) / sum(spectral$d^2)
    rank <- which(energy >= 0.99)[1]
    projection <- matrix(0, ncol(X), rank + 1)
    projection[1, 1] <- 1
    projection[-1, -1] <- spectral$v[, seq_len(rank), drop = FALSE]
    X <- X %*% projection
  }
  list(X = X, basis = basis, scale = scale, projection = projection, energy = energy)
}

rate_features <- function(weather, design) {
  sweep(evaluate_basis(weather, design$basis), 2, design$scale, "/") %*% design$projection
}
