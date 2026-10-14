# Split, rank-normalised R-hat and Geyer initial-positive-sequence ESS.
split_chains <- function(x) {
  n <- floor(nrow(x) / 2)
  cbind(x[seq_len(n), , drop = FALSE], x[(nrow(x) - n + 1):nrow(x), , drop = FALSE])
}

basic_rhat <- function(x) {
  n <- nrow(x); W <- mean(apply(x, 2, var)); B <- n * var(colMeans(x))
  if (W == 0) return(if (B == 0) 1 else Inf)
  sqrt(((n - 1) * W / n + B / n) / W)
}

chain_ess <- function(x) {
  n <- nrow(x); m <- ncol(x)
  W <- mean(apply(x, 2, var)); B <- n * var(colMeans(x))
  vplus <- (n - 1) * W / n + B / n
  if (vplus < 1e-20) return(n * m)
  covariances <- vapply(seq_len(m), function(j) {
    v <- x[, j] - mean(x[, j]); size <- 2^ceiling(log2(2 * n))
    padded <- c(v, rep(0, size - n)); transform <- fft(padded)
    Re(fft(transform * Conj(transform), inverse = TRUE))[seq_len(n)] / size / n
  }, numeric(n))
  rho <- 1 - (W - rowMeans(covariances)) / vplus; rho[1] <- 1
  pairs <- rho[seq(1, n - 1, 2)] + rho[seq(2, n, 2)]
  first_negative <- which(pairs <= 0)[1]
  if (!is.na(first_negative)) pairs <- head(pairs, first_negative - 1)
  if (length(pairs) == 0) return(n * m)
  pairs <- cummin(pairs)
  min(n * m, n * m / max(1, -1 + 2 * sum(pairs)))
}

diagnose_matrix <- function(x) {
  x <- split_chains(x)
  rank_normalise <- function(a) matrix(qnorm((rank(a, ties.method = "average") - 3/8) /
                                             (length(a) + 1/4)), nrow(a), ncol(a))
  ranked <- rank_normalise(x)
  folded <- rank_normalise(abs(x - median(x)))
  limits <- quantile(x, c(.05, .95))
  c(rhat = max(basic_rhat(ranked), basic_rhat(folded)),
    ess_bulk = chain_ess(ranked),
    ess_tail = min(chain_ess(1 * (x <= limits[1])), chain_ess(1 * (x >= limits[2]))))
}

diagnose_fit <- function(fit) {
  rows <- list(); k <- 0L
  add <- function(name, matrices) {
    k <<- k + 1L
    rows[[k]] <<- as.data.table(c(list(parameter = name),
      as.list(diagnose_matrix(do.call(cbind, matrices)))))
  }
  for (j in seq_len(fit$p)) {
    add(paste0("beta[", j, "]"), lapply(fit$chains, function(c) c$beta[, j]))
    if (fit$mode != "independent") add(paste0("lambda0[", j, "]"), lapply(fit$chains, function(c) c$lambda0[, j]))
    if (fit$mode == "hierarchical") add(paste0("lambdaG[", j, "]"), lapply(fit$chains, function(c) c$lambdaG[, j]))
    if (fit$mode != "shared") for (g in seq_along(fit$genotypes))
      add(paste0("a[", fit$genotypes[g], ",", j, "]"), lapply(fit$chains, function(c) c$effective[, j, g]))
    for (l in seq_along(fit$locations))
      add(paste0("site[", fit$locations[l], ",", j, "]"), lapply(fit$chains, function(c) c$site[, j, l]))
  }
  add("sigma2", lapply(fit$chains, function(c) c$sigma2))
  rbindlist(rows)
}

posterior_mean <- function(fit) {
  effective <- Reduce(`+`, lapply(fit$chains, function(c) apply(c$effective, c(2, 3), mean))) / length(fit$chains)
  site <- if (length(fit$locations)) Reduce(`+`, lapply(fit$chains, function(c) apply(c$site, c(2, 3), mean))) /
            length(fit$chains) else matrix(0, fit$p, 0)
  list(effective = effective, site = site,
       beta = Reduce(`+`, lapply(fit$chains, function(c) colMeans(c$beta))) / length(fit$chains))
}

posterior_draws <- function(fit, genotype, location, count = 160L) {
  # Equal deterministic selections from all chains, preserving each joint draw.
  per_chain <- ceiling(count / length(fit$chains))
  selected <- unique(round(seq(1, fit$samples, length.out = per_chain)))
  g <- match(genotype, fit$genotypes); l <- match(location, fit$locations)
  theta <- do.call(rbind, lapply(fit$chains, function(c) {
    a <- c$effective[selected, , g, drop = FALSE]; dim(a) <- c(length(selected), fit$p)
    if (!is.na(l)) {
      b <- c$site[selected, , l, drop = FALSE]; dim(b) <- dim(a); a <- a + b
    }
    a
  }))
  list(theta = theta, sigma2 = unlist(lapply(fit$chains, function(c) c$sigma2[selected])))
}

evaluate_fit <- function(fit, design, intervals, test_rows) {
  means <- posterior_mean(fit)
  increment <- numeric(nrow(intervals))
  groups <- split(seq_len(nrow(intervals)), paste(intervals$genotype, intervals$location, sep = "|"))
  for (ii in groups) {
    g <- match(intervals$genotype[ii[1]], fit$genotypes)
    l <- match(intervals$location[ii[1]], fit$locations)
    theta <- means$effective[, g]
    if (!is.na(l)) theta <- theta + means$site[, l]
    increment[ii] <- drop(design$X[ii, , drop = FALSE] %*% theta)
  }
  predictions <- copy(intervals)
  predictions[, `:=`(predicted_increment = increment, held_out = test_rows)]
  predictions[, `:=`(observed = cumsum(dy), predicted = cumsum(predicted_increment)),
              by = .(environment, genotype)]
  metric <- function(d, label) {
    error <- d$predicted - d$observed
    data.table(scope = label, n = nrow(d), rmse = sqrt(mean(error^2)),
      r2 = 1 - sum(error^2) / sum((d$observed - mean(d$observed))^2),
      correlation = cor(d$observed, d$predicted))
  }
  test <- predictions[held_out == TRUE]
  terminal <- test[, .SD[.N], by = .(environment, genotype)]
  list(predictions = predictions,
       metrics = rbind(metric(test, "trajectory checkpoints"), metric(terminal, "terminal biomass")))
}

ale_curve <- function(weather, design, theta, driver, bins = 10L) {
  boundaries <- unique(as.numeric(quantile(weather[[driver]], seq(0, 1, length.out = bins + 1))))
  assignments <- findInterval(weather[[driver]], boundaries, all.inside = TRUE)
  differences <- matrix(0, length(boundaries) - 1, nrow(theta)); counts <- integer(nrow(differences))
  for (k in seq_len(nrow(differences))) {
    ii <- which(assignments == k); counts[k] <- length(ii)
    lower <- copy(weather[ii]); upper <- copy(lower)
    lower[[driver]] <- rep(boundaries[k], nrow(lower)); upper[[driver]] <- rep(boundaries[k + 1], nrow(upper))
    delta <- rate_features(upper, design) - rate_features(lower, design)
    differences[k, ] <- drop(colMeans(delta) %*% t(theta))
  }
  accumulated <- apply(differences, 2, cumsum)
  # Bin-centre interpolation and empirical-distribution centring.
  at_centres <- accumulated - differences / 2
  at_centres <- sweep(at_centres, 2, drop(crossprod(counts / sum(counts), at_centres)), "-")
  list(x = (head(boundaries, -1) + tail(boundaries, -1)) / 2,
       counts = counts, boundaries = boundaries, draws = at_centres)
}
