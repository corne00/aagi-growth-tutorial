# Gaussian random-effects / ARD model in effective genotype coefficients:
# a_g = beta + b_g;  a_g | beta, lambda_G ~ N(beta, diag(lambda_G)^-1).
# beta | lambda_0 ~ N(0, diag(lambda_0)^-1).
# lambda_0,k ~ Gamma(2, rate=1); lambda_G,k ~ Gamma(2, rate=.02).
# sigma2 ~ InvGamma(2, scale=.1). Priors do NOT depend on sigma2.
# Location coefficients use independent N(0, 1/100) priors in
# orthonormal sum-to-zero contrast coordinates. These precisions are fixed.

sufficient_statistics <- function(X, y, group) {
  lapply(split(seq_along(y), group), function(ii) list(
    A = crossprod(X[ii, , drop = FALSE]), t = drop(crossprod(X[ii, , drop = FALSE], y[ii])),
    yy = sum(y[ii]^2), n = length(ii)))
}

draw_gaussian <- function(Q, h) {
  R <- chol(Q)
  mean <- backsolve(R, forwardsolve(t(R), h))
  drop(mean + backsolve(R, rnorm(length(h))))
}

joint_genotype_update <- function(stats, lambda0, lambdaG, sigma2, offsets = NULL) {
  # Block-Schur Gaussian update: sample beta and every a_g JOINTLY.
  # This avoids the strong beta / b_g dependence of naive alternating updates.
  p <- length(lambda0); G <- length(stats)
  D <- diag(lambdaG, p); S <- diag(lambda0, p); h <- numeric(p)
  inverses <- vector("list", G); targets <- vector("list", G)
  for (g in seq_len(G)) {
    A <- stats[[g]]$A / sigma2
    target <- stats[[g]]$t / sigma2
    if (!is.null(offsets)) target <- target - offsets[[g]] / sigma2
    R <- chol(A + D); V <- chol2inv(R)
    inverses[[g]] <- list(R = R, V = V)
    targets[[g]] <- target
    S <- S + D - D %*% V %*% D
    h <- h + drop(D %*% V %*% target)
  }
  beta <- draw_gaussian((S + t(S)) / 2, h)
  effective <- vapply(seq_len(G), function(g) {
    V <- inverses[[g]]$V; R <- inverses[[g]]$R
    drop(V %*% (targets[[g]] + D %*% beta) + backsolve(R, rnorm(p)))
  }, numeric(p))
  list(beta = beta, effective = effective)
}

fit_gibbs <- function(X, y, genotype, location = NULL, mode = "hierarchical",
                      warmup = 1200L, samples = 1200L, chains = 4L, seed = 20261014L) {
  levels <- sort(unique(genotype)); group <- match(genotype, levels)
  stats <- sufficient_statistics(X, y, group)
  p <- ncol(X); G <- length(stats); n <- length(y)
  site_levels <- if (is.null(location)) character() else sort(unique(location))
  L <- length(site_levels)
  contrasts <- if (L > 1) qr.Q(qr(contr.helmert(L))) else matrix(0, L, 0)
  site_group <- if (L > 0) match(location, site_levels) else NULL
  cross <- if (L > 0) lapply(seq_len(G), function(g) lapply(seq_len(L), function(l) {
    ii <- which(group == g & site_group == l)
    list(A = crossprod(X[ii, , drop = FALSE]), t = drop(crossprod(X[ii, , drop = FALSE], y[ii])))
  })) else NULL
  site_stats <- if (L > 0) sufficient_statistics(X, y, site_group) else NULL
  all_stats <- sufficient_statistics(X, y, rep(1L, n))[[1]]
  output <- vector("list", chains)
  for (chain in seq_len(chains)) {
    set.seed(seed + chain * 1000L)
    lambda0 <- rep(1, p); lambdaG <- rep(100, p); sigma2 <- var(y) * runif(1, .4, 1.4)
    beta <- rnorm(p, 0, .1); effective <- matrix(0, p, G); site <- matrix(0, p, L)
    B <- matrix(NA_real_, samples, p); A <- array(NA_real_, c(samples, p, G))
    C <- array(NA_real_, c(samples, p, L)); LG <- matrix(NA_real_, samples, p)
    L0 <- matrix(NA_real_, samples, p); S2 <- numeric(samples)
    for (iteration in seq_len(warmup + samples)) {
      if (mode == "shared") {
        beta <- draw_gaussian(all_stats$A / sigma2 + diag(lambda0, p), all_stats$t / sigma2)
        effective <- matrix(beta, p, G)
      } else if (mode == "independent") {
        for (g in seq_len(G)) effective[, g] <- draw_gaussian(
          stats[[g]]$A / sigma2 + diag(1, p), stats[[g]]$t / sigma2)
        beta <- rowMeans(effective)
      } else {
        offsets <- if (L > 0) lapply(seq_len(G), function(g)
          Reduce(`+`, lapply(seq_len(L), function(l) drop(cross[[g]][[l]]$A %*% site[, l])))) else NULL
        draw <- joint_genotype_update(stats, lambda0, lambdaG, sigma2, offsets)
        beta <- draw$beta; effective <- draw$effective
        lambdaG <- rgamma(p, shape = 2 + G / 2,
                           rate = .02 + rowSums((effective - beta)^2) / 2)
      }
      lambda0 <- rgamma(p, shape = 2.5, rate = 1 + beta^2 / 2)
      if (L > 1) {
        Q <- diag(100, p * (L - 1)); h <- numeric(p * (L - 1))
        for (l in seq_len(L)) {
          K <- kronecker(matrix(contrasts[l, ], 1), diag(p))
          target <- site_stats[[l]]$t - Reduce(`+`, lapply(seq_len(G), function(g)
            drop(cross[[g]][[l]]$A %*% effective[, g])))
          Q <- Q + crossprod(K, site_stats[[l]]$A %*% K) / sigma2
          h <- h + drop(crossprod(K, target)) / sigma2
        }
        coordinates <- matrix(draw_gaussian(Q, h), p, L - 1)
        site <- coordinates %*% t(contrasts)
      }
      residual_ss <- sum(vapply(seq_len(G), function(g) {
        a <- effective[, g]
        if (L == 0) return(stats[[g]]$yy - 2 * sum(a * stats[[g]]$t) + sum(a * (stats[[g]]$A %*% a)))
        stats[[g]]$yy + sum(vapply(seq_len(L), function(l) {
          theta <- a + site[, l]; s <- cross[[g]][[l]]
          -2 * sum(theta * s$t) + sum(theta * (s$A %*% theta))
        }, numeric(1)))
      }, numeric(1)))
      stopifnot(residual_ss > -1e-6)
      sigma2 <- 1 / rgamma(1, shape = 2 + n / 2, rate = .1 + max(residual_ss, 0) / 2)
      if (iteration > warmup) {
        k <- iteration - warmup; B[k, ] <- beta; A[k, , ] <- effective
        if (L > 0) C[k, , ] <- site
        LG[k, ] <- lambdaG; L0[k, ] <- lambda0; S2[k] <- sigma2
      }
    }
    output[[chain]] <- list(beta = B, effective = A, site = C, lambdaG = LG, lambda0 = L0, sigma2 = S2)
    message("  chain ", chain, " complete: sigma = ", round(sqrt(mean(S2)), 4))
  }
  list(chains = output, genotypes = levels, locations = site_levels,
       mode = mode, warmup = warmup, samples = samples, seed = seed,
       p = p, n = n, priors = list(global = c(shape = 2, rate = 1),
       genotype = c(shape = 2, rate = .02), sigma2 = c(shape = 2, scale = .1)))
}
