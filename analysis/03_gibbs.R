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
  S <- diag(lambda0 + G * lambdaG, p); h <- numeric(p)
  # Elementwise diagonal multiplication replaces two dense matrix products.
  precision_outer <- tcrossprod(lambdaG)
  diagonal_indices <- seq.int(1L, p * p, by = p + 1L)
  inverses <- vector("list", G); targets <- vector("list", G)
  for (g in seq_len(G)) {
    A <- stats[[g]]$A / sigma2
    target <- stats[[g]]$t / sigma2
    if (!is.null(offsets)) target <- target - offsets[, g] / sigma2
    A[diagonal_indices] <- A[diagonal_indices] + lambdaG
    R <- chol(A); V <- chol2inv(R)
    inverses[[g]] <- list(R = R, V = V)
    targets[[g]] <- target
    S <- S - V * precision_outer
    h <- h + lambdaG * drop(V %*% target)
  }
  beta <- draw_gaussian((S + t(S)) / 2, h)
  effective <- vapply(seq_len(G), function(g) {
    V <- inverses[[g]]$V; R <- inverses[[g]]$R
    drop(V %*% (targets[[g]] + lambdaG * beta) + backsolve(R, rnorm(p)))
  }, numeric(p))
  list(beta = beta, effective = effective)
}

joint_location_update <- function(stats, lambda0, lambdaG, sigma2, cross_design,
                                  site_precision, site_target, contrasts) {
  # Jointly sample population, genotype and location coefficients. This is the
  # SAME Gaussian posterior as alternating location updates, with less dependence.
  p <- length(lambda0); G <- length(stats); q <- nrow(site_precision)
  Sbb <- diag(lambda0 + G * lambdaG, p); Sbc <- matrix(0, p, q)
  Scc <- site_precision / sigma2 + diag(100, q)
  hb <- numeric(p); hc <- as.vector(site_target %*% contrasts) / sigma2
  D <- diag(lambdaG, p); parts <- vector("list", G)
  precision_outer <- tcrossprod(lambdaG)
  for (g in seq_len(G)) {
    R <- chol(stats[[g]]$A / sigma2 + D); V <- chol2inv(R)
    target <- stats[[g]]$t / sigma2
    K <- cross_design[[g]] / sigma2; VK <- V %*% K; Vt <- drop(V %*% target)
    Sbb <- Sbb - V * precision_outer
    Sbc <- Sbc + VK * lambdaG
    Scc <- Scc - crossprod(K, VK)
    hb <- hb + lambdaG * Vt; hc <- hc - drop(crossprod(K, Vt))
    parts[[g]] <- list(R = R, V = V, target = target, K = K)
  }
  S <- rbind(cbind(Sbb, Sbc), cbind(t(Sbc), Scc)); h <- c(hb, hc)
  jj <- p + seq_len(q)
  common <- draw_gaussian((S + t(S)) / 2, h)
  effective <- vapply(parts, function(s)
    drop(s$V %*% (s$target + lambdaG * common[1:p] - s$K %*% common[jj]) + backsolve(s$R, rnorm(p))), numeric(p))
  list(beta = common[1:p], effective = effective,
       site = matrix(common[jj], p, ncol(contrasts)) %*% t(contrasts),
       common_precision = S, common_target = h)
}

fit_gibbs <- function(X, y, genotype, location = NULL, mode = "hierarchical",
                      warmup = 1200L, samples = 1200L, chains = 4L, seed = 20261014L,
                      resume_fit = NULL) {
  levels <- sort(unique(genotype)); group <- match(genotype, levels)
  stats <- sufficient_statistics(X, y, group)
  p <- ncol(X); G <- length(stats); n <- length(y)
  site_levels <- if (is.null(location)) character() else sort(unique(location))
  L <- length(site_levels)
  contrasts <- if (L > 1) qr.Q(qr(contr.helmert(L))) else matrix(0, L, 0)
  site_group <- if (L > 0) match(location, site_levels) else NULL
  pair_stats <- if (L > 0) sufficient_statistics(X, y, paste(group, site_group, sep = ":")) else NULL
  cross <- if (L > 0) lapply(seq_len(G), function(g) lapply(seq_len(L), function(l)
    pair_stats[[paste(g, l, sep = ":")]])) else NULL
  # Stack genotype crossproducts per site for BLAS matrix-vector operations.
  cross_stacked <- if (L > 0) lapply(seq_len(L), function(l)
    do.call(rbind, lapply(seq_len(G), function(g) cross[[g]][[l]]$A))) else NULL
  site_stats <- if (L > 0) sufficient_statistics(X, y, site_group) else NULL
  # The location likelihood precision is fixed apart from sigma2. Build it once.
  site_precision <- if (L > 1) Reduce(`+`, lapply(seq_len(L), function(l)
    kronecker(tcrossprod(contrasts[l, ]), site_stats[[l]]$A))) else NULL
  cross_design <- if (L > 1) lapply(seq_len(G), function(g)
    do.call(cbind, lapply(seq_len(L - 1), function(j)
      Reduce(`+`, lapply(seq_len(L), function(l) cross[[g]][[l]]$A * contrasts[l, j]))))) else NULL
  site_target <- if (L > 1) vapply(site_stats, function(s) s$t, numeric(p)) else NULL
  all_stats <- sufficient_statistics(X, y, rep(1L, n))[[1]]
  output <- vector("list", chains)
  for (chain in seq_len(chains)) {
    set.seed(seed + chain * 1000L)
    lambda0 <- rep(1, p); lambdaG <- rep(100, p); sigma2 <- var(y) * runif(1, .4, 1.4)
    beta <- rnorm(p, 0, .1); effective <- matrix(0, p, G); site <- matrix(0, p, L)
    if (!is.null(resume_fit)) {
      previous <- resume_fit$chains[[chain]]; last <- resume_fit$samples
      lambda0 <- previous$lambda0[last, ]; lambdaG <- previous$lambdaG[last, ]
      sigma2 <- previous$sigma2[last]
      if (L > 0) site <- matrix(previous$site[last, , ], p, L)
    }
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
        draw <- if (L > 1) joint_location_update(stats, lambda0, lambdaG, sigma2,
          cross_design, site_precision, site_target, contrasts) else
          joint_genotype_update(stats, lambda0, lambdaG, sigma2)
        beta <- draw$beta; effective <- draw$effective
        if (L > 1) site <- draw$site
        lambdaG <- rgamma(p, shape = 2 + G / 2,
                           rate = .02 + rowSums((effective - beta)^2) / 2)
      }
      lambda0 <- rgamma(p, shape = 2.5, rate = 1 + beta^2 / 2)
      if (L > 1) {
        site_targets <- vapply(seq_len(L), function(l)
          site_stats[[l]]$t - drop(crossprod(cross_stacked[[l]], as.vector(effective))), numeric(p))
      }
      residual_ss <- sum(vapply(seq_len(G), function(g) {
        a <- effective[, g]
        stats[[g]]$yy - 2 * sum(a * stats[[g]]$t) + sum(a * (stats[[g]]$A %*% a))
      }, numeric(1)))
      if (L > 1) residual_ss <- residual_ss + sum(vapply(seq_len(L), function(l)
        -2 * sum(site[, l] * site_targets[, l]) +
          sum(site[, l] * (site_stats[[l]]$A %*% site[, l])), numeric(1)))
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
