# Tests use observed APSIM intervals. No fabricated workshop results.
library(data.table)
for (file in c("01_data.R", "02_design.R", "03_gibbs.R", "04_posterior.R")) source(file.path("analysis", file))
p <- readRDS("results/prepared.rds"); intervals <- p$intervals
training <- intervals$year > 1990
design <- prepare_design(p$weather, intervals, training)
raw_integrals <- integrated_features(p$weather, intervals, design$basis)
stopifnot(max(abs(raw_integrals[, 1] - (intervals$t1 - intervals$t0))) < 1e-10)
# Training transformations cannot depend on held-out weather.
w <- copy(p$weather); w[year <= 1990, temperature := temperature + 100]
b1 <- make_basis(p$weather[year > 1990]); b2 <- make_basis(w[year > 1990])
stopifnot(identical(b1, b2))
# Compare the block-Schur conditional with a dense Gaussian calculation,
# using real sufficient statistics. Test deterministic moments, not RNG equality.
rows <- which(training & intervals$genotype %in% c("g008", "g052"))[1:80]
X <- design$X[rows, ]; y <- intervals$dy[rows]
groups <- match(intervals$genotype[rows], sort(unique(intervals$genotype[rows])))
stats <- sufficient_statistics(X, y, groups); P <- ncol(X); G <- length(stats)
D <- diag(100, P); P0 <- diag(1, P); sigma2 <- .2
Q <- matrix(0, P * (G + 1), P * (G + 1)); h <- numeric(nrow(Q))
Q[1:P, 1:P] <- P0 + G * D
S <- P0; v <- numeric(P)
for (g in seq_len(G)) {
  ii <- P * g + seq_len(P); C <- stats[[g]]$A / sigma2 + D; V <- solve(C)
  Q[ii, ii] <- C; Q[1:P, ii] <- -D; Q[ii, 1:P] <- -D
  h[ii] <- stats[[g]]$t / sigma2
  S <- S + D - D %*% V %*% D; v <- v + drop(D %*% V %*% h[ii])
}
dense_mean <- solve(Q, h); dense_covariance <- solve(Q)
stopifnot(max(abs(solve(S, v) - dense_mean[1:P])) < 1e-7,
          max(abs(solve(S) - dense_covariance[1:P, 1:P])) < 1e-7)
message("Data, integration, leakage and Gaussian conditional checks passed.")
