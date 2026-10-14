library(data.table)
source("analysis/03_gibbs.R")
r <- readRDS("results/location_1983_1990.rds"); prepared <- readRDS("results/prepared.rds")
rows <- which(prepared$intervals$genotype %in% c("g008", "g052") &
              prepared$intervals$location %in% c("Emerald", "Yanco") & prepared$intervals$year == 1997)
X <- r$design$X[rows, ]; y <- prepared$intervals$dy[rows]; meta <- prepared$intervals[rows]
g <- match(meta$genotype, sort(unique(meta$genotype))); l <- match(meta$location, sort(unique(meta$location)))
stats <- sufficient_statistics(X, y, g); site_stats <- sufficient_statistics(X, y, l)
pair <- sufficient_statistics(X, y, paste(g, l, sep = ":")); p <- ncol(X); G <- length(stats); L <- length(site_stats)
contrasts <- qr.Q(qr(contr.helmert(L)))
site_precision <- Reduce(`+`, lapply(1:L, function(l) kronecker(tcrossprod(contrasts[l, ]), site_stats[[l]]$A)))
cross_design <- lapply(1:G, function(g) do.call(cbind, lapply(1:(L-1), function(j)
  Reduce(`+`, lapply(1:L, function(l) pair[[paste(g, l, sep = ":")]]$A * contrasts[l,j])))))
target <- vapply(site_stats, function(s) s$t, numeric(p))
draw <- joint_location_update(stats, rep(1,p), rep(100,p), .2, cross_design, site_precision, target, contrasts)
q <- nrow(site_precision); m <- p+q; Q <- matrix(0,m+G*p,m+G*p); h <- numeric(nrow(Q))
Q[1:p,1:p] <- diag(1+G*100,p); jj <- p+1:q
Q[jj,jj] <- site_precision/.2 + diag(100,q); h[jj] <- as.vector(target%*%contrasts)/.2
for(g in 1:G){
  ii <- m + (g-1)*p + 1:p; link <- rbind(-diag(100,p), t(cross_design[[g]])/.2)
  Q[ii,ii] <- stats[[g]]$A/.2 + diag(100,p); Q[1:m,ii] <- link; Q[ii,1:m] <- t(link)
  h[ii] <- stats[[g]]$t/.2
}
V <- solve(Q); mean <- solve(Q,h)
stopifnot(max(abs(solve(draw$common_precision)-V[1:m,1:m]))<1e-7,
          max(abs(solve(draw$common_precision,draw$common_target)-mean[1:m]))<1e-7)
message("Joint population / genotype / location conditional matches dense posterior moments.")
