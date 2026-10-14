library(data.table)
for(file in c("03_gibbs.R", "04_posterior.R"))source(file.path("analysis",file))
p <- readRDS("results/prepared.rds"); r <- readRDS("results/location_1983_1990.rds")
keep <- p$intervals$year > 1990
fit <- fit_gibbs(r$design$X[keep,],p$intervals$dy[keep],p$intervals$genotype[keep],
  location=p$intervals$location[keep],warmup=600,samples=1800,chains=4)
r$fit <- fit;r$diagnostics <- diagnose_fit(fit)
r$evaluation <- evaluate_fit(fit,r$design,p$intervals,!keep)
print(r$diagnostics[,.(max_rhat=max(rhat),min_bulk_ess=min(ess_bulk),min_tail_ess=min(ess_tail))])
saveRDS(r,"results/location_1983_1990.rds",compress=FALSE)
