# Fast final integrity pass. Recompute selected summaries from saved JOINT draws.
# This also permits completing a long optional fit after the main export.
library(data.table);library(jsonlite)
for(file in c("01_data.R","02_design.R","03_gibbs.R","04_posterior.R"))source(file.path("analysis",file))
args<-commandArgs(trailingOnly=TRUE)
input<-if(length(args))args[1]else "assets/results.js"
text <- paste(readLines(input),collapse="\n")
data <- fromJSON(sub(";$","",sub("^window.WORKSHOP_DATA=","",text)),simplifyVector=FALSE)
prepared <- readRDS("results/prepared.rds"); intervals <- prepared$intervals
raw <- read_apsim("data/SimulatedData124Envs_Biomass.txt")
# Validate exported daily rows against exact genotype/environment identifiers.
for(key in names(data$trajectories)){
  selected <- data$trajectories[[key]]
  selected_genotype <- selected$genotype
  environment_id <- paste(selected$location,selected$year,sep="_")
  d <- raw$observations[.(environment_id,selected_genotype)]
  stopifnot(!anyDuplicated(d$time),max(d$time)==max(unlist(selected$time)))
  data$trajectories[[key]]$daily_time <- d$time
  data$trajectories[[key]]$daily_biomass <- d$biomass
}
rm(raw);gc()
full <- readRDS("results/hierarchical_full.rds"); path <- prepared$weather[environment=="Yanco_2010"]
Phi <- rate_features(path,full$design)
population <- do.call(rbind,lapply(full$fit$chains,function(c)c$beta[unique(round(seq(1,nrow(c$beta),length.out=40))),]))
summaries <- function(values){q<-apply(values,1,quantile,c(.025,.5,.975));list(lower=unname(q[1,]),median=unname(q[2,]),upper=unname(q[3,]))}
for(g in unlist(data$genotypes)){
  draw<-posterior_draws(full$fit,g,"Yanco")
  data$hierarchy$genotypes[[g]]$deviation<-summaries(Phi%*%t(draw$theta-population))
}
rm(full);gc()
for(file in list.files("results",pattern="^(shared|hierarchical)_[0-9].*rds$",full.names=TRUE)){
  r<-readRDS(file)
  for(key in names(data$trajectories)){
    c<-data$trajectories[[key]]
    if(c$model!=r$fit$mode || !c$year %in% r$holdout_years)next
    ii<-which(intervals$genotype==c$genotype & intervals$environment==paste(c$location,c$year,sep="_"))
    draw<-posterior_draws(r$fit,c$genotype,c$location)
    increment<-r$design$X[ii,,drop=FALSE]%*%colMeans(draw$theta)
    c$mean_biomass<-c(0,cumsum(increment))
    environment_id<-paste(c$location,c$year,sep="_")
    weather<-prepared$weather[environment==environment_id]
    rate_weather<-weather[time<=max(unlist(c$daily_time))]
    c$weather<-list(time=weather$time,temperature=weather$temperature,Radn=weather$Radn,Rain=weather$Rain,
      thermal_time=weather$thermal_time,rain14=weather$rain14)
    c$rate<-list(time=rate_weather$time,mean=drop(rate_features(rate_weather,r$design)%*%colMeans(draw$theta)))
    data$trajectories[[key]]<-c
  }
  rm(r);gc()
}
diagnostic_summary<-function(d)list(max_rhat=max(d$rhat),min_bulk_ess=min(d$ess_bulk),min_tail_ess=min(d$ess_tail),
  parameters=nrow(d),passed=max(d$rhat)<1.01&&min(d$ess_bulk)>=400&&min(d$ess_tail)>=400)
for(extension in c("location","svd")){
  r<-readRDS(paste0("results/",extension,"_1983_1990.rds"))
  data[[extension]]<-list(metrics=r$evaluation$metrics,diagnostics=diagnostic_summary(r$diagnostics),
    features=r$fit$p,fold="1983–1990",energy=r$design$energy,chains=length(r$fit$chains),
    warmup=r$fit$warmup,samples_per_chain=r$fit$samples)
  if(extension=="svd")data$svd$full_energy<-r$design$energy
  data$diagnostics[[paste0(extension,"_1983_1990")]]<-diagnostic_summary(r$diagnostics)
  fwrite(r$diagnostics,paste0("results/diagnostics_",extension,"_1983_1990.csv"))
}
stopifnot(all(vapply(data$diagnostics,function(d)isTRUE(d$passed),logical(1))))
data$code<-setNames(lapply(names(data$code),function(name)paste(deparse(get(name),width.cutoff=95),collapse="\n")),names(data$code))
writeLines(paste0("window.WORKSHOP_DATA=",toJSON(data,auto_unbox=TRUE,digits=7,dataframe="rows",null="null"),";"),"assets/results.js")
manifest<-fromJSON("results/manifest.json",simplifyVector=FALSE);manifest$diagnostics<-data$diagnostics
manifest$extensions<-list(location=data$location,svd=data$svd)
manifest$analysis_file_md5<-as.list(tools::md5sum(list.files("analysis",pattern="[.]R$",full.names=TRUE)))
write_json(manifest,"results/manifest.json",auto_unbox=TRUE,pretty=TRUE,digits=9)
message("Final joint-draw summaries and all convergence checks passed.")
