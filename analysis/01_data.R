# All analysis scripts use free, open-source R packages only.
library(data.table)
setDTthreads(2)

read_apsim <- function(path) {
  columns <- c("Env", "loc", "year", "geno", "das", "biomass",
               "MinT", "MaxT", "Radn", "Rain", "vp", "day_length")
  d <- fread(path, select = columns)
  setnames(d, c("Env", "loc", "geno", "das", "biomass"),
              c("environment", "location", "genotype", "time", "biomass"))
  d[, biomass := biomass / 1000] # kg/ha -> t/ha
  stopifnot(all(d[time == 0, biomass] == 0))
  setkey(d, environment, genotype, time)
  stopifnot(!anyDuplicated(d[, .(environment, genotype, time)]),
            all(is.finite(as.matrix(d[, .(time, biomass, MinT, MaxT, Radn, Rain, vp, day_length)]))))
  # Only weather, never crop states or APSIM water-status outputs.
  weather <- unique(d[, .(environment, location, year, time, MinT, MaxT,
                          Radn, Rain, vp, day_length)])
  stopifnot(!anyDuplicated(weather[, .(environment, time)]))
  setorder(weather, environment, time)
  weather[, temperature := (MinT + MaxT) / 2]
  weather[, thermal_time := cumsum(pmax(temperature, 0)), by = environment]
  weather[, rain14 := frollsum(Rain, 14, align = "right", fill = NA_real_), by = environment]
  weather[, rain14 := fifelse(is.na(rain14), cumsum(Rain), rain14), by = environment]
  list(observations = d, weather = weather)
}

observation_intervals <- function(d, spacing = 14L) {
  # Include the native final observation; never extrapolate beyond the record.
  d[, {
    stopifnot(all(diff(time) == 1), time[1] == 0)
    end <- unique(c(seq.int(1L, .N, by = spacing), .N))
    data.table(start_index = head(end, -1), end_index = tail(end, -1),
      t0 = time[head(end, -1)], t1 = time[tail(end, -1)],
      dy = diff(biomass[end]))
  }, by = .(environment, location, year, genotype)]
}
