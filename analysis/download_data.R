# Public CC BY 4.0 dataset, DOI: 10.6084/m9.figshare.24843705.v1.
# Download ONLY the biomass/weather file; no markers are required.
dir.create("data", showWarnings = FALSE)
path <- "data/SimulatedData124Envs_Biomass.txt"
expected <- "b208ae77ec2c9f9e18dc605d62d2812e"
if (!file.exists(path)) download.file("https://ndownloader.figshare.com/files/43720239",
                                    path, mode = "wb", method = "libcurl")
stopifnot(unname(tools::md5sum(path)) == expected)
message("Verified APSIM source file: ", path)
