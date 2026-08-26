# `liang_park_sim` is simulated data, drawn from the published parameters of
# an existing three-class latent profile solution rather than from
# participant responses. The original study data are not redistributable;
# the generator and the full generating-parameter tables that produced this
# file are kept outside the distributed package.

liang_park_sim <- read.csv("data-raw/liang_park_sim/liang_park_sim.csv")
liang_park_sim$boundary <- as.logical(liang_park_sim$boundary)

usethis::use_data(liang_park_sim, overwrite = TRUE)
