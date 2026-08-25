# `liang_ark_sim` is simulated data, drawn from the published parameters of
# an existing three-class latent profile solution rather than from
# participant responses. The original study data are not redistributable;
# the generator and the full generating-parameter tables that produced this
# file are kept outside the distributed package.

liang_ark_sim <- read.csv("data-raw/liang_ark_sim/liang_ark_sim.csv")
liang_ark_sim$boundary <- as.logical(liang_ark_sim$boundary)

usethis::use_data(liang_ark_sim, overwrite = TRUE)
