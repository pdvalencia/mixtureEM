# Zero-inflated Poisson: the growth family (fit_lcga(family = "zip")) and the
# cross-sectional count emission (measurement = "zip").

off <- list(latent = 0, categorical = 0, poisson = 0, variances = 0)

sim_zip <- function(n, Tn, b0, b1, pz, seed) {
  set.seed(seed)
  mu <- exp(b0 + b1 * (seq_len(Tn) - 1))
  y  <- matrix(rpois(n * Tn, rep(mu, each = n)), n, Tn)
  y[matrix(runif(n * Tn) < pz, n, Tn)] <- 0
  y
}

test_that("zero inflation fixed at zero reproduces the Poisson likelihood exactly", {
  y <- sim_zip(150, 4, 0.5, 0.1, 0.3, 1)
  y[3, 2] <- NA
  D <- mixtureEM:::.lcga_design(0:3, 1)

  po <- lcga_model(2, D, "poisson")
  po$parameters <- list(coefs = rbind(c(0.2, 0.1), c(1, -0.05)), dispersion = c(1, 1))
  zi <- lcga_model(2, D, "zip")
  zi$parameters <- c(po$parameters, list(inflation = matrix(-Inf, 2, 1)))
  expect_equal(log_likelihood(zi, y), log_likelihood(po, y), tolerance = 1e-12)

  pc <- poisson_model(2, "poisson_nan")
  pc$parameters$rates <- rbind(c(1, 2, 3, 4), c(.5, .5, 1, 2))
  zc <- poisson_model(2, "zip_nan")
  zc$parameters <- list(rates = pc$parameters$rates, inflation = matrix(0, 2, 4))
  expect_equal(log_likelihood(zc, y), log_likelihood(pc, y), tolerance = 1e-12)
})

test_that("one-class zero-inflated counts: Poisson over-extracts, ZIP does not", {
  y <- sim_zip(400, 5, 1, 0, 0.35, 2)
  bic <- function(fam, k)
    fit_lcga(y, times = 5, n_classes = k, family = fam, n_init = 5,
             random_state = 1, n_cores = 1, bayes_constants = off)$metrics$bic
  expect_lt(bic("poisson", 2), bic("poisson", 1))
  zb <- vapply(1:2, function(k) bic("zip", k), numeric(1))
  expect_lt(zb[1], zb[2])
})

test_that("a ZIP growth fit recovers its classes and reports the inflation", {
  y <- rbind(sim_zip(300, 5, 1.2, 0.1, 0.2, 3), sim_zip(300, 5, -0.5, 0.2, 0.6, 4))
  f <- fit_lcga(y, times = 5, n_classes = 2, family = "zip", n_init = 10,
                random_state = 1, n_cores = 1, bayes_constants = off)
  expect_equal(f$metrics$n_params, 7)
  hi <- which.max(f$growth$coefficients[, 1])
  expect_equal(unname(f$growth$coefficients[hi, ]), c(1.2, 0.1), tolerance = 0.15)
  expect_equal(unname(f$growth$inflation_probability[c(hi, 3 - hi)]), c(0.2, 0.6),
               tolerance = 0.25)
  ge <- growth_estimates(f)
  expect_true(all(is.finite(ge$se[ge$parameter == "inflation_logit"])))
  expect_output(print(f), "ZERO INFLATION")
  expect_identical(mixtureEM:::.emission_family(f), "zero-inflated Poisson")
})

test_that("cross-sectional ZIP fits, with missing cells, and reports the inflation", {
  y <- rbind(sim_zip(250, 4, 1.3, 0, 0.3, 5), sim_zip(250, 4, -0.3, 0, 0.1, 6))
  y[sample(length(y), 40)] <- NA
  f <- fit_mixture(y, n_classes = 2, measurement = "zip", n_init = 5,
                   random_state = 1, n_cores = 1, bayes_constants = off)
  expect_equal(f$metrics$n_params, 17)
  # One item of the low class has no excess zeros in this sample: its
  # inflation sits at 0, is flagged, and only it loses its standard error.
  me <- suppressMessages(measurement_estimates(f))
  infl <- me[me$parameter == "inflation", ]
  expect_equal(sum(infl$boundary), 1L)
  expect_true(all(is.finite(infl$se[!infl$boundary])))
  expect_true(all(is.na(infl$se[infl$boundary])))
  expect_output(measurement_summary(f), "ZERO INFLATION")
})
