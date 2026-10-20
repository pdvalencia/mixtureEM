# Negative binomial (NB2): the growth family (fit_lcga(family = "negbin")) and
# the cross-sectional count emission (measurement = "negbin").

off <- list(latent = 0, categorical = 0, poisson = 0, variances = 0)

sim_nb <- function(n, Tn, b0, b1, alpha, seed) {
  set.seed(seed)
  mu <- exp(b0 + b1 * (seq_len(Tn) - 1))
  matrix(rnbinom(n * Tn, mu = rep(mu, each = n), size = 1 / alpha), n, Tn)
}

test_that("an overdispersion at the Poisson boundary reproduces the Poisson likelihood", {
  y <- sim_nb(150, 4, 0.5, 0.1, 0.5, 1)
  y[3, 2] <- NA
  D <- mixtureEM:::.lcga_design(0:3, 1)

  po <- lcga_model(2, D, "poisson")
  po$parameters <- list(coefs = rbind(c(0.2, 0.1), c(1, -0.05)), dispersion = c(1, 1))
  nb <- lcga_model(2, D, "negbin")
  nb$parameters <- c(po$parameters, list(overdispersion = matrix(1e-13, 2, 4)))
  expect_equal(log_likelihood(nb, y), log_likelihood(po, y), tolerance = 1e-12)
  # Approaching the boundary from inside the NB density, too.
  nb$parameters$overdispersion[] <- 1e-9
  expect_equal(log_likelihood(nb, y), log_likelihood(po, y), tolerance = 1e-6)

  pc <- poisson_model(2, "poisson_nan")
  pc$parameters$rates <- rbind(c(1, 2, 3, 4), c(.5, .5, 1, 2))
  nc <- poisson_model(2, "negbin_nan")
  nc$parameters <- list(rates = pc$parameters$rates,
                        overdispersion = matrix(1e-13, 2, 4))
  expect_equal(log_likelihood(nc, y), log_likelihood(pc, y), tolerance = 1e-12)
})

test_that("one-class overdispersed counts: Poisson over-extracts, NB does not", {
  y <- sim_nb(400, 5, 1, 0, 1, 2)
  bic <- function(fam, k)
    fit_lcga(y, times = 5, n_classes = k, family = fam, n_init = 5,
             random_state = 1, n_cores = 1, bayes_constants = off)$metrics$bic
  expect_lt(bic("poisson", 2), bic("poisson", 1))
  nb <- vapply(1:2, function(k) bic("negbin", k), numeric(1))
  expect_lt(nb[1], nb[2])
})

test_that("an NB growth fit recovers its classes under both alpha structures", {
  y <- rbind(sim_nb(300, 4, 2.3, 0.1, 0.1, 3), sim_nb(300, 4, 0.5, 0.2, 0.1, 4))
  fit <- function(d)
    fit_lcga(y, times = 4, n_classes = 2, family = "negbin", dispersion = d,
             n_init = 5, random_state = 1, n_cores = 1, bayes_constants = off)
  fo <- fit("occasion")
  fc <- fit("class")
  expect_equal(fo$metrics$n_params, 4 + 4 + 1)
  expect_equal(fc$metrics$n_params, 4 + 2 + 1)
  hi <- which.max(fo$growth$coefficients[, 1])
  expect_equal(unname(fo$growth$coefficients[hi, ]), c(2.3, 0.1), tolerance = 0.1)
  # Shared across the classes, or across the occasions.
  expect_equal(fo$growth$overdispersion[1, ], fo$growth$overdispersion[2, ])
  expect_equal(unname(apply(fc$growth$overdispersion, 1, sd)), c(0, 0))
  expect_equal(mean(fo$growth$overdispersion), 0.1, tolerance = 0.3)

  ge <- growth_estimates(fo)
  od <- ge[ge$parameter == "overdispersion", ]
  expect_equal(nrow(od), 2 * 4)
  expect_true(all(is.finite(od$se)))
  expect_output(print(fo), "OVERDISPERSION")
  expect_output(print(fc), "per class")
  expect_identical(mixtureEM:::.emission_family(fo), "negative binomial")
})

test_that("cross-sectional NB fits, with missing cells, and reports alpha", {
  y <- rbind(sim_nb(250, 4, 1.8, 0, 0.3, 5), sim_nb(250, 4, 0.2, 0, 0.3, 6))
  y[sample(length(y), 40)] <- NA
  f <- fit_mixture(y, n_classes = 2, measurement = "negbin", n_init = 5,
                   random_state = 1, n_cores = 1, bayes_constants = off)
  expect_equal(f$metrics$n_params, 2 * 4 + 4 + 1)
  me <- suppressMessages(measurement_estimates(f))
  od <- me[me$parameter == "overdispersion", ]
  expect_equal(nrow(od), 2 * 4)
  expect_true(all(is.finite(od$se)))
  expect_equal(od$estimate[od$class == 1], od$estimate[od$class == 2])
  expect_output(measurement_summary(f), "OVERDISPERSION")
})
