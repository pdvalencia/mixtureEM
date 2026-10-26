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

test_that("the ZIP M-step is the per-record EM update, prior and linear inflation included", {
  # The M-step works on each class's occasion sums. It must land where the
  # per-record form does: the structural share z of every zero, the inflation
  # logit fitted to z with the class weights, the rate to the counts with
  # weights w (1 - z) plus the prior's pseudo-rows, each by glm.fit() at a
  # far tighter stop than the package uses.
  y <- sim_zip(300, 5, 0.6, 0.1, 0.25, 7)
  y[sample(length(y), 60)] <- NA
  n <- nrow(y); Tn <- ncol(y)
  D <- mixtureEM:::.lcga_design(0:4, 1)
  mm <- lcga_model(2, D, "zip", inflation_design = D)
  mm <- init_params(mm, y, NULL, random_state = 3)
  set.seed(4); R <- matrix(runif(n * 2), n, 2); R <- R / rowSums(R)
  got <- m_step(mm, y, R)$parameters

  pr <- mixtureEM:::.lcga_prior_rows(mm, mixtureEM:::.lcga_marginal(y))
  v  <- as.vector(y); ok <- !is.na(v); DL <- D[rep(1:Tn, each = n), ]
  ctl <- list(epsilon = 1e-14, maxit = 100)
  for (k in 1:2) {
    w  <- rep(R[, k], Tn)
    pk <- plogis(as.vector(DL %*% mm$parameters$inflation[k, ]))
    mk <- exp(as.vector(DL %*% mm$parameters$coefs[k, ]))
    z  <- ifelse(ok & v == 0, pk / (pk + (1 - pk) * exp(-mk)), 0)
    gi <- suppressWarnings(glm.fit(DL[ok, ], z[ok], weights = w[ok],
                                   family = binomial(), control = ctl))$coefficients
    gc <- suppressWarnings(glm.fit(rbind(DL[ok, ], pr$D), c(v[ok], pr$y),
                                   weights = c((w * (1 - z))[ok], pr$w),
                                   family = poisson(), control = ctl))$coefficients
    expect_equal(got$inflation[k, ], unname(gi), tolerance = 1e-9)
    expect_equal(got$coefs[k, ], unname(gc), tolerance = 1e-9)
  }
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

test_that("a linear inflation logit is fitted, counted and reported", {
  sim_zip_t <- function(n, b, g, seed) {
    set.seed(seed)
    t  <- 0:4
    mu <- exp(b[1] + b[2] * t)
    p  <- plogis(g[1] + g[2] * t)
    y  <- matrix(rpois(n * 5, rep(mu, each = n)), n, 5)
    y[matrix(runif(n * 5) < rep(p, each = n), n, 5)] <- 0
    y
  }
  y <- rbind(sim_zip_t(400, c(1.3, -0.1), c(-2, 0.6), 7),
             sim_zip_t(400, c(0.6, 0.1), c(0.8, -0.5), 8))
  f0 <- fit_lcga(y, times = 5, n_classes = 2, family = "zip", n_init = 5,
                 random_state = 1, n_cores = 1, bayes_constants = off)
  f1 <- fit_lcga(y, times = 5, n_classes = 2, family = "zip",
                 inflation_degree = 1, n_init = 5, random_state = 1,
                 n_cores = 1, bayes_constants = off)
  # Degree 0 is the constant model it always was.
  expect_null(dim(f0$growth$inflation_probability))
  expect_equal(f1$metrics$n_params, 9)
  expect_gt(f1$metrics$ll, f0$metrics$ll)
  hi <- which.max(f1$growth$coefficients[, 1])
  expect_gt(f1$growth$inflation[hi, 2], 0)
  expect_lt(f1$growth$inflation[3 - hi, 2], 0)
  expect_equal(dim(f1$growth$inflation_probability), c(2L, 5L))
  ge <- growth_estimates(f1)
  gl <- ge[ge$parameter == "inflation_logit", ]
  expect_setequal(unique(gl$item), c("intercept", "linear"))
  expect_true(all(is.finite(gl$se)))
  expect_output(print(f1), "by occasion")
})

test_that("inflation_degree is refused where it cannot apply", {
  y <- sim_zip(100, 4, 0.5, 0.1, 0.3, 9)
  expect_error(fit_lcga(y, times = 4, family = "poisson", inflation_degree = 1),
               "zip")
  expect_error(fit_lcga(y, times = 4, family = "zip", inflation_degree = 2),
               "0 (constant) or 1", fixed = TRUE)
  expect_error(fit_lcga(y[, 1:2], times = 2, degree = 0, family = "zip",
                        inflation_degree = 1), "at least 3 occasions")
})
