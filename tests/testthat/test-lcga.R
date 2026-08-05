# ==============================================================================
# Latent class growth analysis and the weighted-GLM M-step behind it
# ==============================================================================

# ------------------------------------------------------------------------------
# The weighted-GLM M-step
#
# The M-step's whole job is to reproduce a weighted GLM fit, so the reference is
# R's own glm(): if these agree, the M-step maximises the weighted likelihood it
# claims to. Non-integer weights are the case that matters, because a
# responsibility is never 0 or 1.
# ------------------------------------------------------------------------------

test_that("the weighted-GLM M-step reproduces glm() for all three families", {
  set.seed(42)
  n <- 300
  D <- cbind(1, rnorm(n), rnorm(n))
  w <- runif(n, 0.1, 3)

  expect_glm_agreement <- function(family, y) {
    ours <- .wglm_fit(D, y, w, .wglm_family(family))$coefficients
    ref  <- suppressWarnings(
      stats::glm(y ~ D[, 2] + D[, 3], weights = w, family = family))
    expect_equal(ours, unname(stats::coef(ref)), tolerance = 1e-6,
                 label = paste("coefficients,", family))
  }

  expect_glm_agreement("gaussian", as.vector(D %*% c(1, 2, -1)) + rnorm(n))
  expect_glm_agreement("binomial", rbinom(n, 1, plogis(D %*% c(0.2, 1, -0.5))))
  expect_glm_agreement("poisson",  rpois(n, exp(D %*% c(0.5, 0.4, -0.3))))
})

test_that("the gaussian M-step returns the weighted residual variance", {
  set.seed(1)
  n <- 500
  D <- cbind(1, rnorm(n))
  w <- runif(n, 0.5, 2)
  y <- as.vector(D %*% c(2, -1)) + rnorm(n, sd = 1.5)

  fit <- .wglm_fit(D, y, w, .wglm_family("gaussian"))
  mu  <- as.vector(D %*% fit$coefficients)
  expect_equal(fit$dispersion, sum(w * (y - mu)^2) / sum(w), tolerance = 1e-10)
})

test_that("cases with no weight are dropped, and an empty class keeps its start", {
  set.seed(3)
  n <- 200
  D <- cbind(1, rnorm(n))
  y <- rbinom(n, 1, 0.4)
  w <- c(rep(1, 100), rep(0, 100))

  # Zero weights must not merely shrink a case's influence but remove it, which
  # is what makes the same routine usable for FIML masking.
  zeroed <- .wglm_fit(D, y, w, .wglm_family("binomial"))$coefficients
  subset <- .wglm_fit(D[1:100, ], y[1:100], w[1:100],
                      .wglm_family("binomial"))$coefficients
  expect_equal(zeroed, subset, tolerance = 1e-8)

  start <- c(0.3, -0.2)
  empty <- .wglm_fit(D, y, rep(0, n), .wglm_family("binomial"), start = start)
  expect_equal(empty$coefficients, start)
  expect_false(empty$converged)
})

# ------------------------------------------------------------------------------
# The LCGA model
# ------------------------------------------------------------------------------

test_that("fit_lcga() recovers known trajectories", {
  set.seed(7)
  n   <- 800
  cls <- rbinom(n, 1, 0.5)
  eta <- outer(ifelse(cls == 1, 0.0, -0.8), rep(1, 4)) +
         outer(ifelse(cls == 1, 1.1, -0.1), 0:3)
  y   <- matrix(rbinom(n * 4, 1, plogis(eta)), n, 4)

  fit <- fit_lcga(y, times = 4, n_classes = 2, n_init = 10, random_state = 11)

  # Classes come back sorted by size, so identify them by shape rather than
  # position: the generating classes differ mainly in slope.
  co  <- fit$growth$coefficients
  up  <- which.max(co[, 2])
  flat <- setdiff(seq_len(2), up)

  # Tolerance is set from the estimator's own sampling variability rather than
  # from what one seed happens to produce: over 20 replications of this design
  # the coefficients have a standard deviation of 0.07 to 0.15 and no bias
  # larger than a couple of Monte Carlo standard errors, so 0.35 is roughly
  # 2.5 standard deviations and is not seed-dependent.
  expect_equal(co[up,   ], c(0.0, 1.1),  tolerance = 0.35)
  expect_equal(co[flat, ], c(-0.8, -0.1), tolerance = 0.35)
  expect_equal(sum(fit$weights), 1)
  expect_equal(unname(fit$weights[up]), 0.5, tolerance = 0.08)
})

test_that("a degree-d trajectory has d + 1 coefficients per class", {
  set.seed(5)
  y <- matrix(rbinom(400 * 5, 1, 0.5), 400, 5)

  lin  <- fit_lcga(y, times = 5, n_classes = 2, degree = 1, n_init = 2,
                   random_state = 1)
  quad <- fit_lcga(y, times = 5, n_classes = 2, degree = 2, n_init = 2,
                   random_state = 1)

  expect_equal(dim(lin$growth$coefficients),  c(2, 2))
  expect_equal(dim(quad$growth$coefficients), c(2, 3))
  # 2 classes x (d + 1) coefficients, plus the one free class proportion.
  expect_equal(lin$metrics$n_params,  5)
  expect_equal(quad$metrics$n_params, 7)
})

test_that("fit_lcga() rejects inputs it cannot model", {
  y <- matrix(rbinom(200 * 4, 1, 0.5), 200, 4)

  # A degree-3 curve through 4 occasions reproduces the occasion means exactly,
  # leaving nothing for the growth constraint to restrict.
  expect_error(fit_lcga(y, times = 4, degree = 3), "at least 5 occasions")
  expect_error(fit_lcga(y, times = 4, time_scores = c(0, 1, 2)),
               "one value per occasion")
  expect_error(fit_lcga(y, times = 2, n_classes = 2),
               "LCGA models one repeated outcome")
  expect_error(fit_lcga(matrix(rnorm(200 * 4), 200, 4), times = 4),
               "values in \\{0, 1\\}")
})

test_that("unobserved occasions are handled by FIML, not deletion", {
  set.seed(9)
  n   <- 600
  cls <- rbinom(n, 1, 0.5)
  eta <- outer(ifelse(cls == 1, 0.0, -0.9), rep(1, 4)) +
         outer(ifelse(cls == 1, 1.0, 0.0), 0:3)
  y   <- matrix(rbinom(n * 4, 1, plogis(eta)), n, 4)

  y_miss <- y
  y_miss[sample(seq_len(n), 120), 3] <- NA   # wave 3 attrition

  fit <- fit_lcga(y_miss, times = 4, n_classes = 2, n_init = 5,
                  random_state = 4)

  # Every case is still in the analysis: FIML drops the missing cell, not the row.
  expect_equal(nrow(fit$data), n)
  expect_equal(fit$missing_data$n_missing, 120L)
  expect_true(all(is.finite(fit$growth$coefficients)))

  co <- fit$growth$coefficients
  up <- which.max(co[, 2])
  expect_equal(co[up, ], c(0.0, 1.0), tolerance = 0.2)
})

test_that("time_scores set the meaning of the intercept and slope", {
  set.seed(13)
  n   <- 600
  cls <- rbinom(n, 1, 0.5)
  eta <- outer(ifelse(cls == 1, 0.0, -0.8), rep(1, 4)) +
         outer(ifelse(cls == 1, 1.0, 0.0), 0:3)
  y   <- matrix(rbinom(n * 4, 1, plogis(eta)), n, 4)

  a <- fit_lcga(y, times = 4, n_classes = 2, n_init = 5, random_state = 2)
  # Doubling the spacing halves every slope and leaves the fit itself untouched.
  b <- fit_lcga(y, times = 4, n_classes = 2, n_init = 5, random_state = 2,
                time_scores = c(0, 2, 4, 6))

  expect_equal(a$metrics$ll, b$metrics$ll, tolerance = 1e-6)
  expect_equal(a$growth$coefficients[, 1], b$growth$coefficients[, 1],
               tolerance = 1e-4)
  expect_equal(a$growth$coefficients[, 2], 2 * b$growth$coefficients[, 2],
               tolerance = 1e-4)
})

test_that("wide and long input describe the same model", {
  set.seed(17)
  n <- 500
  y <- matrix(rbinom(n * 4, 1, 0.45), n, 4)
  long <- data.frame(id = rep(seq_len(n), 4), t = rep(0:3, each = n),
                     u = as.vector(y))

  wide <- fit_lcga(y, times = 4, n_classes = 2, n_init = 5, random_state = 9)
  lng  <- fit_lcga(long, id = "id", time = "t", item = "u", n_classes = 2,
                   n_init = 5, random_state = 9)

  expect_equal(wide$metrics$ll, lng$metrics$ll, tolerance = 1e-8)
  expect_equal(wide$growth$coefficients, lng$growth$coefficients,
               tolerance = 1e-6)
})

test_that("a distal outcome rides the existing three-step machinery", {
  set.seed(23)
  n   <- 600
  cls <- rbinom(n, 1, 0.5)
  eta <- outer(ifelse(cls == 1, 0.0, -0.8), rep(1, 4)) +
         outer(ifelse(cls == 1, 1.1, -0.1), 0:3)
  y   <- matrix(rbinom(n * 4, 1, plogis(eta)), n, 4)
  distal <- 2 * cls + rnorm(n)

  # The first argument is `indicator`, not `outcome`, precisely so that this
  # call is possible: `outcome` has to stay free for the distal outcome that
  # fit_mixture() consumes. Naming the repeated variable `outcome` silently
  # rebound this argument and the fit failed inside the reshaping code.
  fit <- suppressMessages(
    fit_lcga(y, times = 4, n_classes = 2, outcome = distal,
             n_init = 5, random_state = 1))

  # Trajectory class -> distal outcome is the commonest applied design of all,
  # so the recovered class means matter as much as the trajectories.
  means <- sort(as.vector(fit$sm$parameters$means))
  expect_equal(means, c(0, 2), tolerance = 0.2)
})

test_that("an LCGA fit reaches the rest of the package", {
  set.seed(21)
  n   <- 400
  cls <- rbinom(n, 1, 0.5)
  eta <- outer(ifelse(cls == 1, 0.0, -0.8), rep(1, 4)) +
         outer(ifelse(cls == 1, 1.1, -0.1), 0:3)
  y   <- matrix(rbinom(n * 4, 1, plogis(eta)), n, 4)

  fit <- fit_lcga(y, times = 4, n_classes = 2, n_init = 5, random_state = 1)

  expect_no_error(capture.output(print(fit)))
  expect_no_error(capture.output(summary(fit)))
  expect_no_error(capture.output(measurement_summary(fit)))
  expect_no_error(capture.output(classification_diagnostics(fit)))

  pdf(NULL); on.exit(dev.off())
  expect_no_error(plot(fit))

  # Class enumeration is usually the research question for a trajectory model,
  # so the BLRT's data generator has to understand a growth emission: its class
  # parameters are coefficients, not one value per column, and the trajectory
  # must be evaluated before anything can be drawn.
  gen <- generate_synthetic_data(fit$mm, sample(seq_len(2), n, TRUE), n)
  expect_equal(dim(gen), c(n, 4L))
  expect_true(all(gen %in% c(0, 1)))
})

