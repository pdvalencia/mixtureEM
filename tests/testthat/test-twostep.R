# ==============================================================================
# Two-step estimation (n_steps = 2; Bakk & Kuha, 2018)
# ==============================================================================
#
# The two-step holds the measurement model at its step-1 estimate and
# maximises the full likelihood over the structural parameters, re-scoring
# every case under the joint model at each iteration. Three exact
# consequences are checked here, none of which needs an external target:
#
# (a) the freeze bites: the measurement parameters and the pooled class sizes
#     are bit-identical to the step-1 fit, for every measurement family and
#     with the priors on or off;
# (b) with near-perfect separation the one-step, two-step and ML three-step
#     agree, because a covariate can then move no class probability;
# (c) freeing the measurement block from the two-step solution and continuing
#     EM lands on the one-step maximum -- the freeze is the only difference.
#
# The remaining assertions check that the estimator is not the uncorrected
# third step it replaced, and that add_covariates(steps = 2) is the same fit.

.ts_sim <- function(measurement = "binary", n = 400, seed = 20260921, rho = 0.85) {
  set.seed(seed)
  z   <- rnorm(n)
  cls <- 1L + rbinom(n, 1, plogis(-0.4 + 0.9 * z))
  X <- switch(measurement,
    binary = matrix(rbinom(n * 6, 1, ifelse(rep(cls, 6) == 1L, rho, 1 - rho)),
                    n, 6),
    continuous = matrix(rnorm(n * 6, ifelse(rep(cls, 6) == 1L, 1, -1)), n, 6),
    categorical = matrix(vapply(rep(cls, 6), function(k)
      sample(1:3, 1, prob = if (k == 1L) c(.7, .2, .1) else c(.1, .2, .7)),
      integer(1)), n, 6))
  y <- rnorm(n, mean = ifelse(cls == 1L, 2, 0))
  list(X = X, Z = data.frame(z = z), y = y)
}

.ts_fit <- function(d, measurement, structural, n_steps, correction = "none",
                    bayes_constants = NULL, n_init = 5, ...) {
  args <- list(d$X, n_classes = 2, measurement = measurement,
               n_steps = n_steps, correction = correction, n_init = n_init,
               random_state = 1, bayes_constants = bayes_constants, ...)
  if (measurement == "continuous") args$variances_equal <- FALSE
  if (structural == "predictors") args$predictors <- d$Z else args$outcome <- d$y
  suppressMessages(do.call(fit_mixture, args))
}

.ts_anchor <- function(fit) {
  B <- fit$sm$parameters$beta
  sweep(B, 2, B[nrow(B), ], "-")[1, ]
}

# ------------------------------------------------------------------------------
# (a) the freeze
# ------------------------------------------------------------------------------

test_that("the two-step leaves the measurement model and class sizes at step 1", {
  # Same seed, same step-1 search: the uncorrected third step (n_steps = 3,
  # correction = "none") never touches the measurement block, so it is the
  # reference for what step 1 produced. A grid, not a fixture: every family
  # the engine fits (categorical is the one outside the L-BFGS polish, so it
  # exercises the EM-only path), with the priors on and off, under a class
  # predictor and under a distal outcome.
  off <- list(latent = 0, categorical = 0, variances = 0)
  for (measurement in c("binary", "continuous", "categorical")) {
    d <- .ts_sim(measurement)
    for (priors in list(NULL, off)) for (structural in c("predictors", "outcome")) {
      fit2 <- .ts_fit(d, measurement, structural, n_steps = 2,
                      bayes_constants = priors)
      fit3 <- .ts_fit(d, measurement, structural, n_steps = 3,
                      bayes_constants = priors)
      info <- sprintf("%s / %s / priors %s", measurement, structural,
                      if (is.null(priors)) "on" else "off")
      expect_true(fit2$converged, info = info)
      expect_identical(fit2$weights, fit3$weights, info = info)
      for (nm in names(fit3$mm$parameters))
        expect_identical(fit2$mm$parameters[[nm]], fit3$mm$parameters[[nm]],
                         info = paste(info, nm))
      expect_null(fit2$frozen, info = info)
      expect_equal(fit2$n_steps, 2L, info = info)
    }
  }
})

# ------------------------------------------------------------------------------
# It is not the uncorrected third step
# ------------------------------------------------------------------------------

test_that("the two-step climbs the joint likelihood the naive third step does not", {
  d   <- .ts_sim("binary")
  off <- list(latent = 0, categorical = 0)
  for (priors in list(NULL, off)) {
    fit2 <- .ts_fit(d, "binary", "predictors", n_steps = 2, bayes_constants = priors)
    fit3 <- .ts_fit(d, "binary", "predictors", n_steps = 3, bayes_constants = priors)
    # Different estimators, different coefficients.
    expect_gt(max(abs(.ts_anchor(fit2) - .ts_anchor(fit3))), 1e-3)
    # metrics$ll is the joint log-likelihood at the two-step point; the
    # naive coefficients score lower on the same objective. With the priors
    # off that objective is exactly what EM maximised, so the inequality is
    # a property of the estimator and not of the penalty.
    if (!is.null(priors)) {
      ll_naive <- sum(fit3$sample_weights *
                        e_step(fit3, fit3$data, fit3$Y)$log_prob_norm)
      expect_gt(fit2$metrics$ll, ll_naive)
      ll_self <- sum(fit2$sample_weights *
                       e_step(fit2, fit2$data, fit2$Y)$log_prob_norm)
      expect_equal(fit2$metrics$ll, ll_self, tolerance = 1e-10)
    }
  }
  # A two-step fit reports a step 1 of its own.
  expect_false(is.null(fit2$step1_metrics))
})

# ------------------------------------------------------------------------------
# (b) near-perfect separation
# ------------------------------------------------------------------------------

test_that("with entropy near 1 the one-, two- and three-step agree", {
  d    <- .ts_sim("binary", rho = 0.98)
  fit1 <- .ts_fit(d, "binary", "predictors", n_steps = 1)
  fit2 <- .ts_fit(d, "binary", "predictors", n_steps = 2)
  fit3 <- .ts_fit(d, "binary", "predictors", n_steps = 3, correction = "ML")
  expect_gt(fit2$metrics$entropy, 0.95)
  b1 <- .ts_anchor(fit1); b2 <- .ts_anchor(fit2); b3 <- .ts_anchor(fit3)
  # Classes may be labelled the other way round in a different fit: align on
  # the sign of the slope.
  if (sign(b1["z"]) != sign(b2["z"])) b2 <- -b2
  if (sign(b1["z"]) != sign(b3["z"])) b3 <- -b3
  expect_lt(max(abs(b2 - b1)), 0.1)
  expect_lt(max(abs(b3 - b1)), 0.1)
})

# ------------------------------------------------------------------------------
# (c) the freeze is the only difference from one-step
# ------------------------------------------------------------------------------

test_that("freeing the measurement block from the two-step solution recovers one-step", {
  d    <- .ts_sim("binary")
  fit1 <- .ts_fit(d, "binary", "predictors", n_steps = 1, n_init = 10)
  fit2 <- .ts_fit(d, "binary", "predictors", n_steps = 2)
  # Resume the engine's own loop from the two-step state with nothing frozen:
  # the same EM the one-step fit runs, started where the two-step stopped.
  st <- fit2
  st$frozen <- NULL
  st <- fit_single_init(st, fit2$data, fit2$Y, max_iter = 2000,
                        refine = FALSE, init_state = st)
  ll_freed <- sum(st$sample_weights * st$lower_bound)
  expect_true(st$converged)
  expect_lt(abs(ll_freed - fit1$metrics$ll), 1e-2)
})

# ------------------------------------------------------------------------------
# The add_*() interface
# ------------------------------------------------------------------------------

test_that("add_covariates(steps = 2) is the two-step on the fitted model", {
  d    <- .ts_sim("binary")
  fit2 <- .ts_fit(d, "binary", "predictors", n_steps = 2)
  fit0 <- suppressMessages(fit_mixture(d$X, n_classes = 2, measurement = "binary",
                                       n_init = 5, random_state = 1))
  fita <- suppressMessages(add_covariates(fit0, d$Z, steps = 2))
  expect_equal(fita$n_steps, 2L)
  expect_equal(fita$correction, "none")
  expect_equal(.ts_anchor(fita), .ts_anchor(fit2), tolerance = 1e-6)
  expect_equal(fita$metrics$ll, fit2$metrics$ll, tolerance = 1e-8)
  expect_false(is.null(fita$Y))

  expect_error(add_covariates(fit0, d$Z, steps = 2, correction = "ML"),
               "three-step only")
  expect_error(add_outcome(fit0, d$y, steps = 2, correction = "BCH"),
               "three-step only")
  expect_error(add_covariates(fit0, d$Z, steps = 4), "must be 3")

  # No step-1 term yet: the printed estimator says which matrix it is.
  out <- capture.output(summary(fita))
  expect_true(any(grepl("Q-function Hessian", out, fixed = TRUE)))

  fitb <- suppressMessages(add_outcome(fit0, d$y, steps = 2))
  expect_equal(fitb$n_steps, 2L)
  expect_identical(fitb$weights, fit0$weights)
  expect_identical(fitb$mm$parameters$pis, fit0$mm$parameters$pis)
})
