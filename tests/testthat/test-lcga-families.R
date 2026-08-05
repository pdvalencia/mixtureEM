# ==============================================================================
# LCGA families beyond the binomial: gaussian and Poisson
# ==============================================================================
#
# Checks that need no reference program: the likelihood against closed-form
# densities, corroboration against flexmix on simulated data, and wrapper
# behaviour for the two families.

test_that("the gaussian LCGA likelihood is the normal density along the curve", {
  set.seed(31)
  n <- 200; Tn <- 4; ts <- 0:(Tn - 1)
  Y <- matrix(rnorm(n * Tn, 2, 1.5), n, Tn)

  coefs <- rbind(c(0.5, 0.30), c(2.0, -0.20))
  disp  <- c(0.8, 1.5)

  mm <- lcga_model(2, design = .lcga_design(ts, 1), family = "gaussian")
  mm$parameters$coefs      <- coefs
  mm$parameters$dispersion <- disp

  # Written out the long way on purpose: a class contributes the product over
  # occasions of a normal density centred on its own straight line, with a
  # variance that does not change over occasions.
  expected <- matrix(0, n, 2)
  for (k in 1:2) {
    mu <- coefs[k, 1] + coefs[k, 2] * ts
    for (t in seq_len(Tn))
      expected[, k] <- expected[, k] +
        dnorm(Y[, t], mean = mu[t], sd = sqrt(disp[k]), log = TRUE)
  }

  expect_equal(log_likelihood(mm, Y), expected, tolerance = 1e-12)
})

test_that("the Poisson LCGA likelihood is the Poisson density along the curve", {
  set.seed(37)
  n <- 200; Tn <- 4; ts <- 0:(Tn - 1)
  Y <- matrix(rpois(n * Tn, 2), n, Tn)

  coefs <- rbind(c(0.4, 0.25), c(1.1, -0.10))

  mm <- lcga_model(2, design = .lcga_design(ts, 1), family = "poisson")
  mm$parameters$coefs      <- coefs
  mm$parameters$dispersion <- c(1, 1)

  expected <- matrix(0, n, 2)
  for (k in 1:2) {
    lambda <- exp(coefs[k, 1] + coefs[k, 2] * ts)   # log link
    for (t in seq_len(Tn))
      expected[, k] <- expected[, k] + dpois(Y[, t], lambda = lambda[t], log = TRUE)
  }

  expect_equal(log_likelihood(mm, Y), expected, tolerance = 1e-12)
})

test_that("an unobserved occasion drops out of the likelihood, for every family", {
  set.seed(41)
  n <- 150; Tn <- 4; ts <- 0:(Tn - 1)

  check_fiml <- function(family, Y) {
    Y_miss <- Y
    Y_miss[1:20, 3] <- NA

    mm <- lcga_model(2, design = .lcga_design(ts, 1), family = family)
    mm$parameters$coefs      <- rbind(c(0.3, 0.2), c(0.9, -0.1))
    mm$parameters$dispersion <- c(1.2, 0.7)

    full <- log_likelihood(mm, Y_miss)
    # The same cases with occasion 3 removed from the model entirely must give
    # the same contribution: that is what "the missing cell informs nothing"
    # means, as against imputing it or deleting the case.
    mm3 <- lcga_model(2, design = .lcga_design(ts[-3], 1), family = family)
    mm3$parameters$coefs      <- mm$parameters$coefs
    mm3$parameters$dispersion <- mm$parameters$dispersion
    drop <- log_likelihood(mm3, Y_miss[1:20, -3, drop = FALSE])

    expect_equal(full[1:20, ], drop, tolerance = 1e-12,
                 label = paste("FIML contribution,", family))
  }

  check_fiml("gaussian", matrix(rnorm(n * Tn), n, Tn))
  check_fiml("poisson",  matrix(rpois(n * Tn, 2), n, Tn))
  check_fiml("binomial", matrix(rbinom(n * Tn, 1, 0.5), n, Tn))
})

# ------------------------------------------------------------------------------
# Level 3: corroboration against flexmix, on simulated data
#
# flexmix's `y ~ t | id` is the same model: the `| id` grouping makes the
# posterior a case-level quantity, so each class is a single trajectory fitted
# to whole cases rather than to individual observations. Our estimates maximise
# a *penalised* likelihood — the emission carries alpha/K pseudo-observations
# per occasion at the pooled trajectory — so ours must sit just below flexmix's
# unpenalised maximum, never above it, and the coefficients must agree to within
# what a penalty that small can move them.
#
# These exercise designs the closed-form likelihood checks above do not cover —
# five occasions rather than four, and a binary family on simulated data —
# with an independent open-source implementation as the referee.
# ------------------------------------------------------------------------------

# Returns coefficients as a K x p matrix in our orientation, plus the
# log-likelihood and class proportions.
.fx_lcga <- function(Y, family, k = 2, nrep = 10) {
  n <- nrow(Y); Tn <- ncol(Y)
  d <- data.frame(t = rep(0:(Tn - 1), each = n), id = rep(seq_len(n), Tn))
  # flexmix's binomial branch wants successes and failures, not a 0/1 vector.
  d$y <- if (family == "binomial") cbind(as.vector(Y), 1 - as.vector(Y)) else
    as.vector(Y)

  m <- flexmix::stepFlexmix(
    y ~ t | id, data = d, k = k, nrep = nrep,
    model = flexmix::FLXMRglm(family = family),
    # minprior = 0 stops flexmix from discarding a small class, which would
    # change the model being compared rather than the fit of it.
    control = list(minprior = 0, tolerance = 1e-12, iter.max = 5000),
    verbose = FALSE)

  par <- flexmix::parameters(m)
  # Read off the slot rather than through logLik(): flexmix's method is S4, and
  # dispatch on stats::logLik() only happens once flexmix is attached, which a
  # test using `flexmix::` alone never does.
  list(coefs = t(par[c("coef.(Intercept)", "coef.t"), , drop = FALSE]),
       ll     = as.numeric(m@logLik),
       sigma  = if ("sigma" %in% rownames(par)) par["sigma", ] else NULL,
       prior  = flexmix::prior(m))
}

# Classes are labelled arbitrarily by both packages, so compare them ordered by
# intercept rather than by position.
.by_intercept <- function(m) unname(m[order(m[, 1]), , drop = FALSE])

expect_matches_flexmix <- function(fit, ref, tag, coef_tol = 0.03,
                                   ll_slack = 0.05) {
  expect_equal(.by_intercept(fit$growth$coefficients), .by_intercept(ref$coefs),
               tolerance = coef_tol, label = paste("coefficients,", tag))
  # Penalised, so strictly below their maximum...
  expect_lt(fit$metrics$ll, ref$ll)
  # ...but only just: the penalty is one pseudo-observation split across the
  # classes at each occasion, against n x T real ones.
  expect_lt(ref$ll - fit$metrics$ll, ll_slack)
}

test_that("gaussian LCGA agrees with flexmix", {
  skip_on_cran()
  skip_if_not_installed("flexmix")

  set.seed(2024)
  n <- 600; Tn <- 5; ts <- 0:(Tn - 1)
  cls <- rbinom(n, 1, 0.45)
  mu  <- outer(ifelse(cls == 1, 2.0, 0.0), rep(1, Tn)) +
         outer(ifelse(cls == 1, 0.5, -0.3), ts)
  Y   <- mu + matrix(rnorm(n * Tn, 0, 1), n, Tn)

  fit <- fit_lcga(Y, times = Tn, n_classes = 2, family = "gaussian",
                  n_init = 20, random_state = 5)
  ref <- .fx_lcga(Y, "gaussian")

  expect_matches_flexmix(fit, ref, "gaussian")

  # The residual variance is ours to report, so it is checked too: flexmix
  # parameterises the same quantity as a standard deviation.
  ord <- order(fit$growth$coefficients[, 1])
  expect_equal(fit$growth$residual_variance[ord],
               (ref$sigma[order(ref$coefs[, 1])])^2,
               tolerance = 0.02, ignore_attr = TRUE)
})

test_that("Poisson LCGA agrees with flexmix", {
  skip_on_cran()
  skip_if_not_installed("flexmix")

  set.seed(2024)
  n <- 600; Tn <- 5; ts <- 0:(Tn - 1)
  cls <- rbinom(n, 1, 0.45)
  eta <- outer(ifelse(cls == 1, 1.2, 0.2), rep(1, Tn)) +
         outer(ifelse(cls == 1, 0.2, -0.15), ts)
  Y   <- matrix(rpois(n * Tn, exp(eta)), n, Tn)

  fit <- fit_lcga(Y, times = Tn, n_classes = 2, family = "poisson",
                  n_init = 20, random_state = 5)
  expect_matches_flexmix(fit, .fx_lcga(Y, "poisson"), "poisson")
})

test_that("binomial LCGA agrees with flexmix too", {
  skip_on_cran()
  skip_if_not_installed("flexmix")

  # Running the binary family through the same external check keeps the three
  # families on one footing and would catch a change to the shared emission
  # that spared the other two.
  set.seed(2024)
  n <- 600; Tn <- 5; ts <- 0:(Tn - 1)
  cls <- rbinom(n, 1, 0.45)
  eta <- outer(ifelse(cls == 1, 0.0, -0.8), rep(1, Tn)) +
         outer(ifelse(cls == 1, 1.1, -0.1), ts)
  Y   <- matrix(rbinom(n * Tn, 1, plogis(eta)), n, Tn)

  fit <- fit_lcga(Y, times = Tn, n_classes = 2, family = "binomial",
                  n_init = 20, random_state = 5)
  # A logit slope is the least well determined quantity here, so the tolerance
  # is wider than for the other two families: ours is 1.222 against flexmix's
  # 1.241, and the gap is in the direction the prior pulls.
  expect_matches_flexmix(fit, .fx_lcga(Y, "binomial"), "binomial",
                         coef_tol = 0.05)
})


# ------------------------------------------------------------------------------
# The rest of the package, for the newly available families
# ------------------------------------------------------------------------------

test_that("the gaussian family counts its residual variances as parameters", {
  set.seed(51)
  Y <- matrix(rnorm(300 * 5), 300, 5)

  g <- fit_lcga(Y, times = 5, n_classes = 3, family = "gaussian",
                n_init = 3, random_state = 1)
  # 3 x 2 coefficients + 3 residual variances + 2 free proportions.
  expect_equal(g$metrics$n_params, 11)

  p <- fit_lcga(matrix(rpois(300 * 5, 2), 300, 5), times = 5, n_classes = 3,
                family = "poisson", n_init = 3, random_state = 1)
  # A Poisson has no free dispersion: 3 x 2 coefficients + 2 proportions.
  expect_equal(p$metrics$n_params, 8)
})

test_that("time_scores rescale the slope for the new families as well", {
  set.seed(53)
  n <- 500; Tn <- 4
  cls <- rbinom(n, 1, 0.5)
  Y <- outer(ifelse(cls == 1, 2, 0), rep(1, Tn)) +
       outer(ifelse(cls == 1, 0.6, -0.2), 0:3) +
       matrix(rnorm(n * Tn), n, Tn)

  a <- fit_lcga(Y, times = 4, n_classes = 2, family = "gaussian",
                n_init = 5, random_state = 3)
  b <- fit_lcga(Y, times = 4, n_classes = 2, family = "gaussian",
                n_init = 5, random_state = 3, time_scores = c(0, 2, 4, 6))

  expect_equal(a$metrics$ll, b$metrics$ll, tolerance = 1e-6)
  expect_equal(a$growth$coefficients[, 1], b$growth$coefficients[, 1],
               tolerance = 1e-4)
  expect_equal(a$growth$coefficients[, 2], 2 * b$growth$coefficients[, 2],
               tolerance = 1e-4)
  # Rescaling time is a reparameterisation, so the residual variance is untouched.
  expect_equal(a$growth$residual_variance, b$growth$residual_variance,
               tolerance = 1e-5)
})

test_that("FIML keeps cases with unobserved occasions, for the new families", {
  set.seed(59)
  n <- 600; Tn <- 4
  cls <- rbinom(n, 1, 0.5)
  Y <- outer(ifelse(cls == 1, 2, 0), rep(1, Tn)) +
       outer(ifelse(cls == 1, 0.6, -0.2), 0:3) +
       matrix(rnorm(n * Tn), n, Tn)
  Y[sample(seq_len(n), 120), 3] <- NA

  fit <- fit_lcga(Y, times = 4, n_classes = 2, family = "gaussian",
                  n_init = 5, random_state = 4)

  expect_equal(nrow(fit$data), n)
  expect_equal(fit$missing_data$n_missing, 120L)
  expect_true(all(is.finite(fit$growth$coefficients)))
  expect_true(all(fit$growth$residual_variance > 0))

  up <- which.max(fit$growth$coefficients[, 2])
  expect_equal(fit$growth$coefficients[up, ], c(2, 0.6), tolerance = 0.2)
})

test_that("print, plot and the BLRT generator handle the new families", {
  set.seed(61)
  n <- 300; Tn <- 4
  cls <- rbinom(n, 1, 0.5)
  Yg <- outer(ifelse(cls == 1, 2, 0), rep(1, Tn)) +
        outer(ifelse(cls == 1, 0.6, -0.2), 0:3) +
        matrix(rnorm(n * Tn), n, Tn)
  Yp <- matrix(rpois(n * Tn, exp(outer(ifelse(cls == 1, 1, 0.2), rep(1, Tn)))),
               n, Tn)

  g <- fit_lcga(Yg, times = 4, n_classes = 2, family = "gaussian",
                n_init = 3, random_state = 1)
  p <- fit_lcga(Yp, times = 4, n_classes = 2, family = "poisson",
                n_init = 3, random_state = 1)

  # The residual variance is printed for the gaussian family and for no other,
  # since the other two have no free dispersion and the emission carries a
  # placeholder 1 that would read as an estimate.
  out_g <- paste(capture.output(print(g)), collapse = "\n")
  out_p <- paste(capture.output(print(p)), collapse = "\n")
  expect_match(out_g, "RESIDUAL VARIANCE")
  expect_false(grepl("RESIDUAL VARIANCE", out_p))
  expect_match(out_g, "identity link")
  expect_match(out_p, "log link")

  pdf(NULL); on.exit(dev.off())
  expect_no_error(plot(g))
  expect_no_error(plot(p))

  # Class enumeration is the usual research question for a trajectory model, so
  # the BLRT's generator has to know how to draw from each family's curve.
  classes <- sample(1:2, n, TRUE)
  gen_g <- generate_synthetic_data(g$mm, classes, n)
  gen_p <- generate_synthetic_data(p$mm, classes, n)
  expect_equal(dim(gen_g), c(n, 4L))
  expect_equal(dim(gen_p), c(n, 4L))
  expect_true(all(gen_p >= 0) && all(gen_p == round(gen_p)))
  expect_false(all(gen_g == round(gen_g)))   # continuous, not counts
})

test_that("each family rejects outcomes it cannot model", {
  set.seed(67)
  cont  <- matrix(rnorm(200 * 4), 200, 4)
  count <- matrix(rpois(200 * 4, 2), 200, 4)

  expect_error(fit_lcga(cont, times = 4, family = "poisson"),
               "non-negative counts")
  expect_error(fit_lcga(count + 0.5, times = 4, family = "poisson"),
               "whole-number counts")
  expect_error(fit_lcga(cont, times = 4, family = "binomial"),
               "values in \\{0, 1\\}")

  inf <- cont; inf[1, 1] <- Inf
  expect_error(fit_lcga(inf, times = 4, family = "gaussian"), "finite outcome")

  # Counts and binary data are legitimate gaussian outcomes — badly modelled,
  # perhaps, but the package does not police that — so these must not error.
  expect_no_error(fit_lcga(count, times = 4, family = "gaussian",
                           n_init = 2, random_state = 1))
})
