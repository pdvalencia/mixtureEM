# ==============================================================================
# Covariates on the growth factors in growth mixture models
# ==============================================================================
#
# Checks that need no reference program: the conditional-mean algebra of the
# emission, recovery of known regressions from simulated data, M-step
# monotonicity, and wrapper behaviour (naming, alignment, FIML, sorting).

# ------------------------------------------------------------------------------
# Level 1: the likelihood
# ------------------------------------------------------------------------------

test_that("the conditional mean is Lambda (alpha + Gamma x), case by case", {
  set.seed(211)
  n <- 200
  L <- .lcga_design(0:3, 1)
  Y <- matrix(rnorm(n * 4, 2, 2), n, 4)
  X <- cbind(rnorm(n), runif(n))

  alpha <- rbind(c(0.5, 0.30), c(2.0, -0.20))
  gamma <- list(matrix(c(0.4, 0.1, -0.2, 0.3), 2, 2),
                matrix(c(0.7, -0.3, 0.5, 0.0), 2, 2))
  psi   <- matrix(c(0.8, 0.1, 0.1, 0.3), 2, 2)
  theta <- rbind(c(0.5, 0.6, 0.7, 0.8), c(0.5, 0.6, 0.7, 0.8))

  mm <- structured_normal_model(2, L, growth_covariates = X,
                                growth_covariates_equal = FALSE)
  mm$parameters$alpha <- alpha
  mm$parameters$gamma <- gamma
  mm$parameters$psi   <- list(psi, psi)
  mm$parameters$theta <- theta

  # Written out one case at a time on purpose: the covariance is untouched by
  # the covariates -- they move the mean and nothing else -- so any disagreement
  # localises to the mean rather than to the structured Sigma.
  expected <- matrix(0, n, 2)
  for (k in 1:2) {
    Sig <- L %*% psi %*% t(L) + diag(theta[k, ])
    ch  <- chol(Sig)
    for (i in seq_len(n)) {
      mu <- as.vector(L %*% (alpha[k, ] + gamma[[k]] %*% X[i, ]))
      z  <- backsolve(ch, Y[i, ] - mu, transpose = TRUE)
      expected[i, k] <- -0.5 * (4 * log(2 * pi) + 2 * sum(log(diag(ch))) +
                                  sum(z^2))
    }
  }

  expect_equal(log_likelihood(mm, Y), expected, tolerance = 1e-10)
})

test_that("a covariate with a zero coefficient leaves the model where it was", {
  set.seed(213)
  n <- 150
  L <- .lcga_design(0:3, 1)
  Y <- matrix(rnorm(n * 4, 1, 1.5), n, 4)

  plant <- function(mm) {
    mm$parameters$alpha <- rbind(c(0.5, 0.3), c(2.0, -0.2))
    psi <- matrix(c(0.7, 0.05, 0.05, 0.2), 2, 2)
    mm$parameters$psi   <- list(psi, psi)
    mm$parameters$theta <- rbind(c(0.4, 0.5, 0.6, 0.7), c(0.4, 0.5, 0.6, 0.7))
    mm
  }

  bare <- plant(structured_normal_model(2, L))
  with_x <- plant(structured_normal_model(2, L,
                                          growth_covariates = matrix(rnorm(n))))
  with_x$parameters$gamma <- rep(list(matrix(0, 2, 1)), 2)

  # The unconditional model is the conditional one at Gamma = 0, and the two
  # code paths through .sn_mu() are different, so this pins the second to the
  # first rather than merely to itself.
  expect_equal(log_likelihood(with_x, Y), log_likelihood(bare, Y),
               tolerance = 1e-12)
})

test_that("estimation recovers class-specific growth-factor regressions", {
  set.seed(2027)
  n <- 900; Tn <- 5
  cls <- rbinom(n, 1, 0.5)
  x   <- rnorm(n)
  # The covariate raises the intercept in one class and lowers it in the other,
  # which no single shared coefficient can represent.
  icept <- ifelse(cls == 1, 4, 1) + ifelse(cls == 1, 0.8, -0.8) * x +
    rnorm(n, 0, 0.8)
  slope <- ifelse(cls == 1, 1.0, 0.1) + 0.25 * x + rnorm(n, 0, 0.3)
  Y <- outer(icept, rep(1, Tn)) + outer(slope, 0:(Tn - 1)) +
    matrix(rnorm(n * Tn, 0, 0.6), n, Tn)

  fit <- fit_gmm(Y, times = Tn, n_classes = 2, growth_predictors = x,
                 growth_predictors_equal = FALSE, n_init = 10,
                 random_state = 5)

  up <- which.max(fit$growth$means[, 2])
  expect_equal(unname(fit$growth$means[up, ]), c(4, 1.0), tolerance = 0.15)
  expect_equal(unname(fit$growth$means[-up, ]), c(1, 0.1), tolerance = 0.15)
  expect_equal(unname(as.vector(fit$growth$coefficients[[up]])), c(0.8, 0.25),
               tolerance = 0.12)
  expect_equal(unname(as.vector(fit$growth$coefficients[[-up]])), c(-0.8, 0.25),
               tolerance = 0.12)
})

test_that("a covariate the data do not need estimates near zero", {
  set.seed(2029)
  n <- 600; Tn <- 4
  cls   <- rbinom(n, 1, 0.5)
  icept <- ifelse(cls == 1, 3, 0) + rnorm(n, 0, 1)
  slope <- ifelse(cls == 1, 0.8, -0.2) + rnorm(n, 0, 0.35)
  Y <- outer(icept, rep(1, Tn)) + outer(slope, 0:3) +
    matrix(rnorm(n * Tn, 0, 0.7), n, Tn)
  noise <- rnorm(n)

  fit <- fit_gmm(Y, times = 4, n_classes = 2, growth_predictors = noise,
                 n_init = 5, random_state = 7)
  expect_lt(max(abs(fit$growth$coefficients[[1L]])), 0.12)
})

# ------------------------------------------------------------------------------
# The M-step
# ------------------------------------------------------------------------------

test_that("EM stays monotone with covariates on the growth factors", {
  # The property that says the joint solve is the exact conditional maximiser
  # rather than an approximation to it. It is worth checking directly here for
  # the same reason it was worth checking when the emission was written: the
  # M-step is an ECM whose first stage now solves a *constrained* system, and a
  # constraint imposed by averaging unconstrained solutions -- the obvious
  # shortcut -- would still look convergent while quietly decreasing the
  # likelihood on some steps.
  set.seed(2047)
  n <- 400; Tn <- 5
  cls <- rbinom(n, 1, 0.45)
  x   <- cbind(rnorm(n), rbinom(n, 1, 0.4))
  ic  <- ifelse(cls == 1, 3, 0) + x %*% c(0.6, -0.4) + rnorm(n, 0, 0.9)
  sl  <- ifelse(cls == 1, 0.9, -0.1) + x %*% c(0.25, 0.1) + rnorm(n, 0, 0.3)
  Y   <- outer(as.vector(ic), rep(1, Tn)) +
    outer(as.vector(sl), seq.int(0L, Tn - 1L)) +
    matrix(rnorm(n * Tn, 0, 0.7), n, Tn)
  Y_miss <- Y
  Y_miss[sample(n, 100), Tn] <- NA

  L <- .lcga_design(seq.int(0L, Tn - 1L), 1)

  run <- function(data, ...) {
    mm <- init_params(
      structured_normal_model(2, L, growth_covariates = x, ...),
      data, NULL, random_state = 4)
    w  <- c(0.5, 0.5)
    ll <- numeric(200L)
    for (it in seq_len(200L)) {
      lp  <- sweep(log_likelihood(mm, data), 2, log(w), "+")
      mx  <- apply(lp, 1, max)
      lse <- mx + log(rowSums(exp(lp - mx)))
      ll[it] <- sum(lse)
      resp <- exp(lp - lse)
      w    <- colMeans(resp)
      mm   <- m_step(mm, data, resp)
    }
    ll
  }

  specs <- list(
    "gamma equal"        = list(),
    "gamma free"         = list(growth_covariates_equal = FALSE),
    "psi and theta free" = list(psi = "free", residual_equal = FALSE),
    "random intercept"   = list(random_effects = "intercept"),
    "no random effects"  = list(random_effects = "none"))

  for (nm in names(specs)) {
    ll <- do.call(run, c(list(Y), specs[[nm]]))
    # A tolerance of one part in 1e10 of the log-likelihood, which is far above
    # the ~1e-13 of floating-point noise an exactly monotone run produces and
    # far below any real decrease.
    expect_gt(min(diff(ll)), -abs(ll[length(ll)]) * 1e-10, label = nm)
  }

  # And with an occasion unobserved for a fifth of the sample, since the
  # covariate design enters the normal equations once per missingness pattern.
  ll <- run(Y_miss)
  expect_gt(min(diff(ll)), -abs(ll[length(ll)]) * 1e-10)
})

# ------------------------------------------------------------------------------
# Input handling
# ------------------------------------------------------------------------------

.sim_cond <- function(n = 400, Tn = 4) {
  cls   <- rbinom(n, 1, 0.5)
  x     <- rnorm(n)
  icept <- ifelse(cls == 1, 3, 0) + 0.5 * x + rnorm(n, 0, 0.9)
  slope <- ifelse(cls == 1, 0.8, -0.2) + 0.3 * x + rnorm(n, 0, 0.3)
  list(Y = outer(icept, rep(1, Tn)) + outer(slope, seq.int(0L, Tn - 1L)) +
         matrix(rnorm(n * Tn, 0, 0.7), n, Tn),
       x = x)
}

test_that("growth predictors are named, dummy-coded and length-checked", {
  set.seed(2031)
  d <- .sim_cond(n = 300)
  grp <- factor(sample(c("low", "high"), 300, TRUE), levels = c("low", "high"))

  fit <- fit_gmm(d$Y, times = 4, n_classes = 2,
                 growth_predictors = data.frame(age = d$x, grp = grp),
                 n_init = 3, random_state = 1)
  # A two-level factor becomes one dummy against the first level, as everywhere
  # else in the package.
  expect_equal(colnames(fit$growth$coefficients[[1L]]), c("age", "grp.high"))
  expect_equal(rownames(fit$growth$coefficients[[1L]]), c("intercept", "linear"))

  # A bare vector keeps the name it was written with, so print() is readable.
  age <- d$x
  named <- fit_gmm(d$Y, times = 4, n_classes = 2, growth_predictors = age,
                   n_init = 2, random_state = 1)
  expect_equal(colnames(named$growth$coefficients[[1L]]), "age")

  expect_error(fit_gmm(d$Y, times = 4, n_classes = 2,
                       growth_predictors = d$x[1:10]),
               "one row per case")
})

test_that("a missing growth predictor is refused rather than silently dropped", {
  set.seed(2033)
  d <- .sim_cond(n = 200)
  x <- d$x; x[c(5, 17)] <- NA

  # Deleting these cases here would leave `weights`, `strata` and `cluster`
  # pointing at the wrong rows, so the package refuses and says why.
  expect_error(fit_gmm(d$Y, times = 4, n_classes = 2, growth_predictors = x),
               "must be complete")
  expect_error(fit_gmm(d$Y, times = 4, n_classes = 2, growth_predictors = x),
               "missing occasions are handled by FIML")
})

test_that("covariates stay aligned when an empty case is deleted", {
  set.seed(2035)
  d <- .sim_cond(n = 400)
  Y <- d$Y
  Y[c(3, 200), ] <- NA               # no observed occasion at all

  # fit_mixture_internal() deletes these rows before the emission is built, so
  # the covariate matrix has to lose the same rows or every case after the first
  # deletion carries someone else's covariate. Comparing against the fit on the
  # data with those rows already removed is what makes that visible: a
  # misalignment would move the estimates, not just the bookkeeping.
  fit <- suppressWarnings(
    fit_gmm(Y, times = 4, n_classes = 2, growth_predictors = d$x,
            n_init = 5, random_state = 9))
  ref <- fit_gmm(Y[-c(3, 200), ], times = 4, n_classes = 2,
                 growth_predictors = d$x[-c(3, 200)], n_init = 5,
                 random_state = 9)

  expect_equal(fit$missing_data$n_empty_rows, 2L)
  expect_equal(nrow(fit$data), 398L)
  expect_equal(fit$metrics$ll, ref$metrics$ll, tolerance = 1e-8)
  # ignore_attr because the two calls derive different display names for the
  # covariate (`d$x` against a subsetted expression); the numbers are the point.
  expect_equal(fit$growth$coefficients[[1L]], ref$growth$coefficients[[1L]],
               tolerance = 1e-8, ignore_attr = TRUE)
})

test_that("FIML on the occasions still works with covariates present", {
  set.seed(2037)
  d <- .sim_cond(n = 700, Tn = 5)
  Y <- d$Y
  Y[sample(700, 120), 4] <- NA
  Y[sample(700, 90), 5]  <- NA

  fit <- fit_gmm(Y, times = 5, n_classes = 2, growth_predictors = d$x,
                 n_init = 5, random_state = 3)
  expect_equal(nrow(fit$data), 700L)
  expect_equal(fit$missing_data$n_missing, 210L)
  expect_equal(unname(as.vector(fit$growth$coefficients[[1L]])), c(0.5, 0.3),
               tolerance = 0.12)
})

# ------------------------------------------------------------------------------
# Downstream: sorting, printing, generating
# ------------------------------------------------------------------------------

test_that("class sorting permutes the regressions with everything else", {
  set.seed(2039)
  n <- 700; Tn <- 4
  cls   <- rbinom(n, 1, 0.25)              # deliberately unequal
  x     <- rnorm(n)
  icept <- ifelse(cls == 1, 4, 0) + ifelse(cls == 1, 1.2, -0.4) * x +
    rnorm(n, 0, 0.7)
  slope <- ifelse(cls == 1, 1.0, -0.1) + rnorm(n, 0, 0.3)
  Y <- outer(icept, rep(1, Tn)) + outer(slope, 0:3) +
    matrix(rnorm(n * Tn, 0, 0.6), n, Tn)

  fit <- fit_gmm(Y, times = 4, n_classes = 2, growth_predictors = x,
                 growth_predictors_equal = FALSE, n_init = 10,
                 random_state = 11)

  expect_true(fit$weights[1] >= fit$weights[2])
  # The larger class is the one at zero. If Gamma were left unpermuted its rows
  # would be paired with the other class's intercepts, and the small class would
  # be reported with the large class's covariate effect.
  small <- which.max(fit$growth$means[, 1])
  expect_gt(fit$growth$coefficients[[small]][1, 1], 0.5)
  expect_lt(fit$growth$coefficients[[-small]][1, 1], 0)
})

test_that("print reports intercepts and regressions rather than means", {
  set.seed(2041)
  d <- .sim_cond(n = 300)
  fit <- fit_gmm(d$Y, times = 4, n_classes = 2, growth_predictors = d$x,
                 n_init = 3, random_state = 1)

  out <- paste(capture.output(print(fit)), collapse = "\n")
  # With covariates alpha is an intercept, not a mean, and saying "mean" would
  # be wrong rather than merely terse.
  expect_match(out, "GROWTH FACTOR INTERCEPTS")
  expect_false(grepl("GROWTH FACTOR MEANS", out))
  expect_match(out, "GROWTH FACTORS ON COVARIATES \\(held equal across classes\\)")
  expect_match(out, "at the mean of the covariates")

  bare <- fit_gmm(d$Y, times = 4, n_classes = 2, n_init = 3, random_state = 1)
  expect_match(paste(capture.output(print(bare)), collapse = "\n"),
               "GROWTH FACTOR MEANS")
})

test_that("the reported trajectory is the curve at the covariate mean", {
  set.seed(2043)
  d <- .sim_cond(n = 400)
  fit <- fit_gmm(d$Y, times = 4, n_classes = 2, growth_predictors = d$x,
                 n_init = 3, random_state = 1)

  xbar  <- mean(d$x)
  gamma <- fit$growth$coefficients[[1L]]
  by_hand <- t(fit$growth$design %*%
                 t(fit$growth$means + rep(as.vector(gamma %*% xbar),
                                          each = fit$n_components)))
  expect_equal(fit$growth$fitted, by_hand, tolerance = 1e-10,
               ignore_attr = TRUE)

  pdf(NULL); on.exit(dev.off())
  expect_no_error(plot(fit))
  expect_no_error(plot(fit, observed = "cases"))
})

test_that("the BLRT generator conditions on the observed covariates", {
  set.seed(2045)
  d <- .sim_cond(n = 500)
  fit <- fit_gmm(d$Y, times = 4, n_classes = 2, growth_predictors = d$x,
                 n_init = 5, random_state = 1)

  set.seed(13)
  N   <- nrow(d$Y)
  gen <- generate_synthetic_data(fit$mm, rep(1L, N), N)
  expect_equal(dim(gen), c(N, 4L))

  # Case i is drawn around its *own* mean, so the generated first occasion
  # tracks the covariate. Drawing every case around the class curve instead
  # would give a null distribution for the unconditional model, which is not the
  # one being tested.
  mu1 <- .sn_mu(fit$mm, 1L)
  expect_equal(colMeans(gen - mu1), rep(0, 4), tolerance = 0.12,
               ignore_attr = TRUE)
  expect_gt(cor(gen[, 1], d$x), 0.1)
  # And the residual covariance is still the structured Sigma.
  expect_equal(cov(gen - mu1), .sn_sigma(fit$mm, 1L), tolerance = 0.2,
               ignore_attr = TRUE)

  # A replicate of a different size cannot be generated conditionally, and
  # silently reusing the wrong covariates would be worse than refusing.
  expect_error(generate_synthetic_data(fit$mm, rep(1L, 10L), 10L),
               "one case per original case")
})
