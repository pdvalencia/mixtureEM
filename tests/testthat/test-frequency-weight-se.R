# A frequency weight of w is w identical cases, so a step-3 sandwich on the
# weighted rows must equal the unweighted one on the data expanded to one row
# per case: the meat is sum_i w_i s_i s_i', not the sampling-weight form
# sum_i w_i^2 s_i s_i'.

test_that("frequency-weighted step-3 sandwiches equal the row-expanded fit's", {
  set.seed(11)
  n   <- 150
  x   <- rnorm(n)
  cl  <- 1L + rbinom(n, 1, plogis(-0.85 + x))
  X   <- sapply(1:5, function(j) rbinom(n, 1, c(0.85, 0.15)[cl]))
  z   <- rbinom(n, 1, c(0.3, 0.7)[cl])
  y   <- rnorm(n, c(0, 1)[cl])
  cnt <- sample(1:3, n, TRUE)
  e   <- rep(seq_len(n), cnt)
  off <- list(categorical = 0, latent = 0, variances = 0, poisson = 0)
  q   <- function(e) suppressWarnings(suppressMessages(e))

  ses <- function(rows, w) {
    s1 <- fit_mixture(X[rows, ], n_classes = 2, measurement = "binary",
                      weights = w, weight_type = "frequency", n_init = 3,
                      random_state = 1, bayes_constants = off, n_cores = 1)
    oz <- q(add_outcome(s1, outcome = z[rows], outcome_type = "categorical",
                        correction = "ML", assignment = "modal", se = "robust"))
    oy <- q(add_outcome(s1, outcome = y[rows], outcome_type = "continuous",
                        variances = "class_specific", correction = "ML", assignment = "modal", se = "robust"))
    ox <- q(add_covariates(s1, predictors = cbind(x = x[rows]),
                           correction = "ML", assignment = "modal",
                           se = "robust"))
    m  <- oy$sm$parameters$means[, 1]
    o  <- order(m)
    list(means = m[o],
         z = sqrt(diag(pinv(-oz$sm$parameters$hessian)))[o],
         y = sqrt(diag(oy$sm$parameters$V_pmv))[3:6][c(o, 2 + o)],
         x = sqrt(diag(ox$sm$parameters$V_robust))[1:2])
  }
  f <- ses(seq_len(n), cnt)
  ex <- ses(e, NULL)
  expect_equal(f$means, ex$means, tolerance = 1e-6)
  expect_equal(f$z, ex$z, tolerance = 1e-5)
  expect_equal(f$y, ex$y, tolerance = 1e-5)
  expect_equal(f$x, ex$x, tolerance = 1e-5)
})

test_that("the two-step sandwich and the outer-product step-1 variance count cases", {
  set.seed(12)
  n   <- 150
  x   <- rnorm(n)
  cl  <- 1L + rbinom(n, 1, plogis(-0.85 + x))
  X   <- sapply(1:5, function(j) rbinom(n, 1, c(0.85, 0.15)[cl]))
  cnt <- sample(1:3, n, TRUE)
  e   <- rep(seq_len(n), cnt)
  off <- list(categorical = 0, latent = 0, variances = 0, poisson = 0)
  two <- function(rows, w) suppressWarnings(suppressMessages(fit_mixture(
    X[rows, ], n_classes = 2, measurement = "binary", predictors = cbind(x = x[rows]),
    n_steps = 2, se = "robust", weights = w, weight_type = "frequency",
    n_init = 3, random_state = 1, bayes_constants = off, n_cores = 1)))
  f  <- two(seq_len(n), cnt)
  ex <- two(e, NULL)
  expect_equal(sqrt(diag(f$sm$parameters$V_robust)),
               sqrt(diag(ex$sm$parameters$V_robust)), tolerance = 1e-5)

  # The outer-product step-1 variance, which the corrected SEs take for a
  # large measurement model, on the same two step-1 fits.
  d1 <- function(fit, w) {
    th <- .step1_pack(fit)
    .step1_variance(fit, fit$data, th, .step1_case_scores(fit, fit$data, th, w),
                    w, "outer")$V
  }
  expect_equal(d1(f, cnt), d1(ex, rep(1, length(e))), tolerance = 1e-5)
})

test_that("the LTA three-step BCH sandwich counts frequency-weighted cases", {
  set.seed(13)
  n  <- 300
  z  <- rnorm(n)
  s1 <- 1L + rbinom(n, 1, 0.5)
  s2 <- ifelse(runif(n) < plogis(1 + 0.8 * z), s1, 3L - s1)
  p  <- c(0.85, 0.15)
  X  <- cbind(sapply(1:5, function(j) rbinom(n, 1, p[s1])),
              sapply(1:5, function(j) rbinom(n, 1, p[s2])))
  cnt <- sample(1:3, n, TRUE)
  e   <- rep(seq_len(n), cnt)
  Z   <- data.frame(z = z)
  bch <- function(rows, w) suppressWarnings(suppressMessages(fit_lta(
    X[rows, ], n_statuses = 2, times = 2, measurement = "binary",
    predictors_transition = Z[rows, , drop = FALSE], n_steps = 3,
    correction = "BCH", weights = w, weight_type = "frequency",
    n_init = 3, random_state = 1, n_cores = 1)))
  f  <- bch(seq_len(n), cnt)
  ex <- bch(e, NULL)
  expect_equal(f$loglik, ex$loglik, tolerance = 1e-6)
  expect_equal(sqrt(diag(f$se$vcov)), sqrt(diag(ex$se$vcov)), tolerance = 1e-4)
})
