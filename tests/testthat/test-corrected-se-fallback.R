# se = "corrected" on a measurement model with no step-1 packing reports the
# uncorrected estimator. That substitution is warned about at fit time, and an
# explicit request for the estimator actually reported is not.

.fallback_sim <- function(n = 300, seed = 1) {
  set.seed(seed)
  cls <- sample(1:2, n, TRUE)
  Y   <- sapply(0:4, function(t)
    rnorm(n, ifelse(cls == 1, 1 + 0.2 * t, 3 - 0.3 * t)))
  list(Y = Y, x = rnorm(n, mean = 0.7 * cls), cls = cls)
}

test_that("a growth step 1 warns under se = 'corrected', three-step and two-step", {
  s   <- .fallback_sim()
  fit <- suppressMessages(suppressWarnings(
    fit_gmm(s$Y, times = 5, n_classes = 2, n_init = 2, random_state = 1)))

  expect_warning(f3 <- suppressMessages(add_covariates(fit, s$x)),
                 "not available for a `structured_normal`")
  expect_match(f3$sm$parameters$V_method, "unavailable")

  expect_warning(f2 <- suppressMessages(add_covariates(fit, s$x, steps = 2)),
                 "step 2's observed information")
  expect_match(f2$sm$parameters$V_method, "unavailable")
})

test_that("asking for the estimator reported is silent, and gives the same numbers", {
  s   <- .fallback_sim()
  fit <- suppressMessages(suppressWarnings(
    fit_gmm(s$Y, times = 5, n_classes = 2, n_init = 2, random_state = 1)))

  f_cor <- suppressWarnings(suppressMessages(add_covariates(fit, s$x)))
  expect_no_warning(f_rob <- suppressMessages(
    add_covariates(fit, s$x, se = "robust")))
  expect_equal(f_rob$sm$parameters$V_robust, f_cor$sm$parameters$V_robust)

  f2_cor <- suppressWarnings(suppressMessages(
    add_covariates(fit, s$x, steps = 2)))
  expect_no_warning(f2_hes <- suppressMessages(
    add_covariates(fit, s$x, steps = 2, se = "hessian")))
  expect_equal(f2_hes$sm$parameters$V_robust, f2_cor$sm$parameters$V_robust)
})

test_that("a binary step 1 carries the correction and does not warn", {
  s   <- .fallback_sim()
  X   <- sapply(1:5, function(j)
    rbinom(length(s$cls), 1, ifelse(s$cls == 1, 0.2, 0.8)))
  fit <- suppressMessages(suppressWarnings(
    fit_mixture(X, n_classes = 2, measurement = "binary", n_init = 2,
                random_state = 1)))
  expect_no_warning(f <- suppressMessages(add_covariates(fit, s$x)))
  expect_match(f$sm$parameters$V_method, "Bakk-Oberski-Vermunt")
})
