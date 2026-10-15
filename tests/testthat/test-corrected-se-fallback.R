# se = "corrected" on a measurement model with no step-1 packing reports the
# uncorrected estimator. That substitution is warned about at fit time, and an
# explicit request for the estimator actually reported is not. A growth step 1
# carries the three-step correction; its two-step term is not built yet.

.fallback_sim <- function(n = 300, seed = 1) {
  set.seed(seed)
  cls <- sample(1:2, n, TRUE)
  Y   <- sapply(0:4, function(t)
    rnorm(n, ifelse(cls == 1, 1 + 0.2 * t, 3 - 0.3 * t)))
  list(Y = Y, x = rnorm(n, mean = 0.7 * cls), cls = cls)
}

.se_of <- function(f) { s <- sqrt(diag(f$sm$parameters$V_robust)); s[s > 0] }

test_that("a growth step 1 carries the three-step correction, and widens the SEs", {
  s <- .fallback_sim()
  fits <- list(
    gmm  = suppressMessages(suppressWarnings(
      fit_gmm(s$Y, times = 5, n_classes = 2, n_init = 2, random_state = 1))),
    lcga = suppressMessages(suppressWarnings(
      fit_lcga(s$Y, times = 5, n_classes = 2, family = "gaussian", n_init = 2,
               random_state = 1))))
  for (fit in fits) {
    expect_no_warning(f3 <- suppressMessages(add_covariates(fit, s$x)))
    expect_match(f3$sm$parameters$V_method, "Bakk-Oberski-Vermunt")
    f_rob <- suppressMessages(add_covariates(fit, s$x, se = "robust"))
    expect_true(all(is.finite(.se_of(f3))))
    expect_true(all(.se_of(f3) >= .se_of(f_rob) - 1e-10))
    expect_gt(max(.se_of(f3) / .se_of(f_rob)), 1.0001)
  }
})

test_that("the growth step-1 term does not depend on how step 1 is packed", {
  # The term is J D1 J'; a smooth change of step-1 coordinates moves J and D1
  # by the same Jacobian. Doubling the Psi unit the packing divides by is such
  # a change, so the corrected variance must not move.
  s   <- .fallback_sim()
  fit <- suppressMessages(suppressWarnings(
    fit_gmm(s$Y, times = 5, n_classes = 2, n_init = 2, random_state = 1)))
  f1 <- suppressMessages(add_covariates(fit, s$x))
  local_mocked_bindings(.qn_sn_psi_unit = function(ms) {
    Lr <- ms$design[, ms$r_cols, drop = FALSE]
    2 * tcrossprod(sqrt(mean(ms$parameters$theta) /
                          pmax(apply(Lr^2, 2, max), 1e-12)))
  })
  f2 <- suppressMessages(add_covariates(fit, s$x))
  expect_equal(.se_of(f2), .se_of(f1), tolerance = 1e-4)
})

test_that("a growth step 1 still warns under a two-step se = 'corrected'", {
  s   <- .fallback_sim()
  fit <- suppressMessages(suppressWarnings(
    fit_gmm(s$Y, times = 5, n_classes = 2, n_init = 2, random_state = 1)))
  expect_warning(f2 <- suppressMessages(add_covariates(fit, s$x, steps = 2)),
                 "step 2's observed information")
  expect_match(f2$sm$parameters$V_method, "unavailable")

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
