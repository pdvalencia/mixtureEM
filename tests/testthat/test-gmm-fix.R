# Fixed growth-factor (co)variances in fit_gmm(fix = ).

sim_fix <- function(n = 400, seed = 11) {
  set.seed(seed)
  z <- rbinom(n, 1, 0.4) + 1
  i <- ifelse(z == 1, rnorm(n, 2, 1), 4)
  s <- ifelse(z == 1, rnorm(n, 0.5, 0.5), -0.3)
  outer(i, rep(1, 4)) + outer(s, 0:3) + matrix(rnorm(n * 4, 0, 0.7), n, 4)
}
zero_class <- function(k)
  data.frame(parameter = c("growth_variance", "growth_variance", "growth_covariance"),
             item = c("intercept", "linear", "intercept with linear"),
             class = k, value = 0)

test_that("a class with its growth-factor variances fixed at zero keeps them there", {
  y <- sim_fix()
  f <- fit_gmm(y, times = 4, n_classes = 2, fix = zero_class(2), n_init = 5,
               random_state = 1)
  free <- fit_gmm(y, times = 4, n_classes = 2, n_init = 5, random_state = 1)
  # Psi held equal has three free elements; fixing class 2 leaves class 1's
  # three, so the count is unchanged.
  expect_equal(f$metrics$n_params, free$metrics$n_params)
  k0 <- unique(f$growth$fix$class)
  expect_length(k0, 1L)
  expect_true(all(f$growth$psi[[k0]] == 0))
  expect_true(all(f$growth$psi[[3L - k0]] != 0))
  expect_length(f$growth$boundary, 0L)

  est <- growth_estimates(f, se = "hessian")
  fixed <- est$class == k0 & est$parameter %in% c("growth_variance", "growth_covariance")
  expect_true(all(is.na(est$se[fixed])))
  expect_true(all(is.finite(est$se[!fixed & est$parameter == "growth_mean"])))
})

test_that("fixing a factor's variance at zero in every class is random_effects without it", {
  y <- sim_fix(seed = 12)
  fx <- do.call(rbind, lapply(1:2, function(k)
    data.frame(parameter = c("growth_variance", "growth_covariance"),
               item = c("linear", "intercept with linear"), class = k, value = 0)))
  a <- fit_gmm(y, times = 4, n_classes = 2, fix = fx, n_init = 5, random_state = 2)
  b <- fit_gmm(y, times = 4, n_classes = 2, random_effects = "intercept",
               n_init = 5, random_state = 2)
  expect_lt(abs(a$metrics$ll - b$metrics$ll), 1e-3)
  expect_equal(a$metrics$n_params, b$metrics$n_params)
})

test_that("fix refuses what it cannot hold", {
  y <- sim_fix(n = 100)
  go <- function(fix, ...) fit_gmm(y, times = 4, n_classes = 2, fix = fix,
                                   n_init = 1, ...)
  expect_error(go(data.frame(parameter = "growth_variance", item = "linear",
                             class = 2, value = 0)), "covariance with intercept")
  expect_error(go(data.frame(parameter = "growth_mean", item = "linear",
                             class = 2, value = 0)), "growth_variance and")
  expect_error(go(data.frame(parameter = "growth_variance", item = "slope",
                             class = 1, value = 1)), "not a random growth factor")
  expect_error(go(zero_class(3)), "between 1 and 2")
  expect_error(go(data.frame(parameter = "growth_variance", item = "linear",
                             class = 1, value = -1)), "cannot be negative")
  expect_error(go(zero_class(2), random_effects = "none"), "there are none")
  expect_error(go(list(a = 1)), "data frame")
})

test_that("blrt() refuses a fit with fixed values", {
  y <- sim_fix(n = 150)
  f <- fit_gmm(y, times = 4, n_classes = 2, fix = zero_class(2), n_init = 2,
               random_state = 3)
  expect_error(blrt(k_small = 1, k_large = 2, from_fit = f, n_reps = 2,
                    verbose = FALSE), "fixed values")
})
