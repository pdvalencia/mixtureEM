# The legacy `X` / `n_components` / `Y` / `structural` spelling of fit_mixture()
# has to reach the same model the current spelling does. It is served by its own
# `return()` early in fit_mixture(), and that is exactly the kind of second path
# that drifts: an argument that is a named formal is not carried in `...`, so one
# left out of that call list is dropped silently rather than erroring, and the
# two spellings then fit different models with nothing to show for it.

test_that("the legacy interface fits the same model as the current one", {
  set.seed(21)
  n   <- 200
  cls <- sample(1:2, n, TRUE)
  X   <- matrix(rnorm(n * 4, ifelse(rep(cls, 4) == 1, 1, -1)), n, 4)

  modern <- suppressMessages(
    fit_mixture(X, n_classes = 2, measurement = "continuous",
                n_init = 3, random_state = 2))
  legacy <- suppressMessages(
    fit_mixture(X, n_components = 2, measurement = "continuous",
                n_init = 3, random_state = 2))

  # The homoscedastic default has to be resolved above the bridge: otherwise the
  # legacy call reaches the engine with no value at all and fits free variances.
  expect_true(isTRUE(legacy$mm$variances_equal))
  expect_equal(legacy$metrics$ll, modern$metrics$ll)
  expect_equal(legacy$mm$parameters$covariances,
               modern$mm$parameters$covariances)

  # An explicit value binds on the legacy path too, and gives the other model.
  free <- suppressMessages(
    fit_mixture(X, n_components = 2, measurement = "continuous",
                variances_equal = FALSE, n_init = 3, random_state = 2))
  expect_false(isTRUE(free$mm$variances_equal))
  expect_false(isTRUE(all.equal(free$mm$parameters$covariances[1, ],
                                free$mm$parameters$covariances[2, ])))

  # And the gate that rejects the constraint for indicators it has no meaning
  # for answers to both spellings.
  B <- matrix(rbinom(n * 4, 1, 0.5), n, 4)
  expect_error(fit_mixture(B, n_components = 2, measurement = "binary",
                           variances_equal = TRUE), "no meaning")
})
