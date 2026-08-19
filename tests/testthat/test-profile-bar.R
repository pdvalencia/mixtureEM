# ==============================================================================
# plot(fit, type = "bar") -- the standardized profile bar chart
# ==============================================================================
# The drawing itself is not worth asserting on; the heights behind it are. The
# property that matters is invariance: because the bars are z-scores against
# the observed marginal, the figure is the same whether or not the indicators
# were standardized before fitting.

make_continuous_fit <- function(seed = 2, n = 200) {
  set.seed(seed)
  X <- cbind(a = rnorm(n), b = rnorm(n) + 2, c = rnorm(n))
  fit_mixture(X, n_classes = 2, measurement = "continuous", n_init = 2)
}

test_that("the bar chart draws for a continuous measurement model", {
  fit <- make_continuous_fit()
  path <- tempfile(fileext = ".png")

  png(path)
  expect_invisible(plot(fit, type = "bar"))
  dev.off()

  png(path)
  expect_invisible(plot(fit, type = "bar", scale = "within"))
  dev.off()

  # The profile path is untouched by the new argument.
  png(path)
  expect_invisible(plot(fit))
  dev.off()
})

test_that("heights are invariant to standardizing the indicators first", {
  set.seed(2)
  X  <- cbind(a = rnorm(200), b = rnorm(200) + 2, c = rnorm(200))
  f  <- fit_mixture(X, n_classes = 2, measurement = "continuous", n_init = 2)
  fz <- fit_mixture(scale(X), n_classes = 2, measurement = "continuous",
                    n_init = 2)

  H  <- .profile_bar_heights(f,  "total")
  Hz <- .profile_bar_heights(fz, "total")

  # Loose: the two fits are separate EM runs, so this is agreement up to
  # convergence noise, not up to floating point.
  expect_lt(max(abs(H - Hz)), 1e-2)
})

test_that("total and within differ only by the denominator", {
  fit <- make_continuous_fit()
  Ht  <- .profile_bar_heights(fit, "total")
  Hw  <- .profile_bar_heights(fit, "within")

  # Same centring, so the ratio is constant down each column.
  ratio <- Hw / Ht
  expect_equal(diff(range(apply(ratio, 2, function(r) diff(range(r))))), 0,
               tolerance = 1e-8)
  # Within-class dispersion is the smaller denominator, so the bars are taller.
  expect_true(all(abs(Hw) >= abs(Ht) - 1e-8))
})

test_that("the bar chart refuses a measurement model it cannot standardize", {
  set.seed(3)
  Xb  <- matrix(rbinom(600, 1, 0.4), ncol = 3)
  fit <- fit_mixture(Xb, n_classes = 2, measurement = "binary", n_init = 2)

  expect_error(plot(fit, type = "bar"), "all-continuous")
  # And names the alternative that does work.
  expect_error(plot(fit, type = "bar"), "profile")
})

test_that("an unnamed indicator matrix still standardizes, positionally", {
  set.seed(4)
  X   <- matrix(rnorm(400), ncol = 2)
  fit <- fit_mixture(X, n_classes = 2, measurement = "continuous", n_init = 2)

  expect_silent(H <- .profile_bar_heights(fit, "total"))
  expect_equal(dim(H), c(2L, 2L))
  expect_true(all(is.finite(H)))
})

# ==============================================================================
# plot(fit, type = "line") -- the same z-scores drawn as lines
# ==============================================================================
# The two renderers must never disagree about the numbers, so what is asserted
# here is that the line plot draws, and that it draws the bar chart's heights.

test_that("the line plot draws and reuses the bar chart's heights", {
  fit  <- make_continuous_fit()
  path <- tempfile(fileext = ".png")

  png(path)
  expect_invisible(plot(fit, type = "line"))
  dev.off()

  png(path)
  expect_invisible(plot(fit, type = "line", scale = "within"))
  dev.off()

  # One helper behind both, so `type` cannot change a height.
  expect_equal(.profile_bar_heights(fit, "total", type = "line"),
               .profile_bar_heights(fit, "total", type = "bar"))
})

test_that("the refusal names the type that was actually asked for", {
  set.seed(3)
  Xb  <- matrix(rbinom(600, 1, 0.4), ncol = 3)
  fit <- fit_mixture(Xb, n_classes = 2, measurement = "binary", n_init = 2)

  expect_error(plot(fit, type = "line"), "type = \"line\"", fixed = TRUE)
  expect_error(plot(fit, type = "bar"),  "type = \"bar\"",  fixed = TRUE)
})
