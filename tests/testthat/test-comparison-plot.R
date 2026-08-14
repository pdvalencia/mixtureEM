# ==============================================================================
# plot() for a model-selection sweep
# ==============================================================================
#
# The picture itself is not testable without comparing rasters, so these check
# the two things that are: that the sweep returns something plot() can dispatch
# on, and that drawing the two-panel version leaves the graphics state as it
# found it. The second matters because `entropy = TRUE` sets mfrow, and a plot
# method that leaks that setting silently splits every subsequent plot the user
# draws in the same session.

.cmp_sweep <- function() {
  set.seed(20260813)
  n   <- 200
  cls <- rbinom(n, 1, 0.5)
  X   <- matrix(rbinom(n * 6, 1, ifelse(rep(cls, 6) == 1L, 0.85, 0.15)), n, 6)
  capture.output(
    out <- suppressMessages(suppressWarnings(
      compare_mixtures(X, k_range = 1:3, measurement = "binary", n_init = 3))))
  out
}

test_that("a sweep returns something plot() can dispatch on", {
  out <- .cmp_sweep()
  expect_s3_class(out, "mixture_comparison")

  # Classing it must not change how it is used. Every documented element still
  # comes out of a plain `$`.
  expect_true(is.data.frame(out$fit_table))
  expect_equal(out$fit_table$Classes, 1:3)
  expect_length(out$models, 3L)
  expect_true(is.numeric(out$best_k))
})

test_that("entropy = TRUE leaves the graphics parameters unchanged", {
  out <- .cmp_sweep()

  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  before <- par(no.readonly = TRUE)

  expect_invisible(plot(out, indices = c("BIC", "AIC"), entropy = TRUE))
  expect_equal(par("mfrow"), before$mfrow)
  expect_equal(par("mar"), before$mar)

  # The single-panel path too, and a bad index is refused rather than silently
  # dropped onto an axis it does not belong on.
  plot(out)
  expect_equal(par("mfrow"), before$mfrow)
  expect_error(plot(out, indices = "LL"))
  expect_error(plot(out, indices = "Entropy"))
})
