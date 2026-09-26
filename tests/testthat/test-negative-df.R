# Goodman's count: a categorical model with more free parameters than the
# response-pattern table has free cells is not identified, and says so.
# Four binary items give W - 1 = 15 free cells; K classes cost 4K + K - 1.

.ndf_items <- function(n = 400, seed = 7) {
  set.seed(seed)
  cls <- sample(1:2, n, TRUE)
  X   <- sapply(1:4, function(j) rbinom(n, 1, ifelse(cls == 1, 0.2, 0.8)))
  colnames(X) <- paste0("u", 1:4)
  X
}

test_that("four classes on four binary items warn, three do not", {
  X <- .ndf_items()
  expect_warning(
    f4 <- suppressMessages(fit_mixture(X, n_classes = 4, measurement = "binary",
                                       n_init = 2, random_state = 1)),
    "19 free parameters.*15 free cells \\(df = -4\\)")
  pd <- .pattern_df(f4)
  expect_equal(c(pd$df, pd$cells, pd$n_params), c(-4, 15, 19))

  f3 <- suppressMessages(suppressWarnings(
    fit_mixture(X, n_classes = 3, measurement = "binary", n_init = 2,
                random_state = 1)))
  expect_equal(.pattern_df(f3)$df, 1)
  expect_no_warning(.warn_negative_df(f3))
})

test_that("the count agrees with absolute_fit() on a measurement-only fit", {
  X <- .ndf_items()
  f <- suppressMessages(suppressWarnings(
    fit_mixture(X, n_classes = 2, measurement = "binary", n_init = 2,
                random_state = 1)))
  expect_equal(.pattern_df(f)$df, absolute_fit(f)$df)
})

test_that("continuous indicators and one-step covariate fits are not counted", {
  set.seed(3)
  Y <- matrix(rnorm(600), ncol = 3)
  fc <- suppressMessages(suppressWarnings(
    fit_mixture(Y, n_classes = 2, measurement = "continuous", n_init = 2,
                random_state = 1)))
  expect_null(.pattern_df(fc))

  X  <- .ndf_items()
  fz <- suppressMessages(suppressWarnings(
    fit_mixture(X, n_classes = 2, measurement = "binary",
                predictors = data.frame(z = rnorm(nrow(X))), n_steps = 1,
                n_init = 2, random_state = 1)))
  expect_null(.pattern_df(fz))
})

test_that("compare_mixtures() names every K past the boundary, once", {
  X <- .ndf_items()
  w <- NULL
  withCallingHandlers(
    capture.output(suppressMessages(
      compare_mixtures(X, k_range = 2:4, measurement = "binary", n_init = 2,
                       random_state = 1, n_cores = 1))),
    warning = function(cnd) {
      w <<- c(w, conditionMessage(cnd))
      invokeRestart("muffleWarning")
    })
  hits <- grep("negative df", w, value = TRUE)
  expect_length(hits, 1L)
  expect_match(hits, "The 4-class model has")
})
