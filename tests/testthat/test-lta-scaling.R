# The parameter-vector likelihood and the MLR scaling factor for an LTA.
# `.lta_refine_sim()` and `.ml` come from helper-lta-packing.R, shared with
# test-lta-refine.R.

test_that(".lta_ll_case() reproduces the fit's own log-likelihood", {
  X <- .lta_refine_sim()
  fit <- suppressMessages(suppressWarnings(
    fit_lta(X, n_statuses = 3, times = 3, measurement = "binary",
            smoothing = 0, bayes_constants = .ml, n_init = 3, random_state = 1,
            standard_errors = FALSE)))

  layout <- .lta_par_layout(fit)
  par    <- .lta_par_pack(fit, layout)
  w      <- fit$weights_vec

  # An exact identity, not an approximation: with smoothing = 0 and no
  # categorical prior the fitted log-likelihood is sum(w * ll) with no penalty
  # term, which is what .lta_ll_case() evaluates.
  expect_equal(sum(w * .lta_ll_case(fit, X, par, layout)), fit$loglik,
               tolerance = 1e-8)
})

test_that(".lta_scaling_pieces() returns a well-formed factor on a clean fit", {
  X <- .lta_refine_sim()
  fit <- suppressMessages(suppressWarnings(
    fit_lta(X, n_statuses = 3, times = 3, measurement = "binary",
            smoothing = 0, bayes_constants = .ml, n_init = 3, random_state = 1,
            standard_errors = FALSE)))

  pieces <- .lta_scaling_pieces(.nested_fit_info(fit))
  expect_false(is.null(pieces))
  expect_equal(pieces$p, fit$n_params)
  # Under correct specification and no weights the sandwich collapses towards
  # the information equality, so c sits near 1. The band is wide on purpose:
  # this asserts "not broken", not a target.
  expect_gt(pieces$c, 0.5)
  expect_lt(pieces$c, 2.0)
})

test_that(".lta_scaling_pieces() refuses what it cannot pack", {
  X <- .lta_refine_sim()
  Z <- data.frame(z = stats::rnorm(nrow(X)))

  cov_fit <- suppressMessages(suppressWarnings(
    fit_lta(X, n_statuses = 3, times = 3, measurement = "binary",
            predictors_initial = Z, n_init = 2, random_state = 1,
            standard_errors = FALSE)))
  expect_null(.lta_scaling_pieces(.nested_fit_info(cov_fit)))

  mix_fit <- suppressMessages(suppressWarnings(
    fit_lta(X, n_statuses = 2, times = 3, measurement = "binary",
            n_classes = 2, n_init = 2, random_state = 1,
            standard_errors = FALSE)))
  expect_null(.lta_scaling_pieces(.nested_fit_info(mix_fit)))
})

test_that(".pieces_for() dispatches on the model class", {
  X <- .lta_refine_sim()
  fit <- suppressMessages(suppressWarnings(
    fit_lta(X, n_statuses = 3, times = 3, measurement = "binary",
            smoothing = 0, bayes_constants = .ml, n_init = 2, random_state = 1,
            standard_errors = FALSE)))
  info <- .nested_fit_info(fit)
  expect_equal(.pieces_for(info), .lta_scaling_pieces(info))

  # A cross-sectional mixture must still go to the old path, unchanged.
  lca <- suppressMessages(suppressWarnings(
    fit_mixture(X[, 1:4], n_classes = 2, measurement = "binary",
                n_init = 2, random_state = 1)))
  linfo <- .nested_fit_info(lca)
  expect_equal(.pieces_for(linfo), .scaling_pieces(linfo))
})

test_that("fit_lta() accepts standard_errors = \"robust\", and only the three values", {
  X <- .lta_refine_sim()
  fit <- suppressMessages(suppressWarnings(
    fit_lta(X, n_statuses = 3, times = 3, measurement = "binary",
            smoothing = 0, bayes_constants = .ml, n_init = 2, random_state = 1,
            standard_errors = "robust")))
  expect_true(isTRUE(fit$se$robust))

  expect_error(
    suppressMessages(suppressWarnings(
      fit_lta(X, n_statuses = 3, times = 3, measurement = "binary",
              n_init = 1, random_state = 1, standard_errors = "banana"))),
    "must be TRUE, FALSE"
  )
})
