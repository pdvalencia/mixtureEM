# ==============================================================================
# Boundary detection for categorical item probabilities under a continuous
# random intercept (R/categorical_boundary.R)
# ==============================================================================
#
# Mirrors test-gmm.R's degenerate-fit pattern: fit a small, fast model, then
# construct a forced-degenerate copy by overwriting its parameters directly
# rather than trying to drive EM to the boundary, which is neither fast nor
# reliable to reproduce on demand. See RECORDS.md, "R11" for why the guard
# checks only fit$ri$A / fit$ri$theta and not mm$parameters$pis.

.small_ri_fit <- function() {
  sim <- .lta_ri_sim(n = 300, Tn = 3, J = 4, seed = 7)
  fit_lta(sim$X, n_statuses = 2, times = 3, measurement = "binary",
         random_intercept = "continuous", n_quadrature = 5,
         n_init = 2, n_cores = 1, random_state = 1, max_iter = 100,
         standard_errors = FALSE)
}

test_that("a clean continuous-RI fit is not flagged", {
  fit <- suppressWarnings(.small_ri_fit())
  expect_null(.categorical_boundary(fit, fit$data))
  expect_null(fit$degenerate)
})

test_that("a response probability pinned at the boundary is flagged", {
  fit <- suppressWarnings(.small_ri_fit())
  degenerate <- fit
  degenerate$ri$A[1, 1] <- 20   # p = 1 - 2.1e-9, well past the |logit| > 8 rule
  flagged <- .categorical_boundary(degenerate, degenerate$data)
  expect_false(is.null(flagged))
  expect_equal(nrow(flagged), 1L)
  expect_equal(flagged$kind, "probability")
  expect_gt(flagged$prob, 0.999)
})

test_that("an ordinary logit well inside the boundary is not flagged", {
  fit <- suppressWarnings(.small_ri_fit())
  fit$ri$A[1, 1] <- 3   # p = 0.953, nowhere near |logit| > 8
  expect_null(.categorical_boundary(fit, fit$data))
})

test_that("the warning names the categorical prior remedy and stays short", {
  fit <- suppressWarnings(.small_ri_fit())
  fit$ri$A[1, 1] <- -25
  expect_warning(
    out <- .check_gaussian_degeneracy(fit, fit$data),
    "categorical")
  expect_false(is.null(out$degenerate))
  expect_equal(out$degenerate$kind, "probability")
})

test_that("no random intercept means no categorical boundary check", {
  Xo <- .lta_ordinal_sim(n = 150, Tn = 3, K = 2, cats = c(3L, 3L, 2L), seed = 1)
  fit <- suppressWarnings(fit_lta(Xo, n_statuses = 2, times = 3,
                                  measurement = "ordinal", n_init = 2,
                                  n_cores = 1, random_state = 1, max_iter = 100,
                                  standard_errors = FALSE))
  expect_null(.categorical_boundary(fit, fit$data))
})
