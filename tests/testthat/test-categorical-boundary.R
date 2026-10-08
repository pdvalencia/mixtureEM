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

# A boundary cell under a random intercept is a note, not a warning: on
# simulated data every such fit was usable, and the flag cannot tell a cell
# empty in the population from one empty in the sample. It must not mark the
# fit degenerate either, or lr_test() and the BIC line would disown a fit
# whose likelihood is fine.
test_that("a boundary cell is noted on print(), without a warning or a degenerate flag", {
  fit <- suppressWarnings(.small_ri_fit())
  fit$ri$A[1, 1] <- -25
  expect_no_warning(out <- .check_gaussian_degeneracy(fit, fit$data))
  expect_null(out$degenerate)
  expect_equal(out$ri_boundary$kind, "probability")
  printed <- paste(capture.output(print(out)), collapse = " ")
  expect_match(printed, "at the boundary", fixed = TRUE)
  expect_match(printed, "in status 1", fixed = TRUE)
  expect_match(printed, "bayes_constants = list(categorical = 2)", fixed = TRUE)
  ms <- paste(capture.output(measurement_summary(out)), collapse = " ")
  expect_match(ms, "Read the probability, not the logit", fixed = TRUE)
  expect_no_match(ms, "stronger prior than the default", fixed = TRUE)
})

test_that("a collapsed class variance still warns", {
  set.seed(3)
  X <- cbind(c(rnorm(150, 0), rnorm(150, 4)), c(rnorm(150, 0), rnorm(150, 4)))
  fit <- suppressWarnings(fit_mixture(X, n_classes = 2,
    measurement = "continuous", n_init = 2, random_state = 1))
  fit$mm$parameters$covariances[1, 1] <- 1e-8
  expect_warning(out <- .check_gaussian_degeneracy(fit, X), "variance")
  expect_equal(out$degenerate$kind, "variance")
})

test_that("no random intercept means no categorical boundary check", {
  Xo <- .lta_ordinal_sim(n = 150, Tn = 3, K = 2, cats = c(3L, 3L, 2L), seed = 1)
  fit <- suppressWarnings(fit_lta(Xo, n_statuses = 2, times = 3,
                                  measurement = "ordinal", n_init = 2,
                                  n_cores = 1, random_state = 1, max_iter = 100,
                                  standard_errors = FALSE))
  expect_null(.categorical_boundary(fit, fit$data))
})
