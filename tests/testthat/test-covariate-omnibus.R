# ==============================================================================
# Term grouping and the omnibus Wald test in the covariate output
# ==============================================================================
#
# A covariate with more than two categories enters the model as k-1 dummy
# columns, and the per-class coefficient table therefore contains no test of the
# covariate itself - only of one level against the reference level, one class at
# a time. The omnibus Wald test is the one that answers the question the user
# asked, and it needs to know which columns belong to which original variable.
#
# That grouping is recorded by prepare_covariates() and carried onto the fitted
# model, rather than recovered from the column names afterwards. These tests
# cover both halves: that the grouping survives every path a covariate matrix
# can take into the model, and that the test built on it is the same statistic
# analytical_wald_test() has always computed.
#
# The external anchor for the statistic's value and its (K-1) x (levels-1)
# degrees of freedom lives in test-covariate-se.R, against a published run.

.omni_sim <- function(n = 500, seed = 20260805) {
  set.seed(seed)
  cls <- sample(1:3, n, replace = TRUE)
  p   <- rbind(c(.85, .85, .85, .15, .15, .15),
               c(.15, .15, .15, .85, .85, .85),
               c(.85, .15, .85, .15, .85, .15))
  X <- t(vapply(seq_len(n), function(i) rbinom(6, 1, p[cls[i], ]), numeric(6)))

  # Marital is the three-level case. Age and Age_Decades are the pair that name
  # matching cannot separate: "Age" is a substring of "Age_Decades", so a
  # grep-based lookup sweeps both into one test.
  Z <- data.frame(
    Marital     = factor(sample(c("Married", "Single", "Other"), n, TRUE),
                         levels = c("Married", "Single", "Other")),
    Age         = rnorm(n) + 0.8 * (cls == 3),
    Age_Decades = rnorm(n))
  Z$Marital[cls == 2] <- sample(c("Single", "Other"), sum(cls == 2), TRUE)

  list(X = X, Z = Z, cls = cls)
}

.omni_fit <- function(d, ...) {
  suppressMessages(fit_mixture(
    d$X, n_classes = 3, measurement = "binary", predictors = d$Z,
    n_steps = 3, correction = "ML", n_init = 5, random_state = 1, ...))
}

# ------------------------------------------------------------------------------
# The grouping itself
# ------------------------------------------------------------------------------

test_that("prepare_covariates records which columns came from which variable", {
  d <- .omni_sim(n = 60)
  P <- prepare_covariates(d$Z)

  expect_equal(colnames(P),
               c("Marital.Single", "Marital.Other", "Age", "Age_Decades"))
  expect_equal(attr(P, "covariate_terms"),
               c("Marital", "Marital", "Age", "Age_Decades"))
})

test_that("a plain numeric matrix falls back to one term per column", {
  M <- matrix(rnorm(20), 10, 2, dimnames = list(NULL, c("a", "b")))
  expect_equal(.covariate_terms(M), c("a", "b"))

  # No names at all: positional labels, never NULL, so callers need no guard.
  U <- matrix(rnorm(20), 10, 2)
  expect_equal(.covariate_terms(U), c("V1", "V2"))
})

test_that("the grouping survives cbind, which drops attributes", {
  d <- .omni_sim(n = 60)
  P <- prepare_covariates(d$Z)
  out <- matrix(rnorm(60), ncol = 1L, dimnames = list(NULL, "y"))

  expect_null(attr(cbind(out, P), "covariate_terms"))
  expect_equal(.covariate_terms(.cbind_covariates(out, P)),
               c("y", "Marital", "Marital", "Age", "Age_Decades"))
})

test_that("the grouping reaches the fitted model, aligned with beta", {
  fit <- .omni_fit(.omni_sim())
  terms <- fit$sm$parameters$terms

  expect_equal(length(terms), ncol(fit$sm$parameters$beta))
  expect_equal(terms,
               c("Intercept", "Marital", "Marital", "Age", "Age_Decades"))
})

# ------------------------------------------------------------------------------
# The test built on it
# ------------------------------------------------------------------------------

test_that("a covariate is tested as one term, not one dummy at a time", {
  fit <- .omni_fit(.omni_sim())
  K   <- 3L

  w <- analytical_wald_test(fit, "Marital")
  expect_equal(w$df, (K - 1L) * 2L)      # (K-1) x (levels-1), LG's convention
  expect_equal(w$Covariate, "Marital")

  # Naming one dummy column tests that column alone.
  expect_equal(analytical_wald_test(fit, "Marital.Single")$df, K - 1L)
})

test_that("a covariate whose name is a prefix of another is not swept in", {
  fit <- .omni_fit(.omni_sim())

  # This is the whole point of recording the grouping: grep("Age", ...) matches
  # Age and Age_Decades both, and would silently report a 4 df test of two
  # covariates as a 2 df test of one.
  expect_equal(analytical_wald_test(fit, "Age")$df, 2L)
  expect_equal(analytical_wald_test(fit, "Age_Decades")$df, 2L)
  expect_false(isTRUE(all.equal(analytical_wald_test(fit, "Age")$Wald_Chi2,
                                analytical_wald_test(fit, "Age_Decades")$Wald_Chi2)))
})

test_that("an unknown covariate names the ones that exist", {
  fit <- .omni_fit(.omni_sim())
  expect_error(analytical_wald_test(fit, "Nonesuch"),
               "Marital, Age, Age_Decades")
})

test_that("name matching still works when no grouping was recorded", {
  fit <- .omni_fit(.omni_sim())
  fit$sm$parameters$terms <- NULL     # a model fitted before the grouping existed

  expect_equal(analytical_wald_test(fit, "Age_Decades")$df, 2L)
  # ... and says so when the match is ambiguous, rather than testing two
  # covariates as though they were one.
  expect_warning(analytical_wald_test(fit, "Age"), "more than one variable")
})

# ------------------------------------------------------------------------------
# The printed output
# ------------------------------------------------------------------------------

test_that("summary() prints one omnibus test per covariate", {
  fit <- .omni_fit(.omni_sim())
  out <- capture.output(summary(fit))

  expect_true(any(grepl("OMNIBUS TEST PER COVARIATE", out)))
  expect_true(any(grepl("^  Marital ", out)))
  expect_true(any(grepl("^  Age ", out)))
  expect_true(any(grepl("^  Age_Decades ", out)))
  # The Hauck-Donner caveat travels with the test, since it is what stops a
  # non-significant omnibus from being read as a clean null.
  expect_true(any(grepl("Hauck-Donner", out)))
})

test_that("the printed statistic is the one analytical_wald_test() returns", {
  fit <- .omni_fit(.omni_sim())
  capture.output(tab <- .print_covariate_omnibus(fit, fit$sm, ref_class = 1L))

  expect_equal(tab$term, c("Marital", "Age", "Age_Decades"))
  for (nm in tab$term) {
    w <- analytical_wald_test(fit, nm, ref_class = 1L)
    expect_lt(abs(tab$chi2[tab$term == nm] - w$Wald_Chi2), 1e-3)
    expect_equal(tab$df[tab$term == nm], w$df)
  }
})

test_that("the omnibus test inherits whichever variance the model was fitted with", {
  d <- .omni_sim()
  for (m in c("corrected", "robust", "hessian")) {
    fit <- .omni_fit(d, se = m)
    expect_equal(analytical_wald_test(fit, "Marital")$Method,
                 fit$sm$parameters$V_method)
  }
})

test_that("nothing is printed when every omnibus test duplicates a printed z", {
  # Two classes and single-column covariates: each omnibus statistic is the
  # square of a z already in the coefficient table, so the block is noise.
  set.seed(20260805)
  n   <- 300
  z   <- rnorm(n)
  cls <- 1L + rbinom(n, 1, plogis(-0.4 + 0.9 * z))
  X   <- matrix(rbinom(n * 6, 1, ifelse(rep(cls, 6) == 1L, .85, .15)), n, 6)
  fit <- suppressMessages(fit_mixture(
    X, n_classes = 2, measurement = "binary", predictors = data.frame(z = z),
    n_steps = 3, correction = "ML", n_init = 3, random_state = 1))

  out <- capture.output(summary(fit))
  expect_false(any(grepl("OMNIBUS TEST PER COVARIATE", out)))

  # The statistic is still available on request, and is that square.
  w <- analytical_wald_test(fit, "z")
  expect_equal(w$df, 1L)
})
