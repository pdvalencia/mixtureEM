# Cases missing on every indicator are deleted before estimation, which is
# the conventional treatment. Such a case contributes exactly zero to the log-likelihood but, if
# retained, still inflates n in BIC/SABIC, is assigned modally to the largest
# class, and enters the entropy calculation as a maximally uncertain case.
# These tests pin the deletion down against the same fit on data that never
# contained the empty rows.

make_lca_data <- function(n = 400, seed = 42) {
  set.seed(seed)
  p1 <- c(.85, .85, .80, .20, .15, .20)
  p2 <- c(.15, .20, .15, .80, .85, .80)
  z  <- rbinom(n, 1, 0.5)
  X  <- t(vapply(seq_len(n),
                 function(i) rbinom(6, 1, if (z[i] == 1) p1 else p2),
                 numeric(6)))
  colnames(X) <- paste0("q", 1:6)
  X
}

test_that("all-missing rows are removed and leave every reported quantity intact", {
  X <- make_lca_data()
  X_empty <- rbind(X, matrix(NA_real_, 40, 6, dimnames = list(NULL, colnames(X))))

  clean <- fit_mixture(X, n_classes = 2, measurement = "binary",
                       n_init = 3, random_state = 1)
  padded <- suppressWarnings(
    fit_mixture(X_empty, n_classes = 2, measurement = "binary",
                n_init = 3, random_state = 1))

  # The padded fit must be indistinguishable from the fit on clean data.
  expect_equal(padded$metrics$ll,      clean$metrics$ll)
  expect_equal(padded$metrics$bic,     clean$metrics$bic)
  expect_equal(padded$metrics$sabic,   clean$metrics$sabic)
  expect_equal(padded$metrics$entropy, clean$metrics$entropy)
  expect_equal(padded$weights,         clean$weights)
  expect_equal(padded$mm$parameters$pis, clean$mm$parameters$pis)

  # Sample size and modal class counts describe the analysed cases only.
  expect_equal(nrow(padded$data), nrow(X))
  expect_equal(padded$n_eff, nrow(X))
  expect_equal(table(max.col(exp(padded$log_resp))),
               table(max.col(exp(clean$log_resp))))
})

test_that("removal is reported on the fitted object and to the user", {
  X <- make_lca_data(n = 100)
  X[c(3, 17, 88), ] <- NA

  expect_warning(
    fit <- fit_mixture(X, n_classes = 2, measurement = "binary",
                       n_init = 2, random_state = 1),
    "no observed value on any indicator")

  expect_equal(fit$missing_data$n_empty_rows, 3L)
  expect_equal(fit$missing_data$empty_rows, c(3L, 17L, 88L))
  expect_equal(fit$missing_data$n_input_rows, 100L)
  expect_output(print(fit), "Cases Removed")

  # The FIML summary must describe the analysed cases, not the deleted ones:
  # after deletion these data are complete, so no missingness is reported.
  expect_false(isTRUE(fit$missing_data$any_missing))
})

test_that("row-aligned arguments are subset with the indicators", {
  X <- make_lca_data(n = 100)
  X[c(3, 17, 88), ] <- NA
  keep <- setdiff(seq_len(100), c(3, 17, 88))

  set.seed(9)
  w  <- runif(100, 0.5, 2)
  st <- rep(1:2, each = 50)
  cl <- rep(1:20, each = 5)

  fit <- suppressWarnings(
    fit_mixture(X, n_classes = 2, measurement = "binary",
                weights = w, strata = st, cluster = cl,
                n_init = 2, random_state = 1))

  expect_length(fit$sample_weights, length(keep))
  expect_equal(as.integer(fit$strata),  as.integer(st[keep]))
  expect_equal(as.integer(fit$cluster), as.integer(cl[keep]))
  # Sampling weights are rescaled to the retained sample size, not the input one.
  expect_equal(fit$sample_weights, w[keep] / sum(w[keep]) * length(keep))
})

test_that("a mis-specified row-aligned argument still raises its own error", {
  X <- make_lca_data(n = 100)
  X[3, ] <- NA
  # Truncating to the retained rows would otherwise mask the user's mistake.
  expect_error(
    fit_mixture(X, n_classes = 2, measurement = "binary", weights = runif(50)),
    "one entry per case")
})

test_that("data with no observed value anywhere is refused", {
  X <- matrix(NA_real_, 50, 6, dimnames = list(NULL, paste0("q", 1:6)))
  expect_error(
    fit_mixture(X, n_classes = 2, measurement = "binary"),
    "no data to fit")
})

test_that("fit_lta deletes cases observed at no occasion", {
  # fit_lta has its own EM driver, so the deletion is implemented separately
  # there and needs its own regression test.
  set.seed(3)
  n  <- 600
  s1 <- rbinom(n, 1, .5)
  s2 <- ifelse(s1 == 1, rbinom(n, 1, .8), rbinom(n, 1, .3))
  mk <- function(s) t(vapply(s, function(k)
    rbinom(4, 1, if (k == 1) c(.85, .8, .85, .8) else c(.15, .2, .15, .2)),
    numeric(4)))
  X <- cbind(mk(s1), mk(s2))
  colnames(X) <- c(paste0("i", 1:4, "_t1"), paste0("i", 1:4, "_t2"))
  X_empty <- rbind(X, matrix(NA_real_, 60, 8, dimnames = list(NULL, colnames(X))))

  clean <- fit_lta(X, n_statuses = 2, times = 2, n_init = 4,
                   random_state = 1, standard_errors = FALSE)
  padded <- suppressWarnings(
    fit_lta(X_empty, n_statuses = 2, times = 2, n_init = 4,
            random_state = 1, standard_errors = FALSE))

  expect_equal(padded$metrics$bic,     clean$metrics$bic)
  expect_equal(padded$metrics$entropy, clean$metrics$entropy)
  expect_equal(padded$metrics$n_eff,   clean$metrics$n_eff)
  expect_equal(padded$prevalences,     clean$prevalences)
  expect_equal(padded$missing_data$n_empty_rows, 60L)
  expect_output(print(padded), "Cases Removed")
})

test_that("complete data are untouched by the empty-row check", {
  X <- make_lca_data(n = 200)
  expect_silent(
    fit <- fit_mixture(X, n_classes = 2, measurement = "binary",
                       n_init = 2, random_state = 1))
  expect_equal(fit$missing_data$n_empty_rows, 0L)
  expect_equal(nrow(fit$data), 200L)
})
