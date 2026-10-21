test_that("Imhof's integral reproduces the chi-square it generalises", {
  # One unit weight is a chi-square on one degree of freedom, three of them a
  # chi-square on three. These pin the orientation of the 1/2 + (1/pi) form:
  # get the sign wrong and the tails swap.
  for (q in c(0.5, 1, 3.84, 10)) {
    expect_lt(abs(.imhof(q, 1) - pchisq(q, 1, lower.tail = FALSE)), 1e-6)
    expect_lt(abs(.imhof(q, c(1, 1, 1)) - pchisq(q, 3, lower.tail = FALSE)),
              1e-6)
  }
})

test_that("Imhof's integral handles negative weights", {
  # W is not positive semi-definite here — Vuong (1989, Theorem 3.3) says the
  # eigenvalues are real and possibly negative, so this branch is live on every
  # real call and the two tests above do not exercise it. The difference of two
  # independent chi-squares on one degree of freedom is symmetric about zero,
  # which gives an exact target for nothing.
  expect_lt(abs(.imhof(0, c(1, -1)) - 0.5), 1e-6)
  for (q in c(0.5, 2)) {
    expect_lt(abs(.imhof(q, c(1, -1)) + .imhof(-q, c(1, -1)) - 1), 1e-6)
  }
})

test_that("compare_mixtures() adds finite VLMR columns and leaves the last NA", {
  set.seed(4)
  n <- 200
  cls <- rbinom(n, 1, 0.5)
  X <- vapply(1:6, function(j) rbinom(n, 1, ifelse(cls == 1, 0.8, 0.2)),
              numeric(n))

  res <- suppressMessages(capture.output(
    fit <- compare_mixtures(X, k_range = 1:3, measurement = "binary",
                            n_init = 2, vlmr = "both")))
  tab <- fit$fit_table

  expect_true(all(c("VLMR_LR", "VLMR_p", "VLMR_p_robust") %in% names(tab)))
  # Row i tests K_i against K_{i+1}, so the last row has nothing to test.
  expect_true(all(is.finite(tab$VLMR_LR[1:2])))
  expect_true(is.na(tab$VLMR_LR[3]))
  expect_true(all(tab$VLMR_p[1:2] >= 0 & tab$VLMR_p[1:2] <= 1))
  expect_true(all(tab$VLMR_p_robust[1:2] >= 0 & tab$VLMR_p_robust[1:2] <= 1))

  # The likelihood-ratio statistic is ordinary arithmetic on two
  # log-likelihoods and must agree with the table's own LL column.
  expect_lt(abs(tab$VLMR_LR[1] - 2 * (tab$LL[2] - tab$LL[1])), 1e-4)

  # The moments of the reference distribution are returned but not printed.
  expect_equal(length(fit$vlmr), 3L)
  expect_true(all(is.finite(c(fit$vlmr[[1]]$standard$mean,
                              fit$vlmr[[1]]$standard$sd))))
})

test_that("VLMR is off by default", {
  set.seed(5)
  X <- matrix(rbinom(300, 1, 0.5), nrow = 100)
  res <- capture.output(fit <- compare_mixtures(X, k_range = 1:2,
                                                measurement = "binary",
                                                n_init = 2))
  expect_false(any(grepl("VLMR", names(fit$fit_table))))
  expect_null(fit$vlmr)
})

test_that("vlmr_test() on two growth fits: the adjusted statistic and its p-value", {
  set.seed(51)
  n <- 300
  cls <- rbinom(n, 1, 0.5)
  y <- outer(ifelse(cls == 1, 3, 1), rep(1, 4)) +
    outer(ifelse(cls == 1, 0.8, 0.2), 0:3) + matrix(rnorm(n * 4, 0, 0.7), n, 4)
  off <- list(categorical = 0, latent = 0, variances = 0)
  f <- lapply(1:2, function(k) fit_gmm(y, n_classes = k, times = 4,
    random_effects = "none", n_init = 3, random_state = 1,
    bayes_constants = off, n_cores = 1))
  v <- vlmr_test(f[[1]], f[[2]])
  expect_s3_class(v, "vlmr_test")
  expect_equal(v$df, f[[2]]$metrics$n_params - f[[1]]$metrics$n_params)
  expect_equal(v$statistic, 2 * (f[[2]]$metrics$ll - f[[1]]$metrics$ll),
               tolerance = 1e-6)
  # Lo, Mendell & Rubin's finite-sample adjustment, on the same distribution:
  # a smaller statistic, so never a smaller p-value.
  expect_equal(v$adjusted, v$statistic / (1 + 1 / (v$df * log(n))))
  expect_gte(v$p_adjusted, v$p_value)
  expect_output(print(v), "Adjusted")
  expect_error(vlmr_test(f[[2]], f[[1]]), "reversed")
})

test_that("a one-class growth model is finished like a larger one", {
  # EM's covariance update cannot leave the positive semi-definite cone, so
  # without the finish a one-class fit stops short wherever its maximum does.
  set.seed(52)
  n <- 300
  y <- outer(rnorm(n, 2, 1), rep(1, 4)) + outer(rnorm(n, 0.4, 0.2), 0:3) +
    matrix(rnorm(n * 4, 0, 0.7), n, 4)
  f <- fit_gmm(y, n_classes = 1, times = 4, n_init = 1, random_state = 1,
               bayes_constants = list(categorical = 0, latent = 0, variances = 0),
               n_cores = 1)
  expect_false(is.null(f$qn_finish))
})

test_that("vlmr_test() on two LTAs, one status against two included", {
  set.seed(9)
  n <- 300
  s1 <- rbinom(n, 1, 0.4) + 1
  s2 <- ifelse(runif(n) < 0.85, s1, 3 - s1)
  pr <- rbind(c(.8, .8, .75, .7), c(.2, .25, .2, .3))
  gen <- function(s) t(vapply(s, function(k) rbinom(4, 1, pr[k, ]), numeric(4)))
  y <- cbind(gen(s1), gen(s2))
  f <- lapply(1:3, function(k) fit_lta(y, n_statuses = k, times = 2,
    n_init = 5, random_state = 1, smoothing = 0,
    bayes_constants = list(categorical = 0), n_cores = 1,
    standard_errors = FALSE))
  # One status has no status posterior to score, so its scores are central
  # differences; at the maximum they still sum to zero.
  a <- .vlmr_lta_pieces(f[[1]])
  expect_equal(a$p, f[[1]]$n_params)
  expect_lt(max(abs(colSums(a$s))), 1e-4)
  for (k in 1:2) for (v in c("standard", "robust")) {
    r <- vlmr_test(f[[k]], f[[k + 1]], v)
    expect_equal(r$statistic, 2 * (f[[k + 1]]$metrics$ll - f[[k]]$metrics$ll),
                 tolerance = 1e-6)
    expect_equal(r$df, f[[k + 1]]$n_params - f[[k]]$n_params)
    expect_true(is.finite(r$mean) && r$sd > 0 && r$p_value >= 0 &&
                  r$p_value <= 1)
  }
  expect_output(print(r), "statuses")
  g <- fit_mixture(y[, 1:4], n_classes = 2, measurement = "binary",
                   n_init = 2, random_state = 1,
                   n_cores = 1)
  expect_error(vlmr_test(g, f[[3]]), "another LTA")
})
