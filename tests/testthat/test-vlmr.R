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
