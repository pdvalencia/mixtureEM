# The block of fit indices under print()'s header.
#
# The trap the second assertion exists for: a three-step fit carries both
# step-1 and combined metrics, and printing a step-1 log-likelihood beside a
# step-3 BIC is worse than printing no criteria at all.

test_that("print() shows the same six indices compare_mixtures() tabulates", {
  set.seed(1)
  X   <- matrix(rbinom(600, 1, 0.5), ncol = 6)
  fit <- suppressWarnings(fit_mixture(X, n_classes = 2,
                                      measurement = "binary",
                                      n_init = 5,
                                      random_state = 1))
  out <- capture.output(print(fit))

  expect_true(any(grepl("Parameters", out)))
  expect_true(any(grepl("AIC", out)))
  expect_true(any(grepl("BIC", out)))
  expect_true(any(grepl("SABIC", out)))
  expect_true(any(grepl(sprintf("BIC%s: %.2f", "\\s*", fit$metrics$bic), out)))
  expect_true(any(grepl(sprintf("SABIC%s: %.2f", "\\s*", fit$metrics$sabic),
                        out)))
})

test_that("a three-step fit reads its criteria off the step it reports", {
  set.seed(2)
  n     <- 200
  cl    <- rbinom(n, 1, 0.4)
  items <- sapply(1:5, function(j) rbinom(n, 1, ifelse(cl == 1, 0.8, 0.2)))
  y     <- rnorm(n, 25 + 2 * cl)
  fit0  <- suppressMessages(suppressWarnings(
    fit_mixture(items, n_classes = 2,
                measurement = "binary",
                n_init = 5, random_state = 2)))
  fit3  <- suppressMessages(suppressWarnings(add_outcome(fit0, y)))

  out <- capture.output(print(fit3))
  expect_true(any(grepl("(Step 1)", out, fixed = TRUE)))
  # Every line comes from step 1 - the BIC as much as the log-likelihood.
  expect_true(any(grepl(sprintf("%.2f", fit3$step1_metrics$bic), out,
                        fixed = TRUE)))
  expect_false(any(grepl(sprintf("BIC\\s*\\(Step 1\\)\\s*: %.2f",
                                 fit3$metrics$bic), out)))
})
