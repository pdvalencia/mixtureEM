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

test_that("compare_mixtures() carries the class separation for continuous indicators only", {
  set.seed(2)
  z <- sample(2, 300, TRUE)
  Y <- matrix(rnorm(900, rep(c(0, 1.5)[z], 3)), ncol = 3)
  cm <- NULL
  invisible(capture.output(cm <- suppressWarnings(compare_mixtures(
    Y, k_range = 1:3, measurement = "continuous", n_init = 3, n_cores = 1,
    random_state = 1))))
  sep <- cm$fit_table$Separation
  expect_true(is.na(sep[1]))
  for (k in 2:3) {
    cd <- NULL
    invisible(capture.output(cd <- classification_diagnostics(cm$models[[paste0("K", k)]])))
    expect_equal(sep[k], attr(cd$separation, "min"))
  }
  expect_true(any(grepl("Separation", capture.output(print(cm)))))
  X <- matrix(rbinom(600, 1, 0.5), ncol = 6)
  invisible(capture.output(cb <- suppressWarnings(compare_mixtures(
    X, k_range = 1:2, measurement = "binary", n_init = 2, n_cores = 1, random_state = 1))))
  expect_null(cb$fit_table$Separation)
})

test_that("compare_longitudinal() carries the status separation of a continuous LTA", {
  set.seed(3)
  n  <- 300
  s1 <- sample(2, n, TRUE, c(.6, .4))
  s2 <- ifelse(runif(n) < .8, s1, 3L - s1)
  Y  <- cbind(matrix(rnorm(3 * n, c(0, 1.8)[s1]), n),
              matrix(rnorm(3 * n, c(0, 1.8)[s2]), n))
  cl <- NULL
  invisible(capture.output(cl <- suppressMessages(suppressWarnings(
    compare_longitudinal(Y, k_range = 2:3, times = 2, measurement = "continuous",
                         n_init = 2, n_cores = 1, random_state = 1,
                         standard_errors = FALSE)))))
  for (k in 2:3) {
    f <- cl$models[[paste0("K", k)]]
    # By hand: each occasion's status means against its variances pooled with
    # that occasion's prevalences; the column is the smallest over occasions.
    d <- sapply(1:2, function(t) {
      p  <- f$mm$models[[t]]$parameters
      s2 <- colSums(f$prevalences[t, ] * p$covariances)
      min(utils::combn(k, 2, function(ab)
        sqrt(sum((p$means[ab[1], ] - p$means[ab[2], ])^2 / s2))))
    })
    expect_equal(cl$fit_table$Separation[k - 1], min(d), tolerance = 1e-12)
    cd <- NULL
    invisible(capture.output(cd <- classification_diagnostics(f)))
    expect_equal(attr(cd[[2]]$separation, "min"), d[2], tolerance = 1e-12)
  }
  expect_true(any(grepl("Separation", capture.output(print(cl)))))
})
