# Cases missing a class predictor: listwise deletion (the default) or the
# mean fill of `missing_predictors = "mean"`.

.mp_data <- function(n = 300, seed = 11) {
  set.seed(seed)
  cls <- rbinom(n, 1, 0.4) + 1
  p   <- rbind(c(.85, .8, .8, .2, .2), c(.15, .2, .25, .8, .85))
  X   <- t(sapply(cls, function(k) rbinom(5, 1, p[k, ])))
  x   <- rnorm(n, mean = cls - 1.5)
  x[sample(n, 30)] <- NA
  list(X = X, x = x)
}

test_that("add_covariates() leaves incomplete cases out of step 3 only", {
  d   <- .mp_data()
  fit <- fit_mixture(d$X, n_classes = 2, measurement = "binary",
                     n_init = 5, random_state = 1)
  ok  <- !is.na(d$x)
  for (corr in c("ML", "BCH", "none")) {
    expect_message(
      f <- suppressWarnings(add_covariates(fit, data.frame(x = d$x),
                                           correction = corr)),
      "30 cases missing a predictor were left out")
    expect_identical(f$struct_keep, ok)
    expect_true(all(is.finite(f$sm$parameters$V_robust)))
    # A dropped case keeps its measurement-model posterior.
    expect_equal(f$log_resp[!ok, ], fit$log_resp[!ok, ], tolerance = 1e-8)
  }
  f2 <- suppressMessages(add_covariates(fit, data.frame(x = d$x), steps = 2))
  expect_equal(f2$log_resp[!ok, ], fit$log_resp[!ok, ], tolerance = 1e-8)
  expect_true(all(is.finite(f2$sm$parameters$V_robust)))

  # "mean" keeps every case, and complete predictors never meet the switch.
  fm <- suppressMessages(add_covariates(fit, data.frame(x = d$x),
                                        missing_predictors = "mean"))
  expect_null(fm$struct_keep)
  expect_false(isTRUE(all.equal(fm$sm$parameters$beta,
                                suppressMessages(add_covariates(
                                  fit, data.frame(x = d$x)))$sm$parameters$beta)))
  xc <- ifelse(ok, d$x, 0)
  expect_equal(
    suppressMessages(add_covariates(fit, xc))$sm$parameters$beta,
    suppressMessages(add_covariates(fit, xc, missing_predictors = "mean"))$sm$parameters$beta)
})

test_that("uncorrected step 3 under listwise is the regression on complete cases", {
  d   <- .mp_data()
  fit <- fit_mixture(d$X, n_classes = 2, measurement = "binary",
                     n_init = 5, random_state = 1)
  ok  <- !is.na(d$x)
  f   <- suppressWarnings(suppressMessages(
    add_covariates(fit, data.frame(x = d$x), correction = "none")))
  ref <- .fit_mnl(cbind(1, d$x[ok]), exp(fit$log_resp[ok, ]),
                  augment = .bayes_alpha(f$sm, "latent") > 0,
                  alpha = .bayes_alpha(f$sm, "latent"))
  b   <- f$sm$parameters$beta
  expect_equal(unname(b[1, ] - b[2, ]), ref$beta[1, ] - ref$beta[2, ],
               tolerance = 1e-5)
})

test_that("a one-step fit drops the incomplete cases whole", {
  d  <- .mp_data()
  ok <- !is.na(d$x)
  expect_message(
    f <- fit_mixture(d$X, n_classes = 2, measurement = "binary",
                     predictors = d$x, n_steps = 1, n_init = 5,
                     random_state = 1),
    "left out of the model \\(n = 270 used\\)")
  g <- fit_mixture(d$X[ok, ], n_classes = 2, measurement = "binary",
                   predictors = d$x[ok], n_steps = 1, n_init = 5,
                   random_state = 1)
  expect_equal(f$metrics$ll, g$metrics$ll, tolerance = 1e-8)
  expect_equal(nrow(f$log_resp), 270L)
})
