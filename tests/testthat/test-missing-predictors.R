# Cases missing a class predictor or an outcome covariate: listwise deletion
# (the default) or the mean fill of `missing_predictors = "mean"`.

.mp_data <- function(n = 300, seed = 11) {
  set.seed(seed)
  cls <- rbinom(n, 1, 0.4) + 1
  p   <- rbind(c(.85, .8, .8, .2, .2), c(.15, .2, .25, .8, .85))
  X   <- t(sapply(cls, function(k) rbinom(5, 1, p[k, ])))
  x   <- rnorm(n, mean = cls - 1.5)
  x[sample(n, 30)] <- NA
  list(X = X, x = x, cls = cls)
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

# A distal outcome on top of .mp_data(): continuous (`y`) and binary (`yb`),
# both depending on the class and on `x`, with a few outcomes missing too.
.mp_outcome <- function(d, seed = 12) {
  set.seed(seed)
  n  <- length(d$x)
  xs <- ifelse(is.na(d$x), 0, d$x)
  y  <- rnorm(n, mean = 2 * d$cls + 0.5 * xs)
  yb <- rbinom(n, 1, plogis(-1 + 1.5 * (d$cls - 1) + 0.5 * xs))
  miss_y <- sample(which(!is.na(d$x)), 10)
  y[miss_y] <- NA
  list(y = y, yb = yb, miss_y = miss_y)
}

test_that("add_outcome() leaves cases missing a covariate out of step 3 only", {
  d   <- .mp_data()
  o   <- .mp_outcome(d)
  fit <- fit_mixture(d$X, n_classes = 2, measurement = "binary",
                     n_init = 5, random_state = 1)
  ok  <- !is.na(d$x)
  # The dropped cases weigh nothing in step 3, so their outcome can be
  # anything without moving an estimate.
  y2  <- o$y; y2[!ok] <- 100
  for (corr in c("ML", "BCH", "none")) {
    expect_message(
      f <- suppressWarnings(add_outcome(fit, o$y, covariates = data.frame(x = d$x),
                                        outcome_type = "continuous",
                                        correction = corr)),
      "30 cases missing a covariate were left out of the outcome regression")
    expect_identical(f$struct_keep, ok)
    expect_equal(f$log_resp[!ok, ], fit$log_resp[!ok, ], tolerance = 1e-8)
    g <- suppressWarnings(suppressMessages(
      add_outcome(fit, y2, covariates = data.frame(x = d$x),
                  outcome_type = "continuous", correction = corr)))
    expect_equal(g$sm$parameters$beta_pooled, f$sm$parameters$beta_pooled,
                 tolerance = 1e-8)
    expect_equal(g$metrics$ll, f$metrics$ll, tolerance = 1e-8)
  }

  # A categorical outcome, under ML, with its sandwich over the kept cases.
  fb <- suppressMessages(add_outcome(fit, o$yb, covariates = data.frame(x = d$x),
                                     outcome_type = "categorical",
                                     correction = "ML"))
  yb2 <- o$yb; yb2[!ok] <- 1 - yb2[!ok]
  gb <- suppressMessages(add_outcome(fit, yb2, covariates = data.frame(x = d$x),
                                     outcome_type = "categorical",
                                     correction = "ML"))
  expect_equal(gb$sm$parameters$beta_pooled, fb$sm$parameters$beta_pooled,
               tolerance = 1e-6)
  expect_equal(gb$sm$parameters$hessian, fb$sm$parameters$hessian,
               tolerance = 1e-5)
  expect_true(all(is.finite(fb$sm$parameters$hessian)))

  # A missing outcome alone deletes nothing, and "mean" keeps every case.
  xc <- ifelse(ok, d$x, 0)
  fy <- suppressMessages(add_outcome(fit, o$y, covariates = data.frame(x = xc),
                                     outcome_type = "continuous"))
  expect_null(fy$struct_keep)
  fm <- suppressMessages(add_outcome(fit, o$y, covariates = data.frame(x = d$x),
                                     outcome_type = "continuous",
                                     missing_predictors = "mean"))
  expect_null(fm$struct_keep)
  expect_false(isTRUE(all.equal(fm$sm$parameters$beta_pooled,
                                f$sm$parameters$beta_pooled)))
})

test_that("the joint step-3 model drops incomplete cases and masks a missing outcome", {
  d   <- .mp_data()
  o   <- .mp_outcome(d)
  fit <- fit_mixture(d$X, n_classes = 2, measurement = "binary",
                     n_init = 5, random_state = 1)
  ok  <- !is.na(d$x)
  expect_message(
    f <- suppressWarnings(add_outcome(fit, o$y, predictors = data.frame(x = d$x),
                                      outcome_type = "continuous",
                                      correction = "ML", assignment = "modal")),
    "30 cases missing a predictor or covariate were left out of the step-3 model")
  y2 <- o$y; y2[!ok] <- 100
  g  <- suppressWarnings(suppressMessages(
    add_outcome(fit, y2, predictors = data.frame(x = d$x),
                outcome_type = "continuous", correction = "ML",
                assignment = "modal")))
  expect_equal(g$metrics$ll, f$metrics$ll, tolerance = 1e-8)
  expect_true(all(f$lower_bound[!ok] == 0))

  # A case whose outcome is missing contributes the class regression's term
  # alone, log sum_k P(k | x) P(a | k): its outcome is not filled in.
  r  <- exp(fit$log_resp); A <- get_modal_resp(r)
  C  <- sweep(t(r) %*% A, 1, colSums(r), "/")
  b  <- f$sm$models$predictor$parameters$beta
  m  <- o$miss_y
  eta <- cbind(1, d$x[m]) %*% t(b)
  p   <- exp(eta - apply(eta, 1, max)); p <- p / rowSums(p)
  hand <- log(rowSums(p * t(C[, max.col(A)[m], drop = FALSE])))
  expect_equal(f$lower_bound[m], hand, tolerance = 1e-8)

  fb <- suppressWarnings(suppressMessages(
    add_outcome(fit, o$y, predictors = data.frame(x = d$x),
                outcome_type = "continuous", correction = "BCH")))
  gb <- suppressWarnings(suppressMessages(
    add_outcome(fit, y2, predictors = data.frame(x = d$x),
                outcome_type = "continuous", correction = "BCH")))
  expect_identical(fb$struct_keep, ok)
  expect_equal(gb$sm$models$distal$parameters$means,
               fb$sm$models$distal$parameters$means, tolerance = 1e-8)
  expect_equal(gb$sm$models$predictor$parameters$beta,
               fb$sm$models$predictor$parameters$beta, tolerance = 1e-6)
})

test_that("a one-step outcome fit drops the cases missing a covariate whole", {
  d  <- .mp_data()
  o  <- .mp_outcome(d)
  ok <- !is.na(d$x)
  expect_message(
    f <- fit_mixture(d$X, n_classes = 2, measurement = "binary",
                     outcome = o$y, outcome_covariates = d$x,
                     outcome_type = "continuous", n_steps = 1, n_init = 5,
                     random_state = 1),
    "30 cases missing a covariate were left out of the model \\(n = 270 used\\)")
  g <- suppressMessages(fit_mixture(
    d$X[ok, ], n_classes = 2, measurement = "binary", outcome = o$y[ok],
    outcome_covariates = d$x[ok], outcome_type = "continuous", n_steps = 1,
    n_init = 5, random_state = 1))
  expect_equal(f$metrics$ll, g$metrics$ll, tolerance = 1e-8)
  expect_equal(nrow(f$log_resp), 270L)
})

test_that("fit_lta() drops a case missing a covariate from a one-step fit", {
  s  <- .lta_cov_refine_sim()
  z  <- s$Z$z
  z[c(3, 17, 40, 41, 99)] <- NA
  ok <- !is.na(z)
  args <- list(n_statuses = 2, times = 3, measurement = "binary", n_init = 3,
               random_state = 2, max_iter = 200, n_cores = 1,
               standard_errors = FALSE)
  expect_message(
    f <- suppressWarnings(do.call(fit_lta, c(list(
      s$X, predictors_initial = data.frame(z = z),
      predictors_transition = data.frame(z = z)), args))),
    "5 cases missing a covariate were left out of the model \\(n = 295 used\\)")
  g <- suppressWarnings(do.call(fit_lta, c(list(
    s$X[ok, ], predictors_initial = data.frame(z = z[ok]),
    predictors_transition = data.frame(z = z[ok])), args)))
  expect_equal(f$metrics$ll, g$metrics$ll, tolerance = 1e-8)
  expect_equal(nrow(f$data), 295L)

  # "mean" keeps every case, and so, for now, does the two-step estimator.
  fm <- suppressWarnings(do.call(fit_lta, c(list(
    s$X, predictors_transition = data.frame(z = z),
    missing_predictors = "mean"), args)))
  expect_equal(nrow(fm$data), 300L)
  f2 <- suppressWarnings(do.call(fit_lta, c(list(
    s$X, predictors_transition = data.frame(z = z), n_steps = 2), args)))
  expect_equal(nrow(f2$data), 300L)
})
