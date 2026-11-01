# fit_gmm() with several outcomes measured at the same occasions (a
# multi-process growth mixture), each with its own growth curve, the growth
# factors' covariance shared by the classes, and free residual covariances
# between outcomes at an occasion and across occasions.

.mp_q <- function(e) suppressMessages(suppressWarnings(e))

# Two outcomes x four occasions, two classes; residuals correlated within an
# occasion (0.06) and the first outcome's across occasions 2-3 (0.03).
.mp_data <- function(n = 400, seed = 1, na = FALSE) {
  set.seed(seed)
  ts <- 0:3
  cl <- sample(1:2, n, TRUE, c(.4, .6))
  B  <- rbind(c(2.0, 0.30, 1.5, 0.20), c(1.0, -0.10, 1.0, 0.00))
  P  <- matrix(c(.10, .00, .04, .00,
                 .00, .01, .00, .00,
                 .04, .00, .08, .00,
                 .00, .00, .00, .01), 4)
  eta <- B[cl, ] + matrix(rnorm(n * 4), n) %*% chol(P)
  Th <- diag(.10, 8)
  for (t in 1:4) Th[2 * t - 1, 2 * t] <- Th[2 * t, 2 * t - 1] <- .06
  Th[3, 5] <- Th[5, 3] <- .03
  Y <- matrix(0, n, 8)
  for (t in 1:4) {
    Y[, 2 * t - 1] <- eta[, 1] + ts[t] * eta[, 2]
    Y[, 2 * t]     <- eta[, 3] + ts[t] * eta[, 4]
  }
  Y <- Y + matrix(rnorm(n * 8), n) %*% chol(Th)
  if (na) Y[sample(length(Y), round(.08 * length(Y)))] <- NA
  colnames(Y) <- paste0(c("dep_", "anx_"), rep(1:4, each = 2))
  Y
}
.mp_fit <- function(Y, ...) .mp_q(fit_gmm(Y, times = 4, n_classes = 2,
                                          n_init = 4, random_state = 1,
                                          n_cores = 1, ...))
.mp_pairs <- list("occasion", "dep@T2 with dep@T3")

test_that("several outcomes get one growth curve each, in time-major order", {
  Y <- .mp_data()
  f <- .mp_fit(Y)
  g <- f$growth
  expect_equal(g$n_items, 2L)
  expect_equal(g$items, c("dep", "anx"))
  expect_equal(g$factor_names, c("intercept_dep", "linear_dep",
                                 "intercept_anx", "linear_anx"))
  # Row (t, j) of the design is outcome j's polynomial at occasion t.
  expect_equal(unname(g$design[3, ]), c(1, 1, 0, 0))
  expect_equal(unname(g$design[4, ]), c(0, 0, 1, 1))
  expect_equal(g$random_names, g$factor_names)
  # 2 x 4 means, 10 Psi, 8 residual variances, 1 class logit.
  expect_equal(f$metrics$n_params, 8 + 10 + 8 + 1)
})

test_that("free residual covariances are estimated, counted and reported", {
  Y <- .mp_data()
  f0 <- .mp_fit(Y)
  f1 <- .mp_fit(Y, residual_covariances = .mp_pairs)
  expect_equal(f1$metrics$n_params, f0$metrics$n_params + 5L)
  expect_gt(f1$metrics$ll, f0$metrics$ll)
  rc <- f1$growth$residual_covariance
  expect_equal(colnames(rc), c("dep@T1 with anx@T1", "dep@T2 with anx@T2",
                               "dep@T3 with anx@T3", "dep@T4 with anx@T4",
                               "dep@T2 with dep@T3"))
  expect_equal(nrow(unique(rc)), 1L)          # held equal across classes
  expect_lt(max(abs(rc[1, 1:4] - .06)), .025)
  expect_lt(abs(rc[1, 5] - .03), .025)
  ge <- growth_estimates(f1)
  r  <- ge[ge$parameter == "residual_covariance", ]
  expect_equal(nrow(r), 10L)
  expect_true(all(is.finite(r$se) & r$se > 0))
  out <- capture.output(print(f1))
  expect_true(any(grepl("RESIDUAL COVARIANCE", out)))
  expect_true(any(grepl("Outcomes           : dep, anx", out, fixed = TRUE)))
  pdf(NULL); on.exit(dev.off())
  expect_silent(plot(f1))
})

test_that("missing answers and named random effects leave one-outcome fits alone", {
  Y <- .mp_data(na = TRUE)
  f <- .mp_fit(Y, residual_covariances = .mp_pairs,
               random_effects = c("intercept_dep", "linear_dep", "intercept_anx"))
  expect_equal(f$growth$random_names,
               c("intercept_dep", "linear_dep", "intercept_anx"))
  expect_equal(f$metrics$n_params, 8 + 6 + 8 + 5 + 1)
  expect_true(f$converged)
  # A single outcome given by keyword is fitted as before.
  y1 <- Y[, c(1, 3, 5, 7)]
  a <- .mp_fit(y1)
  b <- .mp_fit(y1, random_effects = c("intercept", "linear"))
  expect_equal(a$metrics$ll, b$metrics$ll, tolerance = 1e-10)
  expect_null(a$growth$residual_covariance)
  # ... and takes a residual covariance across occasions.
  d <- .mp_fit(y1, residual_covariances = "dep@T2 with dep@T3")
  expect_equal(d$metrics$n_params, a$metrics$n_params + 1L)
})

test_that("predict() and blrt() carry the residual covariances", {
  Y <- .mp_data()
  f <- .mp_fit(Y, residual_covariances = .mp_pairs)
  p <- predict(f, newdata = Y[1:20, ], type = "posterior")
  expect_equal(unname(as.matrix(p)), unname(exp(f$log_resp[1:20, ])),
               tolerance = 1e-8)
  spec <- .blrt_growth_spec(f)
  expect_equal(spec$args$theta_pairs, f$growth$residual_pairs)
})

test_that("unsupported combinations are refused with a reason", {
  Y <- .mp_data()
  expect_error(.mp_fit(Y, family = "binomial"), "one repeated outcome")
  expect_error(.mp_fit(Y, psi = "free"), "held equal across classes")
  expect_error(.mp_fit(Y, residual_covariances = .mp_pairs,
                       residual_equal = FALSE), "held equal across classes")
  expect_error(.mp_fit(Y, growth_predictors = rnorm(nrow(Y))),
               "not available yet")
  expect_error(.mp_fit(Y, residual_covariances = "dep@T1 with nope@T1"),
               "two of the data's columns")
  expect_error(.mp_fit(Y[, c(1, 3, 5, 7)], residual_covariances = "occasion"),
               "several outcomes")
  expect_error(.mp_fit(Y, residual_covariances = .mp_pairs,
                       residual = "constant"), "residual variance per occasion")
  expect_error(.mp_fit(Y, random_effects = c("intercept_dep", "slope")),
               "not a growth factor")
})
