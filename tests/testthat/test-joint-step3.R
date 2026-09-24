# add_outcome(predictors = ): class predictors and distal outcomes in one
# bias-corrected (ML or BCH) step-3 model.

.joint_sim <- function(n = 600, seed = 11) {
  set.seed(seed)
  z   <- rbinom(n, 1, 0.5)
  cls <- 1L + rbinom(n, 1, plogis(-0.4 + 1.0 * z))
  p   <- rbind(rep(0.85, 6), rep(0.15, 6))
  X   <- t(sapply(cls, function(k) rbinom(6, 1, p[k, ])))
  y1  <- c(0, 2)[cls] + 0.5 * z + rnorm(n, sd = c(1, 1.6)[cls])
  y2  <- rbinom(n, 1, plogis(c(-1, 0.5)[cls] - 0.7 * z))
  list(X = X, d = data.frame(z = z, y1 = y1, y2 = y2))
}

.joint_fit <- function(s, correction = "ML", ...) {
  fit <- fit_mixture(s$X, n_classes = 2, measurement = "binary", n_init = 5,
                     random_state = 1, n_cores = 1)
  suppressMessages(add_outcome(fit, outcome = s$d[c("y1", "y2")],
                               covariates = s$d["z"], predictors = s$d["z"],
                               correction = correction, ...))
}

test_that("the joint model counts its parameters and reproduces its own log-likelihood", {
  s <- .joint_sim()
  j <- .joint_fit(s, variances = "class_specific", assignment = "modal",
                  se = "hessian")
  # (K-1)*2 class regression + (K intercepts + 1 slope + K variances) + (K + 1)
  expect_equal(n_parameters(j$sm), 10L)

  # The step-3 log-likelihood written out by hand from the fitted parameters.
  fit <- fit_mixture(s$X, n_classes = 2, measurement = "binary", n_init = 5,
                     random_state = 1, n_cores = 1)
  r   <- exp(fit$log_resp)
  a   <- max.col(r, ties.method = "first")
  A   <- diag(2)[a, ]
  C   <- t(r) %*% A
  C   <- C / rowSums(C)
  m   <- j$sm$models
  eta <- cbind(1, s$d$z) %*% t(m$predictor$parameters$beta)
  pk  <- exp(eta) / rowSums(exp(eta))
  b1  <- m$distal$parameters$beta_pooled
  v1  <- as.vector(m$distal$parameters$covariances)
  b2  <- m$distal2$parameters$beta_pooled
  f   <- sapply(1:2, function(k) {
    dnorm(s$d$y1, b1[k] + b1[3] * s$d$z, sqrt(v1[k])) *
      dbinom(s$d$y2, 1, plogis(b2[k] + b2[3] * s$d$z))
  })
  ll <- sum(log(rowSums(pk * f * t(C[, a]))))
  expect_equal(j$metrics$ll, ll, tolerance = 1e-8)
  expect_true(v1[1] != v1[2])

  se <- sqrt(diag(j$sm$parameters$vcov_joint))
  expect_length(se, 10L)
  expect_true(all(is.finite(se) & se > 0))
})

test_that("the joint BCH model maximises its weighted log-likelihood and reports the case sandwich", {
  s  <- .joint_sim()
  jh <- .joint_fit(s, correction = "BCH", variances = "class_specific",
                   assignment = "modal", se = "hessian")
  jr <- .joint_fit(s, correction = "BCH", variances = "class_specific",
                   assignment = "modal", se = "robust")
  expect_equal(n_parameters(jh$sm), 10L)

  # The BCH weights and the weighted log-likelihood written out by hand.
  fit <- fit_mixture(s$X, n_classes = 2, measurement = "binary", n_init = 5,
                     random_state = 1, n_cores = 1)
  r   <- exp(fit$log_resp)
  A   <- diag(2)[max.col(r, ties.method = "first"), ]
  C   <- t(A) %*% r
  C   <- sweep(C, 2, colSums(r), "/")
  dw  <- A %*% t(solve(C))
  m   <- jh$sm$models
  eta <- cbind(1, s$d$z) %*% t(m$predictor$parameters$beta)
  lpk <- eta - log(rowSums(exp(eta)))
  b1  <- m$distal$parameters$beta_pooled
  v1  <- as.vector(m$distal$parameters$covariances)
  b2  <- m$distal2$parameters$beta_pooled
  lf  <- sapply(1:2, function(k) {
    dnorm(s$d$y1, b1[k] + b1[3] * s$d$z, sqrt(v1[k]), log = TRUE) +
      dbinom(s$d$y2, 1, plogis(b2[k] + b2[3] * s$d$z), log = TRUE)
  })
  expect_equal(jh$metrics$ll, sum(dw * (lpk + lf)), tolerance = 1e-8)

  # `se` does not change the BCH variance: it is always the sandwich.
  expect_equal(jh$sm$parameters$vcov_joint, jr$sm$parameters$vcov_joint)
  expect_match(jh$sm$parameters$V_method, "BCH")
  se <- sqrt(diag(jh$sm$parameters$vcov_joint))
  expect_true(all(is.finite(se) & se > 0))
})

test_that("each outcome can be contrasted and summarised", {
  s <- .joint_sim()
  j <- .joint_fit(s)
  cy1 <- outcome_contrasts(j)
  cy2 <- outcome_contrasts(j, outcome = "y2")
  expect_equal(attr(cy1, "outcome"), "y1")
  expect_equal(attr(cy2, "outcome"), "y2")
  expect_equal(attr(cy2, "outcome_type"), "categorical")
  expect_true(all(is.finite(c(cy1$se, cy2$se))))
  expect_error(outcome_contrasts(j, outcome = "nope"), "one of this fit")

  txt <- capture.output(sm <- summary(j))
  expect_true(any(grepl("Distal outcome 2 of 2: y2", txt)))
  expect_named(sm$outcomes, c("y1", "y2"))
})

test_that("the sandwich and the Hessian share the estimates", {
  s  <- .joint_sim()
  jh <- .joint_fit(s, se = "hessian")
  jr <- .joint_fit(s, se = "robust")
  expect_equal(jh$metrics$ll, jr$metrics$ll)
  expect_false(isTRUE(all.equal(jh$sm$parameters$vcov_joint,
                                jr$sm$parameters$vcov_joint)))
})

test_that("unsupported combinations are refused", {
  s   <- .joint_sim(n = 200)
  fit <- fit_mixture(s$X, n_classes = 2, measurement = "binary", n_init = 2,
                     random_state = 1, n_cores = 1)
  expect_error(add_outcome(fit, s$d$y1, predictors = s$d$z,
                           correction = "none"), "not yet available")
  expect_error(add_outcome(fit, s$d$y1, predictors = s$d$z, steps = 2),
               "three-step")
  expect_error(suppressMessages(
    add_outcome(fit, s$d$y1, covariates = s$d["z"], predictors = s$d$z,
                slopes = "class_specific", correction = "ML")),
    "pooled")
  expect_error(suppressMessages(
    add_outcome(fit, s$d$y1, variances = "class_specific")),
    "class_specific")
})

test_that("variances = \"equal\" is the one-variance model", {
  s   <- .joint_sim(n = 300)
  fit <- fit_mixture(s$X, n_classes = 2, measurement = "binary", n_init = 2,
                     random_state = 1, n_cores = 1)
  o <- suppressMessages(add_outcome(fit, s$d$y1, covariates = s$d["z"]))
  v <- o$sm$parameters$covariances
  expect_equal(v[1], v[2])
  expect_equal(n_parameters(o$sm), ncol(o$sm$parameters$beta_pooled) + 1L)
})
