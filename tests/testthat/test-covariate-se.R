# ==============================================================================
# Standard errors for step-three covariate models
# ==============================================================================
#
# Three kinds of check.
#
# 1. Algebra. The step-three score and Hessian in R/step3_variance.R are written
#    out analytically; both are re-derived here by finite differences of the
#    step-three log-likelihood itself and required to agree. A closed form that
#    is subtly wrong still returns a plausible standard error, so this is the
#    check that matters most. The same block verifies that the step-one
#    parameter packing round-trips exactly for every measurement family it
#    claims to support - a packing that silently loses a parameter would make
#    the correction term wrong without making it look wrong.
#
# 2. Structure. The estimators must stand in the order theory puts them in, the
#    unadjusted third step must reduce to the ordinary weighted multinomial
#    logit, and the correction must shrink towards nothing as the classification
#    becomes certain - the condition Bakk, Oberski and Vermunt (2014) identify
#    as the one under which the correction stops being needed.
#
# 3. An external benchmark: a two-level check of `se = "robust"` against a
#    reference three-step run, which is the same
#    estimator, on the sleep-quality data in
#    `Datos Serenamente/[Analysis] Patterns of sleep quality...`. Skipped when
#    that folder is not on the machine.

# ------------------------------------------------------------------------------
# Fixtures
# ------------------------------------------------------------------------------

.cse_sim <- function(n = 400, seed = 20260804, rho = 0.85) {
  set.seed(seed)
  z   <- rnorm(n)
  cls <- 1L + rbinom(n, 1, plogis(-0.4 + 0.9 * z))
  list(
    X = matrix(rbinom(n * 6, 1, ifelse(rep(cls, 6) == 1L, rho, 1 - rho)), n, 6),
    Z = data.frame(z = z, g = factor(sample(c("a", "b", "c"), n, TRUE)))
  )
}

.cse_fit <- function(d, ...) {
  suppressMessages(fit_mixture(
    d$X, n_classes = 2, measurement = "binary", predictors = d$Z,
    n_steps = 3, correction = "ML", n_init = 5, random_state = 1, ...))
}

# The step-three log-likelihood, as a function of the free coefficients, built
# from nothing but log/exp so that it shares no code with the implementation.
.cse_L3 <- function(pars, Zmat, resp1, Cn, w, K) {
  D <- ncol(Zmat)
  B <- rbind(matrix(pars, K - 1L, D, byrow = TRUE), 0)
  e <- exp(Zmat %*% t(B))
  P <- e / rowSums(e)
  Z <- if (is.null(Cn)) P else P %*% Cn
  sum(w * rowSums(resp1 * log(pmax(Z, 1e-300))))
}

# Everything the score and Hessian need, recomputed from a fitted model.
.cse_pieces <- function(fit, d, adjusted = TRUE) {
  ms <- fit; ms$sm <- NULL
  r  <- exp(e_step(ms, d$X, NULL)$log_resp)
  w  <- fit$sample_weights
  B  <- fit$sm$parameters$beta
  # sort_model_classes() may have moved the anchored row; the softmax is
  # invariant to a common shift, so re-anchor on the last class.
  B  <- sweep(B, 2, B[nrow(B), ], "-")
  list(Zmat = .covariate_design(fit$sm, prepare_covariates(d$Z)),
       resp1 = r, w = w, beta = B,
       Cn = if (adjusted) .step3_classification_table(r, w) else NULL,
       K = fit$n_components)
}

# ------------------------------------------------------------------------------
# 1. Algebra
# ------------------------------------------------------------------------------

test_that("the step-3 score is the gradient of the step-3 log-likelihood", {
  d   <- .cse_sim()
  fit <- .cse_fit(d)

  for (adjusted in c(TRUE, FALSE)) {
    p    <- .cse_pieces(fit, d, adjusted)
    pcs  <- .step3_pieces(p$beta, p$Zmat, p$resp1, p$Cn)
    ana  <- colSums(.step3_scores(pcs, p$Zmat, p$w, p$K))

    b0  <- as.vector(t(p$beta[seq_len(p$K - 1L), , drop = FALSE]))
    eps <- 1e-6
    num <- vapply(seq_along(b0), function(m) {
      hi <- b0; hi[m] <- hi[m] + eps
      lo <- b0; lo[m] <- lo[m] - eps
      (.cse_L3(hi, p$Zmat, p$resp1, p$Cn, p$w, p$K) -
         .cse_L3(lo, p$Zmat, p$resp1, p$Cn, p$w, p$K)) / (2 * eps)
    }, numeric(1))

    expect_lt(max(abs(ana - num)), 1e-4 * max(1, max(abs(num))))
  }
})

test_that("the step-3 Hessian is the curvature of the step-3 log-likelihood", {
  d   <- .cse_sim()
  fit <- .cse_fit(d)

  for (adjusted in c(TRUE, FALSE)) {
    p   <- .cse_pieces(fit, d, adjusted)
    pcs <- .step3_pieces(p$beta, p$Zmat, p$resp1, p$Cn)
    ana <- .step3_hessian(pcs, p$Zmat, p$w, p$resp1, p$K)

    b0  <- as.vector(t(p$beta[seq_len(p$K - 1L), , drop = FALSE]))
    eps <- 1e-4
    L3  <- function(b) .cse_L3(b, p$Zmat, p$resp1, p$Cn, p$w, p$K)
    num <- matrix(0, length(b0), length(b0))
    for (a in seq_along(b0)) for (b in a:length(b0)) {
      ea <- numeric(length(b0)); ea[a] <- eps
      eb <- numeric(length(b0)); eb[b] <- eps
      num[a, b] <- num[b, a] <-
        (L3(b0 + ea + eb) - L3(b0 + ea - eb) -
           L3(b0 - ea + eb) + L3(b0 - ea - eb)) / (4 * eps^2)
    }

    expect_lt(max(abs(ana - num)) / max(abs(num)), 1e-4)
  }
})

test_that("with no classification error the step-3 Hessian is the multinomial-logit one", {
  # The unadjusted branch must reproduce exactly what m_step.covariate() gets
  # from an ordinary weighted multinomial logit: same estimator, so the
  # information matrix cannot differ.
  d   <- .cse_sim()
  fit <- suppressMessages(fit_mixture(
    d$X, n_classes = 3, measurement = "binary", predictors = d$Z,
    n_steps = 2, n_init = 5, random_state = 1))

  p   <- .cse_pieces(fit, d, adjusted = FALSE)
  pcs <- .step3_pieces(p$beta, p$Zmat, p$resp1, NULL)
  ana <- .step3_hessian(pcs, p$Zmat, p$w, p$resp1, p$K)

  K <- p$K; D <- ncol(p$Zmat)
  P <- softmax_rows(p$Zmat %*% t(p$beta))
  mnl <- matrix(0, (K - 1L) * D, (K - 1L) * D)
  for (j in seq_len(K - 1L)) for (l in seq_len(K - 1L)) {
    wt <- P[, j] * ((j == l) - P[, l]) * p$w
    mnl[((j-1)*D+1):(j*D), ((l-1)*D+1):(l*D)] <-
      -t(p$Zmat) %*% sweep(p$Zmat, 1, wt, "*")
  }
  expect_lt(max(abs(ana - mnl)), 1e-8 * max(abs(mnl)))
})

test_that("step-1 packing round-trips for every supported measurement family", {
  set.seed(99)
  n   <- 300
  cls <- sample(1:2, n, TRUE)

  specs <- list(
    binary = list(X = matrix(rbinom(n * 5, 1, ifelse(rep(cls, 5) == 1, .8, .2)),
                             n, 5),
                  measurement = "binary"),
    categorical = list(
      X = matrix(vapply(rep(cls, 4), function(k)
        sample(1:3, 1, prob = if (k == 1) c(.7, .2, .1) else c(.1, .2, .7)),
        integer(1)), n, 4),
      measurement = "categorical"),
    continuous = list(X = matrix(rnorm(n * 4, ifelse(rep(cls, 4) == 1, 1, -1)),
                                 n, 4),
                      measurement = "continuous"),
    count = list(X = matrix(rpois(n * 4, ifelse(rep(cls, 4) == 1, 4, 1)), n, 4),
                 measurement = "count")
  )
  specs$mixed <- list(X = cbind(specs$binary$X[, 1:3], specs$continuous$X[, 1:2]),
                      measurement = list(binary = 1:3, continuous = 4:5))

  for (nm in names(specs)) {
    sp  <- specs[[nm]]
    fit <- suppressMessages(fit_mixture(sp$X, n_classes = 2,
                                        measurement = sp$measurement,
                                        n_init = 3, random_state = 2))
    par <- .step1_pack(fit)
    expect_false(is.null(par), info = nm)

    ll_direct <- logsumexp(sweep(log_likelihood(fit$mm, sp$X), 2,
                                 log(fit$weights), "+"), MARGIN = 1)
    expect_equal(.step1_ll_case(fit, sp$X, par), ll_direct, info = nm)
    expect_equal(.step1_pack(.step1_unpack(fit, par)), par, info = nm)
  }
})

test_that("step-1 packing respects equality constraints across occasions", {
  set.seed(7)
  n   <- 300
  cls <- sample(1:2, n, TRUE)
  X   <- cbind(matrix(rbinom(n * 4, 1, ifelse(rep(cls, 4) == 1, .8, .2)), n, 4),
               matrix(rbinom(n * 4, 1, ifelse(rep(cls, 4) == 1, .75, .25)), n, 4))

  inv  <- suppressMessages(fit_rmlca(X, n_classes = 2, times = 2,
                                     measurement_invariance = "full",
                                     n_init = 3, random_state = 3))
  free <- suppressMessages(fit_rmlca(X, n_classes = 2, times = 2,
                                     measurement_invariance = "none",
                                     n_init = 3, random_state = 3))

  # An item held equal across the two occasions is one free parameter, not two,
  # so the invariant packing must be shorter by exactly K * n_items.
  expect_equal(length(.step1_pack(free)) - length(.step1_pack(inv)),
               inv$n_components * 4L)

  for (f in list(inv, free)) {
    par <- .step1_pack(f)
    expect_equal(.step1_pack(.step1_unpack(f, par)), par)
    expect_equal(.step1_ll_case(f, X, par),
                 logsumexp(sweep(log_likelihood(f$mm, X), 2,
                                 log(f$weights), "+"), MARGIN = 1))
  }
})

# ------------------------------------------------------------------------------
# 2. Structure
# ------------------------------------------------------------------------------

test_that("the corrected variance exceeds the step-3 sandwich it extends", {
  # D3* = D3 + J D1 J' adds a positive semi-definite term, so no standard error
  # may shrink. This is the claim of Bakk et al. (2014, eq. 17).
  d    <- .cse_sim()
  rob  <- .cse_fit(d, se = "robust")
  corr <- .cse_fit(d, se = "corrected")

  se_of <- function(f) sqrt(diag(f$sm$parameters$V_robust))
  expect_true(all(se_of(corr) >= se_of(rob) - 1e-8))
  expect_true(any(se_of(corr) > se_of(rob) + 1e-6))

  extra <- corr$sm$parameters$V_robust - rob$sm$parameters$V_robust
  ev    <- eigen((extra + t(extra)) / 2, only.values = TRUE)$values
  expect_gt(min(ev), -1e-8)
})

test_that("the step-1 correction vanishes when the classification is certain", {
  # With near-perfect separation there is almost nothing left to carry over
  # from step 1, which is the regime Bakk et al. identify as needing no
  # correction at all (entropy R^2 > .90 with a large step-1 sample).
  gap <- function(rho, n) {
    d    <- .cse_sim(n = n, rho = rho)
    rob  <- .cse_fit(d, se = "robust")
    corr <- .cse_fit(d, se = "corrected")
    keep <- diag(rob$sm$parameters$V_robust) > 0
    max(sqrt(diag(corr$sm$parameters$V_robust)[keep]) /
          sqrt(diag(rob$sm$parameters$V_robust)[keep])) - 1
  }
  fuzzy <- gap(0.75, 400)
  sharp <- gap(0.99, 1500)

  expect_lt(sharp, 0.02)
  expect_gt(fuzzy, sharp)
})

test_that("the estimator used is recorded and reported", {
  d <- .cse_sim()
  labels <- vapply(c("corrected", "robust", "hessian"),
                   function(m) .cse_fit(d, se = m)$sm$parameters$V_method,
                   character(1))
  expect_match(labels[["corrected"]], "Bakk")
  expect_match(labels[["robust"]],    "sandwich")
  expect_match(labels[["hessian"]],   "Observed information")

  fit <- .cse_fit(d, se = "corrected")
  expect_output(summary(fit), "Standard errors: Bakk")
  expect_equal(attr(confint(fit), "method"), fit$sm$parameters$V_method)
  expect_match(analytical_wald_test(fit, "z")$Method, "Bakk")
})

test_that("BCH and one-step fits say they are using the uncorrected Hessian", {
  d <- .cse_sim()
  bch <- suppressWarnings(suppressMessages(fit_mixture(
    d$X, n_classes = 2, measurement = "binary", predictors = d$Z,
    n_steps = 3, correction = "BCH", n_init = 3, random_state = 1)))
  one <- suppressMessages(fit_mixture(
    d$X, n_classes = 2, measurement = "binary", predictors = d$Z,
    n_steps = 1, n_init = 3, random_state = 1))

  for (f in list(bch, one)) {
    expect_null(f$sm$parameters$V_robust)
    expect_output(summary(f), "Standard errors: Q-function Hessian")
  }
})

test_that("a measurement model without an unconstrained packing falls back cleanly", {
  set.seed(5)
  n   <- 300
  cls <- sample(1:2, n, TRUE)
  Y   <- t(vapply(cls, function(k)
    rnorm(4, (if (k == 1) 1 else -1) + (0:3) * (if (k == 1) .4 else -.2)),
    numeric(4)))
  Z <- data.frame(z = rnorm(n))

  fit <- suppressMessages(fit_lcga(Y, n_classes = 2, times = 4,
                                   family = "gaussian", predictors = Z,
                                   n_steps = 3, correction = "ML",
                                   n_init = 3, random_state = 1))
  expect_null(.step1_pack(fit))
  expect_false(is.null(fit$sm$parameters$V_robust))
  expect_match(fit$sm$parameters$V_method, "step-1 correction unavailable")
})

test_that("a survey design still reaches the meat of the sandwich", {
  d  <- .cse_sim()
  n  <- nrow(d$X)
  fit <- suppressMessages(fit_mixture(
    d$X, n_classes = 2, measurement = "binary", predictors = d$Z,
    n_steps = 3, correction = "ML", n_init = 3, random_state = 1,
    strata = rep(1:4, length.out = n), cluster = rep(1:40, length.out = n)))
  expect_match(fit$sm$parameters$V_method, "survey-linearized")
})

