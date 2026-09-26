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
    n_steps = 3, correction = "none", n_init = 5, random_state = 1))

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
                      measurement = "continuous",
                      variances_equal = FALSE),
    count = list(X = matrix(rpois(n * 4, ifelse(rep(cls, 4) == 1, 4, 1)), n, 4),
                 measurement = "count")
  )
  specs$mixed <- list(X = cbind(specs$binary$X[, 1:3], specs$continuous$X[, 1:2]),
                      measurement = list(binary = 1:3, continuous = 4:5))

  for (nm in names(specs)) {
    sp  <- specs[[nm]]
    fit <- suppressMessages(fit_mixture(sp$X, n_classes = 2,
                                        measurement = sp$measurement,
                                        variances_equal = sp$variances_equal,
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

test_that("step-1 packing respects equal variances across classes", {
  set.seed(11)
  n   <- 300
  cls <- sample(1:2, n, TRUE)
  X   <- matrix(rnorm(n * 4, ifelse(rep(cls, 4) == 1, 1, -1)), n, 4)

  free <- suppressMessages(fit_mixture(X, n_classes = 2, measurement = "continuous",
                                       variances_equal = FALSE,
                                       n_init = 3, random_state = 4))
  inv  <- suppressMessages(fit_mixture(X, n_classes = 2, measurement = "continuous",
                                       variances_equal = TRUE,
                                       n_init = 3, random_state = 4))

  # K variance cells per item collapse to one, so the constrained packing is
  # shorter by (K - 1) * n_items.
  expect_equal(length(.step1_pack(free)) - length(.step1_pack(inv)),
               (inv$n_components - 1L) * 4L)

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

test_that("the Hessian-based standard errors invert the free block, never the padded matrix", {
  # m_step.covariate() stores the K*D Hessian with the anchor class's block
  # padded by a -1e8 diagonal. pinv()'s cutoff is relative to the largest
  # singular value, so inverting the padded matrix zeroes every direction of
  # real curvature below 1e8 * sqrt(eps) = 1.5 -- and reports a variance that
  # is too small, silently. A covariate on a small scale has exactly that kind
  # of curvature, which is what this fixture sets up: the same classes as
  # .cse_sim(), with z in units twenty times smaller.
  d <- .cse_sim()
  d$Z$z <- d$Z$z / 20
  one <- suppressMessages(fit_mixture(
    d$X, n_classes = 2, measurement = "binary", predictors = d$Z,
    n_steps = 1, n_init = 3, random_state = 1))

  H <- one$sm$parameters$hessian
  V <- vcov(one)
  K <- one$n_components
  D <- ncol(one$sm$parameters$beta)
  ref <- attr(V, "ref_class")
  idx <- as.vector(vapply(setdiff(seq_len(K), ref),
                          function(k) ((k - 1L) * D + 1L):(k * D), integer(D)))

  # The fixture has the property the test needs: the free block carries a
  # singular value under the padded cutoff, so the two computations differ.
  expect_lt(min(svd(-H[idx, idx])$d), 1.5)

  # solve() has no relative cutoff, so it is an independent oracle for the
  # free block's inverse.
  expect_equal(matrix(as.numeric(V), nrow(V)), solve(-H[idx, idx]),
               tolerance = 1e-8)

  # And the padded inverse is not what vcov() returns: at least one standard
  # error is understated by half or more the old way.
  old_se <- sqrt(pmax(diag(pinv(-H)[idx, idx]), 0))
  expect_gt(max(1 - old_se / sqrt(diag(V))), 0.5)

  # confint(), analytical_wald_test() and summary()'s table all read the
  # same covariance, so they must agree with vcov() on this fit.
  se <- sqrt(diag(V))
  ci <- suppressWarnings(confint(one, ref_class = ref))
  half <- unlist(lapply(names(ci), function(nm) {
    free <- setdiff(seq_len(K), ref)
    log(ci[[nm]]$Upper[free]) - log(ci[[nm]]$OR[free])
  }))
  expect_equal(sort(half), sort(unname(qnorm(0.975) * se)), tolerance = 1e-8)

  b <- coef(one, exponentiate = FALSE, ref_class = ref)
  free <- setdiff(seq_len(K), ref)
  w <- analytical_wald_test(one, "z", ref_class = ref)
  expect_equal(w$df, K - 1L)
  # analytical_wald_test() rounds the statistic it returns to three decimals.
  expect_equal(w$Wald_Chi2, round((b[free, "z"] / se[["Class 1:z"]])^2, 3))

  s <- capture.output(tab <- summary(one, ref_class = ref)$coefficients)
  expect_equal(sort(tab$se), sort(unname(se)), tolerance = 1e-8)
})

test_that("a fitted covariate model is always anchored on its last class", {
  # .fit_mnl() pins the last class while estimating, but the size sort at the
  # end of fit_mixture() permutes the rows and used to carry the anchor to
  # whatever rank its class's weight happened to land on. The re-anchoring is a
  # reparameterisation, so nothing a user reads from the fit may move.
  d <- .cse_sim()
  for (seed in 1:4) {
    fit <- suppressMessages(fit_mixture(
      d$X, n_classes = 3, measurement = "binary", predictors = d$Z,
      n_steps = 1, n_init = 3, random_state = seed))
    B <- fit$sm$parameters$beta
    expect_equal(which(rowSums(abs(B)) == 0), nrow(B), info = seed)
  }

  # Applied by hand to an arbitrary anchor, the recentring leaves every fitted
  # class probability where it was and gives vcov() the same free-block
  # standard errors as the original parameterisation implies.
  fit  <- suppressMessages(fit_mixture(
    d$X, n_classes = 3, measurement = "binary", predictors = d$Z,
    n_steps = 1, n_init = 3, random_state = 1))
  Zmat <- .covariate_design(fit$sm, prepare_covariates(d$Z))
  P_of <- function(sm) softmax_rows(Zmat %*% t(sm$parameters$beta))
  moved <- .recenter_covariate_beta(fit$sm, 1L)
  expect_equal(which(rowSums(abs(moved$parameters$beta)) == 0), 1L)
  expect_equal(P_of(moved), P_of(fit$sm), tolerance = 1e-12)
  back <- .recenter_covariate_beta(moved, 3L)
  expect_equal(back$parameters$beta, fit$sm$parameters$beta, tolerance = 1e-10)
  # Only the free block is compared: the anchor's padded -1e8 diagonal would
  # swamp a relative comparison of the whole matrix.
  fidx <- seq_len(2L * ncol(fit$sm$parameters$beta))
  expect_equal(back$parameters$hessian[fidx, fidx],
               fit$sm$parameters$hessian[fidx, fidx], tolerance = 1e-6)
})

test_that("a measurement model without an unconstrained packing falls back cleanly", {
  set.seed(5)
  n   <- 300
  cls <- sample(1:2, n, TRUE)
  Y   <- t(vapply(cls, function(k)
    rnorm(4, (if (k == 1) 1 else -1) + (0:3) * (if (k == 1) .4 else -.2)),
    numeric(4)))
  Z <- data.frame(z = rnorm(n))

  expect_warning(
    fit <- suppressMessages(fit_lcga(Y, n_classes = 2, times = 4,
                                     family = "gaussian", predictors = Z,
                                     n_steps = 3, correction = "ML",
                                     n_init = 3, random_state = 1)),
    "not available for a `lcga`")
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


# ------------------------------------------------------------------------------
# The estimates are recoverable on the log scale
# ------------------------------------------------------------------------------
#
# The printed output reports odds ratios, which is the scale these effects are
# published on. What these two check is that the log-scale quantities behind it
# can still be got at, at full precision, by anyone who needs them on that
# scale.

test_that("confint() returns full precision and still prints to three decimals", {
  fit <- .cse_fit(.cse_sim())
  ci  <- confint(fit)

  # The defect this replaced rounded inside the returned object, which put a
  # 0.001 floor under any comparison of these numbers.
  or <- ci$z$OR
  expect_false(isTRUE(all.equal(or, round(or, 3), tolerance = 0)))

  # And the console output is unchanged, because the rounding moved into the
  # print method rather than disappearing.
  out <- paste(capture.output(print(ci)), collapse = "\n")
  expect_match(out, "CONFIDENCE INTERVALS FOR ODDS RATIOS", fixed = TRUE)
  expect_match(out, sprintf("%7.3f", or[2]), fixed = TRUE)
})

test_that("vcov() gives the standard errors and coef() the log-scale estimates", {
  fit <- .cse_fit(.cse_sim())
  K   <- fit$n_components
  D   <- ncol(fit$sm$parameters$beta)

  V  <- vcov(fit)
  se <- sqrt(diag(V))
  expect_equal(dim(V), c((K - 1L) * D, (K - 1L) * D))
  expect_length(se, (K - 1L) * D)
  expect_true(all(is.finite(se)))
  expect_true(all(se > 0))
  expect_equal(attr(V, "method"), fit$sm$parameters$V_method)

  # The anchor class is found rather than assumed, because sort_model_classes()
  # reorders the classes after estimation and can put it in any row. The rows of
  # vcov() must be the classes that are *not* the anchor, in order, and must
  # line up with the coefficients taken against that same reference.
  ref  <- attr(V, "ref_class")
  free <- setdiff(seq_len(K), ref)
  expect_true(all(fit$sm$parameters$beta[ref, ] == 0))
  expect_equal(rownames(V),
               paste(rep(paste("Class", free), each = D),
                     colnames(fit$sm$parameters$beta), sep = ":"))
  b <- coef(fit, exponentiate = FALSE, ref_class = ref)
  expect_true(all(b[ref, ] == 0))

  # The alignment check with teeth: against the same reference, confint()'s
  # half-width on the log scale is z * se, so these standard errors must
  # reproduce the intervals confint() already prints. This only holds if the
  # rows of vcov() and the rows of coef() name the same coefficients.
  ci <- suppressWarnings(confint(fit, ref_class = ref))
  half <- unlist(lapply(names(ci), function(nm)
    log(ci[[nm]]$Upper[free]) - log(ci[[nm]]$OR[free])))
  expect_equal(sort(half), sort(unname(qnorm(0.975) * se)), tolerance = 1e-8)

  # exponentiate = FALSE is exactly the log of the default, and the default is
  # exactly today's behaviour.
  expect_equal(log(coef(fit)), coef(fit, exponentiate = FALSE),
               tolerance = 1e-12)
  expect_equal(coef(fit), exp(coef(fit, exponentiate = FALSE)),
               tolerance = 1e-12)
})
