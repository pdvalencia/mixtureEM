# ==============================================================================
# outcome_contrasts() -- which classes differ on the distal outcome
# ==============================================================================
#
# The contrasts are only as good as the covariance they are formed from, and an
# indexing error there is silent: every number still looks plausible. What
# pins it down is that the joint Wald over the K-1 reference contrasts must
# reproduce the omnibus statistic the package already reports, which it can only
# do if the contrast matrix and the covariance agree about which entry is which
# class. That identity is asserted for both outcome types.

.sim_outcome_fit <- function(seed = 7, n = 400, outcome = c("continuous",
                                                            "categorical")) {
  outcome <- match.arg(outcome)
  set.seed(seed)
  z     <- rbinom(n, 1, 0.5)
  items <- matrix(rbinom(n * 5, 1, ifelse(z == 1, 0.85, 0.15)), nrow = n)
  y     <- if (outcome == "continuous") rnorm(n, mean = 10 + 3 * z)
           else factor(ifelse(runif(n) < ifelse(z == 1, 0.7, 0.2), "hi", "lo"))
  fit <- fit_mixture(items, n_classes = 3, measurement = "binary", n_init = 10)
  suppressMessages(add_outcome(fit, y))
}

# ------------------------------------------------------------------------------
# Continuous outcome
# ------------------------------------------------------------------------------

test_that("all pairs are reported once, as class minus reference", {
  fo <- .sim_outcome_fit()
  oc <- outcome_contrasts(fo)

  expect_s3_class(oc, "outcome_contrasts")
  expect_identical(nrow(oc), 3L)                       # choose(3, 2)
  expect_true(all(oc$class > oc$reference))
  expect_true(all(is.na(oc$category)))

  mu <- as.vector(.distal_submodel(fo$sm, "distal_continuous")$parameters$means)
  expect_equal(oc$estimate,
               c(mu[2] - mu[1], mu[3] - mu[1], mu[3] - mu[2]),
               tolerance = 1e-12)

  # Naming a reference flips the sign of the contrast it appears in.
  r2 <- outcome_contrasts(fo, ref = 2)
  expect_identical(nrow(r2), 2L)
  expect_true(all(r2$reference == 2L))
  expect_equal(r2$estimate[r2$class == 1L], -oc$estimate[oc$class == 2L],
               tolerance = 1e-12)
})

test_that("the standard error carries the covariance between the class means", {
  fo <- .sim_outcome_fit()
  S  <- .distal_submodel(fo$sm, "distal_continuous")$parameters$Sigma_mu
  oc <- outcome_contrasts(fo, ref = 1)

  expect_equal(oc$se[oc$class == 2L],
               sqrt(S[2, 2] + S[1, 1] - 2 * S[1, 2]), tolerance = 1e-12)

  # And is not the naive root-sum-of-squares, which is the mistake the function
  # exists to prevent. The two agree only when the covariance is zero.
  expect_false(isTRUE(all.equal(oc$se[oc$class == 2L],
                                sqrt(S[2, 2] + S[1, 1]))))
  expect_identical(attr(oc, "method"),
                   "sandwich covariance of the class means")
})

test_that("the reference contrasts reproduce the omnibus Wald statistic", {
  fo <- .sim_outcome_fit()
  cs <- .distal_submodel(fo$sm, "distal_continuous")
  mu <- as.vector(cs$parameters$means)
  S  <- cs$parameters$Sigma_mu
  K  <- fo$n_components

  R    <- cbind(-1, diag(K - 1L))
  Wref <- as.numeric(t(R %*% mu) %*% solve(R %*% S %*% t(R)) %*% (R %*% mu))
  omni <- .wald_omnibus_means(mu, as.vector(cs$parameters$ses), K,
                              Sigma_mu = S)

  expect_equal(Wref, omni$stat, tolerance = 1e-8)
})

test_that("summary() prints and returns the contrasts for a continuous outcome", {
  fo  <- .sim_outcome_fit()
  out <- utils::capture.output(s <- summary(fo))

  expect_true(any(grepl("Pairwise class differences", out)))
  expect_true(any(grepl("Class 2 vs 1", out)))
  expect_s3_class(s$outcome$contrasts, "outcome_contrasts")
  expect_equal(s$outcome$contrasts$estimate, outcome_contrasts(fo)$estimate,
               tolerance = 1e-12)
})

# ------------------------------------------------------------------------------
# Categorical outcome
# ------------------------------------------------------------------------------

test_that("a categorical outcome contrasts log odds and reports odds ratios", {
  fo <- .sim_outcome_fit(seed = 8, n = 500, outcome = "categorical")
  oc <- outcome_contrasts(fo)

  expect_true(all(oc$category == 2L))          # one non-reference category
  expect_true(all(c("OR", "OR_lower", "OR_upper") %in% names(oc)))
  expect_equal(oc$OR, exp(oc$estimate), tolerance = 1e-12)

  # summary()'s own reference table is the same contrast, so the two must not
  # disagree; what outcome_contrasts() adds is the pair it never showed.
  s  <- summary(fo)
  r1 <- outcome_contrasts(fo, ref = 1)
  expect_equal(r1$estimate, s$outcome$odds_ratios$estimate, tolerance = 1e-12)
  expect_equal(r1$se,       s$outcome$odds_ratios$se,       tolerance = 1e-12)
  expect_true(any(oc$class == 3L & oc$reference == 2L))
})

test_that("the categorical reference contrasts reproduce their omnibus too", {
  fo <- .sim_outcome_fit(seed = 8, n = 500, outcome = "categorical")
  ps <- .distal_submodel(fo$sm, "distal_pooled")
  b  <- ps$parameters$beta_pooled
  K  <- fo$n_components
  L  <- ncol(b)
  Mm <- nrow(b)

  V <- pinv(-ps$parameters$hessian)
  R <- matrix(0, Mm * (K - 1L), Mm * L)
  r <- 1L
  for (m in seq_len(Mm)) for (k in setdiff(seq_len(K), 1L)) {
    R[r, (m - 1L) * L + k]  <-  1
    R[r, (m - 1L) * L + 1L] <- -1
    r <- r + 1L
  }
  rv <- R %*% as.vector(t(b))
  W  <- as.numeric(t(rv) %*% pinv(R %*% V %*% t(R)) %*% rv)

  omni <- .wald_omnibus_pooled(b, ps$parameters$hessian, K, L - K, 1L)
  expect_equal(W, omni$stat, tolerance = 1e-8)
})

# ------------------------------------------------------------------------------
# Arguments and refusals
# ------------------------------------------------------------------------------

test_that("the multiplicity adjustment is applied across the reported set", {
  fo <- .sim_outcome_fit()
  oc <- outcome_contrasts(fo, adjust = "holm")

  expect_true("p_adj" %in% names(oc))
  expect_equal(oc$p_adj, stats::p.adjust(oc$p, "holm"), tolerance = 1e-12)
  expect_true(all(oc$p_adj >= oc$p - 1e-12))
  expect_false("p_adj" %in% names(outcome_contrasts(fo)))
})

test_that("the confidence level widens the interval", {
  fo   <- .sim_outcome_fit()
  wide <- outcome_contrasts(fo, level = 0.99)
  narr <- outcome_contrasts(fo, level = 0.95)
  expect_true(all(wide$lower < narr$lower))
  expect_true(all(wide$upper > narr$upper))
})

test_that("a fit with no outcome, and a bad reference, are refused", {
  set.seed(9)
  X   <- matrix(rbinom(600, 1, 0.4), ncol = 3)
  fit <- fit_mixture(X, n_classes = 2, measurement = "binary", n_init = 2)

  expect_error(outcome_contrasts(fit), "No distal outcome")

  fo <- .sim_outcome_fit()
  expect_error(outcome_contrasts(fo, ref = 0), "between 1 and 3")
  expect_error(outcome_contrasts(fo, ref = 4), "between 1 and 3")
  expect_error(outcome_contrasts(fo, level = 1), "between 0 and 1")
})

test_that("a class-specific outcome model is refused, and says where to go", {
  set.seed(10)
  n     <- 300
  z     <- rbinom(n, 1, 0.5)
  items <- matrix(rbinom(n * 5, 1, ifelse(z == 1, 0.85, 0.15)), nrow = n)
  y     <- rnorm(n, mean = 10 + 3 * z)
  cov1  <- rnorm(n)
  fit   <- fit_mixture(items, n_classes = 2, measurement = "binary", n_init = 5)
  fo    <- suppressMessages(add_outcome(fit, y, covariates = cbind(cov1 = cov1),
                                        slopes = "class_specific"))

  expect_error(outcome_contrasts(fo), "no covariance between the blocks")
  expect_error(outcome_contrasts(fo), "bootstrap_covariates")
})
