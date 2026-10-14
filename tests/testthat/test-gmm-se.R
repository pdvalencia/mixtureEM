# ==============================================================================
# Standard errors for growth mixture parameters (growth_estimates())
# ==============================================================================
#
# The oracle is the mixture-of-normals log-likelihood written out by hand on
# the reported scale and differentiated by stats::optimHess(), which shares
# nothing with the package's packing, scores or delta map.

.gmm_se_data <- function(n = 400, seed = 1) {
  set.seed(seed)
  cls   <- rbinom(n, 1, 0.5)
  icept <- ifelse(cls == 1, 3, 1) + rnorm(n, 0, 1)
  slope <- ifelse(cls == 1, 1.0, 0.2) + rnorm(n, 0, 0.4)
  outer(icept, rep(1, 4)) + outer(slope, 0:3) + matrix(rnorm(n * 4, 0, 0.7), n, 4)
}

test_that("growth_estimates() reproduces the inverse observed information", {
  y   <- .gmm_se_data()
  fit <- suppressWarnings(fit_gmm(y, times = 4, n_classes = 2, n_init = 5,
                                  random_state = 1))
  ge  <- growth_estimates(fit, se = "hessian")

  L  <- cbind(1, 0:3)
  ll <- function(th) {
    a  <- rbind(th[1:2], th[3:4])
    P  <- matrix(c(th[5], th[6], th[6], th[7]), 2)
    Th <- diag(th[8:11])
    p1 <- plogis(th[12])
    d  <- sapply(1:2, function(k)
      log(c(p1, 1 - p1)[k]) +
        mvtnorm_free_logdens(y, as.vector(L %*% a[k, ]), L %*% P %*% t(L) + Th))
    m <- apply(d, 1, max)
    sum(m + log(rowSums(exp(d - m))))
  }
  mvtnorm_free_logdens <- function(y, mu, S) {
    ch <- chol(S)
    z  <- backsolve(ch, t(y) - mu, transpose = TRUE)
    -sum(log(diag(ch))) - ncol(y) / 2 * log(2 * pi) - colSums(z^2) / 2
  }

  g  <- fit$growth
  th <- c(g$means[1, ], g$means[2, ], g$psi[[1]][c(1, 2, 4)],
          g$residual_variance[1, ], qlogis(fit$weights[1]))
  ref <- sqrt(diag(solve(-stats::optimHess(th, ll))))

  pick <- function(parameter, item, class)
    ge$se[ge$parameter == parameter & ge$item == item & ge$class == class]
  ours <- c(pick("growth_mean", "intercept", 1), pick("growth_mean", "linear", 1),
            pick("growth_mean", "intercept", 2), pick("growth_mean", "linear", 2),
            pick("growth_variance", "intercept", 1),
            pick("growth_covariance", "intercept with linear", 1),
            pick("growth_variance", "linear", 1),
            vapply(g$time_labels, function(t) pick("residual_variance", t, 1), 0),
            pick("class_logit", "log(p1/p2)", 1))
  expect_length(ours, 12L)
  expect_lt(max(abs(ours / ref - 1)), 2e-3)
})

test_that("the three estimators agree on a correctly specified model", {
  y   <- .gmm_se_data(n = 600, seed = 2)
  fit <- suppressWarnings(fit_gmm(y, times = 4, n_classes = 2, n_init = 5,
                                  random_state = 1))
  r <- growth_estimates(fit)
  h <- growth_estimates(fit, se = "hessian")
  o <- growth_estimates(fit, se = "opg")
  expect_equal(r$estimate, h$estimate)
  expect_lt(max(abs(r$se / h$se - 1), na.rm = TRUE), 0.25)
  expect_lt(max(abs(o$se / h$se - 1), na.rm = TRUE), 0.25)
  expect_equal(r$z, r$estimate / r$se)
  V <- attr(r, "vcov")
  expect_equal(V, t(V))
  expect_equal(nrow(V), nrow(r))
  # The estimates are the ones measurement_summary() reports.
  ms <- suppressWarnings(utils::capture.output(m <- measurement_summary(fit)))
  expect_equal(r$estimate[seq_len(nrow(m))], m$estimate)
})

test_that("growth regressions and a free Psi get standard errors", {
  set.seed(3)
  n <- 500
  x <- rnorm(n)
  cls <- rbinom(n, 1, 0.5)
  icept <- ifelse(cls == 1, 3, 1) + 0.5 * x + rnorm(n, 0, 1)
  slope <- ifelse(cls == 1, 1.0, 0.2) + rnorm(n, 0, 0.4)
  y <- outer(icept, rep(1, 4)) + outer(slope, 0:3) + matrix(rnorm(n * 4, 0, 0.7), n, 4)
  fit <- suppressWarnings(fit_gmm(y, times = 4, n_classes = 2, psi = "free",
                                  growth_predictors = data.frame(x = x),
                                  n_init = 5, random_state = 1))
  ge <- growth_estimates(fit)
  reg <- ge[ge$parameter == "growth_regression", ]
  expect_true(nrow(reg) > 0L)
  expect_true(all(is.finite(ge$se)))
  b <- reg$estimate[reg$item == "intercept ON x"][1]
  s <- reg$se[reg$item == "intercept ON x"][1]
  expect_lt(abs(b - 0.5), 4 * s)
})

test_that("growth_estimates() refuses a fit that is not a growth mixture", {
  expect_error(growth_estimates(list()), "fit_gmm")
})
