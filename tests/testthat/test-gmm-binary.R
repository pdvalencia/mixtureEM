# fit_gmm(family = "binomial"): the growth mixture model for a binary outcome,
# its random effects integrated out by Gauss-Hermite quadrature
# (R/gmm_binary.R).

.bg_q <- function(e) suppressMessages(suppressWarnings(e))

.bg_sim <- function(n = 300, seed = 1, ts = 0:3) {
  set.seed(seed)
  cl <- sample(1:2, n, TRUE, c(.45, .55))
  mu <- rbind(c(-1.5, .6), c(1, -.4))
  b  <- cbind(rnorm(n, 0, 1.1), rnorm(n, 0, .3))
  eta <- outer(mu[cl, 1] + b[, 1], rep(1, length(ts))) +
    outer(mu[cl, 2] + b[, 2], ts)
  matrix(rbinom(length(eta), 1, plogis(eta)), n)
}

test_that("one quadrature node is the latent class growth model", {
  Y <- .bg_sim()
  g <- .bg_q(fit_gmm(Y, times = 4, n_classes = 2, family = "binomial",
                     n_quadrature = 1, n_init = 3, random_state = 2, n_cores = 1))
  l <- .bg_q(fit_lcga(Y, times = 4, n_classes = 2, family = "binomial",
                      n_init = 3, random_state = 2, n_cores = 1))
  expect_equal(g$metrics$ll, l$metrics$ll, tolerance = 1e-6)
  expect_equal(g$growth$means, l$growth$coefficients, tolerance = 1e-4,
               ignore_attr = TRUE)
})

test_that("the class likelihood is the integral over the random effects", {
  ms <- binary_growth_model(2, .lcga_design(c(0, 1, 3), 1), n_quadrature = 30)
  ms$parameters$coefs <- rbind(c(-1, .5), c(.8, -.3))
  ms <- .bg_set_psi(ms, .bg_psi_unmap(ms, list(t(chol(matrix(c(1.5, -.2, -.2, .3), 2))))))
  X <- rbind(c(1, 0, 1), c(0, NA, 0), c(1, 1, 1))
  ll <- log_likelihood(ms, X)
  Psi <- ms$parameters$psi[[1]]; Lam <- cbind(1, c(0, 1, 3))
  brute <- function(y, a) {
    o <- !is.na(y)
    f <- function(b1, b2) {
      eta <- Lam %*% (a + c(b1, b2))
      prod(stats::dbinom(y[o], 1, stats::plogis(eta[o])))
    }
    inner <- function(b1) vapply(b1, function(u) stats::integrate(function(v)
      vapply(v, function(w) f(u, w) *
               mvtnorm_free_dens(c(u, w), Psi), 0), -6, 6, rel.tol = 1e-10)$value, 0)
    log(stats::integrate(inner, -8, 8, rel.tol = 1e-10)$value)
  }
  mvtnorm_free_dens <- function(x, S)
    exp(-0.5 * sum(x * solve(S, x))) / (2 * pi * sqrt(det(S)))
  for (i in seq_len(nrow(X))) for (k in 1:2)
    expect_equal(ll[i, k], brute(X[i, ], ms$parameters$coefs[k, ]),
                 tolerance = 1e-5)
  # Nodes moved to any centre and scale are the same integral: the weights
  # carry the change of variable.
  pat <- .bg_patterns(X)
  ms$node_centres <- list(key = pat$key,
    m = rep(list(matrix(c(.4, -.3, 1, .2, -.5, .1), 3)), 2),
    S = rep(list(matrix(c(.8, .1, 0, .9), 3, 4, byrow = TRUE)), 2))
  expect_equal(log_likelihood(ms, X), ll, tolerance = 1e-5)
})

test_that("each psi form counts its parameters", {
  Y <- .bg_sim(n = 200)
  np <- function(psi, re = "intercept_slope")
    .bg_q(fit_gmm(Y, times = 4, n_classes = 2, family = "binomial", psi = psi,
                  random_effects = re, n_quadrature = 5, n_init = 1,
                  random_state = 1, n_cores = 1))$metrics$n_params
  # 4 growth means, 1 class logit, then Psi.
  expect_equal(np("equal"), 5 + 3)
  expect_equal(np("variances"), 5 + 5)
  expect_equal(np("free"), 5 + 6)
  expect_equal(np("equal", "intercept"), 5 + 1)
  expect_equal(np("variances", "intercept"), 5 + 2)
})

test_that("the binary arguments are checked", {
  Y <- .bg_sim(n = 60)
  expect_error(fit_gmm(Y, times = 4, family = "binomial", residual = "constant"),
               "no residual variance")
  expect_error(fit_gmm(Y, times = 4, family = "binomial", random_effects = "none"),
               "fit_lcga")
  expect_error(fit_gmm(Y, times = 4, family = "binomial",
                       growth_predictors = rnorm(60)), "not available yet")
  expect_error(fit_gmm(Y, times = 4, family = "binomial", n_quadrature = 0),
               "n_quadrature")
  expect_error(fit_gmm(Y + 1, times = 4, family = "binomial"), "\\{0, 1\\}")
  expect_error(fit_gmm(cbind(Y, Y[, 1]), times = 5, degree = 2,
                       family = "binomial", random_effects = "all",
                       psi = "variances"), "at most two random effects")
})

test_that("EM climbs the quadrature likelihood at every step", {
  Y <- .bg_sim(n = 200) * 1.0
  for (psi in c("equal", "variances", "free")) {
    # The fixed rule: an adaptive one moves its nodes between iterations.
    ms <- binary_growth_model(2, .lcga_design(0:3, 1), psi = psi, n_quadrature = 7,
                              adaptive = FALSE)
    ms <- init_params(ms, Y, NULL, random_state = 4)
    w <- c(.5, .5); prev <- -Inf
    for (it in 1:15) {
      lj  <- sweep(log_likelihood(ms, Y), 2, log(w), "+")
      tot <- sum(logsumexp(lj, MARGIN = 1))
      expect_gte(tot, prev - 1e-8)
      prev <- tot
      resp <- exp(lj - logsumexp(lj, MARGIN = 1))
      w  <- colMeans(resp)
      ms <- m_step(ms, Y, resp)
    }
  }
})

test_that("the finish's analytic scores match central differences", {
  Y <- .bg_sim(n = 150, ts = c(0, 1, 3, 4))
  Y[sample(length(Y), 40)] <- NA
  for (psi in c("equal", "variances", "free")) {
    f <- .bg_q(fit_gmm(Y, times = 4, time_scores = c(0, 1, 3, 4),
                       family = "binomial", psi = psi, n_quadrature = 5,
                       n_init = 1, random_state = 1, n_cores = 1))
    pr <- .qn_mixture_problem(f, f$data, NULL)
    an <- .qn_mixture_scores(pr$unpack(pr$par), pr$lay, f$data, NULL)$S
    fd <- .qn_fd_scores(pr$case_ll, pr$par, seq_along(pr$par))
    expect_equal(an, fd, tolerance = 1e-6, label = psi)
  }
})

test_that("classes keep their own covariance when they are re-sorted", {
  Y <- .bg_sim(n = 250, seed = 3)
  f <- .bg_q(fit_gmm(Y, times = 4, n_classes = 2, family = "binomial",
                     psi = "free", n_quadrature = 5, n_init = 2,
                     random_state = 1, n_cores = 1))
  expect_true(f$weights[1] >= f$weights[2])
  expect_equal(f$mm$parameters$psi,
               lapply(f$mm$parameters$chol, tcrossprod))
  # The reported log-likelihood is the likelihood of the reported parameters.
  ll <- sum(logsumexp(sweep(log_likelihood(f$mm, f$data), 2, log(f$weights), "+"),
                      MARGIN = 1))
  expect_equal(ll, f$metrics$ll, tolerance = 1e-6)
})

test_that("growth_estimates() reports the binary growth parameters", {
  Y <- .bg_sim(n = 300, seed = 5)
  f <- .bg_q(fit_gmm(Y, times = 4, n_classes = 2, family = "binomial",
                     n_quadrature = 7, n_init = 2, random_state = 1, n_cores = 1))
  ge <- .bg_q(growth_estimates(f, se = "hessian"))
  expect_setequal(intersect(unique(ge$parameter),
                            c("growth_mean", "growth_variance", "growth_covariance",
                              "residual_variance")),
                  c("growth_mean", "growth_variance", "growth_covariance"))
  expect_true(all(is.finite(ge$se[ge$parameter == "growth_mean"])))
  expect_output(print(f), "Gauss-Hermite")
  expect_true(all(f$growth$fitted >= 0 & f$growth$fitted <= 1))
})

test_that("predict() scores the estimation sample as the fit does", {
  Y <- .bg_sim(n = 200, seed = 6)
  Y[sample(length(Y), 30)] <- NA
  f <- .bg_q(fit_gmm(Y, times = 4, n_classes = 2, family = "binomial",
                     n_quadrature = 5, n_init = 2, random_state = 1, n_cores = 1))
  p <- predict(f, newdata = Y, type = "posterior")
  expect_equal(unname(as.matrix(p)), unname(exp(f$log_resp)), tolerance = 1e-8)
})

test_that("replicates are drawn with the random effects", {
  ms <- binary_growth_model(1, .lcga_design(0:3, 1), random_effects = "intercept",
                            n_quadrature = 5)
  ms$parameters$coefs <- matrix(c(0, 0), 1)
  ms <- .bg_set_psi(ms, 2)                        # intercept SD 2
  set.seed(1)
  Y <- generate_synthetic_data(ms, rep(1L, 4000), 4000)
  expect_true(all(Y %in% c(0, 1)))
  # Answers share a case's random intercept: well above zero, which is what
  # occasion-by-occasion draws from the class curve would give.
  expect_gt(stats::cor(Y[, 1], Y[, 2]), 0.3)
})

test_that("blrt() and report_methods() take a binary growth mixture", {
  Y <- .bg_sim(n = 150, seed = 8)
  f <- .bg_q(fit_gmm(Y, times = 4, n_classes = 2, family = "binomial",
                     n_quadrature = 5, n_init = 2, random_state = 1, n_cores = 1))
  b <- .bg_q(blrt(from_fit = f, k_small = 1, k_large = 2, n_reps = 2,
                  n_init_base = 2, n_init_boot = 2, n_cores = 1, verbose = FALSE))
  expect_true(is.finite(b$p_value))
  expect_match(paste(capture.output(report_methods(f)), collapse = " "),
               "binary outcome")
})

test_that("half the starts are drawn around the start pool, inside the edge", {
  ms <- binary_growth_model(3, .lcga_design(0:4, 1), n_quadrature = 5,
                            start_pool = list(rbind(c(7, 0), c(-3, 1), c(30, -9))))
  Y <- .bg_sim(n = 100, seed = 9, ts = 0:4)
  near <- vapply(1:60, function(s) {
    st <- init_params(ms, Y, NULL, random_state = s)
    expect_true(all(abs(st$parameters$coefs) <= .qn_edge))
    abs(st$parameters$coefs[1, 1] - 7) < 2.5
  }, TRUE)
  expect_gt(mean(near), 0.3)
  expect_lt(mean(near), 0.7)
  # Each start is a function of its seed alone.
  expect_identical(init_params(ms, Y, NULL, random_state = 4)$parameters,
                   init_params(ms, Y, NULL, random_state = 4)$parameters)
})

test_that("starts without a pool seed classes from observed response patterns", {
  # Two patterns only, 11000 and 00000: a class seeded from the first has the
  # least-squares line through (A, A, -A, -A, -A), slope -0.6 per unit of
  # intercept whatever A is, and from the second intercept -A, slope 0.
  Y <- rbind(matrix(c(1, 1, 0, 0, 0), 20, 5, byrow = TRUE), matrix(0, 80, 5))
  ms <- binary_growth_model(3, .lcga_design(0:4, 1), n_quadrature = 5)
  rows <- do.call(rbind, lapply(1:40, function(s)
    init_params(ms, Y, NULL, random_state = s)$parameters$coefs))
  step <- abs(rows[, 2] / rows[, 1] + 0.6) < 1e-8 & rows[, 1] > 2 & rows[, 1] < 10
  none <- abs(rows[, 2]) < 1e-8 & rows[, 1] < -2 & rows[, 1] > -10
  expect_gt(sum(step), 5)
  expect_gt(sum(none), 5)
})

test_that("the start pool is the latent class growth fit, and only for a search", {
  Y <- .bg_sim(n = 150, seed = 10)
  D <- .lcga_design(0:3, 1)
  pool <- .bg_q(.bg_start_pool(Y, 2, D, n_init = 3, random_state = 1, n_cores = 1))
  l <- .bg_q(fit_lcga(Y, times = 4, n_classes = 2, family = "binomial",
                      n_init = 3, random_state = 1, n_cores = 1))
  expect_length(pool, 1L)
  expect_equal(pool[[1]], l$growth$coefficients, tolerance = 1e-6,
               ignore_attr = TRUE)
  expect_null(.bg_start_pool(Y, 2, D, n_init = 1))
  expect_null(.bg_start_pool(Y, 2, D, n_init = 5, refine_from = list()))
})
