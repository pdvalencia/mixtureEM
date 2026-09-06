# Covariates on the continuous random intercept (`predictors_random_intercept`,
# internal/parts/part-r8-ri-factor-covariates.md, W8). The factor
# f_i ~ N(x_i'beta, 1) is fitted by reweighting the fixed Gauss-Hermite grid
# rather than shifting it per case; the nine rungs below are the cheap ladder
# that has to be green before the 35-50 minute Dating-data fits in W10-W11 are
# ever run. `.lta_refine_sim()` and `.ml` are in helper-lta-packing.R.

# --- 1. beta = 0 reduces exactly to the no-covariate node prior ------------

test_that("beta = 0 reduces exactly to the no-covariate node prior", {
  X <- .lta_refine_sim(n = 400, K = 2, Tn = 3, J = 4, seed = 11)
  set.seed(11); Z <- matrix(rnorm(400), 400, 1)
  fit <- suppressWarnings(fit_lta(X, n_statuses = 2, times = 3,
    measurement = "binary", random_intercept = "continuous", n_quadrature = 10,
    predictors_random_intercept = Z, n_init = 1, max_iter = 15,
    random_state = 3, standard_errors = FALSE, refine = FALSE))

  # Direct likelihood evaluation, not .lta_em(), so nothing here can be
  # confounded by .lta_ri_sign_normalise()'s own auto-correction.
  st_zero <- fit
  st_zero$ri_beta <- matrix(0, ncol(Z), 1L)
  ll_zero <- sum(fit$weights_vec * .lta_ri_ll_case(st_zero, X))

  st_noZ <- fit
  st_noZ$Z_ri    <- NULL
  st_noZ$ri_beta <- NULL
  ll_noZ <- sum(fit$weights_vec * .lta_ri_ll_case(st_noZ, X))

  expect_equal(ll_zero, ll_noZ, tolerance = 1e-10)
})

# --- 2. n_quadrature = 1 reduces RI-with-covariates to regular LTA ---------

test_that("n_quadrature = 1 reduces RI-with-covariates exactly to regular LTA", {
  X <- .lta_refine_sim(n = 400, K = 2, Tn = 3, J = 4, seed = 11)
  set.seed(11); Z <- matrix(rnorm(400), 400, 1)
  fit_reg <- suppressWarnings(fit_lta(X, n_statuses = 2, times = 3,
    measurement = "binary", smoothing = 0, bayes_constants = .ml,
    n_init = 1, random_state = 3, tol = 1e-12, max_iter = 3000,
    standard_errors = FALSE))
  fit_ri1 <- suppressWarnings(fit_lta(X, n_statuses = 2, times = 3,
    measurement = "binary", smoothing = 0, bayes_constants = .ml,
    n_init = 1, random_state = 3, tol = 1e-12, max_iter = 3000,
    standard_errors = FALSE, random_intercept = "continuous",
    n_quadrature = 1, predictors_random_intercept = Z))
  # At Q = 1 fbar is identically 0 (one node, mass 1), so the WLS closed form
  # returns beta = 0 from the very first M-step and stays there.
  expect_lt(max(abs(fit_ri1$ri_beta)), 1e-8)
  expect_lt(abs(fit_ri1$loglik - fit_reg$loglik), 1e-8)
})

# --- 3. The two RI likelihood paths agree at a nonzero covariate effect ----

test_that("the two RI likelihood paths agree at a nonzero covariate effect", {
  X <- .lta_refine_sim(n = 400, K = 2, Tn = 3, J = 4, seed = 11)
  set.seed(11); Z <- matrix(rnorm(400), 400, 1)
  fit <- suppressWarnings(fit_lta(X, n_statuses = 2, times = 3,
    measurement = "binary", random_intercept = "continuous", n_quadrature = 10,
    predictors_random_intercept = Z, n_init = 1, max_iter = 25,
    random_state = 3, standard_errors = FALSE))
  ll_e    <- .lta_ri_e_step(fit, X, fit$weights_vec)$ll
  ll_case <- .lta_ri_ll_case(fit, X)
  expect_lt(max(abs(ll_e - ll_case)), 1e-10)
})

# --- 4. EM monotonicity on the penalised objective -------------------------

test_that("EM stays monotone on the penalised objective with covariates on the factor", {
  X <- .lta_refine_sim(n = 400, K = 2, Tn = 3, J = 4, seed = 11)
  set.seed(11); Z <- matrix(rnorm(400), 400, 1)
  st <- suppressWarnings(fit_lta(X, n_statuses = 2, times = 3,
    measurement = "binary", random_intercept = "continuous", n_quadrature = 10,
    predictors_random_intercept = Z, n_init = 1, max_iter = 1,
    random_state = 3, standard_errors = FALSE, refine = FALSE))

  obj <- function(s) s$loglik + .lta_log_prior(s, X, 1)
  vals <- numeric(50)
  for (i in seq_len(50)) {
    st <- .lta_em(st, X, max_iter = 1L, tol = 0, alpha = 1)
    vals[i] <- obj(st)
  }
  expect_true(all(diff(vals) > -1e-9))
})

# --- 5. Finite differences on the ri_beta score block, at a non-optimal point

test_that("the ri_beta score block matches finite differences at a perturbed point", {
  X <- .lta_refine_sim(n = 400, K = 2, Tn = 3, J = 4, seed = 11)
  set.seed(11); Z <- matrix(rnorm(400), 400, 1)
  # max_iter kept small on purpose: at the optimum both sides of the check are
  # ~0 and it proves nothing.
  fit <- suppressWarnings(fit_lta(X, n_statuses = 2, times = 3,
    measurement = "binary", random_intercept = "continuous", n_quadrature = 10,
    predictors_random_intercept = Z, n_init = 1, max_iter = 8,
    random_state = 3, standard_errors = FALSE, refine = FALSE))

  layout <- .lta_par_layout(fit)
  par0   <- .lta_par_pack(fit, layout)
  w      <- fit$weights_vec
  sc0    <- .lta_score_matrix(fit, X)

  obj <- function(p) {
    st <- .lta_par_unpack(p, fit, layout)
    sum(w * .lta_score_matrix(st, X)$ll) + .lta_penalty(st, X, layout, 0)$value
  }
  ana <- colSums(sweep(sc0$S, 1, w, "*")) + .lta_penalty(fit, X, layout, 0)$gradient
  eps <- 1e-5
  fd <- vapply(seq_along(par0), function(i) {
    e <- numeric(length(par0)); e[i] <- eps
    (obj(par0 + e) - obj(par0 - e)) / (2 * eps)
  }, numeric(1))
  rel <- abs(ana - fd) / pmax(1, abs(fd))
  expect_lt(max(rel), 1e-4)

  beta_cols <- Find(function(b) identical(b$name, "ri_beta"), sc0$blocks)$cols
  expect_false(is.null(beta_cols))
  expect_lt(max(rel[beta_cols]), 1e-4)
})

# --- 6. Quadrature accuracy of the reweighting identity against integrate() -

test_that("the covariate reweighting identity matches stats::integrate()", {
  # One item, a steep loading (4.6) -- the setting measured in the work order
  # against a deliberately hard response. This is the rung that catches a
  # dropped `- mu^2 / 2`.
  A0 <- 0; lambda <- 4.6
  p_given_f <- function(f) stats::plogis(A0 + lambda * f)
  for (mu in c(0, 1, 2, 3)) {
    exact <- stats::integrate(function(f) p_given_f(f) * stats::dnorm(f, mu, 1),
                              -Inf, Inf)$value
    for (Qn in c(40L, 80L)) {
      gh <- .gauss_hermite(Qn)
      approx <- sum(gh$w * exp(mu * gh$z - mu^2 / 2) * p_given_f(gh$z))
      tol <- if (Qn == 40L) 5e-4 else 5e-5
      expect_lt(abs(approx - exact), tol)
    }
  }
})

# --- 7. Flip invariance: L, Dnode and ri_beta must move together -----------

test_that("flip invariance: L, Dnode and ri_beta must move together", {
  X <- .lta_refine_sim(n = 400, K = 2, Tn = 3, J = 4, seed = 11)
  set.seed(11); Z <- matrix(rnorm(400), 400, 1)
  fit <- suppressWarnings(fit_lta(X, n_statuses = 2, times = 3,
    measurement = "binary", random_intercept = "continuous", n_quadrature = 10,
    predictors_random_intercept = Z, n_init = 1, max_iter = 15,
    random_state = 3, standard_errors = FALSE, refine = FALSE))
  w <- fit$weights_vec

  # Evaluated with .lta_ri_ll_case() directly -- going through .lta_em() would
  # let .lta_ri_sign_normalise() auto-correct a flipped loading right back,
  # which would make this test pass for the wrong reason.
  ll0 <- sum(w * .lta_ri_ll_case(fit, X))

  flipped <- fit
  flipped$ri$L     <- -flipped$ri$L
  flipped$ri$Dnode <- -flipped$ri$Dnode
  flipped$ri_beta  <- -flipped$ri_beta
  ll_flip <- sum(w * .lta_ri_ll_case(flipped, X))
  expect_equal(ll_flip, ll0, tolerance = 1e-10)

  # The second half is what stops the first half from being vacuous: negating
  # only L and Dnode, leaving ri_beta behind, must change the likelihood.
  partial <- fit
  partial$ri$L     <- -partial$ri$L
  partial$ri$Dnode <- -partial$ri$Dnode
  ll_partial <- sum(w * .lta_ri_ll_case(partial, X))
  expect_gt(abs(ll_partial - ll0), 1e-6)
})

# --- 8. Packing invariants on a covariate fit -------------------------------

test_that("packing invariants hold on a random-intercept covariate fit", {
  X <- .lta_refine_sim(n = 400, K = 2, Tn = 3, J = 4, seed = 11)
  set.seed(11); Z <- matrix(rnorm(400), 400, 1)
  fit <- suppressWarnings(fit_lta(X, n_statuses = 2, times = 3,
    measurement = "binary", random_intercept = "continuous", n_quadrature = 10,
    predictors_random_intercept = Z, n_init = 1, max_iter = 15,
    random_state = 3, standard_errors = FALSE))

  layout <- .lta_par_layout(fit)
  sc     <- .lta_score_matrix(fit, X)
  expect_equal(vapply(layout, function(b) b$len, integer(1)),
               vapply(sc$blocks, function(b) length(b$cols), integer(1)))
  par <- .lta_par_pack(fit, layout)
  expect_equal(length(par), ncol(sc$S))
  expect_equal(ncol(sc$S), fit$n_params)

  st2  <- .lta_par_unpack(par, fit, layout)
  par2 <- .lta_par_pack(st2, layout)
  expect_equal(par2, par, tolerance = 1e-12)
  expect_equal(st2$ri_beta, fit$ri_beta, tolerance = 1e-12)

  ll_case <- .lta_ll_case(fit, X, par, layout)
  ll_e    <- .lta_ri_e_step(fit, X, fit$weights_vec)$ll
  expect_equal(ll_case, ll_e, tolerance = 1e-10)
  expect_equal(sum(fit$weights_vec * ll_case), fit$loglik, tolerance = 1e-8)
})

# --- 9. Recovery: one replication of a factor-on-covariates generator ------

test_that("the factor-on-covariates coefficients are recovered", {
  skip_on_cran()
  n <- 3000L; Tn <- 3L; K <- 2L; J <- 5L
  set.seed(202)
  Z <- cbind(x1 = stats::rnorm(n), x2 = stats::rnorm(n))
  beta_true <- c(0.6, -0.4)
  f <- stats::rnorm(n, as.vector(Z %*% beta_true), 1)

  delta <- c(0.5, 0.5)
  tau   <- matrix(c(0.75, 0.25, 0.25, 0.75), K, K, byrow = TRUE)
  A <- matrix(c(1, -1), K, J)
  L <- rep(1.0, J)
  S <- matrix(0L, n, Tn)
  S[, 1] <- sample.int(K, n, replace = TRUE, prob = delta)
  for (t in 2:Tn) for (k in seq_len(K)) {
    idx <- S[, t - 1] == k
    if (any(idx))
      S[idx, t] <- sample.int(K, sum(idx), replace = TRUE, prob = tau[k, ])
  }
  X <- matrix(NA_real_, n, Tn * J)
  for (t in seq_len(Tn)) for (j in seq_len(J)) {
    eta <- A[S[, t], j] + L[j] * f
    X[, (t - 1) * J + j] <- stats::rbinom(n, 1, stats::plogis(eta))
  }
  colnames(X) <- paste0("t", rep(seq_len(Tn), each = J), "_i", rep(seq_len(J), Tn))

  fit <- suppressWarnings(fit_lta(X, n_statuses = K, times = Tn,
    measurement = "binary", random_intercept = "continuous", n_quadrature = 15,
    predictors_random_intercept = Z, smoothing = 0, bayes_constants = .ml,
    n_init = 20, random_state = 9, standard_errors = FALSE))

  s <- sign(fit$ri$L[which.max(abs(fit$ri$L)), 1])
  betahat <- as.vector(s * fit$ri_beta)
  expect_lt(max(abs(betahat - beta_true)), 0.12)
})
