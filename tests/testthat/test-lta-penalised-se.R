# The latent transition model's standard errors are the variance of the
# estimator the fit computes: the inverse curvature of the log-likelihood plus
# the priors (default), and the sandwich H^-1 (B - Hp) H^-1 around it
# (`"robust"`), B the uncentred N / (N - 1) outer product of the case scores.
# Rebuilt here from the fit's own likelihood and prior by differencing the
# log-likelihood itself, which the package does not do.

.pse_sim <- function(n = 300, seed = 5) {
  set.seed(seed)
  S <- matrix(0L, n, 2)
  S[, 1] <- sample.int(3, n, TRUE, prob = c(.08, .32, .60))
  tau <- rbind(c(.85, .10, .05), c(.05, .90, .05), c(.02, .08, .90))
  S[, 2] <- vapply(S[, 1], function(k) sample.int(3, 1, prob = tau[k, ]), 1L)
  rho <- rbind(rep(.85, 6), c(.85, .85, .85, .15, .15, .15), rep(.15, 6))
  X <- matrix(0L, n, 12)
  for (t in 1:2) for (j in 1:6)
    X[, (t - 1) * 6 + j] <- rbinom(n, 1, rho[S[, t], j])
  X
}

.pse_pieces <- function(fit, alpha) {
  X   <- as.matrix(fit$data)
  w   <- fit$weights_vec
  lay <- .lta_par_layout(fit)
  par <- .lta_par_pack(fit, lay)
  ll  <- function(p) sum(w * .lta_ll_case(fit, X, p, lay))
  lp  <- function(p) .lta_log_prior(.lta_par_unpack(p, fit, lay), X, alpha)
  S   <- .lta_score_matrix(fit, X)$S
  n   <- sum(w)
  list(Hd = .fd_hessian_ref(ll, par), Hp = .fd_hessian_ref(lp, par),
       B = n / (n - 1) * crossprod(S, w * S))
}

test_that("the default is the inverse curvature of the penalised objective", {
  fit <- fit_lta(.pse_sim(), n_statuses = 3, times = 2, measurement = "binary",
                 measurement_invariance = "full", n_init = 5, random_state = 1,
                 n_cores = 1)
  expect_equal(fit$se$method, "Observed information")
  pc <- .pse_pieces(fit, 1)
  V  <- solve(-(pc$Hd + pc$Hp))
  # 1e-3: the two Hessians are differenced along different routes.
  expect_lt(max(abs(sqrt(diag(fit$se$vcov)) / sqrt(diag(V)) - 1)), 1e-3)
  # The prior is felt here (a status holding 8 percent): leaving it out moves
  # some standard error by more than the tolerance above.
  expect_gt(max(abs(sqrt(diag(solve(-pc$Hd))) / sqrt(diag(V)) - 1)), 5e-3)
})

test_that("robust is the sandwich around that curvature", {
  fit <- fit_lta(.pse_sim(), n_statuses = 3, times = 2, measurement = "binary",
                 measurement_invariance = "full", n_init = 5, random_state = 1,
                 n_cores = 1, standard_errors = "robust")
  expect_true(fit$se$robust)
  pc <- .pse_pieces(fit, 1)
  Ai <- solve(-(pc$Hd + pc$Hp))
  V  <- Ai %*% (pc$B - pc$Hp) %*% Ai
  expect_lt(max(abs(sqrt(diag(fit$se$vcov)) / sqrt(diag(V)) - 1)), 1e-3)
})

test_that("with every prior off they are the plain observed information", {
  # A larger sample: at n = 300 a transition sits on the boundary without the
  # prior, and the plain information is then singular.
  fit <- fit_lta(.pse_sim(n = 1500), n_statuses = 3, times = 2,
                 measurement = "binary",
                 measurement_invariance = "full", n_init = 5, random_state = 1,
                 n_cores = 1, smoothing = 0,
                 bayes_constants = list(categorical = 0))
  pc <- .pse_pieces(fit, 0)
  expect_equal(max(abs(pc$Hp)), 0)
  expect_lt(max(abs(sqrt(diag(fit$se$vcov)) /
                      sqrt(diag(solve(-pc$Hd))) - 1)), 1e-3)
})
