# The Newton-type finish after EM (R/qn_finish.R). The claim it makes is that
# it lands where EM itself would land if it were run on to a strict absolute
# rule, only without the tens of thousands of iterations; these check exactly
# that, on fits with response probabilities on the 0/1 boundary (the case EM
# crawls on), plus the arithmetic the finish relies on.

.qn_lca_data <- function() {
  set.seed(3)
  n  <- 1500
  cl <- sample(1:3, n, TRUE, c(.40, .35, .25))
  P  <- rbind(c(0,   .9, .8, .7, .2, .3),
              c(.6,  1,  .2, .3, .8, .7),
              c(.3,  .4, .5, .9, .9, .1))
  t(vapply(cl, function(k) rbinom(6, 1, P[k, ]), numeric(6)))
}

.qn_ll <- function(s) sum(s$sample_weights * e_step(s, s$data, s$Y)$log_prob_norm)

test_that("a boundary LCA is finished to the maximum a long EM run reaches", {
  X  <- .qn_lca_data()
  bc <- list(categorical = 0, latent = 0)
  f  <- suppressWarnings(fit_mixture(X, n_classes = 3, measurement = "binary",
          n_init = 5, random_state = 1, bayes_constants = bc))
  f0 <- suppressWarnings(fit_mixture(X, n_classes = 3, measurement = "binary",
          n_init = 5, random_state = 1, bayes_constants = bc, refine = FALSE))
  expect_false(is.null(f$qn_finish))
  expect_true(f$converged)
  expect_gt(f$qn_finish$n_held, 0L)          # the two cells simulated at 0 and 1

  # EM from the unfinished solution, on an absolute rule of 1e-12.
  g <- f0
  g$mm$em_tol <- list(abs = 1e-12, rel = 0)
  g <- fit_single_init(g, g$data, g$Y, max_iter = 2e5, refine = FALSE,
                       init_state = g)
  expect_true(g$converged)
  expect_lt(abs(.qn_ll(f) - .qn_ll(g)), 1e-6)
  expect_gt(.qn_ll(f), .qn_ll(f0))
})

.qn_lta_data <- function() {
  set.seed(7)
  K <- 3; Tn <- 2; J <- 5; n <- 800
  rho <- matrix(c(1,   .80, .85, .75, .90,
                  .50, .50, .50, .50, .50,
                  .10, .20, 0,   .25, .10), nrow = K, byrow = TRUE)
  prev <- rbind(c(.5, .3, .2), c(.3, .4, .3))
  X <- matrix(NA_integer_, n, J * Tn)
  for (t in seq_len(Tn)) {
    s <- sample.int(K, n, TRUE, prob = prev[t, ])
    for (j in seq_len(J)) X[, (t - 1L) * J + j] <- rbinom(n, 1, rho[s, j])
  }
  X
}

test_that("a transition-free LTA is finished to the maximum a long EM run reaches", {
  args <- list(.qn_lta_data(), n_statuses = 3, times = 2,
               measurement = "binary", measurement_invariance = "full",
               n_init = 5, random_state = 1, n_cores = 1,
               .transition_free = TRUE, smoothing = 0,
               bayes_constants = list(categorical = 0))
  f  <- suppressWarnings(do.call(fit_lta, args))
  f0 <- suppressWarnings(do.call(fit_lta, c(args, refine = FALSE)))
  expect_false(is.null(f$qn_finish))
  long <- .lta_em(f0, f0$data, max_iter = 200000L, tol = 1e-15, alpha = 0)
  expect_true(long$converged)
  expect_lt(abs(f$loglik - long$loglik), 1e-6)
  expect_gt(f$loglik, f0$loglik + 1e-3)
})

test_that("the shared-row scores fold to the gradient of the reduced vector", {
  args <- list(.qn_lta_data(), n_statuses = 3, times = 2,
               measurement = "binary", measurement_invariance = "full",
               n_init = 2, random_state = 1, n_cores = 1,
               .transition_free = TRUE, refine = FALSE)
  f   <- suppressWarnings(do.call(fit_lta, args))
  red <- .qn_lta_reduction(f)
  expect_true(red$shared)
  # One free transition vector per occasion, not K rows of it.
  expect_lt(length(red$r0), length(red$src))
  w  <- f$weights_vec
  X  <- as.matrix(f$data)
  sc <- .lta_score_matrix(red$state, X, shared_ok = TRUE)
  ga <- as.vector(rowsum(colSums(sweep(sc$S, 1, w, "*")), red$src))
  gn <- .qn_fd_grad(function(r)
    sum(w * .lta_ll_case(red$state, X, r[red$src], red$layout)), red$r0)
  expect_lt(max(abs(ga - gn)), 1e-5)
  # The gate is unchanged for every other caller.
  expect_null(.lta_score_matrix(f, X))
})

test_that("candidate selection keeps distinct solutions within the band", {
  w   <- rep(1, 4)
  sig <- list(c(0, 0, 0, 0), c(0, 0, 0, 0.004), c(1, 1, 1, 1), c(2, 2, 2, 2))
  # The second is the first solution again; the fourth is outside the band.
  expect_identical(.qn_pick(c(-10, -10.01, -10.02, -11), sig, w), c(1L, 3L))
  expect_identical(.qn_pick(c(-10, -10.01, -10.02, -10.03), sig, w,
                            max_finish = 2L), c(1L, 3L))
})
