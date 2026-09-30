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

test_that("a freely loading random intercept on binary items is finished to the long-EM maximum", {
  sim  <- .lta_ri_sim(n = 400, Tn = 3, J = 4, seed = 11)
  args <- list(sim$X, n_statuses = 2, times = 3, measurement = "binary",
               measurement_invariance = "full", random_intercept = "continuous",
               n_quadrature = 10, n_init = 3, random_state = 1, n_cores = 1,
               smoothing = 0, bayes_constants = list(categorical = 0),
               standard_errors = FALSE)
  f <- suppressWarnings(do.call(fit_lta, args))
  expect_false(is.null(f$qn_finish))
  expect_true(f$qn_finish$converged)
  # EM from the same fit on a strict rule reaches the same maximum.
  long <- .lta_em(f, f$data, max_iter = 200000L, tol = 1e-15, alpha = 0)
  expect_true(long$converged)
  expect_lt(abs(f$loglik - long$loglik), 1e-6)
  expect_equal(as.vector(f$ri$L), as.vector(long$ri$L), tolerance = 1e-4)
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

# Two Gaussian means tied together by a prior with off-diagonal curvature:
# the prior's diagonal, all BHHH uses of it, understates the curvature along
# p1 - p2, and at k = 99 n each BHHH step removes only 1/100 of the error
# there. The maximum is known in closed form. A third, optional coordinate
# the objective ignores makes the exact Hessian singular.
.qn_crawl <- function(k, flat = FALSE) {
  set.seed(5)
  n <- 50; y <- rnorm(n, 1); z <- rnorm(n, -1)
  np <- if (flat) 3L else 2L
  case_ll <- function(p) -0.5 * (y - p[1])^2 - 0.5 * (z - p[2])^2
  scores  <- function(p, idx) cbind(y - p[1], z - p[2], 0)[, idx, drop = FALSE]
  prior   <- function(p) -0.5 * k * (p[1] - p[2])^2
  prior_grad <- function(p) {
    g <- c(-k, k, 0)[seq_len(np)] * (p[1] - p[2])
    list(g = g, h = c(-k, -k, 0)[seq_len(np)])
  }
  A <- matrix(c(n + k, -k, -k, n + k), 2)
  list(case_ll = case_ll, scores = scores, prior = prior,
       prior_grad = prior_grad, w = rep(1, n), par0 = rep(0, np),
       top = solve(A, c(sum(y), sum(z))))
}

test_that("a crawling finish switches to the exact curvature and reaches the maximum", {
  p <- .qn_crawl(k = 99 * 50)
  kind <- rep("free", 2)
  r <- .qn_finish(p$case_ll, p$w, p$par0, p$prior, p$scores, kind = kind,
                  prior_grad = p$prior_grad)
  expect_equal(r$n_hessian, 1L)
  expect_true(r$converged)
  expect_lt(r$iterations, 110L)
  expect_lt(max(abs(r$par - p$top)), 1e-8)
  # Without analytic scores there is no switch, and BHHH is still crawling
  # at the cap.
  b <- .qn_finish(p$case_ll, p$w, p$par0, p$prior, NULL, kind = kind,
                  prior_grad = p$prior_grad)
  expect_equal(b$n_hessian, 0L)
  expect_gt(b$iterations, 500L)
  expect_gt(r$value, b$value)
})

test_that("a finish that converges before the switch point never builds a Hessian", {
  p <- .qn_crawl(k = 0)
  r <- .qn_finish(p$case_ll, p$w, p$par0, p$prior, p$scores,
                  kind = rep("free", 2), prior_grad = p$prior_grad)
  expect_equal(r$n_hessian, 0L)
  expect_lt(r$iterations, 100L)
  expect_lt(max(abs(r$par - p$top)), 1e-8)
})

test_that("a Hessian that is not negative definite is not used", {
  p <- .qn_crawl(k = 99 * 50, flat = TRUE)
  r <- .qn_finish(p$case_ll, p$w, p$par0, p$prior, p$scores,
                  kind = rep("free", 3), prior_grad = p$prior_grad)
  # Tried at iteration 100 and again 100 later, singular both times: BHHH
  # stays in charge throughout and is still climbing at the cap.
  expect_equal(r$n_hessian, 2L)
  expect_gt(r$iterations, 500L)
  expect_identical(r$par[3], 0)
  v0 <- sum(p$case_ll(p$par0)) + p$prior(p$par0)
  expect_gt(r$value, v0)
})
