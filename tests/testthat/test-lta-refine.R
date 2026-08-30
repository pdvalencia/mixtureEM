# The post-EM L-BFGS refinement in fit_lta(). Three properties, in the order
# they matter: the gradient it climbs is the gradient of the objective it
# claims, the objective is the one EM maximises, and the climb never hands back
# a worse fit or a fit that has quietly dropped a constraint.

# `.lta_refine_sim()` and `.ml` are in helper-lta-packing.R, which
# test-lta-scaling.R shares.

test_that("the analytic gradient matches central finite differences", {
  X <- .lta_refine_sim()
  fit <- suppressMessages(suppressWarnings(
    fit_lta(X, n_statuses = 3, times = 3, measurement = "binary",
            smoothing = 0, bayes_constants = .ml, n_init = 2, random_state = 1,
            refine = FALSE, standard_errors = FALSE)))

  layout <- .lta_par_layout(fit)
  par0   <- .lta_par_pack(fit, layout)
  w      <- fit$weights_vec

  # The packing and the score columns must describe the same parameters. This
  # is free and it catches a whole class of layout error on its own.
  expect_equal(length(par0), ncol(.lta_score_matrix(fit, X)$S))

  # Packing and unpacking must be the identity on a fit.
  rt <- .lta_par_unpack(par0, fit, layout)
  expect_equal(rt$delta, fit$delta, tolerance = 1e-10)
  expect_equal(rt$tau, fit$tau, tolerance = 1e-10)

  obj <- function(p) {
    st <- .lta_par_unpack(p, fit, layout)
    sum(w * .lta_score_matrix(st, X)$ll) +
      .lta_penalty(st, X, layout, 0)$value
  }
  ana <- colSums(sweep(.lta_score_matrix(fit, X)$S, 1, w, "*")) +
    .lta_penalty(fit, X, layout, 0)$gradient
  eps <- 1e-5
  fd <- vapply(seq_along(par0), function(i) {
    e <- numeric(length(par0)); e[i] <- eps
    (obj(par0 + e) - obj(par0 - e)) / (2 * eps)
  }, numeric(1))

  expect_lt(max(abs(ana - fd)) / max(1, max(abs(fd))), 1e-5)
})

test_that("the gradient vanishes at an EM fixed point, priors on and off", {
  X <- .lta_refine_sim()
  for (priors in c(FALSE, TRUE)) {
    extra <- if (priors) list() else
      list(smoothing = 0, bayes_constants = .ml)
    fit <- suppressMessages(suppressWarnings(do.call(fit_lta, c(
      list(X, n_statuses = 3, times = 3, measurement = "binary", n_init = 2,
           random_state = 1, tol = 1e-14, max_iter = 20000, refine = FALSE,
           standard_errors = FALSE), extra))))
    layout <- .lta_par_layout(fit)
    g <- colSums(sweep(.lta_score_matrix(fit, X)$S, 1, fit$weights_vec, "*")) +
      .lta_penalty(fit, X, layout, if (priors) 1 else 0)$gradient
    # An EM fixed point is a stationary point of the objective EM maximises.
    # If the refinement's objective were a near neighbour of it rather than the
    # thing itself, this is the check that would fail.
    expect_lt(max(abs(g)), 1e-3)
  }
})

test_that("the refinement never lowers the log-likelihood", {
  X <- .lta_refine_sim()
  for (seed in 1:3) {
    a <- suppressMessages(suppressWarnings(
      fit_lta(X, n_statuses = 3, times = 3, measurement = "binary",
              smoothing = 0, bayes_constants = .ml, n_init = 3,
              random_state = seed, refine = FALSE, standard_errors = FALSE)))
    b <- suppressMessages(suppressWarnings(
      fit_lta(X, n_statuses = 3, times = 3, measurement = "binary",
              smoothing = 0, bayes_constants = .ml, n_init = 3,
              random_state = seed, refine = TRUE, standard_errors = FALSE)))
    expect_gte(b$loglik, a$loglik - 1e-8)
    expect_equal(b$n_params, a$n_params)
  }
})

test_that("the refinement respects measurement invariance", {
  X <- .lta_refine_sim()
  fit <- suppressMessages(suppressWarnings(
    fit_lta(X, n_statuses = 3, times = 3, measurement = "binary",
            measurement_invariance = "full", smoothing = 0,
            bayes_constants = .ml, n_init = 3, random_state = 2,
            refine = TRUE, standard_errors = FALSE)))

  # The constraint is carried by the measurement model, which exists from the
  # moment fit_lta() builds it; `state$longitudinal` is attached only once the
  # search is over, so a refinement reading the invariant items from there sees
  # nothing and frees every item. It then buys log-likelihood by breaking the
  # restriction the model is defined by. These two assertions are what caught
  # exactly that.
  p1 <- fit$mm$models[[1]]$parameters$pis
  for (t in 2:3) expect_equal(fit$mm$models[[t]]$parameters$pis, p1)
  expect_equal(.lta_invariant_items(fit), seq_len(4))
})

test_that("the refinement declines the models whose parameters it cannot pack", {
  X <- .lta_refine_sim()
  Z <- data.frame(z = stats::rnorm(nrow(X)))

  # A covariate on the initial status: delta is a regression, not a simplex.
  cov_fit <- suppressMessages(suppressWarnings(
    fit_lta(X, n_statuses = 3, times = 3, measurement = "binary",
            predictors_initial = Z, n_init = 2, random_state = 1,
            standard_errors = FALSE)))
  expect_false(.lta_scores_full(cov_fit))
  expect_identical(.lta_refine_lbfgs(cov_fit, X), cov_fit)

  # A mixture over chains: every score becomes class-conditional.
  mix_fit <- suppressMessages(suppressWarnings(
    fit_lta(X, n_statuses = 2, times = 3, measurement = "binary",
            n_classes = 2, n_init = 2, random_state = 1,
            standard_errors = FALSE)))
  expect_false(.lta_scores_full(mix_fit))
  expect_identical(.lta_refine_lbfgs(mix_fit, X), mix_fit)
})
