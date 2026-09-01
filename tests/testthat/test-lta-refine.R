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

test_that("the analytic gradient matches finite differences for Gaussian variances", {
  # .lta_gaussian_refine_sim() is in helper-lta-packing.R. Priors on and off,
  # the same two passes the binary fixture above uses, since the log_sd block
  # is the first Gaussian block to carry a prior at all (means never have).
  X <- .lta_gaussian_refine_sim()
  for (priors in c(FALSE, TRUE)) {
    extra <- if (priors) list() else list(smoothing = 0, bayes_constants = .ml)
    fit <- suppressMessages(suppressWarnings(do.call(fit_lta, c(
      list(X, n_statuses = 3, times = 2, measurement = "continuous",
           n_init = 2, random_state = 1, refine = FALSE,
           standard_errors = FALSE), extra))))

    layout <- .lta_par_layout(fit)
    par0   <- .lta_par_pack(fit, layout)
    w      <- fit$weights_vec

    expect_true(any(vapply(layout, function(b) b$kind == "log_sd", logical(1))))
    expect_equal(length(par0), ncol(.lta_score_matrix(fit, X)$S))

    rt <- .lta_par_unpack(par0, fit, layout)
    expect_equal(rt$mm$models[[1]]$parameters$covariances,
                fit$mm$models[[1]]$parameters$covariances, tolerance = 1e-10)

    alpha <- if (priors) 1 else 0
    obj <- function(p) {
      st <- .lta_par_unpack(p, fit, layout)
      sum(w * .lta_score_matrix(st, X)$ll) +
        .lta_penalty(st, X, layout, alpha)$value
    }
    ana <- colSums(sweep(.lta_score_matrix(fit, X)$S, 1, w, "*")) +
      .lta_penalty(fit, X, layout, alpha)$gradient
    eps <- 1e-5
    fd <- vapply(seq_along(par0), function(i) {
      e <- numeric(length(par0)); e[i] <- eps
      (obj(par0 + e) - obj(par0 - e)) / (2 * eps)
    }, numeric(1))

    expect_lt(max(abs(ana - fd)) / max(1, max(abs(fd))), 1e-5)
  }
})

test_that("the refinement moves Gaussian variances and never lowers the log-likelihood", {
  X <- .lta_gaussian_refine_sim()
  a <- suppressMessages(suppressWarnings(
    fit_lta(X, n_statuses = 3, times = 2, measurement = "continuous",
            smoothing = 0, bayes_constants = .ml, n_init = 2, random_state = 1,
            refine = FALSE, standard_errors = FALSE)))
  b <- suppressMessages(suppressWarnings(
    fit_lta(X, n_statuses = 3, times = 2, measurement = "continuous",
            smoothing = 0, bayes_constants = .ml, n_init = 2, random_state = 1,
            refine = TRUE, standard_errors = TRUE)))

  expect_gte(b$loglik, a$loglik - 1e-8)
  expect_equal(b$n_params, a$n_params)
  # The polish is the only thing separating a and b, and it is the first thing
  # in the package that can move a variance for an LTA fit at all.
  expect_false(isTRUE(all.equal(b$mm$models[[1]]$parameters$covariances,
                                a$mm$models[[1]]$parameters$covariances)))

  expect_false(isTRUE(b$se$conditional))
  var_blocks <- Filter(function(bl) grepl("^log_sd", bl$name), b$se$blocks)
  expect_true(length(var_blocks) > 0)
  var_se <- unlist(lapply(var_blocks, function(bl)
    sqrt(diag(b$se$vcov[bl$cols, bl$cols, drop = FALSE]))))
  expect_true(all(is.finite(var_se)))
})

test_that("the refinement's gradient matches finite differences for covariate delta/tau", {
  # delta_beta/tau_beta replace the simplex-valued delta/tau with case-level
  # regressions; the score is the same (observed - predicted) residual as the
  # plain blocks, now weighted by the covariate row instead of read off a
  # fixed probability vector -- .lta_score_matrix()'s "delta_beta"/"tau_beta"
  # branches. Exercised under the default `transition_effects = "by_origin"`.
  #
  # Uses .lta_cov_refine_sim(), not .lta_refine_sim() with a bolted-on random
  # covariate: delta_beta/tau_beta carry no prior of their own (by design),
  # so fitting an unbounded continuous covariate against data it does not
  # actually predict drives the MLE towards separation -- large coefficients
  # whose softmax saturates against the 1e-300 floor, which is genuinely not
  # differentiable there and fails a finite-difference check for a reason
  # that has nothing to do with whether the score formula is correct. A
  # fixture where the covariate has a real, moderate effect stays identified.
  sim <- .lta_cov_refine_sim()
  X <- sim$X; Z <- sim$Z
  fit <- suppressMessages(suppressWarnings(
    fit_lta(X, n_statuses = 2, times = 3, measurement = "binary",
            predictors_initial = Z, predictors_transition = Z,
            smoothing = 0, bayes_constants = .ml, n_init = 2,
            random_state = 1, refine = FALSE, standard_errors = FALSE)))
  expect_true(.lta_scores_full(fit))

  layout <- .lta_par_layout(fit)
  par0   <- .lta_par_pack(fit, layout)
  w      <- fit$weights_vec

  expect_equal(length(par0), ncol(.lta_score_matrix(fit, X)$S))
  expect_equal(length(par0), fit$n_params)

  rt <- .lta_par_unpack(par0, fit, layout)
  expect_equal(rt$delta_beta, fit$delta_beta, tolerance = 1e-10)
  expect_equal(rt$tau_beta, fit$tau_beta, tolerance = 1e-10)

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

test_that("the same gradient check holds under transition_effects = \"common\"", {
  # A structurally different code path: one tau_beta block per matrix, shared
  # across origin statuses via the origin-dummy design (.lta_tau_design()),
  # rather than one block per origin.
  sim <- .lta_cov_refine_sim()
  X <- sim$X; Z <- sim$Z
  fit <- suppressMessages(suppressWarnings(
    fit_lta(X, n_statuses = 2, times = 3, measurement = "binary",
            predictors_transition = Z, transition_effects = "common",
            smoothing = 0, bayes_constants = .ml, n_init = 2,
            random_state = 1, refine = FALSE, standard_errors = FALSE)))
  expect_true(.lta_scores_full(fit))

  layout <- .lta_par_layout(fit)
  par0   <- .lta_par_pack(fit, layout)
  w      <- fit$weights_vec
  expect_equal(length(par0), fit$n_params)

  obj <- function(p) {
    st <- .lta_par_unpack(p, fit, layout)
    sum(w * .lta_score_matrix(st, X)$ll)
  }
  ana <- colSums(sweep(.lta_score_matrix(fit, X)$S, 1, w, "*"))
  eps <- 1e-5
  fd <- vapply(seq_along(par0), function(i) {
    e <- numeric(length(par0)); e[i] <- eps
    (obj(par0 + e) - obj(par0 - e)) / (2 * eps)
  }, numeric(1))

  expect_lt(max(abs(ana - fd)) / max(1, max(abs(fd))), 1e-5)
})

test_that("the refinement never lowers the log-likelihood for a covariate model", {
  sim <- .lta_cov_refine_sim()
  X <- sim$X; Z <- sim$Z
  a <- suppressMessages(suppressWarnings(
    fit_lta(X, n_statuses = 2, times = 3, measurement = "binary",
            predictors_initial = Z, predictors_transition = Z,
            smoothing = 0, bayes_constants = .ml, n_init = 2,
            random_state = 1, refine = FALSE, standard_errors = FALSE)))
  b <- suppressMessages(suppressWarnings(
    fit_lta(X, n_statuses = 2, times = 3, measurement = "binary",
            predictors_initial = Z, predictors_transition = Z,
            smoothing = 0, bayes_constants = .ml, n_init = 2,
            random_state = 1, refine = TRUE, standard_errors = FALSE)))
  expect_true(isTRUE(b$refined_lbfgs))
  expect_gte(b$loglik, a$loglik - 1e-8)
  expect_equal(b$n_params, a$n_params)
})

test_that("the refinement's gradient matches finite differences for a mixture over chains", {
  # Every score is class-conditional here (a class-membership block, plus
  # delta and tau blocks per class), scaled by the posterior class membership
  # per the Fisher identity -- .lta_score_matrix()'s C > 1 branch. This is the
  # same check test one above runs for a single chain, on a model where the
  # class-mixing block and the per-class delta/tau blocks are all exercised
  # together.
  X <- .lta_refine_sim()
  fit <- suppressMessages(suppressWarnings(
    fit_lta(X, n_statuses = 2, times = 3, measurement = "binary",
            n_classes = 2, smoothing = 0, bayes_constants = .ml, n_init = 2,
            random_state = 1, refine = FALSE, standard_errors = FALSE)))
  expect_true(.lta_scores_full(fit))

  layout <- .lta_par_layout(fit)
  par0   <- .lta_par_pack(fit, layout)
  w      <- fit$weights_vec

  expect_equal(length(par0), ncol(.lta_score_matrix(fit, X)$S))
  expect_equal(length(par0), fit$n_params)

  rt <- .lta_par_unpack(par0, fit, layout)
  expect_equal(rt$delta_c, fit$delta_c, tolerance = 1e-10)
  expect_equal(rt$tau_c, fit$tau_c, tolerance = 1e-10)
  expect_equal(rt$class_weights, fit$class_weights, tolerance = 1e-10)

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

test_that("the refinement improves a mixture-over-chains fit and matches EM at a stationary point", {
  X <- .lta_refine_sim()
  for (seed in 1:3) {
    a <- suppressMessages(suppressWarnings(
      fit_lta(X, n_statuses = 2, times = 3, measurement = "binary",
              n_classes = 2, smoothing = 0, bayes_constants = .ml,
              n_init = 3, random_state = seed, refine = FALSE,
              standard_errors = FALSE)))
    b <- suppressMessages(suppressWarnings(
      fit_lta(X, n_statuses = 2, times = 3, measurement = "binary",
              n_classes = 2, smoothing = 0, bayes_constants = .ml,
              n_init = 3, random_state = seed, refine = TRUE,
              standard_errors = FALSE)))
    expect_gte(b$loglik, a$loglik - 1e-8)
    expect_equal(b$n_params, a$n_params)
  }

  fit <- suppressMessages(suppressWarnings(
    fit_lta(X, n_statuses = 2, times = 3, measurement = "binary",
            n_classes = 2, n_init = 2, random_state = 1, tol = 1e-14,
            max_iter = 20000, refine = FALSE, standard_errors = FALSE)))
  layout <- .lta_par_layout(fit)
  g <- colSums(sweep(.lta_score_matrix(fit, X)$S, 1, fit$weights_vec, "*")) +
    .lta_penalty(fit, X, layout, 1)$gradient
  expect_lt(max(abs(g)), 1e-2)
})

test_that("a refined fit reports the iterations EM actually ran", {
  # The refinement re-runs .lta_em() at max_iter = 0 to write its posteriors
  # back at the polished parameters, and .lta_em() returns `converged` and
  # `n_iter` re-initialised because no iteration ran. Restoring only the first
  # of the two left print() reporting "Converged: TRUE (in 0 iterations)" for
  # every fit the polish improved.
  X <- .lta_refine_sim()
  skeleton <- suppressMessages(suppressWarnings(
    fit_lta(X, n_statuses = 3, times = 3, measurement = "binary",
            smoothing = 0, bayes_constants = .ml, n_init = 1, max_iter = 1L,
            refine = FALSE, standard_errors = FALSE, order_by_size = FALSE)))

  # Stopped short on purpose, so the polish has somewhere to climb.
  set.seed(99)
  loose <- .lta_em(.lta_random_start(skeleton, X), X, max_iter = 8L,
                   tol = 1e-8, alpha = 0)
  polished <- .lta_refine_lbfgs(loose, X, alpha = 0)

  expect_true(isTRUE(polished$refined_lbfgs))
  expect_gt(polished$loglik, loose$loglik)
  expect_identical(polished$n_iter, loose$n_iter)
  expect_identical(polished$converged, loose$converged)
})
