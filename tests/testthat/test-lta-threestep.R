# ==============================================================================
# The bias-adjusted three-step estimator for latent transition models, steps
# one and two: .lta_threestep_step12() and the transition-free measurement
# model it fits.
#
# Step 1 is the measurement model equated across occasions with NO transition
# structure and each occasion's status prevalences free. Everything checked
# here follows from that one sentence:
#
#   (a) the constraint bit. Every origin row of an occasion's transition
#       matrix is the same vector, and the parameter count is one free
#       prevalence vector per occasion rather than K rows of them. Nothing
#       else in the file means anything if this one fails.
#   (b) the constraint is the model it claims to be. With no transitions the
#       joint log-likelihood FACTORISES: it equals the sum of the T
#       single-occasion latent class log-likelihoods computed by hand from the
#       fitted parameters. That is the check that step 1 is one model and not
#       T of them glued together, and it is exact rather than approximate.
#   (c) the posteriors are unsmoothed. Occasion t's posterior is a function of
#       occasion t's items and occasion t's base rates alone, which is what
#       makes the assigned status a lone indicator of the true status. A
#       transition-bearing fit's posteriors are not, and the two differ.
#   (d) step 2's tables have the orientation and the normalisation the fixed
#       logits are read off -- rows the assigned status, COLUMNS the true one,
#       columns summing to one.
#   (e) the guard rails.
#
# The fixture is a four-occasion binary panel with well separated items and
# prevalences that move a long way across occasions, because a classification
# error matrix is a function of the base rate and a fixture whose base rate
# stands still cannot tell a per-occasion matrix from a pooled one.
# ==============================================================================

.ts3_cache <- new.env(parent = emptyenv())

.ts3_sim <- function() {
  if (!is.null(.ts3_cache$X)) return(.ts3_cache$X)
  set.seed(7)
  K <- 3; Tn <- 4; J <- 5; n <- 800
  rho <- matrix(c(.90, .80, .85, .75, .90,
                  .50, .50, .50, .50, .50,
                  .10, .20, .15, .25, .10), nrow = K, byrow = TRUE)
  prev <- rbind(c(.6, .3, .1), c(.4, .4, .2), c(.3, .4, .3), c(.2, .3, .5))
  X <- matrix(NA_integer_, n, J * Tn)
  for (t in seq_len(Tn)) {
    s <- sample.int(K, n, TRUE, prob = prev[t, ])
    for (j in seq_len(J))
      X[, (t - 1L) * J + j] <- rbinom(n, 1, rho[s, j])
  }
  .ts3_cache$X <- X
  X
}

.ts3_fit <- function() {
  if (!is.null(.ts3_cache$fit)) return(.ts3_cache$fit)
  .ts3_cache$fit <- fit_lta(
    .ts3_sim(), n_statuses = 3, times = 4, measurement = "binary",
    measurement_invariance = "full", n_init = 8, random_state = 1,
    n_cores = 1, .transition_free = TRUE)
  .ts3_cache$fit
}

# ------------------------------------------------------------------------------
# (a) the constraint bit
# ------------------------------------------------------------------------------

test_that("the transition-free fit ties every origin row to one vector", {
  fit <- .ts3_fit()
  for (t in seq_along(fit$tau)) {
    m <- fit$tau[[t]]
    expect_equal(max(abs(sweep(m, 2, m[1L, ]))), 0, tolerance = 0)
  }
})

test_that("each occasion costs one prevalence vector, not K rows of them", {
  fit <- .ts3_fit()
  # 5 items x 3 statuses equated across occasions, plus 4 occasions of two
  # free prevalences. The unconstrained model of the same shape counts 35.
  expect_equal(fit$n_params, 5L * 3L + 4L * 2L)
})

# ------------------------------------------------------------------------------
# (b) the constraint is the model it claims to be
# ------------------------------------------------------------------------------

test_that("the joint log-likelihood factorises over the occasions", {
  fit <- .ts3_fit()
  X   <- .ts3_sim()
  K   <- 3L; J <- 5L; Tn <- 4L
  pis <- fit$mm$models[[1L]]$parameters$pis
  # Occasion 1's prevalences are `delta`; the rest are any row of the matrix
  # that ends at that occasion, the rows all being equal by (a).
  prev <- rbind(fit$delta,
                t(vapply(fit$tau, function(m) m[1L, ], numeric(K))))

  lls <- vapply(seq_len(Tn), function(t) {
    Xt  <- X[, (t - 1L) * J + seq_len(J), drop = FALSE]
    lik <- matrix(0, nrow(Xt), K)
    for (k in seq_len(K))
      lik[, k] <- prev[t, k] *
        exp(Xt %*% log(pis[k, ]) + (1 - Xt) %*% log(1 - pis[k, ]))
    sum(log(rowSums(lik)))
  }, numeric(1))

  expect_lt(abs(sum(lls) - fit$loglik), 1e-8)
})

# ------------------------------------------------------------------------------
# (c) the posteriors are unsmoothed
# ------------------------------------------------------------------------------

test_that("an occasion's posterior reads that occasion's items alone", {
  fit <- .ts3_fit()
  X   <- .ts3_sim()
  K   <- 3L; J <- 5L
  pis <- fit$mm$models[[1L]]$parameters$pis
  prev <- rbind(fit$delta,
                t(vapply(fit$tau, function(m) m[1L, ], numeric(K))))

  for (t in seq_len(4L)) {
    Xt  <- X[, (t - 1L) * J + seq_len(J), drop = FALSE]
    lik <- matrix(0, nrow(Xt), K)
    for (k in seq_len(K))
      lik[, k] <- prev[t, k] *
        exp(Xt %*% log(pis[k, ]) + (1 - Xt) %*% log(1 - pis[k, ]))
    own <- lik / rowSums(lik)
    expect_equal(unname(fit$gamma[[t]]), unname(own), tolerance = 1e-10)
  }
})

test_that("a transition-bearing fit's posteriors are NOT the same object", {
  # The point of dropping the transitions: with them in place occasion 2's
  # posterior borrows from occasions 1 and 3 through the chain, and the
  # assigned status stops being a lone indicator of the true status.
  free <- .ts3_fit()
  with_tr <- fit_lta(.ts3_sim(), n_statuses = 3, times = 4,
                     measurement = "binary", measurement_invariance = "full",
                     n_init = 8, random_state = 1, n_cores = 1,
                     standard_errors = FALSE)
  expect_gt(max(abs(free$gamma[[2L]] - with_tr$gamma[[2L]])), 1e-3)
})

# ------------------------------------------------------------------------------
# (d) step 2's tables
# ------------------------------------------------------------------------------

.ts3_step12 <- function(assignment = "modal") {
  key <- paste0("s12_", assignment)
  if (!is.null(.ts3_cache[[key]])) return(.ts3_cache[[key]])
  X <- .ts3_sim()
  .ts3_cache[[key]] <- .lta_threestep_step12(
    match.call(fit_lta, quote(fit_lta(
      X, n_statuses = 3, times = 4, measurement = "binary",
      measurement_invariance = "full", n_init = 8, random_state = 1,
      n_cores = 1))),
    environment(), assignment = assignment)
  .ts3_cache[[key]]
}

test_that("step 2 returns one error matrix per occasion, columns summing to 1", {
  res <- .ts3_step12()
  expect_length(res$D, 4L)
  for (t in seq_len(4L)) {
    expect_equal(dim(res$D[[t]]), c(3L, 3L))
    expect_equal(unname(colSums(res$D[[t]])), rep(1, 3L), tolerance = 1e-12)
    expect_true(all(res$D[[t]] > 0))
  }
})

test_that("the matrices differ across occasions when the base rate moves", {
  # A pooled matrix would make these identical. The fixture's prevalences run
  # from .6/.3/.1 to .2/.3/.5, so they must not be.
  res <- .ts3_step12()
  expect_gt(max(abs(res$D[[1L]] - res$D[[4L]])), 1e-2)
})

test_that("the fixed logits are the log-ratio of D against its last row", {
  res <- .ts3_step12()
  for (t in seq_len(4L)) {
    expect_equal(unname(res$logits[[t]]),
                 unname(log(sweep(res$D[[t]], 2, res$D[[t]][3L, ], "/"))),
                 tolerance = 1e-12)
    expect_equal(unname(res$logits[[t]][3L, ]), rep(0, 3L), tolerance = 1e-12)
  }
})

test_that("step 1 is returned whole, with its prevalences by occasion", {
  res <- .ts3_step12()
  expect_s3_class(res$step1, "lta_model")
  expect_equal(res$step1$loglik, .ts3_fit()$loglik, tolerance = 1e-10)
  expect_equal(dim(res$prevalences), c(4L, 3L))
  expect_equal(unname(rowSums(res$prevalences)), rep(1, 4L), tolerance = 1e-10)
  expect_equal(unname(res$prevalences[1L, ]), unname(res$step1$delta),
               tolerance = 1e-12)
  expect_equal(dim(res$modal), c(nrow(.ts3_sim()), 4L))
  expect_true(all(res$modal %in% 1:3))
})

test_that("proportional assignment gives a softer table than modal", {
  m <- .ts3_step12("modal")
  p <- .ts3_step12("proportional")
  expect_equal(p$assignment, "proportional")
  # The diagonal of a proportional table is the posterior mass a status keeps
  # rather than the count of cases sent to it, so it cannot be larger.
  for (t in seq_len(4L))
    expect_true(all(diag(p$D[[t]]) <= diag(m$D[[t]]) + 1e-12))
})

# ------------------------------------------------------------------------------
# (e) the guard rails
# ------------------------------------------------------------------------------

test_that("step 1 refuses a donor fit alongside its own search", {
  X <- .ts3_sim()
  expect_error(
    .lta_threestep_step12(
      match.call(fit_lta, quote(fit_lta(
        X, n_statuses = 3, times = 4, measurement = "binary",
        refine_from = donor))),
      environment()),
    "refine_from")
})

test_that("the transition-free model refuses what it has not been graded in", {
  X <- .ts3_sim()
  base <- function(...) fit_lta(X, n_statuses = 3, times = 4,
                                measurement = "binary", n_init = 2,
                                random_state = 1, n_cores = 1,
                                .transition_free = TRUE, ...)
  z <- data.frame(v = rnorm(nrow(X)))
  expect_error(base(predictors_initial = z), "does not accept")
  expect_error(base(predictors_transition = z), "does not accept")
  expect_error(base(n_classes = 2), "does not accept")
  expect_error(base(random_intercept = "continuous"), "does not accept")
  expect_error(base(forbidden_transitions = matrix(c(0, 0, 0, 1, 0, 0, 0, 0, 0), 3, 3, byrow = TRUE)), "does not accept")
})

test_that("the transition-free fit reports no standard errors", {
  # The packed parameter vector still describes K free transition rows, so a
  # Jacobian read off it would describe a model that was not fitted. The fit
  # reports nothing rather than something wrong; carrying step-1 uncertainty
  # into the three-step's standard errors is a separate piece of work.
  fit <- fit_lta(.ts3_sim(), n_statuses = 3, times = 4, measurement = "binary",
                 measurement_invariance = "full", n_init = 2, random_state = 1,
                 n_cores = 1, .transition_free = TRUE, standard_errors = TRUE)
  expect_null(fit$se)
})

test_that("nothing changes for an ordinary fit", {
  # Every branch the flag adds is guarded by it, so a fit that does not set it
  # is the fit it was before, to the last bit.
  fit <- fit_lta(.ts3_sim(), n_statuses = 3, times = 4, measurement = "binary",
                 measurement_invariance = "full", n_init = 4, random_state = 3,
                 n_cores = 1, standard_errors = FALSE)
  expect_false(isTRUE(fit$tau_independent))
  expect_equal(fit$n_params, 5L * 3L + 2L + 3L * 6L)
  expect_gt(max(abs(sweep(fit$tau[[1L]], 2, fit$tau[[1L]][1L, ]))), 1e-3)
})
