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

test_that("the three-step refuses a tie it has no classes to apply", {
  X <- .ts3_sim()
  expect_error(
    .lta_threestep_step12(
      match.call(fit_lta, quote(fit_lta(
        X, n_statuses = 3, times = 4, measurement = "binary",
        tie_initial_status = TRUE))),
      environment()),
    "applies only to a fit with several latent classes")
})

test_that("the transition-free model refuses what it has not been graded in", {
  X <- .ts3_sim()
  base <- function(...) fit_lta(X, n_statuses = 3, times = 4,
                                measurement = "binary", n_init = 2,
                                random_state = 1, n_cores = 1,
                                .transition_free = TRUE, ...)
  # Covariates are accepted since Part 52 Phase E: a step 1 with
  # `predictors_items` regresses each occasion's status on them
  # (test-lta-dif-items.R).
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

# ==============================================================================
# W3 -- step three: the structural model on the reduced data, with the
# classification error held fixed.
#
# The reduction is the claim being tested. After step 2 the items are gone and
# what is left is one assigned-status variable per occasion plus a K x K table
# per occasion saying how unreliable it is. Step 3 is an ordinary latent
# transition model on exactly that, with the measurement block frozen.
#
#   (f) the emission really is D, in the right orientation, and really is
#       frozen. This is the one mistake in the item that is silent and fatal:
#       planted backwards the fit converges, looks sane and is wrong.
#   (g) the fixed emission costs no parameters, because it was not estimated.
#   (h) the correction is applied in the right DIRECTION. With the error
#       matrices set to the identity the estimator must reduce to the naive
#       classify-analyse baseline -- the raw cross-tabulation of the assigned
#       labels -- exactly. With the real matrices it must move away from that
#       baseline and towards step 1's own prevalences.
#   (i) the proportional reduction's arithmetic.
#   (j) the guard rails.
# ==============================================================================

.ts3_s12 <- function(assignment = "modal") {
  key <- paste0("s12_", assignment)
  if (!is.null(.ts3_cache[[key]])) return(.ts3_cache[[key]])
  X <- .ts3_sim()
  .ts3_cache[[key]] <- mixtureEM:::.lta_threestep_step12(
    match.call(fit_lta, quote(fit_lta(
      X, n_statuses = 3, times = 4, measurement = "binary",
      measurement_invariance = "full", n_init = 8, random_state = 1,
      n_cores = 1, standard_errors = FALSE))),
    environment(), assignment = assignment)
  .ts3_cache[[key]]
}

.ts3_step3 <- function(s12) {
  X <- .ts3_sim()
  mixtureEM:::.lta_threestep_step3(
    s12,
    match.call(fit_lta, quote(fit_lta(
      X, n_statuses = 3, times = 4, n_init = 6, random_state = 1,
      n_cores = 1, smoothing = 0, standard_errors = FALSE))),
    environment())
}

# ------------------------------------------------------------------------------
# (f) the emission is D transposed, and it does not move
# ------------------------------------------------------------------------------

test_that("step 3's emission is the classification error, transposed", {
  # D has rows the ASSIGNED status and columns the TRUE one; an emission has
  # rows the latent status and columns the response. So the planted matrix is
  # t(D), and this assertion is the only thing standing between a typo and a
  # converged fit at the wrong answer.
  s12 <- .ts3_s12()
  fit <- .ts3_step3(s12)

  for (t in seq_len(4L)) {
    planted <- fit$mm$models[[t]]$parameters$pis
    expect_equal(max(abs(planted - t(s12$D[[t]]))), 0, tolerance = 0)
    # Rows are the true status, so they are what sums to one here -- D's
    # columns are what sum to one. If the two ever agreed the table would be
    # symmetric and the orientation would be untestable; it is not.
    expect_lt(max(abs(rowSums(planted) - 1)), 1e-12)
  }
  expect_gt(max(abs(s12$D[[1L]] - t(s12$D[[1L]]))), 1e-3)
})

test_that("frozen means frozen: the emission is the one that went in", {
  # Not a restatement of the test above. That one checks what was planted;
  # this one checks that EM left it alone, which is the whole of what
  # freezing the measurement block buys -- and it is checked against the input
  # rather than against the fit's own copy of it.
  s12  <- .ts3_s12()
  fit  <- .ts3_step3(s12)
  want <- lapply(s12$D, t)
  for (t in seq_len(4L))
    expect_equal(unname(fit$mm$models[[t]]$parameters$pis), unname(want[[t]]),
                 tolerance = 0)
})

# ------------------------------------------------------------------------------
# (g) a constant is not a parameter
# ------------------------------------------------------------------------------

test_that("the fixed emission costs nothing in the parameter count", {
  # Two free initial-status probabilities plus three occasions of a 3x3
  # transition table. The twelve emission numbers were handed to this fit, not
  # estimated by it, so they are not in the count and not in any criterion
  # built on it.
  fit <- .ts3_step3(.ts3_s12())
  expect_equal(fit$n_params, 2L + 3L * 3L * 2L)
  expect_true(isTRUE(fit$mm_fixed))
})

test_that("the two-step's parameter count is untouched by that rule", {
  # The two-step also freezes its measurement block, but there the block WAS
  # estimated -- in step one, on the same data -- so it is counted. The
  # predicate that excludes a fixed emission must not reach this fit.
  X <- .ts3_sim()
  z <- data.frame(v = as.numeric(X[, 1L]))
  fit <- fit_lta(X, n_statuses = 3, times = 4, measurement = "binary",
                 measurement_invariance = "full", n_init = 2, random_state = 1,
                 n_cores = 1, n_steps = 2, predictors_initial = z,
                 standard_errors = FALSE)
  expect_false(isTRUE(fit$mm_fixed))
  expect_gt(fit$n_params, 5L * 3L)
})

# ------------------------------------------------------------------------------
# (h) the correction points the right way
# ------------------------------------------------------------------------------

test_that("an identity error matrix reduces step 3 to classify-analyse", {
  # The naive baseline: believe the labels. With D the identity the model says
  # the assigned status IS the true one, so the maximum-likelihood initial
  # status distribution is just the proportion of cases assigned to each status
  # at occasion 1, and each transition row the corresponding cross-tabulation.
  # Exact, free, and it fixes the DIRECTION of the correction: anything that
  # inverted the table rather than applying it would fail here.
  s12  <- .ts3_s12()
  flat <- s12
  flat$D <- rep(list(diag(3L)), 4L)
  fit <- .ts3_step3(flat)

  W <- s12$modal
  expect_equal(unname(fit$delta), unname(prop.table(tabulate(W[, 1L], 3L))),
               tolerance = 1e-6)
  for (t in seq_len(3L)) {
    tab <- table(factor(W[, t], levels = 1:3), factor(W[, t + 1L], levels = 1:3))
    expect_equal(unname(fit$tau[[t]]), unname(prop.table(as.matrix(tab), 1)),
                 tolerance = 1e-6)
  }
})

test_that("the real error matrices move the estimate towards step 1's own", {
  # The point of the whole estimator. The naive baseline is biased because the
  # labels are wrong some of the time; correcting for how often they are wrong
  # must move the initial-status distribution back towards the one step 1
  # estimated from the items themselves.
  s12  <- .ts3_s12()
  flat <- s12
  flat$D <- rep(list(diag(3L)), 4L)

  target    <- s12$prevalences[1L, ]
  corrected <- .ts3_step3(s12)$delta
  naive     <- .ts3_step3(flat)$delta

  expect_lt(max(abs(corrected - target)), max(abs(naive - target)))
  expect_lt(max(abs(corrected - target)), 0.01)
})

# ------------------------------------------------------------------------------
# (i) the proportional reduction
# ------------------------------------------------------------------------------

test_that("proportional spreads each case over the grid without losing any", {
  # A case is not assigned to one status but to every combination of statuses
  # across the occasions, weighted by the product of its posteriors. The grid
  # is K^T rows however large the sample, and the weights must still add up to
  # the number of cases -- each case's posteriors multiply out to one over it.
  s12 <- .ts3_s12("proportional")
  red <- mixtureEM:::.lta_threestep_reduce(s12)

  expect_lte(nrow(red$W), 3L^4L)
  expect_equal(sum(red$weights), sum(s12$step1$weights_vec), tolerance = 1e-8)
  expect_true(all(red$weights > 0))
  expect_equal(ncol(red$W), 4L)

  # Modal, by contrast, is one row per case.
  expect_equal(nrow(mixtureEM:::.lta_threestep_reduce(.ts3_s12())$W),
               nrow(.ts3_sim()))
})

test_that("the grid is capped rather than allowed to explode", {
  s12 <- .ts3_s12("proportional")
  big <- s12
  big$step1$n_statuses <- 12L           # 12^4 = 20736 cells
  expect_error(mixtureEM:::.lta_threestep_reduce(big), "too many to enumerate")
})

test_that("proportional and modal are different estimators on the same data", {
  # A proportional table is softer than a modal one, so it applies less
  # correction; the two must not come out identical, or one of the rules is
  # not reaching the fit. Read on the transitions: with nothing predicting it,
  # the corrected initial-status distribution is step 1's own under either
  # rule (the error table is built from step 1's posteriors, so the first
  # occasion's margin is reproduced exactly), and taken to the maximum the two
  # agree to 1e-9 there; they differed by 1e-6 only while EM stopped short.
  a <- .ts3_step3(.ts3_s12())$tau
  b <- .ts3_step3(.ts3_s12("proportional"))$tau
  expect_gt(max(abs(unlist(a) - unlist(b))), 1e-3)
})

# ------------------------------------------------------------------------------
# (j) the guard rails
# ------------------------------------------------------------------------------

test_that("the fixed-emission model refuses what it has not been reasoned about", {
  W <- .ts3_s12()$modal
  D <- lapply(.ts3_s12()$D, t)
  base <- function(...) fit_lta(W, n_statuses = 3, times = 4,
                                measurement = "categorical",
                                measurement_invariance = "none",
                                n_init = 2, random_state = 1, n_cores = 1,
                                standard_errors = FALSE,
                                .fixed_emission = D, ...)
  expect_error(base(n_classes = 2), "does not accept")
  expect_error(base(random_intercept = "continuous"), "does not accept")
  expect_error(base(group = rep(1:2, length.out = nrow(W))), "does not accept")
  expect_error(base(.transition_free = TRUE), "does not accept")
  # A structural regression is NOT refused: fitting one is what step 3 is for.
  expect_s3_class(base(predictors_initial =
                         data.frame(v = as.numeric(W[, 1L] == 1L))),
                  "lta_model")
})

test_that("a mis-shaped error matrix is an error, not a fit", {
  # Every one of these is a constant the fit will never move, so a mistake in
  # it cannot show up as a failure to converge. It shows up as a plausible fit
  # at the wrong number, which is why these are errors.
  W <- .ts3_s12()$modal
  D <- lapply(.ts3_s12()$D, t)
  base <- function(dd) fit_lta(W, n_statuses = 3, times = 4,
                               measurement = "categorical",
                               measurement_invariance = "none",
                               n_init = 2, random_state = 1, n_cores = 1,
                               standard_errors = FALSE, .fixed_emission = dd)
  expect_error(base(D[1:3]), "one per occasion")
  expect_error(base(lapply(D, function(m) m[1:2, , drop = FALSE])),
               "ROWS sum to 1")
  expect_error(base(lapply(D, function(m) m * 2)), "ROWS sum to 1")
  # Handed the table without transposing it: its columns sum to one, not its
  # rows, so this is caught rather than fitted. It is only caught because the
  # table is asymmetric -- which is exactly when getting it wrong matters.
  expect_error(base(.ts3_s12()$D), "ROWS sum to 1")
})

test_that("step 3 keeps step 1's status labels", {
  # Sorting by prevalence relabels statuses from most to least common. Step 3's
  # status 2 has to be step 1's status 2, because the fixed table indexes them
  # and cannot be renumbered underneath it.
  s12 <- .ts3_s12()
  fit <- .ts3_step3(s12)
  expect_lt(max(abs(fit$delta - s12$prevalences[1L, ])), 0.01)
})

test_that("step 3 reports standard errors, with the measurement model known", {
  # The uncorrected step-3 standard errors: the measurement model is treated as
  # known, which here it genuinely is -- step 2 computed it. Propagating step
  # one's own uncertainty into them is a separate piece of work.
  s12 <- .ts3_s12()
  X   <- .ts3_sim()
  fit <- mixtureEM:::.lta_threestep_step3(
    s12,
    match.call(fit_lta, quote(fit_lta(
      X, n_statuses = 3, times = 4, n_init = 6, random_state = 1,
      n_cores = 1, smoothing = 0, standard_errors = TRUE))),
    environment())
  expect_true(fit$se$conditional)
  expect_true(all(is.finite(fit$se$prob_se$delta)))
  expect_true(all(fit$se$prob_se$delta > 0))
})

test_that("nothing changes for a fit that plants no emission", {
  # Every branch the flag adds is guarded by it.
  fit <- fit_lta(.ts3_sim(), n_statuses = 3, times = 4, measurement = "binary",
                 measurement_invariance = "full", n_init = 4, random_state = 3,
                 n_cores = 1, standard_errors = FALSE)
  expect_false(isTRUE(fit$mm_fixed))
  expect_equal(fit$n_params, 5L * 3L + 2L + 3L * 6L)
})

# ==============================================================================
# W4 -- covariates on the step-3 model.
#
# Step 3 is an ordinary fit_lta() call, so `predictors_initial` and
# `predictors_transition` reach the regression without anything being built for
# them -- under modal assignment, where the reduced data still has one row per
# case. Under proportional it does not: the sample is collapsed onto the K^T
# grid of status combinations, and a covariate breaks the exchangeability that
# collapse assumes. The remedy is to collapse WITHIN each distinct covariate
# pattern, and the claim being tested is that this is exact rather than an
# approximation.
#
#   (k) the pattern expansion is the case-by-case expansion, summed. Proved by
#       running the same code with one pattern per case and aggregating, rather
#       than by a second hand-written copy of the arithmetic.
#   (l) the covariate reaches the fit, under both rules, and is counted.
# ==============================================================================

# A covariate that genuinely predicts the occasion-1 status, so the regression
# has something to find: occasion 1's item total, dichotomised. Read off the
# fixture rather than drawn, so it costs no second seed.
.ts3_cov <- function() as.integer(rowSums(.ts3_sim()[, 1:5]) >= 3L)

.ts3_step3_cov <- function(s12, se = FALSE) {
  X <- .ts3_sim(); z <- .ts3_cov()
  mixtureEM:::.lta_threestep_step3(
    s12,
    match.call(fit_lta, quote(fit_lta(
      X, n_statuses = 3, times = 4, n_init = 4, random_state = 1,
      n_cores = 1, smoothing = 0, standard_errors = se,
      predictors_initial = z, predictors_transition = z))),
    environment())
}

# ------------------------------------------------------------------------------
# (k) the pattern expansion
# ------------------------------------------------------------------------------

test_that("one pattern per case IS the case-by-case expansion", {
  # The whole design in one assertion. Handing the reduction a covariate that
  # separates every case makes each case its own pattern, which is the
  # expansion carried out case by case -- the thing the collapse is claimed to
  # equal. Summing those weights back over the grid must reproduce the
  # unconditional reduction exactly, and summing them within a two-level
  # covariate's patterns must reproduce that covariate's reduction exactly.
  s12 <- .ts3_s12("proportional")
  n   <- nrow(.ts3_sim())
  z   <- .ts3_cov()

  none <- mixtureEM:::.lta_threestep_reduce(s12)
  each <- mixtureEM:::.lta_threestep_reduce(s12, matrix(seq_len(n), ncol = 1))
  two  <- mixtureEM:::.lta_threestep_reduce(s12, matrix(z, ncol = 1))

  # Every case, every combination it has any weight on.
  expect_gt(nrow(each$W), nrow(two$W))
  expect_equal(sum(each$weights), sum(none$weights), tolerance = 1e-10)

  cell <- function(r) apply(r$W, 1, paste, collapse = "-")
  agg  <- function(r, by) tapply(r$weights, by, sum)

  expect_equal(as.vector(agg(each, cell(each))[cell(none)]), none$weights,
               tolerance = 1e-10)
  expect_equal(
    as.vector(agg(each, paste(cell(each), z[each$Z[, 1]]))[
      paste(cell(two), two$Z[, 1])]),
    two$weights, tolerance = 1e-10)
})

test_that("a constant covariate reproduces the unconditional grid", {
  # One pattern is the unconditional case, so the covariate-aware branch must
  # be bit-for-bit the branch that shipped before it.
  s12 <- .ts3_s12("proportional")
  none <- mixtureEM:::.lta_threestep_reduce(s12)
  one  <- mixtureEM:::.lta_threestep_reduce(
    s12, matrix(1, nrow(.ts3_sim()), 1))

  expect_equal(unname(one$W), unname(none$W))
  expect_equal(one$weights, none$weights)
  expect_null(none$Z)
})

test_that("the expansion carries one covariate row per data row", {
  s12 <- .ts3_s12("proportional")
  red <- mixtureEM:::.lta_threestep_reduce(s12, matrix(.ts3_cov(), ncol = 1))

  expect_equal(nrow(red$Z), nrow(red$W))
  expect_lte(nrow(red$W), 2L * 3L^4L)
  expect_setequal(unique(red$Z[, 1]), c(0, 1))
  expect_equal(sum(red$weights), nrow(.ts3_sim()), tolerance = 1e-8)
  expect_true(all(red$weights > 0))

  # Modal hands its covariates straight back: its rows are still cases.
  modal <- mixtureEM:::.lta_threestep_reduce(
    .ts3_s12(), matrix(.ts3_cov(), ncol = 1))
  expect_equal(nrow(modal$Z), nrow(.ts3_sim()))
  expect_equal(modal$Z[, 1], .ts3_cov())
})

# ------------------------------------------------------------------------------
# (l) the covariate reaches the fit
# ------------------------------------------------------------------------------

test_that("a covariate reaches step 3's initial-status and transition models", {
  fit <- .ts3_step3_cov(.ts3_s12())

  expect_false(is.null(fit$delta_beta))
  expect_equal(dim(fit$delta_beta), c(3L, 2L))       # intercept + one slope
  expect_length(fit$tau_beta, 3L)
  # Intercept, two origin dummies, one slope, under transition_effects="common".
  expect_true(all(vapply(fit$tau_beta, ncol, 1L) == 4L))
  expect_true(fit$converged)

  # 2 x 2 for the initial status, 2 x 4 for each of three transitions, and the
  # fixed emission still costs nothing.
  expect_equal(fit$n_params, 2L * 2L + 3L * (2L * 4L))
  expect_true(isTRUE(fit$mm_fixed))
})

test_that("both assignment rules carry the covariate, and disagree", {
  a <- .ts3_step3_cov(.ts3_s12())
  b <- .ts3_step3_cov(.ts3_s12("proportional"))

  expect_equal(b$n_params, a$n_params)
  expect_true(b$converged)
  expect_gt(max(abs(a$delta_beta - b$delta_beta)), 1e-6)
})

test_that("step 3 reports standard errors for the covariate block", {
  fit <- .ts3_step3_cov(.ts3_s12(), se = TRUE)
  bn  <- vapply(fit$se$blocks, function(x) x$name, "")
  se  <- sqrt(diag(fit$se$vcov))

  expect_true(fit$se$conditional)
  expect_equal(nrow(fit$se$vcov), fit$n_params)
  expect_true("delta_beta" %in% bn)
  cols <- fit$se$blocks[[which(bn == "delta_beta")]]$cols
  expect_true(all(is.finite(se[cols])) && all(se[cols] > 0))
})

# ------------------------------------------------------------------------------
# (m) the pattern start
# ------------------------------------------------------------------------------
# Step 3's emission is known, so a covariate reaches the likelihood only through
# its pattern, and `transition_effects = "by_origin"` on one binary covariate is
# saturated in its two patterns: it IS the unconditional model fitted within
# each. That identity is what the pattern start is built on, so it is tested
# directly, and with `smoothing = 0` it holds only if no pseudo-observation is
# left in the covariate M-step either.

.ts3_groups <- function(s12) {
  z   <- .ts3_cov()
  red <- mixtureEM:::.lta_threestep_reduce(s12, cbind(z, z))
  lapply(0:1, function(v) {
    i <- which(red$Z[, 1] == v)
    fit_lta(red$W[i, , drop = FALSE], n_statuses = 3, times = 4,
            measurement = "categorical", measurement_invariance = "none",
            weights = red$weights[i], weight_type = "frequency",
            .fixed_emission = lapply(s12$D, t), n_init = 4, random_state = 1,
            n_cores = 1, smoothing = 0, standard_errors = FALSE)
  })
}

test_that("by_origin on a binary covariate is its two patterns fitted apart", {
  s12 <- .ts3_s12()
  X <- .ts3_sim(); z <- .ts3_cov()
  fit <- suppressWarnings(mixtureEM:::.lta_threestep_step3(
    s12,
    match.call(fit_lta, quote(fit_lta(
      X, n_statuses = 3, times = 4, n_init = 4, random_state = 1,
      n_cores = 1, smoothing = 0, standard_errors = FALSE,
      predictors_initial = z, predictors_transition = z,
      transition_effects = "by_origin"))),
    environment()))
  groups <- .ts3_groups(s12)
  expect_equal(fit$loglik, sum(vapply(groups, `[[`, 0, "loglik")),
               tolerance = 1e-4)
})

test_that("the pattern start reproduces each pattern's own fit exactly", {
  s12 <- .ts3_s12()
  X <- .ts3_sim(); z <- .ts3_cov()
  cl  <- match.call(fit_lta, quote(fit_lta(
    X, n_statuses = 3, times = 4, n_init = 4, random_state = 1, n_cores = 1,
    smoothing = 0, predictors_initial = z, predictors_transition = z,
    transition_effects = "by_origin")))
  red <- mixtureEM:::.lta_threestep_reduce(s12, cbind(z, z))
  st  <- mixtureEM:::.lta_threestep_pattern_start(s12, red, 1L, 1L, cl,
                                                  environment())
  groups <- .ts3_groups(s12)
  for (v in 0:1) {
    g <- groups[[v + 1L]]
    expect_equal(as.vector(softmax_rows(cbind(1, v) %*% t(st$delta_beta))),
                 g$delta_c[[1]], tolerance = 1e-8)
    # A row nobody in the pattern occupies contributes nothing to its
    # likelihood and is arbitrary in both fits, so only occupied rows compare.
    occ <- g$delta_c[[1]]
    for (t in 1:3) {
      M <- t(vapply(1:3, function(k)
        as.vector(softmax_rows(cbind(1, v) %*% t(st$tau_beta[[t]][[k]]))),
        numeric(3)))
      live <- occ > 1e-8
      expect_equal(M[live, ], g$tau[[t]][live, ], tolerance = 1e-8,
                   ignore_attr = TRUE)
      occ <- as.vector(occ %*% g$tau[[t]])
    }
  }
})

test_that("there is no pattern start without a few covariate patterns", {
  s12 <- .ts3_s12()
  X <- .ts3_sim()
  cl  <- match.call(fit_lta, quote(fit_lta(X, n_statuses = 3, times = 4)))
  red <- mixtureEM:::.lta_threestep_reduce(s12)
  expect_null(mixtureEM:::.lta_threestep_pattern_start(s12, red, 0L, 0L, cl,
                                                       environment()))
  # A continuous covariate: a pattern per case, far more than the cap.
  zc  <- seq_len(nrow(X)) / nrow(X)
  red <- mixtureEM:::.lta_threestep_reduce(s12, cbind(zc))
  expect_null(mixtureEM:::.lta_threestep_pattern_start(s12, red, 1L, 0L, cl,
                                                       environment()))
})

# ==============================================================================
# The user surface: fit_lta(n_steps = 3).
#
#   (k) the public call is the internal pipeline, and nothing else: same step
#       1, same error matrices, same step 3.
#   (l) `correction = "none"` is the naive classify-analyse baseline, which
#       has a closed form under modal assignment.
#   (m) the refusals, each reached before step 1 spends any time.
#   (n) `measurement_invariance = "none"`: step 1's labels are matched across
#       occasions, since its likelihood cannot tell them apart.
# ==============================================================================

.ts3_public <- function(...) {
  fit_lta(.ts3_sim(), n_statuses = 3, times = 4, measurement = "binary",
          n_init = 6, random_state = 1, n_cores = 1, smoothing = 0,
          standard_errors = FALSE, n_steps = 3, ...)
}

test_that("the public three-step is the internal pipeline", {
  X   <- .ts3_sim()
  cl  <- match.call(fit_lta, quote(fit_lta(
    X, n_statuses = 3, times = 4, measurement = "binary", n_init = 6,
    random_state = 1, n_cores = 1, smoothing = 0, standard_errors = FALSE)))
  s12 <- mixtureEM:::.lta_threestep_step12(cl, environment(),
                                           assignment = "modal")
  ref <- mixtureEM:::.lta_threestep_step3(s12, cl, environment())
  fit <- .ts3_public(assignment = "modal")

  expect_identical(fit$n_steps, 3L)
  expect_equal(fit$loglik, ref$loglik, tolerance = 1e-10)
  expect_equal(fit$delta, ref$delta, tolerance = 1e-10)
  expect_equal(fit$threestep$classification_error, s12$D)
  expect_equal(fit$step1$loglik, s12$step1$loglik)
  expect_identical(fit$threestep$correction, "ML")
  expect_identical(fit$threestep$assignment, "modal")
})

test_that("proportional is the default assignment", {
  expect_identical(.ts3_public()$threestep$assignment, "proportional")
})

test_that("correction = \"none\" is the naive cross-tabulation of the labels", {
  fit <- .ts3_public(assignment = "modal", correction = "none")
  for (t in seq_len(4L))
    expect_equal(unname(fit$threestep$classification_error[[t]]), diag(3L))

  W <- fit$threestep$modal
  expect_equal(unname(fit$delta), prop.table(tabulate(W[, 1L], 3L)),
               tolerance = 1e-6)
  for (t in seq_len(3L)) {
    tab <- table(factor(W[, t], levels = 1:3), factor(W[, t + 1L], levels = 1:3))
    expect_equal(unname(fit$tau[[t]]), unname(prop.table(as.matrix(tab), 1)),
                 tolerance = 1e-6)
  }
})

test_that("the three-step refuses what it is not defined for, before fitting", {
  X <- .ts3_sim()
  z <- rep(0:1, length.out = nrow(X))
  base <- list(X, n_statuses = 3, times = 4, measurement = "binary",
               n_init = 1, n_cores = 1, n_steps = 3)
  refused <- list(list(random_intercept = "continuous"),
                  list(n_classes = 2),
                  list(mover_stayer = TRUE),
                  list(group = z),
                  list(predictors_random_intercept = z))
  for (a in refused)
    expect_error(do.call(fit_lta, c(base, a)), "`n_steps = 3` does not support")
  # BCH: modal assignment only, and not yet with a distal or a design.
  expect_error(do.call(fit_lta, c(base, list(correction = "BCH",
                                             assignment = "proportional"))),
               "needs `assignment = \"modal\"`")
  expect_error(do.call(fit_lta, c(base, list(correction = "BCH",
                                             distal = data.frame(y = z)))),
               "does not yet support distal")
  expect_error(do.call(fit_lta, c(base, list(correction = "BCH", strata = z))),
               "does not yet support strata or cluster")
  # A design is carried under modal assignment only (W6b, below).
  expect_error(do.call(fit_lta, c(base, list(strata = z))),
               "only with `assignment = \"modal\"`")
  expect_error(do.call(fit_lta, c(base, list(cluster = seq_along(z),
                                             assignment = "proportional"))),
               "only with `assignment = \"modal\"`")
})

test_that("a non-invariant step 1 has its labels matched across occasions", {
  # Scramble one occasion's labels by hand -- which leaves the likelihood
  # exactly where it was -- and the matching must undo it.
  s1 <- fit_lta(.ts3_sim(), n_statuses = 3, times = 4, measurement = "binary",
                measurement_invariance = "none", n_init = 4, random_state = 1,
                n_cores = 1, standard_errors = FALSE, .transition_free = TRUE)
  s1 <- mixtureEM:::.lta_threestep_align(s1)
  pis <- lapply(s1$mm$models, function(m) m$parameters$pis)

  p  <- c(2L, 3L, 1L)
  sc <- s1
  sc$mm$models[[3]] <- mixtureEM:::.permute_emission_classes(sc$mm$models[[3]], p)
  sc$gamma[[3]] <- sc$gamma[[3]][, p]
  sc$tau[[2]]   <- sc$tau[[2]][, p]
  sc$tau[[3]]   <- sc$tau[[3]][p, ]
  sc$tau_c[[1]] <- sc$tau
  back <- mixtureEM:::.lta_threestep_align(sc)

  expect_identical(back$alignment[[3]], order(p))
  expect_equal(back$mm$models[[3]]$parameters$pis, pis[[3]])
  expect_equal(back$gamma[[3]], s1$gamma[[3]])
  expect_equal(back$prevalences, s1$prevalences)
  # The fixture's statuses are high, middle and low at every occasion, so once
  # matched they rank the same way everywhere. Not a closeness bound: the
  # middle status answers every item at .5 and is loosely estimated when each
  # occasion is fitted on its own, so its profile wanders by up to .5.
  for (t in 2:4)
    expect_identical(order(rowMeans(pis[[t]])), order(rowMeans(pis[[1]])))
})

test_that("measurement_invariance = \"none\" runs end to end", {
  fit <- .ts3_public(measurement_invariance = "none", assignment = "modal")
  expect_identical(fit$n_steps, 3L)
  expect_length(fit$threestep$alignment, 4L)
  expect_true(fit$converged)
  # With free time-varying transitions step 3's likelihood does not depend on
  # the labels, so the matched and the unmatched fits are the same model.
  expect_equal(fit$n_params, 20L)
})

# ==============================================================================
# W6 -- standard errors that know step 1 was estimated.
#
# The corrected variance is V2 + V2 H21' S11 H21 V2. No external program
# prints it for this model, so each of its pieces is checked on its own terms:
#
#   S11   against the Hessian of a transition-free likelihood written from
#         scratch here, which shares no code with the package's packing;
#   H21   through what it means rather than how it is computed: V2 H21' is the
#         slope of step 3's estimate in step 1's parameters, so step 3 is
#         refitted at a nudged step 1 and the two slopes compared;
#   V     is V2 plus a positive semi-definite term, so no standard error can
#         shrink, and with nothing to propagate it is not applied at all.
# ==============================================================================

.ts3_w6 <- function(assignment = "modal", correction = "ML") {
  key <- paste("w6", assignment, correction)
  if (!is.null(.ts3_cache[[key]])) return(.ts3_cache[[key]])
  .ts3_cache[[key]] <- fit_lta(
    .ts3_sim(), n_statuses = 3, times = 4, measurement = "binary",
    n_init = 6, random_state = 1, n_cores = 1, smoothing = 0, n_steps = 3,
    assignment = assignment, correction = correction)
}

.ts3_w6_s12 <- function(fit) {
  list(step1 = fit$step1, assignment = fit$threestep$assignment,
       correction = fit$threestep$correction,
       zero_floor = fit$threestep$zero_floor, modal = fit$threestep$modal,
       D = fit$threestep$classification_error)
}

test_that("step 1's variance is the transition-free model's own", {
  fit <- .ts3_w6()
  s1  <- fit$step1
  X   <- .ts3_sim()
  K <- 3; J <- 5; Tn <- 4
  # Item logits (status fastest, item slowest), then each occasion's
  # prevalence logits against the last status.
  th0 <- c(qlogis(s1$mm$models[[1]]$parameters$pis),
           unlist(lapply(seq_len(Tn), function(t) {
             p <- s1$prevalences[t, ]
             log(p[1:2] / p[3])
           })))
  ll <- function(th) {
    rho <- matrix(plogis(th[seq_len(K * J)]), K, J)
    sum(vapply(seq_len(Tn), function(t) {
      e  <- exp(c(th[K * J + (t - 1) * 2 + 1:2], 0))
      Xt <- X[, (t - 1) * J + seq_len(J)]
      lk <- vapply(seq_len(K), function(k)
        e[k] / sum(e) * exp(Xt %*% log(rho[k, ]) + (1 - Xt) %*% log(1 - rho[k, ])),
        numeric(nrow(X)))
      sum(log(rowSums(lk)))
    }, numeric(1)))
  }
  expect_equal(ll(th0), s1$loglik, tolerance = 1e-8)
  se_scratch <- sqrt(diag(solve(-optimHess(th0, ll))))
  # The package orders the prevalences first.
  se_ours <- sqrt(diag(fit$se$step1_vcov))
  expect_length(se_ours, 23L)
  expect_lt(max(abs(se_ours / c(tail(se_scratch, 8), head(se_scratch, 15)) - 1)),
            5e-3)
})

test_that("the correction only ever adds, and says it has been applied", {
  fit <- .ts3_w6()
  V   <- diag(fit$se$vcov)
  V2  <- diag(fit$se$threestep_V2)
  expect_true(isTRUE(fit$se$threestep))
  expect_false(fit$se$conditional)
  expect_true(all(V >= V2 - 1e-12))
  expect_gt(max(V / V2), 1.1)
  expect_true(all(is.finite(unlist(fit$se$prob_se))))
  # At step 1's own point, what the variance reads off step 1 is exactly what
  # step 3 was fitted with.
  pt <- mixtureEM:::.lta_threestep_vcov(fit, .ts3_w6_s12(fit), parts = TRUE)
  side <- pt$step1_side(pt$r0)
  for (t in 1:4)
    expect_equal(side$emission[[t]], t(fit$threestep$classification_error[[t]]),
                 tolerance = 1e-10)
})

test_that("the propagation term is the slope of step 3's estimate in step 1's", {
  fit <- .ts3_w6()
  pt  <- mixtureEM:::.lta_threestep_vcov(fit, .ts3_w6_s12(fit), parts = TRUE)
  slope <- pt$V2 %*% t(pt$H12)          # d(step 3 estimate) / d(step 1 point)
  a   <- which.max(colSums(slope^2))
  eps <- 1e-3
  refit <- function(s) {
    r <- pt$r0
    r[a] <- r[a] + s * eps
    side <- pt$step1_side(r)
    f <- fit_lta(fit$data, n_statuses = 3, times = 4,
                 measurement = "categorical", measurement_invariance = "none",
                 weights = side$weights, weight_type = "frequency",
                 .fixed_emission = side$emission, n_init = 6, random_state = 1,
                 n_cores = 1, smoothing = 0, tol = 1e-12, max_iter = 20000,
                 standard_errors = FALSE)
    mixtureEM:::.lta_par_pack(f, pt$lay3)
  }
  num <- (refit(1) - refit(-1)) / (2 * eps)
  expect_lt(max(abs(num - slope[, a])) / max(abs(slope[, a])), 0.02)
})

test_that("proportional assignment propagates through the weights too", {
  fit <- .ts3_w6("proportional")
  expect_true(isTRUE(fit$se$threestep))
  expect_true(all(diag(fit$se$vcov) >= diag(fit$se$threestep_V2) - 1e-12))
})

test_that("modal with no correction has nothing to propagate", {
  # Every table is the identity and every weight a case weight, so step 1's
  # parameters never reach step 3 and the variance is left as it was.
  fit <- .ts3_w6("modal", "none")
  expect_null(fit$se$threestep)
  expect_true(fit$se$conditional)
})

# ==============================================================================
# W6b -- strata and cluster under modal assignment: the stacked sandwich.
#
# Graded against printed reference numbers in the validation suite; here, two
# properties of the estimator itself, on fits with the item prior off so that
# doubling the data doubles the likelihood exactly.
#   (o) every case its own PSU, one stratum: the design variance is the robust
#       one, which agrees with the model-based variance to sampling error;
#   (p) every case duplicated inside its own PSU: a copy adds no information,
#       so every standard error is where it was, to the finite-difference
#       noise of step 1's Hessian. That noise is ~0.02 per element at
#       |LL| ~ 1e4 (rounding over h^2 = 1e-10) and moves these SEs by up to
#       1.2%. The 3.4e-3 once recorded here measured two bitwise-identical
#       EM iterates, whose rounding errors were the same; the Newton-type
#       finish leaves the two fits 5e-11 apart, which is enough to show it.
# ==============================================================================

.ts3_w6b <- function(X, ...) {
  fit_lta(X, n_statuses = 3, times = 4, measurement = "binary", n_init = 6,
          random_state = 1, n_cores = 1, smoothing = 0, n_steps = 3,
          assignment = "modal", tol = 1e-10, max_iter = 20000,
          bayes_constants = list(categorical = 0), ...)
}

test_that("with every case its own PSU the design variance is the model's, roughly", {
  X  <- .ts3_sim()
  m  <- .ts3_w6b(X)
  d  <- .ts3_w6b(X, cluster = seq_len(nrow(X)))
  expect_true(d$has_survey_design)
  expect_match(d$se$method, "survey-linearized sandwich over both steps")
  expect_equal(d$delta, m$delta, tolerance = 1e-8)
  r <- sqrt(diag(d$se$vcov) / diag(m$se$vcov))
  expect_gt(min(r), 0.7)
  expect_lt(max(r), 1.5)
  .ts3_cache$w6b_own <- d
})

test_that("duplicating every case inside its own PSU changes no standard error", {
  X  <- .ts3_sim()
  n  <- nrow(X)
  d1 <- .ts3_cache$w6b_own %||% .ts3_w6b(X, cluster = seq_len(n))
  d2 <- .ts3_w6b(X[rep(seq_len(n), each = 2), ],
                 cluster = rep(seq_len(n), each = 2))
  expect_equal(d2$delta, d1$delta, tolerance = 1e-8)
  expect_lt(max(abs(sqrt(diag(d2$se$vcov) / diag(d1$se$vcov)) - 1)), 0.02)
})

# ==============================================================================
# The BCH correction: product weights, step 3 fitted to the weighted paths.
# Graded in the validation suite; here, that it recovers a simulated chain the
# uncorrected estimate gets wrong, and that its weights are what the
# construction says.
# ==============================================================================

.ts3_chain <- function() {
  if (!is.null(.ts3_cache$chain)) return(.ts3_cache$chain)
  set.seed(11)
  n <- 2000; J <- 5
  tau <- rbind(c(.8, .2), c(.2, .8))
  s1 <- sample.int(2, n, TRUE, prob = c(.6, .4))
  s2 <- vapply(s1, function(s) sample.int(2, 1, prob = tau[s, ]), 1L)
  # Items weak enough that one case in six or so is misclassified at each
  # occasion, which is what attenuates the uncorrected transitions.
  rho <- rbind(rep(.75, J), rep(.25, J))
  X <- cbind(matrix(rbinom(n * J, 1, rho[s1, ]), n),
             matrix(rbinom(n * J, 1, rho[s2, ]), n))
  fit <- function(correction)
    fit_lta(X, n_statuses = 2, times = 2, measurement = "binary",
            n_init = 4, random_state = 1, n_cores = 1, n_steps = 3,
            correction = correction, assignment = "modal")
  .ts3_cache$chain <- list(tau = tau, bch = fit("BCH"), none = fit("none"),
                           ml = fit("ML"))
  .ts3_cache$chain
}

test_that("BCH recovers the transitions the uncorrected estimate attenuates", {
  # Measured: stay probabilities .842/.823 BCH, .841/.822 ML, .727/.671
  # uncorrected, against .8/.8.
  ch  <- .ts3_chain()
  bch <- diag(ch$bch$tau[[1]])
  raw <- diag(ch$none$tau[[1]])
  expect_lt(max(abs(bch - diag(ch$tau))), 0.06)
  expect_gt(min(abs(raw - diag(ch$tau))), 0.06)
  # Two corrections built on different principles land on the same answer.
  expect_lt(max(abs(ch$bch$tau[[1]] - ch$ml$tau[[1]])), 0.005)
  expect_equal(ch$bch$threestep$correction, "BCH")
  expect_equal(ch$bch$threestep$assignment, "modal")
})

test_that("each case's path weights multiply out and sum to one", {
  f <- .ts3_chain()$bch
  W <- f$threestep$bch_weights
  expect_equal(dim(W), c(nrow(f$threestep$modal), 4L))
  expect_equal(unname(rowSums(W)), rep(1, nrow(W)), tolerance = 1e-10)
  expect_true(any(W < 0))
  # Column (s1, s2), occasion 1 fastest, is occasion 1's weight times
  # occasion 2's.
  B <- lapply(1:2, function(t)
    solve(t(f$threestep$classification_error[[t]]))[f$threestep$modal[, t], ])
  expect_equal(W[, 3], B[[1]][, 1] * B[[2]][, 2], tolerance = 1e-12)
})

test_that("BCH reports the sandwich, larger than the uncorrected Hessian", {
  ch <- .ts3_chain()
  expect_match(ch$bch$se$method, "Case-clustered sandwich of the BCH-weighted")
  expect_true(all(is.finite(diag(ch$bch$se$vcov))))
  expect_true(all(diag(ch$bch$se$vcov) > diag(ch$none$se$vcov)))
})
