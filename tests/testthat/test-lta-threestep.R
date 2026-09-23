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
  # not reaching the fit.
  a <- .ts3_step3(.ts3_s12())$delta
  b <- .ts3_step3(.ts3_s12("proportional"))$delta
  expect_gt(max(abs(a - b)), 1e-6)
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
