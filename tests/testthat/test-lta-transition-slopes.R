# ==============================================================================
# transition_invariance = "slopes": transition slopes shared across occasions,
# transition intercepts free per occasion.
#
# The whole feature is Tn - 2 occasion-contrast columns appended to the
# transition design, so the checks here are mostly structural and need no
# external number. Two of them are load-bearing:
#
#   * The Tn = 2 identity. With two occasions there are no contrasts to add, so
#     "slopes" and "full" must be the SAME model, and the fits must agree
#     bit for bit rather than merely closely. That is what says the appended
#     block is genuinely zero-width, and therefore that nothing already shipped
#     can have moved.
#   * The per-occasion transition matrices must actually differ.
#     .lta_par_unpack() used to broadcast occasion 1's matrix over the rest
#     whenever the transitions were homogeneous, which is harmless when every
#     occasion is identical by construction and silently destroys occasions
#     2..Tn-1 once they are not.
#
# The fixture is .lta_cov_refine_sim() from helper-lta-packing.R -- a covariate
# with a real, moderate effect on both the initial status and the transitions.
# A noise covariate would drive the unpenalised coefficients towards separation
# and fail the finite-difference check for a reason that has nothing to do with
# the score formula; the same reasoning is written out at length in
# test-lta-refine.R.
#
# On the fit budget. Every check here except the nesting one is structural --
# a parameter count, a design column, a round trip -- and none of them needs a
# converged optimum, so they run at a hard `max_iter` cap and share their fits
# through the memo below. The nesting check is the exception and pays for
# convergence: measured on this fixture, a 150-iteration cap puts "none" BELOW
# "slopes", which is an artefact of stopping early and not a violation of
# anything. Capping everything cost 12.3 minutes' worth of the public suite;
# this arrangement costs about three.
# ==============================================================================

.slp_cache <- new.env(parent = emptyenv())

.slp_sim <- function(Tn) {
  key <- paste0("sim", Tn)
  if (is.null(.slp_cache[[key]])) .slp_cache[[key]] <- .lta_cov_refine_sim(Tn = Tn)
  .slp_cache[[key]]
}

# `max_iter = 150` is the structural default; pass 1000 where the optimum
# itself is the thing under test.
.slp_fit <- function(ti, Tn, max_iter = 150L, covariates = TRUE) {
  key <- paste(ti, Tn, max_iter, covariates, sep = "/")
  if (!is.null(.slp_cache[[key]])) return(.slp_cache[[key]])
  sim <- .slp_sim(Tn)
  args <- list(sim$X, n_statuses = 2, times = Tn, measurement = "binary",
               transition_invariance = ti, n_init = 3, n_cores = 1,
               random_state = 5, max_iter = max_iter, refine = FALSE,
               standard_errors = FALSE)
  if (covariates)
    args <- c(args, list(predictors_initial = sim$Z,
                         predictors_transition = sim$Z))
  .slp_cache[[key]] <- suppressMessages(suppressWarnings(do.call(fit_lta, args)))
  .slp_cache[[key]]
}

test_that("with two occasions \"slopes\" and \"full\" are the same model", {
  f_full <- .slp_fit("full", 2L)
  f_slp  <- .slp_fit("slopes", 2L)

  # identical(), not expect_equal(): there is nothing here to be approximately
  # right about. A difference of any size means the design gained a column.
  expect_true(identical(f_full$loglik, f_slp$loglik))
  expect_identical(f_full$n_params, f_slp$n_params)
  expect_identical(ncol(.lta_tau_design(f_slp, 1L, 1L)),
                   ncol(.lta_tau_design(f_full, 1L, 1L)))
})

test_that("\"slopes\" costs exactly (K - 1) * (Tn - 2) parameters over \"full\"", {
  for (Tn in c(3L, 4L)) {
    mi <- if (Tn == 3L) 1000L else 150L   # the Tn = 3 trio is shared with the
    f_full <- .slp_fit("full", Tn, mi)    # nesting test, which needs 1000
    f_slp  <- .slp_fit("slopes", Tn, mi)
    expect_identical(f_slp$n_params, f_full$n_params + (2L - 1L) * (Tn - 2L))
    # The packed vector and the score matrix have to agree with the count.
    expect_equal(length(.lta_par_pack(f_slp, .lta_par_layout(f_slp))),
                 f_slp$n_params)
    expect_equal(ncol(.lta_score_matrix(f_slp, .slp_sim(Tn)$X)$S),
                 f_slp$n_params)
  }
})

test_that("\"slopes\" is nested between \"full\" and \"none\"", {
  f_full <- .slp_fit("full", 3L, 1000L)
  f_slp  <- .slp_fit("slopes", 3L, 1000L)
  f_none <- .slp_fit("none", 3L, 1000L)

  expect_lte(f_full$loglik, f_slp$loglik + 1e-6)
  expect_lte(f_slp$loglik,  f_none$loglik + 1e-6)

  # Both restrictions are nested, so lr_test() tests them.
  t1 <- lr_test(f_full, f_slp)
  t2 <- lr_test(f_slp, f_none)
  expect_true(is.finite(t1$df) && t1$df > 0)
  expect_true(is.finite(t2$df) && t2$df > 0)
  expect_equal(t1$df, f_slp$n_params - f_full$n_params)
  expect_equal(t2$df, f_none$n_params - f_slp$n_params)
})

test_that("the occasion contrasts are on the right rows", {
  fit <- .slp_fit("slopes", 4L)

  # Tn = 4 gives contrasts for occasions 2 and 3; occasion 1 is the reference
  # and is all zero in both columns.
  d1 <- .lta_tau_design(fit, 1L, 1L)
  d2 <- .lta_tau_design(fit, 1L, 2L)
  d3 <- .lta_tau_design(fit, 1L, 3L)
  expect_true(all(c("occ:2", "occ:3") %in% colnames(d1)))
  expect_true(all(d1[, "occ:2"] == 0) && all(d1[, "occ:3"] == 0))
  expect_true(all(d2[, "occ:2"] == 1) && all(d2[, "occ:3"] == 0))
  expect_true(all(d3[, "occ:2"] == 0) && all(d3[, "occ:3"] == 1))

  # The origin dummies and the covariate must be untouched by the append.
  expect_identical(colnames(d1)[seq_len(3L)], c("Intercept", "from:1", "z"))
})

test_that("the per-occasion transition matrices genuinely differ", {
  # The guard on .lta_par_unpack()'s broadcast. Remove
  # `&& is.null(state$tau_beta)` from that condition and this fails: every
  # occasion's matrix becomes occasion 1's, while the coefficients that
  # produced them stay different, so the fit carries two disagreeing
  # descriptions of itself.
  fit <- .slp_fit("slopes", 3L, 1000L)
  expect_false(isTRUE(all.equal(fit$tau[[1]], fit$tau[[2]])))

  # And the round trip through the packed vector must preserve that.
  layout <- .lta_par_layout(fit)
  rt <- .lta_par_unpack(.lta_par_pack(fit, layout), fit, layout)
  expect_equal(rt$tau_beta, fit$tau_beta, tolerance = 1e-10)
  expect_equal(rt$tau[[1]], fit$tau[[1]], tolerance = 1e-10)
  expect_equal(rt$tau[[2]], fit$tau[[2]], tolerance = 1e-10)
  expect_false(isTRUE(all.equal(rt$tau[[1]], rt$tau[[2]])))

  # A covariate-free "full" fit on the same data must still pool them -- the
  # guard must not have switched the broadcast off for the case it serves.
  f_full <- .slp_fit("full", 3L, 150L, covariates = FALSE)
  expect_true(isTRUE(all.equal(f_full$tau[[1]], f_full$tau[[2]])))
})

test_that("the score matrix agrees with finite differences under \"slopes\"", {
  sim <- .slp_sim(3L)
  X <- sim$X
  fit <- suppressMessages(suppressWarnings(
    fit_lta(X, n_statuses = 2, times = 3, measurement = "binary",
            transition_invariance = "slopes",
            predictors_initial = sim$Z, predictors_transition = sim$Z,
            smoothing = 0, bayes_constants = .ml, n_init = 2, max_iter = 150,
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
    sum(w * .lta_score_matrix(st, X)$ll) + .lta_penalty(st, X, layout, 0)$value
  }
  eps <- 1e-5
  fd_at <- function(p) vapply(seq_along(p), function(i) {
    e <- numeric(length(p)); e[i] <- eps
    (obj(p + e) - obj(p - e)) / (2 * eps)
  }, numeric(1))
  ana_at <- function(st) colSums(sweep(.lta_score_matrix(st, X)$S, 1, w, "*")) +
    .lta_penalty(st, X, layout, 0)$gradient

  # Checked at a perturbed point as well, so a wrong occasion column cannot
  # hide behind a gradient that is near zero in every direction anyway.
  set.seed(2)
  par1 <- par0 + stats::runif(length(par0), -0.05, 0.05)
  for (p in list(par0, par1)) {
    fd  <- fd_at(p)
    ana <- ana_at(.lta_par_unpack(p, fit, layout))
    expect_lt(max(abs(ana - fd)) / max(1, max(abs(fd))), 1e-5)
  }
})

test_that("\"slopes\" refuses the two models it cannot mean", {
  sim <- .slp_sim(3L)

  # Nothing predicts the transitions, so there are no slopes to share.
  expect_error(
    fit_lta(sim$X, n_statuses = 2, times = 3, measurement = "binary",
            transition_invariance = "slopes", n_init = 1, n_cores = 1,
            random_state = 5, standard_errors = FALSE),
    "no slopes to share")

  # Under "by_origin" the restriction would bind a different set of
  # coefficients, so it is refused rather than guessed at.
  expect_error(
    fit_lta(sim$X, n_statuses = 2, times = 3, measurement = "binary",
            transition_invariance = "slopes",
            predictors_transition = sim$Z, transition_effects = "by_origin",
            n_init = 1, n_cores = 1, random_state = 5,
            standard_errors = FALSE),
    "transition_effects")
})
