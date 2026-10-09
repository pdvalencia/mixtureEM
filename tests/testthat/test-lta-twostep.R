# ==============================================================================
# fit_lta(n_steps = 2) -- the two-step estimator for latent transition models.
#
# Step 1 fits the measurement model alone. Step 2 holds it there and maximises
# the initial-status and transition regressions on the full likelihood, the
# E-step still running on the joint model. The estimator is therefore the
# one-step estimator with one block held fixed, and every check here follows
# from that sentence alone:
#
#   (a) the freeze bit. The item probabilities after step 2 are step 1's, to
#       machine precision. This is the check that the M-step actually skipped
#       the measurement update rather than converging back to something near
#       it; nothing else in the file means anything if this one fails.
#   (b) the freeze does something. The one-step fit's measurement block is
#       NOT step 1's -- the covariate reshapes the statuses when it is allowed
#       to, which is the whole reason the estimator exists.
#   (c) the freeze costs likelihood, and unfreezing recovers exactly that
#       cost: the two-step sits below the one-step, and continuing the
#       two-step solution with the measurement block free lands on the
#       one-step maximum.
#   (d) the argument's guard rails.
#
# TWO fixtures, for a reason worth not rediscovering. `.lta_cov_refine_sim()`
# (helper-lta-packing.R) has deliberately weak measurement separation -- item
# probabilities drawn uniform on (0.25, 0.75) -- which is what makes it a good
# test of the covariate score blocks and a bad one for anything that compares
# two optima. Measured 2026-09-22: on that fixture a covariate fit has not
# converged after 1000 iterations, one-step and continuation alike, and the
# two are 1.54 apart purely because they are still climbing at different
# rates. Checks (a) and (b) are exact or one-sided and run there; check (c)
# compares two optima and needs a fixture where the optimum is reachable, so
# it runs on the separated one below, where both fits converge.
#
# `order_by_size = FALSE` throughout. Comparing step 1's item table with step
# 2's requires the two to be in the same status order, and sorting by size is
# free to put them in different ones once the initial-status distribution has
# been replaced by a regression.
# ==============================================================================

.ts_cache <- new.env(parent = emptyenv())

.ts_sim <- function() {
  if (is.null(.ts_cache$sim)) .ts_cache$sim <- .lta_cov_refine_sim(Tn = 3)
  .ts_cache$sim
}

# The separated fixture: the same generating mechanism, with the measurement
# model pinned well away from 0.5 so EM converges in tens of iterations
# instead of thousands.
.ts_sep <- function() {
  if (!is.null(.ts_cache$sep)) return(.ts_cache$sep)
  set.seed(404)
  n <- 600L; Tn <- 3L; J <- 4L
  z  <- stats::rnorm(n)
  s  <- stats::rbinom(n, 1, stats::plogis(0.3 + 0.8 * z)) + 1L
  S  <- matrix(0L, n, Tn); S[, 1] <- s
  for (t in 2:Tn) {
    p_stay <- stats::plogis(0.4 + ifelse(S[, t - 1] == 1, 0.6, -0.6) * z)
    S[, t] <- ifelse(stats::runif(n) < p_stay, S[, t - 1], 3L - S[, t - 1])
  }
  rho <- matrix(c(rep(0.12, J), rep(0.88, J)), 2L, J, byrow = TRUE)
  X <- do.call(cbind, lapply(seq_len(Tn), function(t)
    vapply(seq_len(J), function(j) stats::rbinom(n, 1, rho[S[, t], j]),
           numeric(n))))
  colnames(X) <- paste0("t", rep(seq_len(Tn), each = J), "_i",
                        rep(seq_len(J), Tn))
  .ts_cache$sep <- list(X = X, Z = data.frame(z = z))
  .ts_cache$sep
}

.ts_fit <- function(key, sim, ...) {
  if (!is.null(.ts_cache[[key]])) return(.ts_cache[[key]])
  args <- list(sim$X, n_statuses = 2, times = 3, measurement = "binary",
               n_cores = 1, random_state = 5, refine = FALSE,
               standard_errors = FALSE, order_by_size = FALSE)
  # modifyList rather than c(): the standard-error block below overrides one of
  # the defaults above, and a duplicated argument name is an error in do.call.
  .ts_cache[[key]] <- suppressMessages(suppressWarnings(
    do.call(fit_lta, utils::modifyList(args, list(...)))))
  .ts_cache[[key]]
}

# Weak fixture. The one-step fit is capped: the checks that use it need its
# measurement block, not its optimum.
.ts_two <- function()
  .ts_fit("two", .ts_sim(), predictors_initial = .ts_sim()$Z,
          predictors_transition = .ts_sim()$Z, n_steps = 2, n_init = 6)
.ts_one <- function()
  .ts_fit("one", .ts_sim(), predictors_initial = .ts_sim()$Z,
          predictors_transition = .ts_sim()$Z, n_init = 6, max_iter = 300)
.ts_plain <- function() .ts_fit("plain", .ts_sim(), n_init = 6)

# Separated fixture.
.ts_two_s <- function()
  .ts_fit("two_s", .ts_sep(), predictors_initial = .ts_sep()$Z,
          predictors_transition = .ts_sep()$Z, n_steps = 2, n_init = 4)
.ts_one_s <- function()
  .ts_fit("one_s", .ts_sep(), predictors_initial = .ts_sep()$Z,
          predictors_transition = .ts_sep()$Z, n_init = 4)

# --- (a) the freeze bit -------------------------------------------------------

test_that("n_steps = 2 leaves the measurement model exactly where step 1 put it", {
  two <- .ts_two()

  expect_equal(two$n_steps, 2L)
  expect_s3_class(two$step1, "lta_model")
  # Step 1 is the same model with the predictors dropped, so it has the plain
  # initial-status and transition tables and no regression blocks at all.
  expect_null(two$step1$delta_beta)
  expect_null(two$step1$tau_beta)
  expect_false(is.null(two$delta_beta))
  expect_false(is.null(two$tau_beta))

  # The comparison that matters: every item probability, at every occasion,
  # identical to the step-1 fit's. Not "close" -- identical, because the
  # M-step never touched them.
  p2 <- lapply(two$mm$models, function(m) m$parameters$pis)
  p1 <- lapply(two$step1$mm$models, function(m) m$parameters$pis)
  expect_equal(p2, p1, tolerance = 0)

  # And the fit is genuinely a different fit from its own step 1: the
  # structural side moved, so step 2 did work.
  expect_gt(two$loglik, two$step1$loglik)
  expect_true(two$converged)
})

# --- (b) the freeze does something --------------------------------------------

test_that("the frozen fit is not simply the one-step fit", {
  two <- .ts_two()
  one <- .ts_one()

  # If the freeze did nothing, the two measurement tables would agree.
  d <- max(abs(one$mm$models[[1]]$parameters$pis -
                 two$mm$models[[1]]$parameters$pis))
  expect_gt(d, 1e-6)

  # The coefficients are real numbers of the right shape. `tau_beta` is one
  # matrix per interval, so it is unlisted before it is checked.
  expect_true(all(is.finite(two$delta_beta)))
  expect_true(all(is.finite(unlist(two$tau_beta))))
  expect_equal(dim(two$delta_beta), dim(one$delta_beta))
  expect_equal(lengths(two$tau_beta), lengths(one$tau_beta))
})

# --- (c) the cost of the freeze, and that unfreezing recovers it --------------

test_that("the two-step sits below the one-step, and unfreezing closes the gap", {
  two <- .ts_two_s()
  one <- .ts_one_s()
  # The claim is about two optima, so both have to be at one.
  expect_true(two$converged)
  expect_true(one$converged)

  # Below: the one-step maximises the same objective over strictly more
  # parameters, so it cannot be beaten by the restricted fit.
  expect_lt(two$loglik, one$loglik + 1e-6)
  # Above: step 2 maximises over the structural block, so it must improve on
  # the measurement-only fit it started from.
  expect_gt(two$loglik, two$step1$loglik)

  # The same covariate model, started at the two-step solution, with nothing
  # frozen. If the freeze is the only difference between the estimators, this
  # has to close the gap the first assertion just measured.
  free <- suppressMessages(suppressWarnings(
    fit_lta(.ts_sep()$X, n_statuses = 2, times = 3, measurement = "binary",
            predictors_initial = .ts_sep()$Z,
            predictors_transition = .ts_sep()$Z,
            n_cores = 1, random_state = 5, refine = FALSE,
            standard_errors = FALSE, order_by_size = FALSE,
            refine_from = two)))

  expect_true(free$converged)
  expect_gt(free$loglik, two$loglik)
  # Stated as a fraction of the gap being closed, which is what the claim is
  # about and is free of the fixture's scale. Measured 2026-09-22: the freeze
  # costs 0.0278 of log-likelihood on this fixture -- small, which is itself
  # the theory's prediction when the statuses are well separated -- and the
  # continuation recovers all but 3.2e-4 of it, 1.2 percent. That residual is
  # the stopping rule and not a second optimum: both fits converged, and two
  # runs approaching one maximum from different directions each stop when
  # their own objective stops moving. The bound is 5 percent, so it has real
  # margin without being a statement about the fixture's scale.
  gap <- one$loglik - two$loglik
  expect_gt(gap, 1e-3)
  expect_lt(abs(free$loglik - one$loglik), 0.05 * gap)
  # And the measurement block is free again: it moved off step 1's values.
  d <- max(abs(free$mm$models[[1]]$parameters$pis -
                 two$step1$mm$models[[1]]$parameters$pis))
  expect_gt(d, 1e-6)
})

# --- (d) guard rails ----------------------------------------------------------

test_that("n_steps validates its arguments", {
  sim  <- .ts_sim()
  base <- list(sim$X, n_statuses = 2, times = 3, measurement = "binary",
               n_init = 1, n_cores = 1, random_state = 5, refine = FALSE,
               standard_errors = FALSE)

  # Nothing to put in step two.
  expect_error(do.call(fit_lta, c(base, list(n_steps = 2))),
               "needs a structural model")
  expect_error(do.call(fit_lta, c(base, list(n_steps = 0))),
               "must be 1, 2 or 3")
  # The three-step's own arguments mean nothing to the two-step.
  expect_error(do.call(fit_lta, c(base, list(n_steps = 2, assignment = "modal"))),
               "apply only with `n_steps = 3`")
  # Two ways of handing over a starting point, one of which is the estimator.
  expect_error(
    suppressWarnings(fit_lta(sim$X, n_statuses = 2, times = 3,
                             measurement = "binary",
                             predictors_initial = sim$Z, n_steps = 2,
                             n_cores = 1, refine = FALSE,
                             standard_errors = FALSE,
                             refine_from = .ts_plain())),
    "has nothing to hand over")
})

test_that("an ordinary fit is untouched by the freeze machinery", {
  # The guard in .lta_em() and the one in .lta_refine_lbfgs() both key on a
  # field that is NULL on every path but the two-step's; this is the check
  # that a plain fit reports no step-1 object and carries no freeze.
  plain <- .ts_plain()
  expect_equal(plain$n_steps, 1L)
  expect_null(plain$step1)
  expect_null(plain$frozen)
})

# --- (e) the standard errors --------------------------------------------------
#
# Step 2 holds the measurement block at an ESTIMATE, not at a known constant,
# and the uncertainty in that estimate belongs in the standard errors of the
# coefficients built on top of it. The pseudo-maximum-likelihood variance
# (Bakk & Kuha, 2018, eq. 5) is V = V2 + V1: V2 the inverse observed
# information of the full likelihood in the structural coefficients, V1 the
# step-one sampling variance carried across by the cross-curvature between the
# two blocks. Two checks, and neither needs anything outside the package:
#
#   V2 against an independent numerical Hessian of a from-scratch joint
#   likelihood, which shares no code with the differencing that produced it;
#   and V1 against what it is for -- it is most of the variance when the
#   statuses are poorly separated and nearly none of it when they are not.

.ts_se <- function(key, sim, ...)
  .ts_fit(key, sim, ..., standard_errors = TRUE)

.ts_two_se <- function()
  .ts_se("two_se", .ts_sim(), predictors_initial = .ts_sim()$Z,
         predictors_transition = .ts_sim()$Z, n_steps = 2, n_init = 6)
.ts_two_s_se <- function()
  .ts_se("two_s_se", .ts_sep(), predictors_initial = .ts_sep()$Z,
         predictors_transition = .ts_sep()$Z, n_steps = 2, n_init = 4)

# sqrt of the diagonals of V2 and of V, structural block only, in packing
# order, plus the vector of standard errors an independent Hessian gives.
.ts_se_pieces <- function(fit, X) {
  layout <- mixtureEM:::.lta_par_layout(fit)
  par    <- mixtureEM:::.lta_par_pack(fit, layout)
  idx    <- mixtureEM:::.lta_par_split(layout)$structural
  w      <- fit$weights_vec
  ll     <- function(v) {
    p <- par; p[idx] <- v
    sum(w * mixtureEM:::.lta_ll_case(fit, X, p, layout))
  }
  list(sd2   = sqrt(diag(fit$se$twostep_V2)),
       sd    = sqrt(diag(fit$se$vcov[idx, idx, drop = FALSE])),
       sd_oh = sqrt(diag(solve(-stats::optimHess(par[idx], ll)))))
}

test_that("the two-step variance is V2 plus the step-one term", {
  fit <- .ts_two_se()
  expect_true(isTRUE(fit$se$twostep))
  expect_match(fit$se$method, "Bakk and Kuha")

  p <- .ts_se_pieces(fit, .ts_sim()$X)

  # V2 is the step-two-only variance, and optimHess() reaches it by a
  # different route: a Hessian of the whole structural vector at once against
  # the block-by-block central differences .lta_twostep_information() takes.
  # On this weakly separated fixture the likelihood is flat enough that two
  # finite-difference Hessians agree only to about 1e-3, depending on exactly
  # where EM stopped. Worst coordinate, measured 2026-10-08, was 1.3e-4 at the
  # default tol and 6e-5 to 3e-4 from tol = 1e-8 to 1e-14. On macOS arm64 it
  # was 1.75e-3. A wrong information formula is wrong at order one, so 1e-2
  # still separates the two; the separated fixture below holds 1e-3.
  expect_lt(max(abs(p$sd2 - p$sd_oh) / p$sd_oh), 1e-2)

  # V1 is positive semi-definite, so every standard error grows.
  expect_true(all(p$sd >= p$sd2 - 1e-10))

  # And it is not a rounding correction. This fixture's statuses are weakly
  # separated on purpose (item probabilities uniform on 0.25-0.75), which is
  # the regime Bakk & Kuha's Table 2 is about: measured 2026-09-22, the
  # step-two-only standard error of the initial-status intercept is 0.18 of
  # the whole. A ratio of 1 everywhere would mean the cross-curvature block is
  # wrong, not that it is small.
  expect_lt(min(p$sd2 / p$sd), 0.5)
})

test_that("the step-one term vanishes as the statuses separate", {
  fit <- .ts_two_s_se()
  p   <- .ts_se_pieces(fit, .ts_sep()$X)

  expect_lt(max(abs(p$sd2 - p$sd_oh) / p$sd_oh), 1e-3)
  # Same estimator, same formula, statuses well separated: with almost no
  # classification uncertainty left there is almost nothing for step one to
  # contribute. Measured 2026-09-22: the worst coordinate is 0.986.
  expect_gt(min(p$sd2 / p$sd), 0.95)
  expect_true(all(p$sd >= p$sd2 - 1e-10))
})

test_that("an ordinary fit reports no two-step variance", {
  plain <- .ts_fit("plain_se", .ts_sep(), n_init = 2, standard_errors = TRUE)
  expect_false(isTRUE(plain$se$twostep))
  expect_null(plain$se$twostep_V2)
  # The measurement block is still described, so the fit did not lose its
  # ordinary standard errors on the way past the new branch.
  expect_true(length(plain$se$prob_se) > 0L)
})
