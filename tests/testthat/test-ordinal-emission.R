# ==============================================================================
# The "ordinal" measurement family (roadmap ### 14.18, plain half)
# ==============================================================================
#
# This is the emission-level and fit_lta()-plumbing half of ### 14.18: the
# random-intercept ordinal M-step is not built yet, so every check here uses
# `random_intercept = "none"`. The one identity test that needs no reference
# program at this stage is W6(b): at equal category counts, "ordinal" and
# "categorical" are the same model and must match bit-for-bit, prior
# included.

# ------------------------------------------------------------------------------
# one_hot_ragged()
# ------------------------------------------------------------------------------

test_that("one_hot_ragged rejects non-integer and out-of-range codes", {
  X <- matrix(c(1, 2, 1.5, 2, 1, 3), nrow = 3)
  expect_error(one_hot_ragged(X, c(2L, 3L)), "non-integer")

  X2 <- matrix(c(1, 2, 3, 1, 1, 2), nrow = 3)
  expect_error(one_hot_ragged(X2, c(2L, 3L)), "1\\.\\.2")
})

test_that("one_hot_ragged gives an all-zero row for NA", {
  X <- matrix(c(1, 2, NA, NA, 2, 3), nrow = 3, byrow = TRUE)
  oh <- one_hot_ragged(X, c(2L, 3L))
  expect_equal(oh[2, ], rep(0, 5))
  expect_equal(unname(rowSums(oh)[-2]), c(2, 2))
})

# ------------------------------------------------------------------------------
# n_parameters.ordinal and the ragged .per_item_nparams() branch
# ------------------------------------------------------------------------------

test_that("n_parameters.ordinal counts K*(cats-1) summed over items", {
  m <- ordinal_model(n_components = 5, cats = c(3L, 3L, 2L))
  expect_equal(n_parameters(m), 5L * (2L + 2L + 1L))
})

test_that(".per_item_nparams gives each ragged item its own cost", {
  m <- ordinal_model(n_components = 5, cats = c(3L, 3L, 2L))
  expect_equal(.per_item_nparams(m, 3L), 5L * c(2L, 2L, 1L))
})

# ------------------------------------------------------------------------------
# fit_lta() plumbing
# ------------------------------------------------------------------------------

test_that("an empty interior category is refused with a message naming the item", {
  set.seed(101)
  n <- 100
  # Item 1 never observes category 2 at either occasion.
  it1 <- sample(c(1L, 3L), n, replace = TRUE)
  it2 <- sample(1:3, n, replace = TRUE)
  X <- cbind(it1, it2, sample(c(1L, 3L), n, replace = TRUE), it2)
  expect_error(
    fit_lta(X, n_statuses = 2, times = 2, measurement = "ordinal",
            n_init = 1, random_state = 1, standard_errors = FALSE),
    "empty interior category")
})

test_that("random_intercept + ordinal is accepted (### 14.18, W4)", {
  set.seed(102)
  X <- matrix(sample(1:3, 200 * 4, replace = TRUE), 200, 4)
  expect_no_error(
    fit_lta(X, n_statuses = 2, times = 2, measurement = "ordinal",
            random_intercept = "continuous", n_quadrature = 3,
            n_init = 1, random_state = 1, standard_errors = FALSE))
})

test_that("random_intercept accepts the bernoulli alias (guard alias bug fixed)", {
  set.seed(103)
  X <- matrix(rbinom(200 * 4, 1, 0.5), 200, 4)
  expect_no_error(
    fit_lta(X, n_statuses = 2, times = 2, measurement = "bernoulli",
            random_intercept = "continuous", n_quadrature = 3,
            n_init = 1, random_state = 1, standard_errors = FALSE))
})

# ------------------------------------------------------------------------------
# W6(b) -- the lambda = 0 saturation identity: at equal category counts and no
# random intercept, "ordinal" and "categorical" are the same model. Column
# layout, priors and the RNG draws consumed on the way to init_params() are
# all identical, so this should hold to floating-point equality, not just to
# the 1e-8 the roadmap asks for.
# ------------------------------------------------------------------------------

test_that("ordinal and categorical are the same model at equal category counts", {
  # The saturation identity at the emission level, isolated from fit_lta()'s
  # random search (which can land two runs on different local optima even
  # when the underlying M-step is byte-for-byte the same): with identical
  # starting `pis`, one M-step and the log-likelihood must agree exactly,
  # prior included, because m_step.ordinal is m_step.multinoulli's formula
  # over the same column layout when cats == max_val everywhere.
  set.seed(21)
  n <- 300
  cats_cat <- categorical_model(2, type = "multinoulli", max_val = 3L)
  cats_ord <- ordinal_model(2, cats = c(3L, 3L))

  X <- matrix(sample(1:3, n * 2, replace = TRUE), n, 2)
  resp <- matrix(runif(n * 2), n, 2)
  resp <- resp / rowSums(resp)

  pis0 <- matrix(runif(2 * 6), 2, 6)
  for (j in 1:2) {
    cols <- ((j - 1) * 3 + 1):(j * 3)
    pis0[, cols] <- pis0[, cols] / rowSums(pis0[, cols])
  }
  cats_cat$parameters$pis <- pis0
  cats_ord$parameters$pis <- pis0

  ll_cat0 <- log_likelihood(cats_cat, X)
  ll_ord0 <- log_likelihood(cats_ord, X)
  expect_identical(ll_ord0, ll_cat0)

  fit_cat <- m_step(cats_cat, X, resp)
  fit_ord <- m_step(cats_ord, X, resp)
  expect_identical(fit_ord$parameters$pis, fit_cat$parameters$pis)
  expect_identical(n_parameters(fit_ord), n_parameters(fit_cat))

  expect_identical(log_likelihood(fit_ord, X), log_likelihood(fit_cat, X))
})

test_that("fit_lta() gives the same single-start fit under ordinal and categorical", {
  # End-to-end plumbing check, at n_init = 1 so there is exactly one candidate
  # and no restart-ranking to diverge on; unlike the previous test this
  # exercises .longitudinal_measurement_spec(), time_blocks_model() and the
  # whole EM driver, not just the emission in isolation.
  set.seed(23)
  n <- 200
  cls <- sample(1:2, n, replace = TRUE)
  draw <- function() apply(matrix(c(.7, .2, .1, .1, .2, .7), 2, 3,
                                  byrow = TRUE)[cls, ], 1,
                           function(p) sample(1:3, 1, prob = p))
  X <- cbind(draw(), draw(), draw(), draw())

  fit_cat <- fit_lta(X, n_statuses = 2, times = 2, measurement = "categorical",
                     n_init = 1, random_state = 9, standard_errors = FALSE)
  fit_ord <- fit_lta(X, n_statuses = 2, times = 2, measurement = "ordinal",
                     n_init = 1, random_state = 9, standard_errors = FALSE)

  # Not bit-identical: two different (if mathematically equivalent) code
  # paths accumulate floating-point rounding differently over the EM
  # iterations, which can shift the iteration the loose stopping rule fires
  # on by one. The exact identity is already proven per-iteration, at the
  # emission level, in the test above; this is the looser end-to-end check
  # that the plumbing (measurement spec, time_blocks_model, the EM driver)
  # wires the new family through correctly.
  expect_equal(fit_ord$loglik, fit_cat$loglik, tolerance = 1e-5)
  expect_equal(fit_ord$n_params, fit_cat$n_params)
})

# ------------------------------------------------------------------------------
# W5 -- parameter counts for the ragged 3/3/2 item block, K = 5, T = 3,
# stationary transitions: the article's own Table 7 shape (### 14.18.8), all
# six rows now: regular, continuous RI, binary RI, and their three
# mover-stayer crosses.
# ------------------------------------------------------------------------------

test_that("ordinal RI parameter counts match the W5 table (### 14.18.8)", {
  set.seed(201)
  n <- 150
  X <- cbind(sample(1:3, n, TRUE), sample(1:3, n, TRUE), sample(1:2, n, TRUE),
            sample(1:3, n, TRUE), sample(1:3, n, TRUE), sample(1:2, n, TRUE),
            sample(1:3, n, TRUE), sample(1:3, n, TRUE), sample(1:2, n, TRUE))

  fit_reg <- fit_lta(X, n_statuses = 5, times = 3, measurement = "ordinal",
                     transition_invariance = "full",
                     n_init = 1, random_state = 1, standard_errors = FALSE)
  expect_equal(fit_reg$n_params, 49L)

  fit_con <- fit_lta(X, n_statuses = 5, times = 3, measurement = "ordinal",
                     transition_invariance = "full",
                     random_intercept = "continuous", n_quadrature = 5,
                     n_init = 1, random_state = 1, standard_errors = FALSE)
  expect_equal(fit_con$n_params, 52L)

  fit_bin <- fit_lta(X, n_statuses = 5, times = 3, measurement = "ordinal",
                     transition_invariance = "full",
                     random_intercept = "binary", n_ri = 2,
                     n_init = 1, random_state = 1, standard_errors = FALSE)
  expect_equal(fit_bin$n_params, 53L)

  # The table's structural count (24, shared with the non-mover-stayer rows,
  # plus one for the stayer mixing weight) assumes ONE initial-status vector
  # shared across the mover and stayer classes, not two -- so these three need
  # tie_initial_status = TRUE, the same option test-lta-ri.R's mover-stayer
  # fits use for the same reason.
  fit_ms <- suppressWarnings(fit_lta(X, n_statuses = 5, times = 3,
                     measurement = "ordinal", transition_invariance = "full",
                     mover_stayer = TRUE, tie_initial_status = TRUE,
                     n_init = 1, random_state = 1, standard_errors = FALSE))
  expect_equal(fit_ms$n_params, 50L)

  fit_ms_con <- suppressWarnings(fit_lta(X, n_statuses = 5, times = 3,
                     measurement = "ordinal", transition_invariance = "full",
                     mover_stayer = TRUE, tie_initial_status = TRUE,
                     random_intercept = "continuous", n_quadrature = 5,
                     n_init = 1, random_state = 1, standard_errors = FALSE))
  expect_equal(fit_ms_con$n_params, 53L)

  fit_ms_bin <- suppressWarnings(fit_lta(X, n_statuses = 5, times = 3,
                     measurement = "ordinal", transition_invariance = "full",
                     mover_stayer = TRUE, tie_initial_status = TRUE,
                     random_intercept = "binary", n_ri = 2,
                     n_init = 1, random_state = 1, standard_errors = FALSE))
  expect_equal(fit_ms_bin$n_params, 54L)
})

# ------------------------------------------------------------------------------
# W6(a) -- binary-as-ordinal regression under a random intercept: a 2-category
# ordinal item is algebraically Bernoulli (### 14.18.2), so refitting the same
# data under "binary" and under "ordinal" (1/2-coded) must give the same
# log-likelihood and parameter count. Needs no reference program.
# ------------------------------------------------------------------------------

test_that("binary-as-ordinal RI regression: identical loglik and n_params (W6a)", {
  sim  <- .lta_ri_sim(n = 300, Tn = 3, J = 5, seed = 99)
  X01  <- sim$X
  X12  <- X01 + 1L   # ordinal codes must be 1-based

  fit_bin <- fit_lta(X01, n_statuses = 2, times = 3, measurement = "binary",
                     random_intercept = "continuous", n_quadrature = 5,
                     n_init = 3, random_state = 1, standard_errors = FALSE)
  fit_ord <- fit_lta(X12, n_statuses = 2, times = 3, measurement = "ordinal",
                     random_intercept = "continuous", n_quadrature = 5,
                     n_init = 3, random_state = 1, standard_errors = FALSE)

  expect_equal(fit_ord$loglik, fit_bin$loglik, tolerance = 1e-6)
  expect_equal(fit_ord$n_params, fit_bin$n_params)
})

# ------------------------------------------------------------------------------
# EM monotonicity for the new RI ordinal M-step (### 14.18, W4): both ECM
# cycles (thresholds, then the loading) must only ever increase the
# penalised objective. A short run on a small fixture is enough to catch a
# sign error in the Newton step, which is the failure mode this exists for;
# the roadmap's own 200-iteration check on the full Dating data is W7-sized
# and deferred with the rest of the expensive tail.
# ------------------------------------------------------------------------------

test_that("the RI ordinal M-step is EM-monotone", {
  # No per-iteration trace is exposed, so this truncates the SAME seeded run
  # at increasing iteration counts: since random_state pins the starting
  # point, the log-likelihood at max_iter = m is one point on the one EM
  # trajectory that run follows, and that sequence must never decrease.
  # `refine = FALSE`: ordinal joined `.lta_scores_full()`'s whitelist in
  # ### 14.18's W9, so the default L-BFGS polish now runs after every
  # truncated fit here too -- and each truncation's *own* polish run adds its
  # own small numerical noise, which breaks monotonicity across truncations
  # (a cross-run comparison the polish makes no claim about) even though the
  # M-step this test exists to check is, and stays, exactly monotone.
  set.seed(303)
  n <- 150
  X <- cbind(sample(1:3, n, TRUE), sample(1:2, n, TRUE),
            sample(1:3, n, TRUE), sample(1:2, n, TRUE))
  lls <- vapply(1:8, function(m_it) {
    # Non-convergence at a tiny max_iter is expected and not the point of
    # this test -- only the trajectory's monotonicity is.
    suppressWarnings(fit_lta(
      X, n_statuses = 2, times = 2, measurement = "ordinal",
      random_intercept = "continuous", n_quadrature = 5,
      n_init = 1, random_state = 5, standard_errors = FALSE,
      max_iter = m_it, refine = FALSE))$loglik
  }, numeric(1))
  expect_true(all(diff(lls) > -1e-6))
})

test_that("ordinal handles a ragged 3/3/2 category block without error", {
  set.seed(22)
  n <- 200
  X <- cbind(sample(1:3, n, replace = TRUE), sample(1:3, n, replace = TRUE),
            sample(1:2, n, replace = TRUE), sample(1:3, n, replace = TRUE),
            sample(1:3, n, replace = TRUE), sample(1:2, n, replace = TRUE))
  fit <- fit_lta(X, n_statuses = 2, times = 2, measurement = "ordinal",
                n_init = 2, random_state = 6, standard_errors = FALSE)
  # Measurement block alone (the ragged part being exercised here): the
  # structural (delta/transition) count is not asserted, since that is W5's
  # job and belongs to a later session.
  expect_equal(n_parameters(fit$mm), 2L * (2L + 2L + 1L))
  expect_true(is.finite(fit$loglik))
})

# ------------------------------------------------------------------------------
# .lta_ri_log_prior()'s ordinal arm (R/lta_ri.R), added so bayes_constants()'s
# default categorical = 1 does not crash the very first ordinal-RI EM
# iteration, had no dedicated test of its own -- the smoke tests above only
# exercise it incidentally, by fitting to convergence at defaults. This
# recomputes its value independently: base-R plogis instead of
# .ordinal_cat_probs(), and a hand-rolled weighted marginal instead of
# .lta_ri_item_marginal_ordinal(). alpha = 0 isolates the categorical term
# from .lta_log_prior_dt()'s delta/tau contribution.
# ------------------------------------------------------------------------------

test_that(".lta_ri_log_prior()'s ordinal arm matches an independent computation", {
  set.seed(55)
  n <- 80
  cls  <- sample(1:2, n, replace = TRUE)
  draw <- function() apply(matrix(c(.75, .25, .2, .8), 2, 2, byrow = TRUE)[cls, ],
                           1, function(p) sample(1:2, 1, prob = p))
  X <- cbind(draw(), draw(), draw())   # one item, T = 3, S = 2 (binary-collapsed)

  fit <- suppressWarnings(fit_lta(X, n_statuses = 2, times = 3,
    measurement = "ordinal", random_intercept = "continuous", n_quadrature = 3,
    n_init = 1, max_iter = 3, random_state = 1, standard_errors = FALSE))

  val <- .lta_ri_log_prior(fit, X, alpha = 0)

  K  <- fit$n_statuses
  Q  <- length(fit$ri$mass)
  Tn <- fit$n_times
  a_cat     <- .bayes_alpha(fit$mm$models[[1]], "categorical")
  prior_obs <- Tn * a_cat / (K * Q)
  w  <- fit$weights_vec
  xs <- as.vector(X)
  ws <- rep(w, times = Tn)
  m1 <- sum(ws[xs == 1]) / sum(ws)
  m2 <- sum(ws[xs == 2]) / sum(ws)
  theta <- fit$ri$theta[, 1]                    # one item, S = 2 -> one column
  expected <- 0
  for (q in seq_len(Q)) {
    shift <- sum(fit$ri$L[1, ] * fit$ri$Dnode[q, ])
    p2 <- plogis(theta + shift)
    p1 <- 1 - p2
    expected <- expected + prior_obs *
      sum(m1 * log(pmax(p1, 1e-300)) + m2 * log(pmax(p2, 1e-300)))
  }
  expect_equal(val, expected, tolerance = 1e-8)

  # And the fact this test exists for: a fresh ordinal-RI fit at the
  # bayes_constants() DEFAULT (categorical = 1, not 0) must not crash or
  # produce a non-finite penalised objective on its first iteration.
  fit_default <- suppressWarnings(fit_lta(X, n_statuses = 2, times = 3,
    measurement = "ordinal", random_intercept = "continuous", n_quadrature = 3,
    n_init = 1, max_iter = 1, random_state = 1, standard_errors = FALSE))
  expect_true(is.finite(.lta_ri_log_prior(fit_default, X, alpha = 1)))
})
