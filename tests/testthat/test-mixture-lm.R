# ==============================================================================
# Mixture latent Markov models and the mover-stayer restriction
# ==============================================================================
#
# A mixture latent Markov model puts latent classes above the chain: each class
# has its own initial distribution and its own transition matrices, and they
# share one measurement model. The mover-stayer model is the restriction in
# which the last class's transition matrix is the identity - a group with zero
# probability of change (Vermunt, "Mover-Stayer Models"; Blumen, Kogan &
# McCarthy, 1955).

# ------------------------------------------------------------------------------
# Structural checks, which need no reference program
# ------------------------------------------------------------------------------

test_that("one class reproduces the single-chain model exactly", {
  set.seed(21)
  X <- matrix(rbinom(300 * 6, 1, 0.5), ncol = 6)
  args <- list(X, n_statuses = 2, times = 3, measurement = "binary",
               n_init = 5, random_state = 2, standard_errors = FALSE)

  plain <- do.call(fit_lta, args)
  one   <- do.call(fit_lta, c(args, n_classes = 1))

  expect_equal(plain$loglik, one$loglik)
  expect_equal(plain$delta, one$delta)
  expect_equal(plain$tau, one$tau)
  expect_equal(plain$n_params, one$n_params)
  # And the single-chain shapes are untouched: a vector and a list of matrices,
  # not a one-row matrix and a list of lists.
  expect_true(is.null(dim(plain$delta)))
  expect_true(is.matrix(plain$tau[[1]]))
  expect_null(plain$class_weights)
})

test_that("the stayer restriction costs the parameters it should", {
  set.seed(22)
  X <- matrix(rbinom(400 * 12, 1, 0.5), ncol = 12)
  args <- list(X, n_statuses = 2, times = 3, measurement = "binary",
               n_init = 3, random_state = 6, standard_errors = FALSE,
               max_iter = 300, tol = 1e-8)

  free <- do.call(fit_lta, c(args, n_classes = 2))
  ms   <- do.call(fit_lta, c(args, mover_stayer = TRUE))

  # Free: 1 class weight + 2 x 1 initial + 2 classes x 2 matrices x 2 rows x 1
  # + 4 items x 2 statuses. Mover-stayer drops the stayer's four rows.
  expect_equal(free$n_params, ms$n_params + 4L)
  expect_equal(ms$n_classes, 2L)
  for (m in ms$tau[[2]]) expect_identical(m, diag(2))
})

test_that("a mixture needs at least three occasions", {
  X <- matrix(rbinom(100 * 4, 1, 0.5), ncol = 4)
  expect_error(
    fit_lta(X, n_statuses = 2, times = 2, measurement = "binary",
            mover_stayer = TRUE),
    "at least three occasions")
})

test_that("covariates are refused rather than silently ignored", {
  X <- matrix(rbinom(200 * 6, 1, 0.5), ncol = 6)
  expect_error(
    fit_lta(X, n_statuses = 2, times = 3, measurement = "binary",
            n_classes = 2, predictors_initial = rnorm(200)),
    "not yet available")
})

test_that("standard errors are declined, not reported wrongly", {
  set.seed(23)
  X <- matrix(rbinom(300 * 9, 1, 0.5), ncol = 9)
  fit <- fit_lta(X, n_statuses = 2, times = 3, measurement = "binary",
                 mover_stayer = TRUE, n_init = 3, random_state = 8,
                 max_iter = 300, tol = 1e-8)
  expect_null(fit$se)
  expect_output(print(fit), "standard errors are not available")
})

test_that("the extractors and plots take a class", {
  skip_on_cran()
  # Simulated mover-stayer data: half the cases never change status, the rest
  # move with probability 1/4 per interval; a single strong binary indicator.
  set.seed(24)
  n      <- 500
  stayer <- rbinom(n, 1, 0.5)
  s      <- rbinom(n, 1, 0.6)
  X      <- matrix(0L, n, 5)
  for (t in 1:5) {
    if (t > 1) {
      move <- rbinom(n, 1, 0.25) == 1 & stayer == 0
      s    <- ifelse(move, 1L - s, s)
    }
    X[, t] <- rbinom(n, 1, ifelse(s == 1, 0.9, 0.1))
  }
  fit <- fit_lta(X, n_statuses = 2, times = 5, measurement = "binary",
                 mover_stayer = TRUE, n_init = 3, random_state = 11,
                 standard_errors = FALSE, max_iter = 400, tol = 1e-8)

  tm <- transition_matrix(fit)
  expect_length(tm, 2L)
  expect_length(transition_matrix(fit, class = 1), 4L)
  expect_true(is.matrix(transition_matrix(fit, occasion = 1, class = 1)))

  # The whole-sample prevalences are the class-weighted average of the
  # per-class ones, which is what makes them comparable with the observed
  # proportions.
  overall <- status_prevalences(fit)
  expect_null(attr(overall, "by_class"))    # not printed as a stray attribute
  parts <- lapply(1:2, function(c) status_prevalences(fit, class = c))
  expect_equal(overall,
               parts[[1]] * fit$class_weights[1] +
                 parts[[2]] * fit$class_weights[2],
               tolerance = 1e-10)

  # A stayer's status prevalences do not move, which is the model's whole
  # claim about that class.
  stayer <- status_prevalences(fit, class = 2)
  expect_lt(max(abs(sweep(stayer, 2, stayer[1, ]))), 1e-10)

  pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_silent(plot(fit, "prevalence"))
  expect_silent(plot(fit, "transitions", class = 1))
})

# ------------------------------------------------------------------------------
# The degenerate optimum, recorded on purpose
# ------------------------------------------------------------------------------

test_that("the collapse detector finds classes describing the same chain", {
  # Tested on a constructed state rather than by hoping a fit lands exactly on
  # a degenerate optimum: EM on noise separates the classes a little, so an
  # integration test of this would be asserting a coin flip. What has to be
  # right is the comparison itself.
  same <- matrix(c(.8, .2, .3, .7), 2, 2, byrow = TRUE)
  diff <- matrix(c(.5, .5, .5, .5), 2, 2, byrow = TRUE)

  state <- list(n_classes = 2L,
                delta_c = list(c(.6, .4), c(.6, .4)),
                tau_c   = list(list(same), list(same)))
  expect_equal(.lta_collapsed_classes(state), matrix(c(1L, 2L), 1, 2))

  # A difference in the transitions alone is enough to tell them apart...
  state$tau_c[[2]] <- list(diff)
  expect_null(.lta_collapsed_classes(state))

  # ...and so is a difference in the initial distribution alone.
  state$tau_c[[2]] <- list(same)
  state$delta_c[[2]] <- c(.4, .6)
  expect_null(.lta_collapsed_classes(state))

  # One class is never "collapsed" with itself.
  expect_null(.lta_collapsed_classes(list(n_classes = 1L)))
})

test_that("a fit records whether its classes collapsed", {
  set.seed(24)
  X <- matrix(rbinom(300 * 9, 1, 0.5), ncol = 9)
  fit <- fit_lta(X, n_statuses = 2, times = 3, measurement = "binary",
                 n_classes = 2, n_init = 2, random_state = 31,
                 standard_errors = FALSE, max_iter = 200, tol = 1e-8)
  # Whether EM collapses these classes is not guaranteed - it is noise - but
  # the verdict stored on the object must be the detector's own, so that a fit
  # read back months later still carries it.
  expect_identical(fit$collapsed_classes, .lta_collapsed_classes(fit))
  expect_equal(fit$n_classes, 2L)
})
