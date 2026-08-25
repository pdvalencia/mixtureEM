# ==============================================================================
# The assignment rule behind the three-step corrections, and the class-membership
# prior at step three
# ==============================================================================
#
# Two things are guarded here, and both are about *not* moving:
#
#   1. `assignment = "proportional"` is the default and must reproduce the
#      results the package produced before the argument existed. The anchors
#      below were recorded from that earlier version on this simulated data.
#   2. `bayes_constants = list(latent = 0)` must turn the step-three prior off
#      exactly, coefficients and standard errors alike. That is what makes the
#      changed default opt-out-able.

.sim_assignment_data <- function(n = 250, seed = 42) {
  set.seed(seed)
  cl    <- rbinom(n, 1, 0.4)
  items <- sapply(1:5, function(j) rbinom(n, 1, ifelse(cl == 1, 0.8, 0.2)))
  colnames(items) <- paste0("item", 1:5)
  list(items = items,
       covs  = data.frame(age = rnorm(n) + 0.8 * cl, z = rnorm(n) - 0.5 * cl),
       bmi   = rnorm(n, 25 + 2 * cl))
}

# Free standard errors of a covariate model, in the layout .fit_mnl() stores.
.free_cov_se <- function(fit) {
  V <- fit$sm$parameters$V_robust
  if (is.null(V)) V <- solve(-fit$sm$parameters$hessian)
  s <- sqrt(abs(diag(V)))
  s[s > 1e-8]
}

# ------------------------------------------------------------------------------
# Part 1 - the assignment rule
# ------------------------------------------------------------------------------

test_that("proportional assignment reproduces the pre-argument results", {
  d    <- .sim_assignment_data()
  fit0 <- suppressMessages(suppressWarnings(
    fit_mixture(d$items, n_classes = 2,
                measurement = "binary",
                n_init = 5, random_state = 9)))

  # Recorded before the `assignment` argument was threaded through fit_bch().
  # Proportional assignment leaves A equal to the posteriors, so the correction
  # is algebraically untouched and these must hold to every digit that matters.
  fb <- suppressMessages(suppressWarnings(add_outcome(fit0, d$bmi)))
  expect_equal(as.vector(fb$sm$parameters$means),
               c(24.8576510753, 27.0387838023), tolerance = 1e-10)

  # Naming the default explicitly must change nothing.
  fb2 <- suppressMessages(suppressWarnings(
    add_outcome(fit0, d$bmi, assignment = "proportional")))
  expect_equal(fb$sm$parameters$means, fb2$sm$parameters$means,
               tolerance = 1e-10)
  expect_identical(fb2$assignment, "proportional")
})

test_that("with a classification table at the identity the two rules agree", {
  # Three exact response patterns, so every case is assigned with probability
  # one. The classification table is then the identity, both corrections invert
  # nothing, and the two rules must agree with each other and with no
  # correction at all. The patterns are deterministic on purpose: with any
  # noise at all a handful of cases out of 300 land between two classes, and
  # the table stops being the identity for reasons that say nothing about the
  # assignment rule.
  set.seed(4)
  n  <- 300
  cl <- sample(1:3, n, replace = TRUE)
  items <- rbind(rep(0, 8), rep(1, 8), rep(c(1, 0), 4))[cl, ]
  colnames(items) <- paste0("item", 1:8)
  y <- rnorm(n, c(0, 5, 10)[cl], 0.5)

  fit0 <- suppressMessages(suppressWarnings(
    fit_mixture(items, n_classes = 3,
                measurement = "binary",
                n_init = 10, random_state = 3)))
  resp <- exp(fit0$log_resp)
  skip_if_not(max(abs(resp - get_modal_resp(resp))) < 1e-6,
              "the simulated classes did not separate cleanly enough")

  m <- function(a) as.vector(suppressMessages(suppressWarnings(
    add_outcome(fit0, y, correction = "BCH", assignment = a)
  )$sm$parameters$means))
  none <- as.vector(suppressMessages(suppressWarnings(
    add_outcome(fit0, y, correction = "none"))$sm$parameters$means))

  expect_equal(m("proportional"), m("modal"), tolerance = 1e-6)
  expect_equal(m("modal"), none, tolerance = 1e-6)
})

test_that("modal assignment still corrects when the classes overlap", {
  # The test above cannot see whether the modal correction does anything: it
  # chooses data whose classification table *is* the identity, so both rules
  # agree however the table was built. This one is the complement, and it is
  # the guard that matters. A modal assignment matrix has a single 1 per row,
  # so crossing it with itself gives a diagonal - and a classification table
  # that comes out as the identity by construction rather than because the
  # classes are well separated. The correction would then invert nothing and
  # hand back the naive modal means, silently, on every dataset.
  d    <- .sim_assignment_data()
  fit0 <- suppressMessages(suppressWarnings(
    fit_mixture(d$items, n_classes = 2,
                measurement = "binary",
                n_init = 5, random_state = 9)))

  resp  <- exp(fit0$log_resp)
  modal <- get_modal_resp(resp)

  # The premise: these classes genuinely overlap, so there is something for the
  # correction to correct. Without this the assertions below would pass
  # vacuously on data that happened to separate.
  expect_gt(max(abs(resp - modal)), 0.05)

  bch_modal <- as.vector(suppressMessages(suppressWarnings(
    add_outcome(fit0, d$bmi, correction = "BCH", assignment = "modal")
  )$sm$parameters$means))

  # What the naive modal means are: assign each case to its most likely class
  # and average the outcome within class, with no correction at all. This is
  # exactly what a no-op correction returns.
  naive <- colSums(modal * d$bmi) / colSums(modal)

  expect_gt(max(abs(bch_modal - naive)), 0.01)

  # And the two rules must part company here, where the table is not the
  # identity - the whole reason the argument exists.
  bch_prop <- as.vector(suppressMessages(suppressWarnings(
    add_outcome(fit0, d$bmi, correction = "BCH", assignment = "proportional")
  )$sm$parameters$means))
  expect_gt(max(abs(bch_modal - bch_prop)), 0.01)
})

# ------------------------------------------------------------------------------
# Part 2 - the step-three prior on the class probabilities
# ------------------------------------------------------------------------------

test_that("bayes_constants$latent = 0 reproduces the unpenalised step three", {
  d <- .sim_assignment_data()
  # Both the measurement model and the covariate model run unpenalised here, so
  # these are the numbers the package produced before the prior reached step
  # three. Recorded from that version.
  fit0 <- suppressMessages(suppressWarnings(
    fit_mixture(d$items, n_classes = 2,
                measurement = "binary",
                n_init = 5, random_state = 9,
                bayes_constants = list(latent = 0))))
  fc <- suppressMessages(suppressWarnings(
    add_covariates(fit0, d$covs, correction = "ML")))

  # Row 1 is the free class; row 2 is the anchor and is zero by construction.
  expect_equal(as.vector(fc$sm$parameters$beta[1, ]),
               c(0.745542992496, -0.673660063212, 0.544359886248),
               tolerance = 1e-6)
  expect_equal(.free_cov_se(fc),
               c(0.197280928482, 0.175671366223, 0.170233042872),
               tolerance = 1e-4)
})

test_that("the prior shrinks both the coefficients and their standard errors", {
  # One step-one solution, two step-three fits, so nothing but the prior differs
  # between them. Shrinkage in both is the documented direction of the effect,
  # and a sign error is what would break it.
  d    <- .sim_assignment_data()
  fit0 <- suppressMessages(suppressWarnings(
    fit_mixture(d$items, n_classes = 2,
                measurement = "binary",
                n_init = 5, random_state = 9)))
  fit_off <- fit0
  fit_off$bayes_constants <- list(latent = 0)

  on  <- suppressMessages(suppressWarnings(
    add_covariates(fit0, d$covs, correction = "ML")))
  off <- suppressMessages(suppressWarnings(
    add_covariates(fit_off, d$covs, correction = "ML")))

  b_on  <- as.vector(on$sm$parameters$beta)
  b_off <- as.vector(off$sm$parameters$beta)
  free  <- abs(b_off) > 1e-8
  expect_true(all(abs(b_on[free]) < abs(b_off[free])))
  expect_true(all(.free_cov_se(on) < .free_cov_se(off)))
})
