# ==============================================================================
# The EM stopping rule for emissions L-BFGS does not polish
# ==============================================================================
#
# The package's default stopping rule is loose on purpose: EM only has to reach
# the neighbourhood of the optimum, because refine_lbfgs() then climbs to the
# penalised maximum. That reasoning holds only for the emissions the refinement
# actually covers. For the rest — counts, polytomous items, mixed measurement
# models, block models over those — where EM stops *is* the answer, and the loose
# rule stopped them after five to ten iterations, tens of log-likelihood units
# short, with rates off by 2.8 and response probabilities off by 0.56.
#
# fit_single_init() now derives the rule from the same whitelist the refinement
# reads, so the two cannot drift apart again. These tests pin that down: the
# whitelist decides, the unpolished models really do converge, and the polished
# ones are untouched.

test_that(".is_refinable agrees with what refine_lbfgs() will actually polish", {
  # The flat emissions, both complete-data and FIML variants.
  expect_true(.is_refinable(categorical_model(2, type = "bernoulli")))
  expect_true(.is_refinable(categorical_model(2, type = "bernoulli_nan")))
  expect_true(.is_refinable(gaussian_model(2, type = "gaussian_diag")))
  expect_true(.is_refinable(gaussian_model(2, type = "gaussian_unit_nan")))

  expect_false(.is_refinable(categorical_model(2, type = "multinoulli", max_val = 3)))
  expect_false(.is_refinable(poisson_model(2, type = "poisson")))
  expect_false(.is_refinable(poisson_model(2, type = "poisson_nan")))
  expect_false(.is_refinable(lcga_model(2, design = .lcga_design(0:3, 1))))

  # A mixed measurement model is never refinable: its sub-models are
  # heterogeneous, so there is no single parameter matrix to pack.
  nested <- build_emission(list(a = list(model = "bernoulli", n_columns = 2),
                                b = list(model = "poisson",   n_columns = 2)),
                           n_components = 2)
  expect_false(.is_refinable(nested))

  # A block model inherits refinability from its sub-model, because the flat
  # view maps the blocks onto one wide matrix that L-BFGS can optimise directly.
  expect_true(.is_refinable(
    time_blocks_model(2, n_items = 3, n_times = 2, sub_model = "bernoulli")))
  expect_false(.is_refinable(
    time_blocks_model(2, n_items = 3, n_times = 2, sub_model = "poisson")))
})

test_that("an unpolished emission is given the tighter rule and reaches the optimum", {
  set.seed(4242)
  n <- 800
  cls <- sample(1:3, n, TRUE)
  rate <- rbind(c(0.4, 0.5, 3.0, 2.5, 0.6, 0.3),
                c(3.0, 2.8, 0.4, 0.5, 2.0, 2.4),
                c(1.2, 1.4, 1.3, 1.1, 1.0, 1.5))
  X <- matrix(rpois(n * 6, rate[cls, ]), n, 6)

  fit <- fit_mixture(X, n_classes = 3, measurement = "count",
                     n_init = 5, random_state = 7)

  # Under the old rule this stopped after 6 iterations, 7.0 of log-likelihood
  # short of the optimum. The count is the check that matters: no stopping rule
  # tight enough to be trusted here converges in single digits.
  expect_gt(fit$n_iter, 30)
  expect_true(fit$converged)

  # Rates are recovered to the accuracy the data allow. Under the old rule the
  # middle class came back at (1.43, 1.72, 1.38, 1.24, ...) against a truth of
  # (1.2, 1.4, 1.3, 1.1, ...) — the sort of error that changes how a class is
  # described in a paper.
  est <- fit$mm$parameters$rates
  ord <- order(est[, 1])
  expect_equal(est[ord, ], rate[order(rate[, 1]), ], tolerance = 0.15,
               ignore_attr = TRUE)
})

test_that("polytomous and mixed measurement models converge too", {
  skip_on_cran()

  set.seed(4242)
  n <- 800
  cls <- sample(1:3, n, TRUE)

  probs <- list(rbind(c(.7,.2,.1), c(.1,.2,.7), c(.3,.4,.3)),
                rbind(c(.6,.3,.1), c(.2,.2,.6), c(.3,.3,.4)),
                rbind(c(.8,.1,.1), c(.1,.3,.6), c(.4,.3,.3)),
                rbind(c(.5,.3,.2), c(.2,.3,.5), c(.3,.4,.3)),
                rbind(c(.7,.2,.1), c(.1,.4,.5), c(.3,.3,.4)))
  Xp <- sapply(probs, function(P)
    apply(P[cls, ], 1, function(p) sample(1:3, 1, prob = p)))

  # Polytomous EM is the slowest of the unpolished emissions: it needed several
  # hundred iterations here, against the seven the old rule allowed it.
  poly <- fit_mixture(Xp, n_classes = 3, measurement = "categorical",
                      n_init = 5, random_state = 7)
  expect_gt(poly$n_iter, 50)
  expect_true(poly$converged)

  Xb <- matrix(rbinom(n * 4, 1, c(.2, .8, .5)[cls]), n, 4)
  Xg <- matrix(rnorm(n * 3, c(-1, 1, 0)[cls]), n, 3)
  Xc <- matrix(rpois(n * 3, rbind(c(0.4, 0.5, 3.0),
                                  c(3.0, 2.8, 0.4),
                                  c(1.2, 1.4, 1.3))[cls, ]), n, 3)

  # A mixed model is unpolished even when every one of its parts would be
  # polished on its own, so it needs the tighter rule just as much.
  mixed <- fit_mixture(cbind(Xb, Xg, Xc), n_classes = 3,
                       measurement = list(binary = 1:4, continuous = 5:7,
                                          count = 8:10),
                       n_init = 5, random_state = 7)
  expect_gt(mixed$n_iter, 20)
  expect_true(mixed$converged)

  gm <- mixed$mm$models
  cont <- gm[[which(vapply(gm, function(s) !is.null(s$parameters$means),
                           logical(1)))]]$parameters$means
  expect_equal(sort(cont[, 1]), c(-1, 0, 1), tolerance = 0.15,
               ignore_attr = TRUE)
})

test_that("failing to converge within max_iter is reported, not swallowed", {
  set.seed(71)
  n <- 400
  cls <- sample(1:3, n, TRUE)
  X <- sapply(1:5, function(j)
    apply(rbind(c(.7, .2, .1), c(.1, .2, .7), c(.3, .4, .3))[cls, ], 1,
          function(p) sample(1:3, 1, prob = p)))

  # Under the loose rule these models "converged" in a handful of iterations, so
  # hitting the cap was near-impossible and print() was the only place the flag
  # appeared. Now that they run the iterations they need, a cap can genuinely
  # bite, and returning whatever iterate EM was on without saying so would be
  # the worst of the three outcomes.
  expect_warning(
    fit <- fit_mixture(X, n_classes = 3, measurement = "categorical",
                       n_init = 1, max_iter = 3, random_state = 1),
    "did not converge")
  expect_false(fit$converged)

  # And it stays quiet when the model does converge.
  expect_no_warning(
    fit_mixture(X, n_classes = 2, measurement = "categorical",
                n_init = 1, max_iter = 5000, random_state = 1))
})

test_that("polished emissions keep their loose rule and their previous answers", {
  set.seed(4242)
  n <- 800
  cls <- sample(1:3, n, TRUE)
  # Well-separated profiles, so that recovering them is a statement about the
  # estimator rather than about how much this particular design can identify.
  pis <- rbind(c(.05, .05, .05, .95, .95, .95),
               c(.95, .95, .95, .05, .05, .05),
               c(.95, .05, .95, .05, .95, .05))
  X <- matrix(rbinom(n * 6, 1, pis[cls, ]), n, 6)

  fit <- fit_mixture(X, n_classes = 3, measurement = "binary",
                     n_init = 5, random_state = 7)

  # EM stops early here by design and L-BFGS finishes the job, so a small
  # iteration count is the correct behaviour, not a symptom. If this starts
  # failing, the tighter rule has leaked onto the refined path and every binary
  # and continuous fit in the package just got slower for nothing.
  expect_lt(fit$n_iter, 60)

  # Classes are labelled arbitrarily, and two of these three share a value on
  # the first item, so they are matched by whole profile rather than sorted on
  # any one column.
  est <- fit$mm$parameters$pis
  idx <- apply(pis, 1, function(true_row)
    which.min(colSums((t(est) - true_row)^2)))
  expect_equal(length(unique(idx)), 3L)
  expect_equal(est[idx, ], pis, tolerance = 0.05, ignore_attr = TRUE)
})
