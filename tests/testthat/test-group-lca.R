# Multiple-group LCA (`group=`/`group_effects=` on fit_mixture()), Collins &
# Lanza (2010) sec. 5.7-5.12. Validated against synthetic data with a known
# per-group generating model: two groups with different item-response
# probabilities AND different class prevalences.

.make_group_data <- function(seed = 7, n = 800) {
  set.seed(seed)
  J <- 4L; K <- 2L
  grp  <- factor(sample(c("A", "B"), n, replace = TRUE))
  prev <- list(A = c(.7, .3), B = c(.3, .7))
  pis  <- list(
    A = list(c1 = rep(.1, J), c2 = rep(.9, J)),
    B = list(c1 = rep(.2, J), c2 = rep(.8, J))
  )
  X <- matrix(NA_real_, n, J)
  for (i in seq_len(n)) {
    g  <- as.character(grp[i])
    cl <- sample(1:K, 1, prob = prev[[g]])
    p  <- if (cl == 1) pis[[g]]$c1 else pis[[g]]$c2
    X[i, ] <- rbinom(J, 1, p)
  }
  colnames(X) <- paste0("item", 1:J)
  list(X = X, grp = grp, J = J, K = K)
}

test_that("group_effects='prevalence' matches predictors= (Collins & Lanza sec. 6.10.2)", {
  d <- .make_group_data()
  fit_group <- fit_mixture(d$X, n_classes = d$K, measurement = "binary",
                           group = d$grp, group_effects = "prevalence",
                           n_steps = 1, n_init = 5, random_state = 1)
  fit_pred  <- fit_mixture(d$X, n_classes = d$K, measurement = "binary",
                           predictors = d$grp,
                           n_steps = 1, n_init = 5, random_state = 1)
  expect_equal(fit_group$metrics$ll, fit_pred$metrics$ll, tolerance = 1e-6)
  expect_equal(fit_group$metrics$n_params, fit_pred$metrics$n_params)
})

test_that("n_params accounting matches hand-derived formulas for each group_effects mode", {
  d <- .make_group_data()
  J <- d$J; K <- d$K; G <- 2L

  fit_none <- fit_mixture(d$X, n_classes = K, measurement = "binary",
                          group = d$grp, group_effects = "none",
                          n_steps = 1, n_init = 3, random_state = 1)
  expect_equal(fit_none$metrics$n_params, K * J + (K - 1))

  fit_prev <- fit_mixture(d$X, n_classes = K, measurement = "binary",
                          group = d$grp, group_effects = "prevalence",
                          n_steps = 1, n_init = 3, random_state = 1)
  expect_equal(fit_prev$metrics$n_params, K * J + (K - 1) * G)

  fit_meas <- fit_mixture(d$X, n_classes = K, measurement = "binary",
                          group = d$grp, group_effects = "measurement",
                          n_steps = 1, n_init = 3, random_state = 1)
  expect_equal(fit_meas$metrics$n_params, K * J * G + (K - 1))

  fit_both <- fit_mixture(d$X, n_classes = K, measurement = "binary",
                          group = d$grp, group_effects = "both",
                          n_steps = 1, n_init = 3, random_state = 1)
  expect_equal(fit_both$metrics$n_params, K * J * G + (K - 1) * G)

  fit_partial <- fit_mixture(d$X, n_classes = K, measurement = "binary",
                             group = d$grp, group_effects = "both",
                             group_invariant_items = 1:2,
                             n_steps = 1, n_init = 3, random_state = 1)
  n_inv <- 2L
  expect_equal(fit_partial$metrics$n_params,
               K * (n_inv + (J - n_inv) * G) + (K - 1) * G)
})

# ------------------------------------------------------------------------------
# `group_invariant_items` with a mixed (list) `measurement` -- Part 17.4 Item A.
# Previously refused outright ("Partial invariance is not defined for a mixed
# measurement block"); .item_param_cols()/.copy_item_params() assumed a flat
# sub-model and n_parameters.blocks() averaged a per-item cost that a mixed
# block does not actually have uniformly (see R/time_blocks.R). Two binary
# items and one trichotomous item, one binary item held equal across groups
# along with the trichotomous one, the other binary item left free.
# ------------------------------------------------------------------------------

.make_mixed_group_data <- function(seed = 13, n = 600) {
  set.seed(seed)
  K <- 2L
  grp  <- factor(sample(c("A", "B"), n, replace = TRUE))
  cl   <- sample(1:K, n, replace = TRUE)
  y1   <- rbinom(n, 1, ifelse(cl == 1, .2, .8))
  y2   <- rbinom(n, 1, ifelse(cl == 1, .3, .7))
  z1   <- sample(1:3, n, replace = TRUE,
                prob = c(1, 1, 1)) # marginal draw, class-independence not needed
  z1   <- ifelse(cl == 1, pmin(z1, 2), pmax(z1, 2))
  data.frame(y1 = y1, y2 = y2, z1 = z1, grp = grp)
}

test_that("n_parameters accounts for a mixed-measurement partial-invariance fit", {
  d <- .make_mixed_group_data()
  K <- 2L; G <- 2L
  fit <- fit_mixture(
    d[, c("y1", "y2", "z1")], n_classes = K,
    measurement = list(binary = c("y1", "y2"), categorical = "z1"),
    group = d$grp, group_effects = "both",
    group_invariant_items = c("y2", "z1"),
    n_steps = 1, n_init = 3, random_state = 1)

  # Per-item cost: K per binary item, K * (M - 1) per trichotomous item.
  # y2 and z1 held equal across the G blocks (once each); y1 free (G times).
  n_meas <- (K * 1L) + (K * 2L) + (K * G)
  expect_equal(fit$metrics$n_params, n_meas + (K - 1) * G)
})

test_that("lr_test() rejects both false null hypotheses with correct df", {
  d <- .make_group_data()
  K <- d$K; J <- d$J; G <- 2L

  fit_both <- fit_mixture(d$X, n_classes = K, measurement = "binary",
                          group = d$grp, group_effects = "both",
                          n_steps = 1, n_init = 10, random_state = 1)
  fit_prev <- fit_mixture(d$X, n_classes = K, measurement = "binary",
                          group = d$grp, group_effects = "prevalence",
                          n_steps = 1, n_init = 10, random_state = 1)
  fit_none <- fit_mixture(d$X, n_classes = K, measurement = "binary",
                          group = d$grp, group_effects = "none",
                          n_steps = 1, n_init = 10, random_state = 1)

  # Measurement-invariance test (sec. 5.8): item-response probabilities
  # genuinely differ by group in the generating model, so this must reject.
  inv_test <- lr_test(fit_prev, fit_both)
  expect_equal(inv_test$df, K * J * (G - 1))
  expect_lt(inv_test$p_value, 0.001)

  # Prevalence-equivalence test (sec. 5.11): prevalences genuinely differ by
  # group too, so this must also reject.
  prev_test <- lr_test(fit_none, fit_prev)
  expect_equal(prev_test$df, K - 1)
  expect_lt(prev_test$p_value, 0.001)
})

test_that("group_effects='both' recovers each group's true class-membership probabilities", {
  d <- .make_group_data()
  fit <- fit_mixture(d$X, n_classes = d$K, measurement = "binary",
                     group = d$grp, group_effects = "both",
                     n_steps = 1, n_init = 10, random_state = 1)

  Z      <- cbind(Intercept = 1, fit$group_info$design)
  beta   <- fit$sm$parameters$beta
  logits <- Z %*% t(beta)
  probs  <- exp(logits) / rowSums(exp(logits))

  # Match each group's fitted classes to the truth by item-response pattern
  # (class labels are not guaranteed to align across group blocks; see the
  # group_effects docs), then compare recovered prevalences to the truth.
  is_high_endorsing <- function(block) {
    which.max(rowMeans(block$parameters$pis))
  }
  block_names <- names(fit$mm$models)
  truth <- list(A = c(low = .7, high = .3), B = c(low = .3, high = .7))

  for (i in seq_along(block_names)) {
    g <- substring(block_names[i], 2, 2)  # "G1" -> "1" is the factor level index
    level <- fit$group_info$levels[as.integer(g)]
    hi <- is_high_endorsing(fit$mm$models[[i]])
    row1 <- which(as.character(d$grp) == level)[1]
    p_hi <- unname(probs[row1, hi])
    expect_equal(p_hi, unname(truth[[level]]["high"]), tolerance = 0.1)
  }
})

# ------------------------------------------------------------------------------
# Warm-starting the configural search
# ------------------------------------------------------------------------------
#
# A group-varying measurement model ties its per-group blocks together only
# through the class labels, and random starting values give each block its own
# arbitrary labelling, so the search can settle below the restriction nested
# inside it -- which is impossible at the optimum and makes the invariance test
# a lower bound. `group_effects = "both"`/`"measurement"` therefore runs one
# extra restart built from fitted values rather than drawn: each group's own
# fit, with its classes permuted to match the pooled solution. See
# R/group_blocks.R.

test_that("a group-varying measurement model never scores below its own restriction", {
  d <- .make_group_data()

  restricted <- fit_mixture(d$X, n_classes = d$K, measurement = "binary",
                            group = d$grp, group_effects = "prevalence",
                            n_steps = 1, n_init = 5, random_state = 1)
  full       <- fit_mixture(d$X, n_classes = d$K, measurement = "binary",
                            group = d$grp, group_effects = "both",
                            n_steps = 1, n_init = 5, random_state = 1)

  # The nesting the invariance LRT rests on. Without the warm start this held
  # only by luck of the draw.
  expect_gt(full$metrics$ll, restricted$metrics$ll)

  test <- lr_test(restricted, full)
  expect_gt(test$statistic, 0)
  expect_equal(test$df,
               full$metrics$n_params - restricted$metrics$n_params)
})

test_that("the warm start recovers the per-group generating parameters", {
  # The data are generated with genuinely different item probabilities by group
  # (.1/.9 in A, .2/.8 in B), which is the case the configural model exists to
  # detect and the one a mis-aligned search blurs into a single pooled answer.
  d   <- .make_group_data()
  fit <- fit_mixture(d$X, n_classes = d$K, measurement = "binary",
                     group = d$grp, group_effects = "both",
                     n_steps = 1, n_init = 5, random_state = 1)

  pis <- lapply(fit$mm$models, function(m) m$parameters$pis)
  expect_length(pis, 2L)

  # Order the classes by their first item so the comparison does not depend on
  # which class the sorting happened to put first.
  low  <- vapply(pis, function(p) min(p[, 1]), numeric(1))
  high <- vapply(pis, function(p) max(p[, 1]), numeric(1))

  expect_lt(max(abs(low  - c(.1, .2))), 0.08)
  expect_lt(max(abs(high - c(.9, .8))), 0.08)
})

test_that("the warm start is skipped, not fatal, when a group is too small", {
  d <- .make_group_data()
  # One group with fewer cases than 2 * n_classes: its own fit is not attempted
  # and the block falls back to the pooled parameters.
  grp <- as.character(d$grp)
  grp[seq_len(3L)] <- "C"
  grp[-seq_len(3L)] <- "A"
  grp <- factor(grp)

  expect_no_error(
    fit <- fit_mixture(d$X, n_classes = d$K, measurement = "binary",
                       group = grp, group_effects = "measurement",
                       n_steps = 1, n_init = 3, random_state = 1))
  expect_equal(fit$mm$n_blocks, 2L)
})

test_that("longitudinal_lrt() still works but is deprecated", {
  d <- .make_group_data()
  fit_prev <- fit_mixture(d$X, n_classes = d$K, measurement = "binary",
                          group = d$grp, group_effects = "prevalence",
                          n_steps = 1, n_init = 3, random_state = 1)
  fit_both <- fit_mixture(d$X, n_classes = d$K, measurement = "binary",
                          group = d$grp, group_effects = "both",
                          n_steps = 1, n_init = 3, random_state = 1)

  expect_warning(old <- longitudinal_lrt(fit_prev, fit_both), "deprecated")
  new <- lr_test(fit_prev, fit_both)
  expect_s3_class(old, "lr_test")
  expect_equal(old$statistic, new$statistic)
  expect_equal(old$df, new$df)
})

test_that("a configural fit is not worse than the sum of its separable parts", {
  # With measurement AND prevalences free by group, nothing is shared across
  # groups, so the joint log-likelihood is exactly the sum of each group's own
  # K-class LCA. A joint fit scoring below that sum is a search or estimator
  # failure, not a modelling result. This caught the L-BFGS refinement running
  # on a model whose class priors come from a regression it cannot represent.
  set.seed(21)
  mk <- function(n, p, mix) {
    z <- sample(seq_len(nrow(p)), n, TRUE, prob = mix)
    matrix(rbinom(n * ncol(p), 1, p[z, ]), n, ncol(p))
  }
  pA <- matrix(c(.9, .85, .8, .15,  .1, .2, .15, .9), 2, 4, byrow = TRUE)
  pB <- matrix(c(.7, .2, .9, .25,   .2, .8, .1, .75), 2, 4, byrow = TRUE)
  XA <- mk(300, pA, c(.8, .2))
  XB <- mk(300, pB, c(.25, .75))
  X  <- rbind(XA, XB)
  colnames(X) <- paste0("i", seq_len(4))
  g  <- factor(rep(c("A", "B"), each = 300))

  per_group <- sum(vapply(list(XA, XB), function(z)
    fit_mixture(z, n_classes = 2, measurement = "binary", n_init = 10,
                random_state = 5)$metrics$ll, numeric(1)))
  joint <- fit_mixture(X, n_classes = 2, measurement = "binary", group = g,
                       group_effects = "both", n_steps = 1, n_init = 10,
                       random_state = 5)

  # Tolerance absorbs EM's own stopping rule on each of the three fits, not a
  # difference in what they converge to.
  expect_gt(joint$metrics$ll, per_group - 0.01)
})

test_that("a group-only fit says what the 3-step default costs", {
  d <- .make_group_data(n = 300)

  # `group` reaches the structural engine through the prevalence effect, so the
  # user gets a 3-step fit without asking for one. The message has to name the
  # consequence: metrics$ll is then on the structural model's scale.
  expect_message(
    fit <- suppressWarnings(fit_mixture(d$X, n_classes = d$K,
                                        measurement = "binary", group = d$grp,
                                        group_effects = "prevalence",
                                        n_init = 3, random_state = 1)),
    "n_steps = 1")

  # An explicit n_steps is left alone.
  expect_no_message(
    suppressWarnings(fit_mixture(d$X, n_classes = d$K, measurement = "binary",
                                 group = d$grp, group_effects = "measurement",
                                 n_steps = 1, n_init = 3, random_state = 1)))
})

# ------------------------------------------------------------------------------
# `group_prevalence_equal`: freezing one or more classes' prevalence to a
# single shared value across groups, on the probability scale directly
# (Collins & Lanza sec. 5.11-5.12's restricted models). See
# internal/ROADMAP.md Part 17.4 Item B and R/group_prevalence.R.
# ------------------------------------------------------------------------------

# Three groups, three classes, with class 2 genuinely pinned at the same
# prevalence (.30) in every group while classes 1 and 3 trade off against each
# other -- the one thing `group_effects = "prevalence"`'s ordinary regression
# route cannot represent (a zero coefficient pins a class to the *reference
# group's* share, not to one shared across every group).
.make_pinned_group_data <- function(seed = 11, n = 1500) {
  set.seed(seed)
  J <- 6L; K <- 3L
  grp  <- factor(sample(c("A", "B", "C"), n, replace = TRUE))
  prev <- list(A = c(.55, .30, .15), B = c(.45, .30, .25), C = c(.35, .30, .35))
  pis  <- rbind(rep(.10, J), rep(.50, J), rep(.90, J))
  X <- matrix(NA_real_, n, J)
  trueclass <- integer(n)
  for (i in seq_len(n)) {
    g  <- as.character(grp[i])
    cl <- sample(1:K, 1, prob = prev[[g]])
    trueclass[i] <- cl
    X[i, ] <- rbinom(J, 1, pis[cl, ])
  }
  colnames(X) <- paste0("item", 1:J)
  list(X = X, grp = grp, J = J, K = K, G = 3L, trueclass = trueclass, prev = prev)
}

test_that("group_prevalence_equal is numerically equivalent to the regression route when unconstrained", {
  d <- .make_group_data()
  fit_reg <- fit_mixture(d$X, n_classes = d$K, measurement = "binary",
                         group = d$grp, group_effects = "prevalence",
                         n_steps = 1, n_init = 10, random_state = 1)
  fit_gp  <- fit_mixture(d$X, n_classes = d$K, measurement = "binary",
                         group = d$grp, group_effects = "prevalence",
                         group_prevalence_equal = integer(0),
                         n_steps = 1, n_init = 10, random_state = 1)
  expect_equal(fit_gp$metrics$ll, fit_reg$metrics$ll, tolerance = 1e-4)
  expect_equal(fit_gp$metrics$n_params, fit_reg$metrics$n_params)
})

test_that("n_parameters.group_prevalence matches the hand-derived formula, including the |S| = K special case", {
  d <- .make_pinned_group_data()
  K <- d$K; G <- d$G

  fit_none <- fit_mixture(d$X, n_classes = K, measurement = "binary",
                          group = d$grp, group_effects = "prevalence",
                          group_prevalence_equal = integer(0),
                          n_steps = 1, n_init = 5, random_state = 1)
  fit_one  <- fit_mixture(d$X, n_classes = K, measurement = "binary",
                          group = d$grp, group_effects = "prevalence",
                          group_prevalence_equal = 2,
                          n_steps = 1, n_init = 5, random_state = 1)
  fit_all  <- fit_mixture(d$X, n_classes = K, measurement = "binary",
                          group = d$grp, group_effects = "prevalence",
                          group_prevalence_equal = TRUE,
                          n_steps = 1, n_init = 5, random_state = 1)

  n_meas <- K * d$J
  expect_equal(fit_none$metrics$n_params, n_meas + G * (K - 1))
  expect_equal(fit_one$metrics$n_params,  n_meas + G * (K - 1) - 1 * (G - 1))
  expect_equal(fit_all$metrics$n_params,  n_meas + (K - 1))

  # D_ALL = the pooled (no-group) model: freezing every class is the same
  # restriction as ignoring the grouping entirely.
  fit_pooled <- fit_mixture(d$X, n_classes = K, measurement = "binary",
                            n_steps = 1, n_init = 5, random_state = 1)
  expect_equal(fit_all$metrics$n_params, fit_pooled$metrics$n_params)
  expect_equal(fit_all$metrics$ll, fit_pooled$metrics$ll, tolerance = 1e-3)
})

test_that("group_prevalence_equal can recover a class genuinely pinned across groups", {
  # Class labels are not identified independently of the measurement model
  # (R/group_blocks.R's own warm-start comment makes the same point for the
  # measurement side): a thorough search is free to relabel classes so that
  # whichever real class is cheapest to freeze ends up in the requested slot,
  # for *any* requested slot. That is a property of mixture-model label
  # switching in general, not of this emission, so this test does not assume
  # "freeze slot k" names the same real class across three separate fits.
  # What it does assert is the substantive claim the M-step formula has to
  # get right: freezing *some* class can cost nothing (the data really do
  # have one), the LRT df matches the hand-derived formula, and the fit that
  # finds the free lunch is the one whose frozen share actually lands on the
  # true generating value.
  d <- .make_pinned_group_data()
  K <- d$K; G <- d$G

  fit_free <- fit_mixture(d$X, n_classes = K, measurement = "binary",
                          group = d$grp, group_effects = "prevalence",
                          group_prevalence_equal = integer(0),
                          n_steps = 1, n_init = 20, random_state = 1)

  fits_frozen <- lapply(1:K, function(k)
    suppressWarnings(fit_mixture(d$X, n_classes = K, measurement = "binary",
                                 group = d$grp, group_effects = "prevalence",
                                 group_prevalence_equal = k,
                                 n_steps = 1, n_init = 20, random_state = 1)))

  tests <- lapply(fits_frozen, lr_test, full = fit_free)
  # Every frozen-class test compares nested models over the same `G - 1`
  # df, whichever class ends up in the requested slot.
  for (t in tests) expect_equal(t$df, G - 1)

  p_values <- vapply(tests, function(t) t$p_value, numeric(1))
  best <- which.max(p_values)
  expect_gt(p_values[best], 0.05)

  # The frozen share the best-fitting restriction settled on should sit near
  # the true pinned value (.30 in every group), not merely be internally
  # self-consistent.
  best_fit     <- fits_frozen[[best]]
  frozen_share <- best_fit$sm$parameters$gamma[1, best_fit$sm$frozen]
  expect_lt(abs(frozen_share - 0.30), 0.05)
})

test_that("group_prevalence_equal rejects an out-of-range class index", {
  d <- .make_group_data()
  expect_error(
    fit_mixture(d$X, n_classes = d$K, measurement = "binary",
               group = d$grp, group_effects = "prevalence",
               group_prevalence_equal = d$K + 1,
               n_steps = 1, n_init = 2),
    "group_prevalence_equal"
  )
})

test_that("group_prevalence_equal requires a prevalence effect and forbids predictors", {
  d <- .make_group_data()
  expect_error(
    fit_mixture(d$X, n_classes = d$K, measurement = "binary",
               group = d$grp, group_effects = "measurement",
               group_prevalence_equal = 1, n_init = 2),
    "group_prevalence_equal"
  )
  expect_error(
    fit_mixture(d$X, n_classes = d$K, measurement = "binary",
               group = d$grp, group_effects = "prevalence",
               predictors = rnorm(nrow(d$X)),
               group_prevalence_equal = 1, n_init = 2),
    "group_prevalence_equal"
  )
})

test_that("class_sizes() reports each group's own class sizes for a group-prevalence fit", {
  d <- .make_pinned_group_data(n = 900)
  fit <- fit_mixture(d$X, n_classes = d$K, measurement = "binary",
                     group = d$grp, group_effects = "prevalence",
                     group_prevalence_equal = integer(0),
                     n_steps = 1, n_init = 10, random_state = 1)

  by_group <- attr(class_sizes(fit), "by_group")
  expect_equal(nrow(by_group), d$K * d$G)
  expect_equal(sort(unique(as.character(by_group$group))), sort(levels(d$grp)))

  # Each group's per-class proportions sum to one.
  totals <- tapply(by_group$proportion, by_group$group, sum)
  expect_equal(as.numeric(totals), rep(1, d$G), tolerance = 1e-6)
})

test_that("print() points to ll_knownclass only for a group= fit", {
  d <- .make_group_data()
  fit_group <- fit_mixture(d$X, n_classes = d$K, measurement = "binary",
                           group = d$grp, group_effects = "prevalence",
                           n_steps = 1, n_init = 5, random_state = 1)
  out_group <- capture.output(print(fit_group))
  expect_true(any(grepl("ll_knownclass", out_group, fixed = TRUE)))

  fit_plain <- suppressWarnings(fit_mixture(d$X, n_classes = d$K,
                                            measurement = "binary",
                                            n_init = 5, random_state = 1))
  out_plain <- capture.output(print(fit_plain))
  expect_false(any(grepl("ll_knownclass", out_plain, fixed = TRUE)))
})
