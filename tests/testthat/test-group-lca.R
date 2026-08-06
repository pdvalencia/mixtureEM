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

test_that("longitudinal_lrt() rejects both false null hypotheses with correct df", {
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
  inv_test <- longitudinal_lrt(fit_prev, fit_both)
  expect_equal(inv_test$df, K * J * (G - 1))
  expect_lt(inv_test$p_value, 0.001)

  # Prevalence-equivalence test (sec. 5.11): prevalences genuinely differ by
  # group too, so this must also reject.
  prev_test <- longitudinal_lrt(fit_none, fit_prev)
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

  test <- longitudinal_lrt(restricted, full)
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
