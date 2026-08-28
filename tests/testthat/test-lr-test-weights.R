# ------------------------------------------------------------------------------
# lr_test() under sampling weights and a complex survey design.
#
# Ground truth this file checks against: unweighted c is ~1 (the sanity check
# that the trace formula is assembled right), a survey design changes the
# meat matrix (and so c), the scaled statistic is invariant to the weight
# scale, and a structural (covariate/group) fit refuses rather than silently
# returning an uncorrected number.
# ------------------------------------------------------------------------------

.lrw_sim <- function(n = 800, seed = 1) {
  set.seed(seed)
  p1 <- c(.85, .80, .75, .20, .15, .10)
  p2 <- c(.15, .20, .25, .80, .85, .90)
  cls <- sample(1:2, n, TRUE, prob = c(.5, .5))
  X <- t(vapply(cls, function(k)
    rbinom(6, 1, if (k == 1) p1 else p2), integer(6)))
  colnames(X) <- paste0("y", 1:6)
  list(X = X, n = n)
}

.lrw_pair <- function(d, ...) {
  fit2 <- suppressMessages(suppressWarnings(fit_mixture(
    d$X, n_classes = 2, measurement = "binary",
    n_init = 5, random_state = 1, ...)))
  fit3 <- suppressMessages(suppressWarnings(fit_mixture(
    d$X, n_classes = 3, measurement = "binary",
    n_init = 5, random_state = 1, ...)))
  list(restricted = fit2, full = fit3)
}

test_that("two nested unweighted fits still return the naive statistic", {
  d <- .lrw_sim()
  fp <- .lrw_pair(d)
  t <- lr_test(fp$restricted, fp$full)
  expect_false(isTRUE(t$scaled))
  expect_equal(t$statistic, t$statistic_raw)
  expect_equal(t$statistic, -2 * (fp$restricted$metrics$ll - fp$full$metrics$ll))
})

test_that("frequency weights do not trip the scaling-correction predicate", {
  d <- .lrw_sim()
  w <- sample(1:5, d$n, TRUE)
  fp <- .lrw_pair(d, weights = w, weight_type = "frequency")
  t <- lr_test(fp$restricted, fp$full)
  expect_false(isTRUE(t$scaled))
})

test_that("sampling weights trigger the scaled statistic, and unweighted c is near 1", {
  d <- .lrw_sim()
  w <- runif(d$n, 0.3, 3)
  fp <- .lrw_pair(d, weights = w, weight_type = "sampling")
  t <- lr_test(fp$restricted, fp$full)
  expect_true(isTRUE(t$scaled))
  expect_false(is.na(t$scaling_factor))
  expect_equal(t$statistic, t$statistic_raw / t$scaling_factor)

  # The unweighted sanity check the formula has to pass.
  info <- .nested_fit_info(fp$restricted)
  pieces_unweighted <- .scaling_pieces(.nested_fit_info(
    suppressMessages(suppressWarnings(fit_mixture(
      d$X, n_classes = 2, measurement = "binary",
      n_init = 5, random_state = 1)))))
  expect_lt(abs(pieces_unweighted$c - 1), 0.2)
})

test_that("a survey design changes the scaling factor, not just the weights", {
  d <- .lrw_sim()
  w <- runif(d$n, 0.3, 3)
  strata  <- rep(1:4, length.out = d$n)
  cluster <- rep(1:40, length.out = d$n)

  fit_w      <- suppressMessages(suppressWarnings(fit_mixture(
    d$X, n_classes = 2, measurement = "binary", n_init = 5, random_state = 1,
    weights = w)))
  fit_design <- suppressMessages(suppressWarnings(fit_mixture(
    d$X, n_classes = 2, measurement = "binary", n_init = 5, random_state = 1,
    weights = w, strata = strata, cluster = cluster)))

  c_w      <- .scaling_pieces(.nested_fit_info(fit_w))$c
  c_design <- .scaling_pieces(.nested_fit_info(fit_design))$c
  expect_false(isTRUE(all.equal(c_w, c_design)))
})

test_that("the scaled statistic is invariant to the weight scale", {
  d <- .lrw_sim()
  w <- runif(d$n, 0.3, 3)

  # The exact half, and the mechanism the invariance actually rests on:
  # `.resolve_weights()` rescales sampling weights to sum to n, so declaring
  # w or 5 * w hands the estimator the same vector. This is deterministic on
  # every platform and is what a real normalization bug would break.
  r1 <- .resolve_weights(w,     d$n, "sampling")$weights
  r5 <- .resolve_weights(5 * w, d$n, "sampling")$weights
  expect_lt(max(abs(r1 - r5)), 1e-12)

  # The end-to-end half. The two resolved vectors agree only to floating
  # point, not bit for bit, and EM amplifies that: the two runs can settle on
  # solutions separated by more than the perturbation, by a margin that
  # differs by platform (measured up to ~0.6% relative on macOS/aarch64 where
  # this file previously failed at 1e-4). The band below is therefore
  # empirical -- it is loose enough to absorb that amplification and still far
  # tighter than the factor-of-5 error a missing rescale would produce.
  fp1 <- .lrw_pair(d, weights = w, weight_type = "sampling")
  fp2 <- .lrw_pair(d, weights = 5 * w, weight_type = "sampling")
  t1 <- lr_test(fp1$restricted, fp1$full)
  t2 <- lr_test(fp2$restricted, fp2$full)
  expect_lt(abs(t1$statistic - t2$statistic), 2e-2 * abs(t1$statistic))
})

test_that("a weighted structural-model pair at n_steps = 1 returns the scaled statistic", {
  d <- .lrw_sim(n = 600)
  grp <- sample(c("A", "B"), d$n, TRUE)
  w   <- runif(d$n, 0.3, 3)

  fit_prev <- suppressMessages(suppressWarnings(fit_mixture(
    d$X, n_classes = 2, measurement = "binary", group = grp,
    group_effects = "prevalence", n_steps = 1, n_init = 5, random_state = 1,
    weights = w)))
  fit_both <- suppressMessages(suppressWarnings(fit_mixture(
    d$X, n_classes = 2, measurement = "binary", group = grp,
    group_effects = "both", n_steps = 1, n_init = 5, random_state = 1,
    weights = w)))

  # The packed step-one log-likelihood must reproduce metrics$ll exactly: that
  # is the check that the structural block, not just the measurement one, made
  # it into the packed vector.
  par <- .joint_pack(fit_prev)
  expect_equal(sum(fit_prev$sample_weights *
                     .joint_ll_case(fit_prev, fit_prev$data, par)),
              fit_prev$metrics$ll)

  t <- lr_test(fit_prev, fit_both)
  expect_true(isTRUE(t$scaled))
  expect_false(is.na(t$scaling_factor))
  expect_equal(t$statistic, t$statistic_raw / t$scaling_factor)

  # The unweighted sanity check, same as the plain-LCA case above.
  fit_prev0 <- suppressMessages(suppressWarnings(fit_mixture(
    d$X, n_classes = 2, measurement = "binary", group = grp,
    group_effects = "prevalence", n_steps = 1, n_init = 5, random_state = 1)))
  pieces0 <- .scaling_pieces(.nested_fit_info(fit_prev0))
  expect_lt(abs(pieces0$c - 1), 0.2)
})

test_that("an unsupported structural family (group_prevalence) still refuses to pack", {
  # .step1_pack_sm()/.joint_pack() cover the "covariate" family (predictors,
  # group_effects = "prevalence"/"both") and nothing else: group_prevalence
  # (R/group_prevalence.R) parameterises a constrained simplex per group
  # rather than an unconstrained beta, so it has no representation here and
  # must fall back to NULL (Option A) rather than silently mis-packing.
  ms <- list(n_steps = 1, Y = matrix(1L, 4, 1),
            mm = list(), sm = group_prevalence_model(2, 2))
  expect_null(.step1_pack_sm(ms$sm))
  expect_null(.joint_pack(ms))
})

test_that("a weighted structural-model pair at the default (3-step) n_steps refuses", {
  d <- .lrw_sim(n = 600)
  grp <- sample(c("A", "B"), d$n, TRUE)
  w   <- runif(d$n, 0.3, 3)

  fit_prev <- suppressMessages(suppressWarnings(fit_mixture(
    d$X, n_classes = 2, measurement = "binary", group = grp,
    group_effects = "prevalence", n_init = 5, random_state = 1, weights = w)))
  fit_both <- suppressMessages(suppressWarnings(fit_mixture(
    d$X, n_classes = 2, measurement = "binary", group = grp,
    group_effects = "both", n_init = 5, random_state = 1, weights = w)))

  expect_error(lr_test(fit_prev, fit_both), "step-3")
})

test_that("an unweighted structural-model pair is unaffected", {
  d <- .lrw_sim(n = 600)
  grp <- sample(c("A", "B"), d$n, TRUE)

  fit_prev <- suppressMessages(suppressWarnings(fit_mixture(
    d$X, n_classes = 2, measurement = "binary", group = grp,
    group_effects = "prevalence", n_steps = 1, n_init = 5, random_state = 1)))
  fit_both <- suppressMessages(suppressWarnings(fit_mixture(
    d$X, n_classes = 2, measurement = "binary", group = grp,
    group_effects = "both", n_steps = 1, n_init = 5, random_state = 1)))

  t <- lr_test(fit_prev, fit_both)
  expect_false(isTRUE(t$scaled))
  expect_gt(t$df, 0)
})

test_that("longitudinal_lrt() inherits the scaled statistic without its own code path", {
  d <- .lrw_sim()
  w <- runif(d$n, 0.3, 3)
  fp <- .lrw_pair(d, weights = w, weight_type = "sampling")
  new <- lr_test(fp$restricted, fp$full)
  expect_warning(old <- longitudinal_lrt(fp$restricted, fp$full), "deprecated")
  expect_equal(old$statistic, new$statistic)
  expect_equal(old$scaled, new$scaled)
})

# ------------------------------------------------------------------------------
# The anchor row is not necessarily row K.
#
# `.step1_pack_sm()` used to slice rows 1..K-1 out of `beta` and call the
# remainder the anchor, which is only correct when the all-zero row is the last
# one. `sort_model_classes()` reorders classes after estimation and leaves the
# zero row wherever it lands, so on a real fit it is often not last -- and every
# other test in this file happens to use a fit where it is, which is why the
# defect survived. These two tests fix a not-last anchor deliberately.
# ------------------------------------------------------------------------------

test_that(".step1_pack_sm() round-trips a beta whose zero row is not last", {
  set.seed(4)
  D <- 3L
  Z <- cbind(1, matrix(rnorm(20 * (D - 1L)), 20, D - 1L))

  # One model, written three ways: the softmax is invariant to adding a
  # constant vector to every row, so these differ only in which row is zero.
  base <- matrix(c(0.4, -0.7, 0.9,
                   -1.1, 0.3, 0.2,
                   0.0, 0.0, 0.0), nrow = 3L, byrow = TRUE)
  variants <- lapply(1:3, function(k) sweep(base, 2, base[k, ], "-"))
  expect_equal(which(rowSums(abs(variants[[2]])) == 0), 2L)

  probs <- lapply(variants, function(B) {
    sm <- covariate_model(3L)
    sm$parameters$beta <- B
    par <- .step1_pack_sm(sm)
    expect_length(par, (3L - 1L) * D)
    back <- .step1_unpack_sm(sm, par)$parameters$beta
    # The reconstruction restores the package's own last-row-zero convention,
    # not the input matrix -- equality of the fitted probabilities is the
    # invariant that matters, and it must hold for every variant.
    expect_equal(which(rowSums(abs(back)) == 0), 3L)
    softmax_rows(Z %*% t(back))
  })

  expect_equal(probs[[1]], softmax_rows(Z %*% t(base)))
  expect_equal(probs[[2]], probs[[1]])
  expect_equal(probs[[3]], probs[[1]])
})

test_that("the packed log-likelihood survives an anchor row that is not last", {
  d <- .lrw_sim(n = 600)
  grp <- sample(c("A", "B"), d$n, TRUE)
  w   <- runif(d$n, 0.3, 3)

  fit <- suppressMessages(suppressWarnings(fit_mixture(
    d$X, n_classes = 3, measurement = "binary", group = grp,
    group_effects = "prevalence", n_steps = 1, n_init = 5, random_state = 1,
    weights = w)))

  # Recentre beta on row 1. This changes no probability the model implies and
  # so leaves metrics$ll exactly what it was; all it does is move the all-zero
  # anchor row off the last position, which is the state a sorted fit can
  # arrive at on its own.
  B <- fit$sm$parameters$beta
  expect_gt(nrow(B), 2L)
  moved <- fit
  moved$sm$parameters$beta <- sweep(B, 2, B[1L, ], "-")
  expect_equal(which(rowSums(abs(moved$sm$parameters$beta)) == 0), 1L)

  # The assertion the packing bug failed by hundreds of log-likelihood units:
  # the packed vector has to describe the fit itself, not a different model.
  par <- .joint_pack(moved)
  expect_false(is.null(par))
  expect_equal(sum(moved$sample_weights *
                     .joint_ll_case(moved, moved$data, par)),
              fit$metrics$ll)

  # And the scaling factor computed at that point is the same one the
  # last-row-anchored spelling of the identical model gives.
  c_last  <- .scaling_pieces(.nested_fit_info(fit))$c
  c_moved <- .scaling_pieces(.nested_fit_info(moved))$c
  expect_lt(abs(c_last - c_moved), 1e-6 * abs(c_last))
})
