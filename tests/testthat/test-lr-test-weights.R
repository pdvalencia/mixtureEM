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
  fp1 <- .lrw_pair(d, weights = w, weight_type = "sampling")
  fp2 <- .lrw_pair(d, weights = 5 * w, weight_type = "sampling")
  t1 <- lr_test(fp1$restricted, fp1$full)
  t2 <- lr_test(fp2$restricted, fp2$full)
  # The invariance is exact given identical fitted parameters; two independent
  # EM refits at different weight scales land at slightly different BFGS
  # numerical paths (documented elsewhere as 8th-significant-figure drift), so
  # the tolerance allows for that noise without allowing a real normalization
  # bug to pass.
  expect_lt(abs(t1$statistic - t2$statistic), 1e-4 * abs(t1$statistic))
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
