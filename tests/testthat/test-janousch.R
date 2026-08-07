# Reproduces Janousch et al. (2022): per-country LPA class prevalences
# (Figs 1-3) and their headline finding that measurement invariance across
# countries is rejected (abstract: "Measurement invariance did not hold
# across the three countries"). The continuous indicators are missing for many
# cases; that is handled by full-information ML, which `fit_mixture()` selects
# on its own, rather than by listwise deletion.

items <- c("anxiety", "depression", "personal_competence", "social_competence",
           "structured_style", "social_resources", "family_cohesion")

test_that("per-country LPA prevalences resemble Janousch et al. (2022) Figs 1-3", {
  skip_on_cran()  # several countries x many random starts is slow for CRAN

  # Switzerland: 3-profile solution, 22.1/42.9/34.9%.
  set.seed(1)
  fit_ch <- fit_mixture(janousch[janousch$country == "Switzerland", items],
                         n_classes = 3, measurement = "continuous",
                         n_init = 30, max_iter = 2000)
  expect_true(fit_ch$converged)
  expect_equal(sort(fit_ch$weights, decreasing = TRUE),
               sort(c(.221, .429, .349), decreasing = TRUE), tolerance = 0.03)

  # Germany: 4-profile solution, 15.7/44.2/27.3/12.7%.
  set.seed(1)
  fit_de <- fit_mixture(janousch[janousch$country == "Germany", items],
                         n_classes = 4, measurement = "continuous",
                         n_init = 30, max_iter = 2000)
  expect_true(fit_de$converged)
  expect_equal(sort(fit_de$weights, decreasing = TRUE),
               sort(c(.157, .442, .273, .127), decreasing = TRUE), tolerance = 0.03)

  # Neither per-country fit is degenerate. This is the other half of the
  # default's job: the prior has to be weak enough to leave a healthy
  # small-sample solution where the paper found it. Germany is the binding
  # case, at n = 342 with four classes.
  expect_null(fit_ch$degenerate)
  expect_null(fit_de$degenerate)
})

test_that("a collapsed class variance is detected and named", {
  skip_on_cran()

  # Pooling the countries and freeing only the prevalences leaves one class
  # free to latch onto the cases sitting at the ceiling of `social_resources`,
  # driving that class's variance on that item towards zero. At the default
  # prior strength this is a *finite* spurious optimum rather than a runaway
  # one, which is what makes it detectable: the variance lands orders of
  # magnitude below the item's own spread with the class mean pinned at the
  # top of the scale.
  expect_warning(
    fit <- fit_mixture(janousch[items], n_classes = 4,
                       measurement = "continuous",
                       group = janousch$country, group_effects = "prevalence",
                       n_steps = 1, n_init = 50, max_iter = 2000,
                       random_state = 11),
    "collapsed towards zero")

  expect_false(is.null(fit$degenerate))
  expect_true("social_resources" %in% fit$degenerate$item)
  # The signature, not merely the flag: a variance far below the item's
  # marginal, and a class mean at a data boundary.
  expect_lt(min(fit$degenerate$ratio), 0.01)
  expect_true(any(fit$degenerate$pinned))

  # And the remedy the warning recommends actually works. The warning's starting
  # point is one artificial observation per class -- 4 here -- and it says to
  # raise it if the warning persists. On this data it does persist at 4, so 5 is
  # the value that clears it; the rule is a place to start plus a check, not a
  # number that always works first time.
  fit5 <- fit_mixture(janousch[items], n_classes = 4,
                      measurement = "continuous",
                      group = janousch$country, group_effects = "prevalence",
                      n_steps = 1, n_init = 50, max_iter = 2000,
                      random_state = 11,
                      bayes_constants = list(variances = 5))
  expect_null(fit5$degenerate)
  # `"prevalence"` pools the measurement model, so the variances live on a flat
  # emission rather than in per-group blocks; the accessor covers both shapes.
  expect_gt(min(.gaussian_variance_cells(fit5$mm)$variance), 0.01)
})

test_that("measurement invariance across countries is rejected, as in the paper", {
  skip_on_cran()

  # Collins & Lanza's measurement-invariance test (sec. 5.8) compares a model
  # in which the item parameters are held equal across groups against one in
  # which they are free. Here that is `"prevalence"` (measurement pooled, each
  # country its own class prevalences) against `"both"` (each country its own
  # measurement model as well) -- NOT `"measurement"`, which frees the item
  # parameters and pools the prevalences, and so tests the opposite
  # restriction.
  #
  # Both fits use `variances = 5`, the strength that clears the collapse on this
  # data (see the test above), because the default leaves it degenerate. They
  # must share a prior for the comparison to be a test at all: testing two
  # solutions estimated under different constraints, or two degenerate ones
  # against each other, produces a number but not a test.
  fit_invariant <- fit_mixture(janousch[items], n_classes = 4,
                               measurement = "continuous",
                               group = janousch$country,
                               group_effects = "prevalence",
                               n_steps = 1, n_init = 50, max_iter = 2000,
                               random_state = 11,
                               bayes_constants = list(variances = 5))
  fit_configural <- fit_mixture(janousch[items], n_classes = 4,
                                measurement = "continuous",
                                group = janousch$country,
                                group_effects = "both",
                                n_steps = 1, n_init = 50, max_iter = 2000,
                                random_state = 11,
                                bayes_constants = list(variances = 5))

  expect_null(fit_invariant$degenerate)
  expect_null(fit_configural$degenerate)

  test <- lr_test(fit_invariant, fit_configural)
  expect_equal(test$df, 112)
  expect_true(test$statistic > 0)
  expect_lt(test$p_value, 0.05)
})

test_that("lr_test() refuses to interpret a degenerate fit", {
  skip_on_cran()

  suppressWarnings(
    fit_invariant <- fit_mixture(janousch[items], n_classes = 4,
                                 measurement = "continuous",
                                 group = janousch$country,
                                 group_effects = "prevalence",
                                 n_steps = 1, n_init = 50, max_iter = 2000,
                                 random_state = 11))
  fit_configural <- fit_mixture(janousch[items], n_classes = 4,
                                measurement = "continuous",
                                group = janousch$country,
                                group_effects = "both",
                                n_steps = 1, n_init = 50, max_iter = 2000,
                                random_state = 11,
                                bayes_constants = list(variances = 5))

  # The restricted fit is degenerate, so the test is meaningless in either
  # direction and must say so rather than printing a statistic.
  expect_warning(lr_test(fit_invariant, fit_configural),
                 "collapsed class variance")
})
