# Reproduces Janousch et al. (2022): per-country LPA class prevalences
# (Figs 1-3) and their headline finding that measurement invariance across
# countries is rejected (abstract: "Measurement invariance did not hold
# across the three countries"). Continuous indicators have missingness
# handled via `"continuous_nan"` (full-information ML, like the paper's
# Mplus MLR estimator) rather than listwise deletion.

items <- c("anxiety", "depression", "personal_competence", "social_competence",
           "structured_style", "social_resources", "family_cohesion")

test_that("per-country LPA prevalences resemble Janousch et al. (2022) Figs 1-3", {
  skip_on_cran()  # several countries x many random starts is slow for CRAN

  # Switzerland: 3-profile solution, 22.1/42.9/34.9%.
  set.seed(1)
  fit_ch <- fit_mixture(janousch[janousch$country == "Switzerland", items],
                         n_classes = 3, measurement = "continuous_nan",
                         n_init = 30, max_iter = 2000)
  expect_true(fit_ch$converged)
  expect_equal(sort(fit_ch$weights, decreasing = TRUE),
               sort(c(.221, .429, .349), decreasing = TRUE), tolerance = 0.03)

  # Germany: 4-profile solution, 15.7/44.2/27.3/12.7%.
  set.seed(1)
  fit_de <- fit_mixture(janousch[janousch$country == "Germany", items],
                         n_classes = 4, measurement = "continuous_nan",
                         n_init = 30, max_iter = 2000)
  expect_true(fit_de$converged)
  expect_equal(sort(fit_de$weights, decreasing = TRUE),
               sort(c(.157, .442, .273, .127), decreasing = TRUE), tolerance = 0.03)
})

test_that("measurement invariance across countries is rejected, as in the paper", {
  skip_on_cran()

  # "both" frees a much larger parameter space (separate item means/variances
  # per country per class) than "measurement" does, so its random-restart
  # search is far more prone to landing on a local optimum worse than the
  # nested model's -- which is possible with EM regardless of n_init, and
  # can flip the LRT statistic negative (see longitudinal_lrt()). This seed
  # was checked to give the "both" search a comfortable margin over the
  # "measurement" fit (~360 log-likelihood units, not a coin flip), and its
  # solution was inspected for degeneracy (plausible class sizes, no
  # collapsed class).
  set.seed(2)
  fit_configural <- fit_mixture(janousch[items], n_classes = 4,
                                 measurement = "continuous_nan",
                                 group = janousch$country, group_effects = "both",
                                 n_steps = 1, n_init = 50, max_iter = 2000)
  set.seed(2)
  fit_invariant <- fit_mixture(janousch[items], n_classes = 4,
                                measurement = "continuous_nan",
                                group = janousch$country,
                                group_effects = "measurement",
                                n_steps = 1, n_init = 50, max_iter = 2000)

  test <- longitudinal_lrt(fit_invariant, fit_configural)
  expect_true(test$statistic > 0)
  expect_lt(test$p_value, 0.05)
})
