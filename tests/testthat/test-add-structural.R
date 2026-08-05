# add_covariates() / add_outcome(): stepwise analyses on a fitted model.
#
# The core contract is exact equivalence with the one-call fit_mixture()
# three-step path: given the same seed, running steps 2-3 on the stored
# step-1 solution must reproduce the same user-visible results without
# re-estimating the measurement model.

.sim_step3_data <- function(n = 250, seed = 42) {
  set.seed(seed)
  cl    <- rbinom(n, 1, 0.4)
  items <- sapply(1:5, function(j) rbinom(n, 1, ifelse(cl == 1, 0.8, 0.2)))
  colnames(items) <- paste0("item", 1:5)
  covs <- data.frame(
    age  = rnorm(n) + 0.8 * cl,
    sexo = factor(ifelse(rbinom(n, 1, stats::plogis(cl - 0.5)) == 1, "F", "M"))
  )
  bmi <- rnorm(n, 25 + 2 * cl)
  grp <- factor(ifelse(rbinom(n, 1, stats::plogis(0.6 * cl)) == 1, "a", "b"))
  list(items = items, covs = covs, bmi = bmi, grp = grp)
}

test_that("add_covariates reproduces the one-call 3-step ML fit exactly", {
  d <- .sim_step3_data()

  fitA <- suppressMessages(suppressWarnings(
    fit_mixture(d$items, n_classes = 2, predictors = d$covs,
                n_steps = 3, correction = "ML", n_init = 5, random_state = 9)))
  fit0 <- suppressMessages(suppressWarnings(
    fit_mixture(d$items, n_classes = 2, n_init = 5, random_state = 9)))
  fitB <- suppressMessages(add_covariates(fit0, d$covs))

  # Same solution, same class order.
  expect_equal(fitA$weights, fitB$weights, tolerance = 1e-10)
  expect_equal(fitA$step1_metrics$ll, fitB$step1_metrics$ll, tolerance = 1e-10)
  expect_equal(fitA$metrics$ll, fitB$metrics$ll, tolerance = 1e-8)

  # The raw beta parameterization may anchor a different class (the one-call
  # path applies the correction before class sorting), but every user-visible
  # quantity must match: odds ratios, and the full printed summary.
  expect_equal(coef(fitA), coef(fitB), tolerance = 1e-10)
  expect_identical(capture.output(suppressMessages(summary(fitA))),
                   capture.output(suppressMessages(summary(fitB))))

  # Term grouping for the omnibus Wald test survives the add path.
  expect_identical(fitA$sm$parameters$terms, fitB$sm$parameters$terms)

  # Metadata reflects the stepwise estimation.
  expect_identical(fitB$n_steps, 3L)
  expect_identical(fitB$correction, "ML")
})

test_that("add_outcome reproduces one-call BCH and ML outcome fits", {
  d <- .sim_step3_data()

  fit0 <- suppressMessages(suppressWarnings(
    fit_mixture(d$items, n_classes = 2, n_init = 5, random_state = 9)))

  # Continuous outcome, BCH (the auto default).
  fitC <- suppressMessages(suppressWarnings(
    fit_mixture(d$items, n_classes = 2, outcome = d$bmi,
                n_steps = 3, correction = "BCH", n_init = 5, random_state = 9)))
  fitD <- suppressMessages(add_outcome(fit0, d$bmi))
  expect_identical(fitD$correction, "BCH")
  expect_equal(fitC$sm$parameters$means, fitD$sm$parameters$means,
               tolerance = 1e-10)
  expect_identical(capture.output(suppressMessages(summary(fitC))),
                   capture.output(suppressMessages(summary(fitD))))

  # Categorical outcome, ML (the auto default).
  fitE <- suppressMessages(suppressWarnings(
    fit_mixture(d$items, n_classes = 2, outcome = d$grp,
                n_steps = 3, correction = "ML", n_init = 5, random_state = 9)))
  fitF <- suppressMessages(suppressWarnings(add_outcome(fit0, d$grp)))
  expect_identical(fitF$correction, "ML")
  expect_identical(capture.output(suppressMessages(summary(fitE))),
                   capture.output(suppressMessages(summary(fitF))))
})

test_that("a conditional fit can be reused: the structural model is replaced", {
  d <- .sim_step3_data()

  fit0  <- suppressMessages(suppressWarnings(
    fit_mixture(d$items, n_classes = 2, n_init = 5, random_state = 9)))
  fitB  <- suppressMessages(add_covariates(fit0, d$covs))

  # Reference: outcome added to the clean unconditional fit.
  fitD  <- suppressMessages(add_outcome(fit0, d$bmi))
  # Same outcome added to the conditional fit: must announce the replacement
  # and produce the same result (the ML correction contaminates $log_resp,
  # so this exercises the posterior-restoration path).
  expect_message(fitG <- suppressWarnings(add_outcome(fitB, d$bmi)),
                 "Replacing the existing structural model")
  expect_equal(fitD$sm$parameters$means, fitG$sm$parameters$means,
               tolerance = 1e-10)
  expect_equal(fitD$metrics$ll, fitG$metrics$ll, tolerance = 1e-8)
})

test_that("covariate rows are aligned with cases removed at fitting time", {
  d <- .sim_step3_data(n = 200)
  items <- d$items
  items[c(3, 117), ] <- NA   # two cases with no observed indicator

  fit0 <- suppressWarnings(suppressMessages(
    fit_mixture(items, n_classes = 2, n_init = 3, random_state = 5)))
  expect_identical(fit0$missing_data$n_empty_rows, 2L)

  # Full-length covariates: matching rows dropped, with a message.
  expect_message(
    fit_full <- suppressWarnings(add_covariates(fit0, d$covs,
                                                correction = "ML")),
    "matching rows")
  # Already-aligned covariates: accepted silently.
  covs_kept <- d$covs[-c(3, 117), , drop = FALSE]
  fit_kept  <- suppressMessages(suppressWarnings(
    add_covariates(fit0, covs_kept, correction = "ML")))
  expect_equal(coef(fit_full), coef(fit_kept), tolerance = 1e-10)

  # Factor term grouping survives the row subsetting.
  expect_true("sexo" %in% fit_full$sm$parameters$terms)

  # Any other length is an error that reports both counts.
  expect_error(add_covariates(fit0, d$covs[1:50, ]), "50 rows")
})

test_that("add_covariates and add_outcome reject unusable fits", {
  d <- .sim_step3_data(n = 150)

  # Not a fitted model.
  expect_error(add_covariates(list(), d$covs), "must be a fitted model")

  # LTA models have their own covariate interface.
  fake_lta <- structure(list(), class = "lta_model")
  expect_error(add_covariates(fake_lta, d$covs), "fit_lta")

  # One-step conditional fits cannot serve as a step-1 solution.
  fit1 <- suppressMessages(suppressWarnings(
    fit_mixture(d$items, n_classes = 2, predictors = d$covs,
                n_steps = 1, n_init = 3, random_state = 5)))
  expect_error(add_outcome(fit1, d$bmi), "one step")

  # Group-as-predictor fits must go through fit_mixture(group=, predictors=).
  fitg <- suppressMessages(suppressWarnings(
    fit_mixture(d$items, n_classes = 2, group = d$grp,
                group_effects = "prevalence", n_steps = 3, n_init = 3,
                random_state = 5)))
  expect_error(add_covariates(fitg, d$covs), "group")

  # Missing structural data.
  fit0 <- suppressMessages(suppressWarnings(
    fit_mixture(d$items, n_classes = 2, n_init = 3, random_state = 5)))
  expect_error(add_covariates(fit0), "`predictors` is required")
  expect_error(add_outcome(fit0), "`outcome` is required")
})
