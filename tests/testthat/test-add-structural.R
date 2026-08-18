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
    fit_mixture(d$items, n_classes = 2,
                measurement = "binary",
                predictors = d$covs,
                n_steps = 3, correction = "ML", n_init = 5, random_state = 9)))
  fit0 <- suppressMessages(suppressWarnings(
    fit_mixture(d$items, n_classes = 2,
                measurement = "binary",
                n_init = 5, random_state = 9)))
  fitB <- suppressMessages(add_covariates(fit0, d$covs))

  # Same solution, same class order.
  expect_equal(fitA$weights, fitB$weights, tolerance = 1e-10)
  expect_equal(fitA$step1_metrics$ll, fitB$step1_metrics$ll, tolerance = 1e-10)
  expect_equal(fitA$metrics$ll, fitB$metrics$ll, tolerance = 1e-8)

  # The raw beta parameterization may anchor a different class (the one-call
  # path applies the correction before class sorting), but every user-visible
  # quantity must match: odds ratios, and the full printed summary. The two
  # paths re-run the step-3 correction independently, so their linear algebra
  # is not guaranteed bit-identical across BLAS/LAPACK implementations; 1e-6
  # is tight enough to catch a real discrepancy while tolerating that noise.
  expect_equal(coef(fitA), coef(fitB), tolerance = 1e-6)

  # The summary is checked on the numbers it is built from, rather than on its
  # rendered text. Comparing the printed output byte for byte contradicted the
  # tolerance granted one line earlier: a difference far below 1e-6 still lands
  # on a rounding boundary every so often and flips a printed digit (0.614
  # against 0.613), which says nothing about whether the two paths agree.
  #
  # The tolerance is looser here than for coef() because the summary carries the
  # standard errors, and those come from inverting the step-3 information: the
  # inversion amplifies the BLAS-level difference between the two paths by
  # roughly an order of magnitude (observed ~5e-6 relative on macOS/aarch64).
  # 1e-4 still catches any discrepancy that would change a reported conclusion.
  sumA <- suppressMessages(capture.output(rowsA <- summary(fitA)))
  sumB <- suppressMessages(capture.output(rowsB <- summary(fitB)))
  expect_equal(rowsA, rowsB, tolerance = 1e-4)
  # The rendered summaries must still have the same shape and labels; only the
  # rounded digits are allowed to differ.
  expect_identical(length(sumA), length(sumB))
  expect_identical(gsub("[0-9]", "", sumA), gsub("[0-9]", "", sumB))

  # Term grouping for the omnibus Wald test survives the add path.
  expect_identical(fitA$sm$parameters$terms, fitB$sm$parameters$terms)

  # Metadata reflects the stepwise estimation.
  expect_identical(fitB$n_steps, 3L)
  expect_identical(fitB$correction, "ML")
})

test_that("add_outcome reproduces one-call BCH and ML outcome fits", {
  d <- .sim_step3_data()

  fit0 <- suppressMessages(suppressWarnings(
    fit_mixture(d$items, n_classes = 2,
                measurement = "binary",
                n_init = 5, random_state = 9)))

  # Continuous outcome, BCH (the auto default).
  fitC <- suppressMessages(suppressWarnings(
    fit_mixture(d$items, n_classes = 2,
                measurement = "binary",
                outcome = d$bmi,
                n_steps = 3, correction = "BCH", n_init = 5, random_state = 9)))
  fitD <- suppressMessages(add_outcome(fit0, d$bmi))
  expect_identical(fitD$correction, "BCH")
  expect_equal(fitC$sm$parameters$means, fitD$sm$parameters$means,
               tolerance = 1e-10)
  expect_identical(capture.output(suppressMessages(summary(fitC))),
                   capture.output(suppressMessages(summary(fitD))))

  # Categorical outcome, ML (the auto default).
  fitE <- suppressMessages(suppressWarnings(
    fit_mixture(d$items, n_classes = 2,
                measurement = "binary",
                outcome = d$grp,
                n_steps = 3, correction = "ML", n_init = 5, random_state = 9)))
  fitF <- suppressMessages(suppressWarnings(add_outcome(fit0, d$grp)))
  expect_identical(fitF$correction, "ML")
  expect_identical(capture.output(suppressMessages(summary(fitE))),
                   capture.output(suppressMessages(summary(fitF))))
})

test_that("a conditional fit can be reused: the structural model is replaced", {
  d <- .sim_step3_data()

  fit0  <- suppressMessages(suppressWarnings(
    fit_mixture(d$items, n_classes = 2,
                measurement = "binary",
                n_init = 5, random_state = 9)))
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
    fit_mixture(items, n_classes = 2,
                measurement = "binary",
                n_init = 3, random_state = 5)))
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
    fit_mixture(d$items, n_classes = 2,
                measurement = "binary",
                predictors = d$covs,
                n_steps = 1, n_init = 3, random_state = 5)))
  expect_error(add_outcome(fit1, d$bmi), "one step")

  # Group-as-predictor fits must go through fit_mixture(group=, predictors=).
  fitg <- suppressMessages(suppressWarnings(
    fit_mixture(d$items, n_classes = 2,
                measurement = "binary",
                group = d$grp,
                group_effects = "prevalence", n_steps = 3, n_init = 3,
                random_state = 5)))
  expect_error(add_covariates(fitg, d$covs), "group")

  # Missing structural data.
  fit0 <- suppressMessages(suppressWarnings(
    fit_mixture(d$items, n_classes = 2,
                measurement = "binary",
                n_init = 3, random_state = 5)))
  expect_error(add_covariates(fit0), "`predictors` is required")
  expect_error(add_outcome(fit0), "`outcome` is required")
})

test_that("the categorical Bayes constant bounds a distal logit on a separated class", {
  # A class in which nobody gives a response has no finite unpenalised logit for
  # it. The prior is what keeps the estimate on the scale; without it the
  # intercept simply runs off towards minus infinity.
  set.seed(3)
  K    <- 3L
  cls  <- rep(seq_len(K), length.out = 300)
  resp <- matrix(0, length(cls), K)
  resp[cbind(seq_along(cls), cls)] <- 1
  y    <- ifelse(cls == 1L, 1L, sample(1:2, length(cls), TRUE))
  Xd   <- matrix(y, ncol = 1)

  fit_alpha <- function(a) {
    st <- distal_categorical_model(K, max_iter = 5000, tol = 1e-10)
    st$bayes_constants <- list(latent = 1, categorical = a, poisson = 1,
                               variances = 1)
    st <- init_params(st, Xd, resp, random_state = 1)
    as.vector(m_step(st, Xd, resp)$parameters$beta_pooled)
  }

  b <- fit_alpha(1)
  expect_true(all(is.finite(b)))
  expect_lt(max(abs(b)), 25)

  # The constant reaches this engine in the same pseudo-observations form the
  # categorical measurement M-step uses, so it must reproduce that closed form:
  #   p_km = (n_km + (alpha/K) * q_m) / (n_k + alpha/K)
  prior_obs <- 1 / K
  closed    <- (as.vector(table(cls, y)[, 2]) + prior_obs * mean(y == 2)) /
    (as.vector(table(cls)) + prior_obs)
  expect_equal(exp(b) / (1 + exp(b)), closed, tolerance = 1e-8)

  # Switching it off recovers plain maximum likelihood, which on this data has
  # no finite maximiser for the separated class.
  expect_gt(abs(fit_alpha(0)[1]), abs(b[1]))
})

test_that("outcome takes one column, and says so when given more", {
  d <- .sim_step3_data(n = 150)
  two <- data.frame(bmi = d$bmi, other = rnorm(150))

  # The old failure was a coercion error naming doubles or xtfrm, which said
  # nothing about the cause.
  expect_error(
    suppressMessages(fit_mixture(d$items, n_classes = 2,
                                 measurement = "binary",
                                 outcome = two,
                                 n_init = 3, random_state = 5)),
    "single distal outcome")
  expect_error(
    suppressMessages(fit_mixture(d$items, n_classes = 2,
                                 measurement = "binary",
                                 outcome = as.matrix(two),
                                 n_init = 3, random_state = 5)),
    "single distal outcome")

  # A single column wrapped in a data frame is unambiguous: accepted, with the
  # column name kept as the outcome label.
  spec <- .build_outcome_spec(data.frame(bmi = d$bmi), NULL, "categorical",
                              "pooled", NULL)
  expect_identical(colnames(spec$Y), "bmi")
  expect_identical(ncol(spec$Y), 1L)
})

# The `data =` form is a way of naming columns and nothing else: it must land on
# the same numbers as extracting the columns by hand.

test_that("the formula form and the hand-extracted form agree exactly", {
  d  <- .sim_step3_data()
  df <- data.frame(d$covs, bmi = d$bmi)

  fit0 <- suppressMessages(suppressWarnings(
    fit_mixture(d$items, n_classes = 2,
                measurement = "binary",
                n_init = 5, random_state = 9)))

  covA <- suppressMessages(add_covariates(fit0, df[, c("age", "sexo")]))
  covB <- suppressMessages(add_covariates(fit0, ~ age + sexo, data = df))
  covC <- suppressMessages(add_covariates(fit0, c("age", "sexo"), data = df))
  expect_equal(coef(covA), coef(covB), tolerance = 1e-12)
  expect_equal(coef(covA), coef(covC), tolerance = 1e-12)
  expect_equal(covA$metrics$ll, covB$metrics$ll, tolerance = 1e-12)

  outA <- suppressMessages(add_outcome(fit0, df$bmi))
  outB <- suppressMessages(add_outcome(fit0, ~ bmi, data = df))
  expect_equal(outA$sm$parameters$means, outB$sm$parameters$means,
               tolerance = 1e-12)
  expect_equal(outA$metrics$ll, outB$metrics$ll, tolerance = 1e-12)
})

test_that("a column missing from `data` is named in the error", {
  d  <- .sim_step3_data()
  df <- data.frame(d$covs, bmi = d$bmi)
  fit0 <- suppressMessages(suppressWarnings(
    fit_mixture(d$items, n_classes = 2,
                measurement = "binary",
                n_init = 5, random_state = 9)))

  expect_error(add_covariates(fit0, ~ age + income, data = df), "income")
  expect_error(add_covariates(fit0, c("age", "income"), data = df), "income")
  expect_error(add_outcome(fit0, ~ weight, data = df), "weight")
  expect_error(add_outcome(fit0, ~ bmi), "`data` must be supplied")
})

test_that("add_outcome's formula must name exactly one outcome", {
  d  <- .sim_step3_data()
  df <- data.frame(d$covs, bmi = d$bmi)
  fit0 <- suppressMessages(suppressWarnings(
    fit_mixture(d$items, n_classes = 2,
                measurement = "binary",
                n_init = 5, random_state = 9)))

  expect_error(add_outcome(fit0, ~ bmi + age, data = df),
               "exactly one distal outcome")
})
