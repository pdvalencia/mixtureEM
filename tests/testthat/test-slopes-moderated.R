# `slopes` naming a subset of covariates: class moderates some covariates
# (class-specific slopes) while others are only adjusted for (pooled slopes).
#
# The core algebraic claim is that the joint design .distal_U() builds is
# exactly the expanded weighted-least-squares problem it is meant to be, so
# the primary check is against lm.wfit() on that same design -- not against
# any reference program's output.

.sim_moderation_data <- function(n = 400, seed = 11) {
  set.seed(seed)
  cl    <- sample(1:3, n, TRUE, c(0.5, 0.3, 0.2))
  items <- sapply(1:6, function(j) {
    p <- c(0.85, 0.5, 0.15)[cl]
    rbinom(n, 1, p)
  })
  colnames(items) <- paste0("item", 1:6)

  loc1 <- rbinom(n, 1, 0.3)
  loc2 <- rbinom(n, 1, 0.2)
  age  <- rnorm(n, 30, 5)
  sex  <- rbinom(n, 1, 0.5)

  mod_effect <- c(0, 8, -5)[cl] * loc1 + c(0, -3, 6)[cl] * loc2
  y <- 10 + 2 * (cl == 2) + 4 * (cl == 3) + mod_effect +
    0.3 * age + 1.5 * sex + rnorm(n, sd = 2)

  list(items = items, loc1 = loc1, loc2 = loc2, age = age, sex = sex, y = y)
}

test_that("moderated slopes match lm.wfit() on the same expanded design", {
  d <- .sim_moderation_data()
  covs <- data.frame(loc1 = d$loc1, loc2 = d$loc2, age = d$age, sex = d$sex)

  fit0 <- suppressMessages(suppressWarnings(
    fit_mixture(d$items, n_classes = 3, measurement = "binary",
               n_init = 5, random_state = 3)))

  # correction = "none" uses the frozen measurement-only posteriors directly
  # as weights, with no BCH classification-error adjustment, so the weights
  # used by m_step() are exactly exp(fit0$log_resp) and can be reproduced
  # without re-deriving the correction.
  fitM <- suppressMessages(add_outcome(
    fit0, d$y, covariates = covs, slopes = c("loc1", "loc2"),
    correction = "none"))

  expect_s3_class(fitM$sm, "distal_continuous_pooled")
  expect_identical(fitM$sm$moderated, c(1L, 2L))

  K   <- fit0$n_components
  Z   <- as.matrix(covs)
  N   <- nrow(Z)
  U   <- .distal_U(Z, K, N, mod = c(1L, 2L))
  W   <- as.vector(exp(fit0$log_resp))
  Yf  <- rep(d$y, K)

  ref <- stats::lm.wfit(x = U, y = Yf, w = W)$coefficients
  got <- as.vector(fitM$sm$parameters$beta_pooled)

  expect_equal(got, unname(ref), tolerance = 1e-6)
})

test_that("Npar and the design layout are L = K + D_pool + K*D_mod", {
  d <- .sim_moderation_data()
  covs <- data.frame(loc1 = d$loc1, loc2 = d$loc2, age = d$age, sex = d$sex)

  fit0 <- suppressMessages(suppressWarnings(
    fit_mixture(d$items, n_classes = 3, measurement = "binary",
               n_init = 5, random_state = 3)))
  fitM <- suppressMessages(add_outcome(
    fit0, d$y, covariates = covs, slopes = c("loc1", "loc2"),
    correction = "none"))

  K <- fit0$n_components  # 3
  # D_pool = 2 (age, sex), D_mod = 2 (loc1, loc2)
  expect_equal(ncol(fitM$sm$parameters$beta_pooled), K + 2L + K * 2L)
  expect_equal(n_parameters(fitM$sm), K + 2L + K * 2L + 1L)

  nm <- colnames(fitM$sm$parameters$beta_pooled)
  expect_identical(nm[1:K], paste0("Class_", 1:K))
  expect_true(all(c("age", "sex") %in% nm))
  expect_true(all(paste0("loc1:Class", 1:K) %in% nm))
  expect_true(all(paste0("loc2:Class", 1:K) %in% nm))
})

test_that("naming every covariate reproduces slopes = 'class_specific'", {
  d <- .sim_moderation_data()
  covs <- data.frame(loc1 = d$loc1, loc2 = d$loc2, age = d$age, sex = d$sex)

  fit0 <- suppressMessages(suppressWarnings(
    fit_mixture(d$items, n_classes = 3, measurement = "binary",
               n_init = 5, random_state = 3)))

  fit_cs <- suppressMessages(add_outcome(
    fit0, d$y, covariates = covs, slopes = "class_specific",
    correction = "none"))
  fit_all_mod <- suppressMessages(add_outcome(
    fit0, d$y, covariates = covs,
    slopes = c("loc1", "loc2", "age", "sex"), correction = "none"))

  K <- fit0$n_components
  D <- 4L

  # fit_cs$sm$parameters$betas is K x (1 + D): intercept, then the 4 slopes.
  # fit_all_mod's beta_pooled is [intercepts, moderated blocks in covariate
  # order] since D_pool = 0. Reassemble class k's [intercept, slopes] vector
  # from each and compare.
  betas_cs <- fit_cs$sm$parameters$betas
  theta    <- as.vector(fit_all_mod$sm$parameters$beta_pooled)
  for (k in seq_len(K)) {
    est_mod <- c(theta[k],
                vapply(seq_len(D), function(j) theta[K + (j - 1L) * K + k],
                       numeric(1)))
    expect_equal(est_mod, unname(betas_cs[k, ]), tolerance = 1e-6)
  }
})

test_that("a subset `slopes` on a categorical outcome errors clearly", {
  d <- .sim_moderation_data()
  grp <- factor(sample(c("a", "b"), length(d$y), TRUE))
  covs <- data.frame(loc1 = d$loc1, age = d$age)

  fit0 <- suppressMessages(suppressWarnings(
    fit_mixture(d$items, n_classes = 3, measurement = "binary",
               n_init = 3, random_state = 3)))

  expect_error(
    suppressMessages(add_outcome(fit0, grp, covariates = covs,
                                 slopes = c("loc1"))),
    '"pooled"')
})

test_that("an unresolved `slopes` term is named in the error", {
  d <- .sim_moderation_data()
  covs <- data.frame(loc1 = d$loc1, age = d$age)

  fit0 <- suppressMessages(suppressWarnings(
    fit_mixture(d$items, n_classes = 3, measurement = "binary",
               n_init = 3, random_state = 3)))

  expect_error(
    suppressMessages(add_outcome(fit0, d$y, covariates = covs,
                                 slopes = c("loc1", "not_a_covariate"))),
    "not_a_covariate")
})

test_that("a one-sided formula names the same terms as a character vector", {
  d <- .sim_moderation_data()
  covs <- data.frame(loc1 = d$loc1, loc2 = d$loc2, age = d$age, sex = d$sex)

  fit0 <- suppressMessages(suppressWarnings(
    fit_mixture(d$items, n_classes = 3, measurement = "binary",
               n_init = 3, random_state = 3)))

  fitA <- suppressMessages(add_outcome(
    fit0, d$y, covariates = covs, slopes = c("loc1", "loc2"),
    correction = "none"))
  fitB <- suppressMessages(add_outcome(
    fit0, d$y, covariates = covs, slopes = ~ loc1 + loc2,
    correction = "none"))

  expect_identical(fitA$sm$moderated, fitB$sm$moderated)
  expect_equal(fitA$sm$parameters$beta_pooled, fitB$sm$parameters$beta_pooled,
               tolerance = 1e-10)
})
