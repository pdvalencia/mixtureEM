# Direct covariate effects on binary indicators of a cross-sectional model
# (`predictors_items` on fit_mixture(), R/categorical_dif.R). The external grade
# lives with the validation tests; these are the two checks that need no
# outside number, plus the arithmetic and the refusals.

.make_dif_data <- function(seed = 11, n = 600) {
  set.seed(seed)
  J <- 4L
  z  <- rbinom(n, 1, 0.5)
  cl <- 1L + rbinom(n, 1, plogis(-0.4 + 0.8 * z))
  th <- rbind(c(-1.5, -1.2, -1.8, -1.0), c(1.4, 1.1, 1.6, 1.2))
  X  <- matrix(NA_real_, n, J)
  for (i in seq_len(n)) {
    eta <- th[cl[i], ] + c(0.9 * z[i], 0, -0.7 * z[i], 0)
    X[i, ] <- rbinom(J, 1, plogis(eta))
  }
  colnames(X) <- paste0("item", 1:J)
  list(X = X, d = data.frame(z = z))
}

test_that("with every slope at zero the DIF emission is the plain one", {
  s  <- .make_dif_data()
  X  <- s$X
  X[3, 2] <- NA
  mm <- init_params(build_emission("bernoulli_nan", n_components = 2), X, NULL,
                    random_state = 1)
  spec <- .dif_spec(list(item1 = ~ z, item3 = ~ z), "item3", s$d, X, "binary")
  md <- .dif_emission(mm, spec)
  md$parameters$dif_slopes <- matrix(0, 2, 4)
  expect_equal(log_likelihood(md, X), log_likelihood(mm, X), tolerance = 1e-12)
  # 2 x 4 thresholds, one uniform slope, two class-specific ones.
  expect_identical(as.integer(n_parameters(md)), 11L)
})

test_that("a binary covariate on every item, by class, is the two-group model", {
  s   <- .make_dif_data()
  off <- list(categorical = 0, latent = 0)
  all_items <- setNames(rep(list(~ z), 4), colnames(s$X))
  fit_dif <- fit_mixture(s$X, n_classes = 2, measurement = "binary",
                         predictors = ~ z, data = s$d,
                         predictors_items = all_items,
                         predictors_items_by_class = colnames(s$X),
                         n_init = 5, random_state = 1, n_cores = 1,
                         bayes_constants = off)
  fit_grp <- suppressMessages(fit_mixture(
    s$X, n_classes = 2, measurement = "binary", group = factor(s$d$z),
    group_effects = "both", n_steps = 1, n_init = 5, random_state = 1,
    n_cores = 1, bayes_constants = off))
  expect_lt(abs(fit_dif$metrics$ll - fit_grp$metrics$ll), 1e-4)
  expect_identical(fit_dif$metrics$n_params, fit_grp$metrics$n_params)
})

test_that("a uniform slope is recovered and reported", {
  s   <- .make_dif_data(n = 1500)
  fit <- fit_mixture(s$X, n_classes = 2, measurement = "binary",
                     predictors = ~ z, data = s$d,
                     predictors_items = list(item1 = ~ z, item3 = ~ z),
                     n_init = 5, random_state = 1, n_cores = 1)
  de <- dif_effects(fit)
  expect_identical(de$item, c("item1", "item3"))
  expect_identical(de$class, c("all", "all"))
  expect_lt(abs(de$estimate[1] - 0.9), 3 * de$se[1])
  expect_lt(abs(de$estimate[2] + 0.7), 3 * de$se[2])
  expect_true(all(de$se > 0 & de$se < 1))
})

test_that("predictors_items refuses what it cannot fit", {
  s <- .make_dif_data(n = 100)
  expect_error(fit_mixture(s$X, n_classes = 2, measurement = "binary",
                           data = s$d, predictors_items = list(nope = ~ z)),
               "not among the indicators")
  expect_error(suppressMessages(fit_mixture(
    s$X, n_classes = 2, measurement = "binary", data = s$d,
    predictors_items = list(item1 = ~ z), predictors = ~ z, n_steps = 3)),
    "n_steps = 1")
  expect_error(fit_mixture(s$X, n_classes = 2, measurement = "continuous",
                           data = s$d, predictors_items = list(item1 = ~ z)),
               "binary indicators only")
  d2 <- data.frame(z = rnorm(100))
  expect_error(fit_mixture(s$X, n_classes = 2, measurement = "binary",
                           data = d2, predictors_items = list(item1 = ~ z)),
               "distinct covariate patterns")
})
