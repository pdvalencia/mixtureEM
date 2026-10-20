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
               "binary and ordinal indicators")
  d2 <- data.frame(z = rnorm(100))
  expect_error(fit_mixture(s$X, n_classes = 2, measurement = "binary",
                           data = d2, predictors_items = list(item1 = ~ z)),
               "distinct covariate patterns")
})

# Ordinal items (R/ordinal_dif.R).
.make_ordinal_dif_data <- function(seed = 21, n = 1200) {
  set.seed(seed)
  z  <- rbinom(n, 1, 0.5)
  cl <- 1L + rbinom(n, 1, plogis(-0.3 + 0.6 * z))
  # P(Y >= s), s = 2, 3, on the logit scale, per class.
  th <- list(c(-0.5, -2.0), c(1.8, 0.2))
  X  <- matrix(NA_real_, n, 4)
  for (i in seq_len(n)) for (j in 1:4) {
    sh  <- if (j == 1) 0.8 * z[i] else 0
    cum <- plogis(th[[cl[i]]] + sh + 0.2 * (j - 2))
    X[i, j] <- 1 + sum(runif(1) < cum)
  }
  colnames(X) <- paste0("item", 1:4)
  list(X = X, d = data.frame(z = z))
}

test_that("ordinal DIF at zero slopes is the plain ordinal emission", {
  s  <- .make_ordinal_dif_data(n = 300)
  X  <- s$X
  X[5, 3] <- NA
  mm <- init_params(build_emission("ordinal", n_components = 2, cats = rep(3L, 4)),
                    X, NULL, random_state = 1)
  spec <- .dif_spec(list(item1 = ~ z, item2 = ~ z), "item2", s$d, X, "ordinal")
  md <- .dif_emission(mm, spec)
  expect_s3_class(md, "ordinal_dif")
  md$parameters$dif_slopes <- matrix(0, 2, 4)
  expect_equal(log_likelihood(md, X), log_likelihood(mm, X), tolerance = 1e-8)
  # 2 x 4 x 2 thresholds, one uniform slope, two class-specific ones.
  expect_identical(as.integer(n_parameters(md)), 19L)
  # The packing round-trips.
  md$parameters$dif_slopes[, 1] <- 0.4
  expect_equal(.ordinal_dif_unpack(md, .ordinal_dif_pack(md))$parameters,
               md$parameters, tolerance = 1e-8)
})

test_that("a two-category ordinal item with DIF is the binary DIF model", {
  s   <- .make_dif_data()
  off <- list(categorical = 0, latent = 0)
  args <- list(n_classes = 2, predictors = ~ z, data = s$d,
               predictors_items = list(item1 = ~ z, item3 = ~ z),
               n_init = 5, random_state = 1, n_cores = 1, bayes_constants = off)
  fb <- do.call(fit_mixture, c(list(s$X, measurement = "binary"), args))
  fo <- do.call(fit_mixture, c(list(s$X + 1, measurement = "ordinal",
                                    cats = rep(2L, 4)), args))
  expect_lt(abs(fb$metrics$ll - fo$metrics$ll), 1e-4)
  expect_identical(fb$metrics$n_params, fo$metrics$n_params)
  expect_equal(dif_effects(fo, se = FALSE)$estimate,
               dif_effects(fb, se = FALSE)$estimate, tolerance = 1e-3)
})

test_that("an ordinal uniform slope is recovered and reported", {
  s   <- .make_ordinal_dif_data()
  fit <- fit_mixture(s$X, n_classes = 2, measurement = "ordinal",
                     cats = rep(3L, 4), predictors = ~ z, data = s$d,
                     predictors_items = list(item1 = ~ z),
                     n_init = 5, random_state = 1, n_cores = 1)
  de <- dif_effects(fit)
  expect_identical(de$class, "all")
  expect_lt(abs(de$estimate - 0.8), 3 * de$se)
  expect_true(de$se > 0 && de$se < 1)
  expect_output(print(fit), "Item Covariates")
})
