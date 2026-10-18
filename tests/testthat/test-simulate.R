# simulate(): data sets drawn from a fitted model's measurement model.
#
# The checks are the stats::simulate() seed contract, and that a large draw
# reproduces the model-implied moments: the class proportions, and each
# indicator's marginal mean sum_k pi_k mu_jk. Tolerances are about four
# Monte Carlo standard errors at n = 20000.

.s_fit <- function(measurement, X, ...)
  suppressWarnings(suppressMessages(
    fit_mixture(X, n_classes = 2, measurement = measurement, n_init = 3,
                random_state = 1, n_cores = 1, ...)))

test_that("a seed reproduces the draw and leaves the caller's stream alone", {
  set.seed(1)
  X   <- matrix(rbinom(600, 1, rep(c(.2, .8), each = 300)), ncol = 6,
                dimnames = list(NULL, paste0("u", 1:6)))
  fit <- .s_fit("binary", X)
  set.seed(99)
  before <- .Random.seed
  a <- simulate(fit, nsim = 2, seed = 7, n = 50)
  expect_identical(.Random.seed, before)
  b <- simulate(fit, nsim = 2, seed = 7, n = 50)
  expect_identical(a, b)
  expect_identical(attr(a, "seed")[[1]], 7)
  expect_length(a, 2L)
  expect_named(a[[1]], paste0("u", 1:6))
  expect_identical(nrow(a[[1]]), 50L)
  expect_length(attr(a[[1]], "classes"), 50L)
  expect_identical(nrow(simulate(fit)[[1]]), nrow(X))
  # Without a seed the draw continues the caller's stream.
  set.seed(5); c1 <- simulate(fit, n = 10)
  set.seed(5); c2 <- simulate(fit, n = 10)
  expect_identical(c1[[1]], c2[[1]])
})

test_that("a large draw reproduces the model-implied moments", {
  set.seed(2)
  n <- 400
  z <- sample(2, n, TRUE, c(.6, .4))
  d <- data.frame(b1 = rbinom(n, 1, c(.2, .8)[z]),
                  p1 = sapply(z, function(k) sample(1:3, 1, prob = if (k == 1) c(.6, .3, .1) else c(.1, .3, .6))),
                  c1 = rnorm(n, c(0, 2)[z]),
                  k1 = rpois(n, c(1, 4)[z]))
  fit <- .s_fit(list(binary = 1, categorical = 2, continuous = 3, count = 4), d)
  s   <- simulate(fit, seed = 11, n = 20000)[[1]]
  cl  <- attr(s, "classes")
  expect_lt(max(abs(tabulate(cl, 2) / 20000 - fit$weights)), 0.015)

  w  <- fit$weights
  ms <- fit$mm$models
  pick <- function(type) ms[[which(vapply(ms, function(m) inherits(m, type),
                                          logical(1)))[1]]]$parameters
  expect_lt(abs(mean(s$b1) - sum(w * pick("bernoulli")$pis[, 1])), 0.015)
  expect_lt(abs(mean(s$c1) - sum(w * pick("gaussian_diag")$means[, 1])), 0.04)
  expect_lt(abs(mean(s$k1) - sum(w * pick("poisson")$rates[, 1])), 0.06)
  pp <- pick("multinoulli")$pis
  expect_lt(max(abs(tabulate(s$p1, 3) / 20000 - colSums(w * pp[, 1:3]))),
            0.015)
  expect_true(all(s$b1 %in% 0:1) && all(s$p1 %in% 1:3))
})

test_that("a model with predictors draws classes given the observed predictors", {
  set.seed(3)
  n <- 300
  x <- rnorm(n)
  z <- ifelse(runif(n) < plogis(1.5 * x), 1L, 2L)
  X <- data.frame(y1 = rnorm(n, c(0, 2)[z]), y2 = rnorm(n, c(0, 2)[z]))
  fit <- .s_fit("continuous", X, predictors = data.frame(x = x), n_steps = 1)
  expect_error(simulate(fit, n = 10), "one case per original case")
  s <- simulate(fit, nsim = 50, seed = 4)
  # Over many draws a case's class frequency is its fitted membership
  # probability, so the draws follow the predictor.
  freq1 <- rowMeans(sapply(s, function(d) attr(d, "classes") == 1L))
  p1    <- exp(mixtureEM:::log_likelihood(fit$sm, fit$Y))[, 1]
  expect_gt(cor(freq1, p1), 0.9)
})

test_that("growth models simulate their fitted trajectories", {
  set.seed(1)
  n   <- 400
  cls <- rbinom(n, 1, 0.5)
  eta <- outer(ifelse(cls == 1, 0.0, -0.8), rep(1, 4)) +
    outer(ifelse(cls == 1, 1.1, -0.1), 0:3)
  y   <- matrix(rbinom(n * 4, 1, plogis(eta)), n, 4)
  lc  <- suppressWarnings(suppressMessages(
    fit_lcga(y, times = 4, n_classes = 2, n_init = 3, random_state = 1,
             n_cores = 1)))
  s   <- simulate(lc, seed = 1, n = 20000)[[1]]
  implied <- colSums(lc$weights * mixtureEM:::.lcga_fitted(lc$mm))
  expect_identical(dim(s), c(20000L, 4L))
  expect_lt(max(abs(colMeans(s) - implied)), 0.015)

  icept <- ifelse(cls == 1, 3, 1) + rnorm(n, 0, 1)
  slope <- ifelse(cls == 1, 1.0, 0.2) + rnorm(n, 0, 0.4)
  g <- outer(icept, rep(1, 4)) + outer(slope, 0:3) +
    matrix(rnorm(n * 4, 0, 0.7), n, 4)
  gm <- suppressWarnings(suppressMessages(
    fit_gmm(g, times = 4, n_classes = 2, n_init = 3, random_state = 1,
            n_cores = 1)))
  s  <- simulate(gm, seed = 1, n = 20000)[[1]]
  implied <- colSums(gm$weights * mixtureEM:::.sn_fitted(gm$mm))
  expect_lt(max(abs(colMeans(s) - implied)), 0.06)
  # Within-class occasions are correlated through the random effects.
  expect_gt(cor(s[[1]], s[[4]]), 0.3)
})

test_that("models simulate() does not cover are refused", {
  set.seed(22)
  X <- matrix(rbinom(150 * 3, 1, 0.5), ncol = 3)
  lta <- suppressWarnings(
    fit_lta(X, n_statuses = 2, times = 3, measurement = "binary", n_init = 2,
            random_state = 2, standard_errors = FALSE))
  expect_error(simulate(lta), "latent transition")
})
