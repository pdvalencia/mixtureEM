# Poisson ("count") measurement model: each item is Poisson within class with
# its own rate. These tests check the emission itself (likelihood, M-step,
# FIML variant) and its integration with the machinery that has to know about
# every parameter matrix: class sorting, mixed measurement models, BLRT data
# generation, and the profile plot.

make_count_data <- function(n = 1500, seed = 7) {
  set.seed(seed)
  rates <- rbind(c(0.5, 0.8, 1.0, 6.0, 5.0),
                 c(6.0, 5.5, 5.0, 0.6, 0.9),
                 c(2.5, 2.0, 3.0, 2.5, 2.0))
  w <- c(.45, .35, .20)
  z <- sample(1:3, n, TRUE, w)
  X <- t(vapply(seq_len(n), function(i) rpois(5, rates[z[i], ]), numeric(5)))
  colnames(X) <- paste0("ev", 1:5)
  list(X = X, z = z, rates = rates, weights = w)
}

# Match fitted classes to generating classes by nearest rate profile, so the
# comparison does not depend on the label order the fit happens to produce.
match_classes <- function(est, truth)
  apply(est, 1, function(r) which.min(colSums((t(truth) - r)^2)))

test_that("Poisson LCA recovers known rates and class sizes", {
  d <- make_count_data()
  fit <- fit_mixture(d$X, n_classes = 3, measurement = "count",
                     n_init = 10, random_state = 11)

  expect_true(fit$converged)
  expect_s3_class(fit$mm, "poisson")
  # 3 classes x 5 rates + 2 free class weights
  expect_equal(fit$metrics$n_params, 17)

  est  <- fit$mm$parameters$rates
  perm <- match_classes(est, d$rates)
  expect_setequal(perm, 1:3)                      # one fitted class per true class
  expect_equal(unname(est), d$rates[perm, ], tolerance = 0.1)
  expect_equal(fit$weights, d$weights[perm], tolerance = 0.03)
})

test_that("the reported log-likelihood is the Poisson mixture likelihood", {
  d <- make_count_data(n = 600, seed = 3)
  fit <- fit_mixture(d$X, n_classes = 2, measurement = "count",
                     n_init = 5, random_state = 2)

  rates <- fit$mm$parameters$rates
  n     <- nrow(d$X)
  dens  <- vapply(seq_len(2), function(k)
    fit$weights[k] * exp(rowSums(dpois(
      d$X, matrix(rates[k, ], n, ncol(d$X), byrow = TRUE), log = TRUE))),
    numeric(n))

  expect_equal(fit$metrics$ll, sum(log(rowSums(dens))))
})

test_that("the FIML variant handles missing counts", {
  d <- make_count_data(n = 1200, seed = 5)
  X <- d$X
  set.seed(21)
  X[matrix(runif(length(X)) < 0.15, nrow = nrow(X))] <- NA

  fit <- suppressWarnings(
    fit_mixture(X, n_classes = 3, measurement = "count",
                n_init = 10, random_state = 11))

  expect_s3_class(fit$mm, "poisson_nan")
  est <- fit$mm$parameters$rates
  expect_equal(unname(est), d$rates[match_classes(est, d$rates), ],
               tolerance = 0.15)
})

test_that("counts combine with other indicator types in a mixed model", {
  d <- make_count_data(n = 1000, seed = 8)
  set.seed(8)
  binary <- matrix(rbinom(1000 * 3, 1, ifelse(d$z == 1, .85, .15)), ncol = 3)
  X <- cbind(binary, d$X[, 1:2])
  colnames(X) <- c("b1", "b2", "b3", "c1", "c2")

  fit <- fit_mixture(X, n_classes = 2,
                     measurement = list(binary = 1:3, count = 4:5),
                     n_init = 5, random_state = 2)

  expect_true(fit$converged)
  expect_s3_class(fit$mm$models$binary, "bernoulli")
  expect_s3_class(fit$mm$models$count,  "poisson")
  # 2 x 3 probabilities + 2 x 2 rates + 1 free class weight
  expect_equal(fit$metrics$n_params, 11)
})

test_that("class sorting keeps rates aligned with class weights", {
  d <- make_count_data(n = 800, seed = 4)
  fit <- fit_mixture(d$X, n_classes = 3, measurement = "count",
                     n_init = 8, random_state = 6, order_by_size = TRUE)

  expect_true(all(diff(fit$weights) <= 0))
  # Each sorted class's rates must still be the profile its own posterior
  # implies: reconstruct the rates from the responsibilities and compare.
  resp  <- exp(fit$log_resp)
  implied <- (t(resp) %*% d$X) / colSums(resp)
  expect_equal(unname(fit$mm$parameters$rates), unname(implied), tolerance = 0.05)
})

test_that("non-count indicators are rejected with a usable message", {
  d <- make_count_data(n = 100, seed = 2)
  expect_error(
    fit_mixture(replace(d$X, 1, -1), n_classes = 2, measurement = "count"),
    "non-negative integer counts")
  expect_error(
    fit_mixture(replace(d$X, 1, 2.5), n_classes = 2, measurement = "count"),
    "non-negative integer counts")
})

test_that("counts carry through the repeated-measures engine", {
  set.seed(15)
  n <- 700
  z <- rbinom(n, 1, .5)
  rate_of <- function(k, t)
    if (k == 1) c(1, 1.2) * (1 + t) else c(5, 4.5) / (1 + 0.4 * t)
  X <- do.call(cbind, lapply(0:2, function(t)
    t(vapply(seq_len(n), function(i) rpois(2, rate_of(z[i] + 1, t)), numeric(2)))))
  colnames(X) <- paste0(rep(c("ev1", "ev2"), 3), "_t", rep(1:3, each = 2))

  fit <- fit_rmlca(X, n_classes = 2, times = 3, measurement = "count",
                   n_init = 5, random_state = 1)

  expect_true(fit$converged)
  expect_s3_class(fit$mm$models[[1]], "poisson")
  # Trajectories must be reported as rates, not silently dropped to all-NA.
  traj <- .rmlca_trajectories(fit$mm)
  expect_equal(traj$kind, "Event rate")
  expect_false(anyNA(traj$values))
  expect_equal(dim(traj$values), c(2L, 3L, 2L))
})

test_that("a count block inside a mixed model upgrades to FIML when it has NAs", {
  d <- make_count_data(n = 700, seed = 15)
  X <- cbind(matrix(rbinom(700 * 2, 1, ifelse(d$z == 1, .8, .2)), ncol = 2),
             d$X[, 1:2])
  colnames(X) <- c("b1", "b2", "c1", "c2")
  set.seed(4)
  X[cbind(sample(700, 100), sample(3:4, 100, TRUE))] <- NA

  fit <- suppressWarnings(
    fit_mixture(X, n_classes = 2, measurement = list(binary = 1:2, count = 3:4),
                n_init = 5, random_state = 3))

  expect_s3_class(fit$mm$models$count,  "poisson_nan")
  expect_s3_class(fit$mm$models$binary, "bernoulli")
  expect_true(fit$converged)
})

test_that("BLRT and the profile plot support count models", {
  d <- make_count_data(n = 600, seed = 12)

  # n_reps = 10 is a deliberate speed fixture: this test checks that the count
  # family reaches blrt() at all, not the resolution of the p-value, so the
  # granularity note blrt() emits at that setting is expected here.
  test <- suppressWarnings(blrt(d$X, k_small = 1, k_large = 2,
               measurement = "count",
               n_reps = 10, n_init_base = 3, n_init_boot = 2, verbose = FALSE))
  expect_true(test$obs_diff > 0)
  expect_length(test$null_dist, 10L)
  expect_true(test$p_value >= 0 && test$p_value <= 1)

  fit <- fit_mixture(d$X, n_classes = 2, measurement = "count",
                     n_init = 5, random_state = 1)
  profile <- .prepare_profile_data(fit$mm, fit$data)
  expect_equal(dim(profile), c(2L, 5L))
  expect_true(all(profile >= 0 & profile <= 1))
})
