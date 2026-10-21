# skew_kurtosis_test(): sample skewness and kurtosis of a fitted mixture's
# outcomes against their distribution in data simulated from the fit.
#
# The sample statistics are checked against identities they must satisfy
# (Mardia's measures of one outcome are the square of its skewness and its
# kurtosis, both on n; every measure is invariant to an affine change of
# scale). The reference distribution is checked on a one-class fit, where it
# is known: a normal sample's skewness has mean 0 and standard deviation about
# sqrt(6 / n), its excess kurtosis sd about sqrt(24 / n), and Mardia's
# skewness mean about p(p + 1)(p + 2) / n.

.sk_data <- function(n, seed = 1) {
  set.seed(seed)
  cls <- rbinom(n, 1, 0.2)
  outer(ifelse(cls == 1, 3, 1), rep(1, 4)) +
    outer(ifelse(cls == 1, 0.8, 0.1), 0:3) + matrix(rnorm(n * 4, 0, 0.7), n, 4)
}
.sk_fit <- function(y, K)
  suppressWarnings(suppressMessages(
    fit_gmm(y, n_classes = K, times = 4, time_scores = 0:3,
            random_effects = "none", n_init = 3, random_state = 1,
            n_cores = 1)))

test_that("the sample statistics satisfy their defining identities", {
  y <- .sk_data(150)
  colnames(y) <- paste0("y", 1:4)
  s <- .sk_stats(y)
  x <- y[, 1]; d <- x - mean(x); n <- length(x)
  m2 <- mean(d^2)
  b1 <- s$value[s$level == "bivariate" & s$statistic == "skewness"]
  expect_equal(.sk_stats(cbind(a = x, b = x + rnorm(n)))$value[3],
               mean(d^3) / (sum(d^2) / (n - 1))^1.5)
  # One outcome: b_{1,1} = (m3 / m2^1.5)^2 and b_{2,1} = m4 / m2^2.
  one <- function(z) {
    dz <- z - mean(z); S <- mean(dz^2); D <- outer(dz, dz) / S
    c(mean(D^3), mean(diag(D)^2))
  }
  expect_equal(one(x), c((mean(d^3) / m2^1.5)^2, mean(d^4) / m2^2))
  # Affine invariance of every statistic.
  z <- sweep(y %*% diag(c(2, 1, 0.5, 3)), 2, c(5, -1, 0, 2), "+")
  colnames(z) <- colnames(y)
  expect_equal(.sk_stats(z)$value, s$value)
  expect_equal(nrow(s), 2 + 8 + 12)
  expect_true(all(b1 >= 0))
})

test_that("a one-class fit's reference distribution is the normal one", {
  y   <- .sk_data(400, seed = 2)
  fit <- .sk_fit(y, 1)
  r   <- skew_kurtosis_test(fit, nsim = 200, seed = 1)
  u   <- r$univariate
  sk  <- u$statistic == "skewness"
  expect_true(all(abs(u$mean[sk]) < 4 * sqrt(6 / 400) / sqrt(200)))
  expect_true(all(abs(u$sd[sk] / sqrt(6 / 400) - 1) < 0.2))
  expect_true(all(abs(u$sd[!sk] / sqrt(24 / 400) - 1) < 0.3))
  expect_equal(r$multivariate$mean[1], 4 * 5 * 6 / 400, tolerance = 0.15)
  # The sample is a two-class mixture, so its shape is rejected.
  expect_lt(r$multivariate$p_value[1], 0.05)
})

test_that("the test keeps the seed contract and refuses what it cannot test", {
  y   <- .sk_data(120, seed = 3)
  fit <- .sk_fit(y, 2)
  set.seed(9)
  before <- get(".Random.seed", envir = globalenv())
  a <- skew_kurtosis_test(fit, nsim = 5, seed = 4)
  expect_identical(get(".Random.seed", envir = globalenv()), before)
  expect_identical(a, skew_kurtosis_test(fit, nsim = 5, seed = 4))
  expect_s3_class(a, "skew_kurtosis_test")
  expect_output(print(a), "Multivariate")
  expect_error(skew_kurtosis_test(fit, nsim = 1), "at least 2")
  X <- matrix(rbinom(300, 1, 0.5), ncol = 3)
  bin <- suppressWarnings(fit_mixture(X, n_classes = 2, measurement = "binary",
                                      n_init = 2, random_state = 1))
  expect_error(skew_kurtosis_test(bin), "normal outcomes")
})
