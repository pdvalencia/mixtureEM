# Binary indicators arrive in whatever coding the source data used. The
# Bernoulli emission is arithmetic on 0/1, so anything else has to be converted
# before it reaches the likelihood -- and until it was, a 1/2-coded item in a
# mixed measurement list reached it unchecked and produced a wrong answer with
# no error at all.

.recode_fixture <- function(n = 300) {
  set.seed(3)
  p <- matrix(c(.9, .85, .8, .15,  .1, .2, .15, .9), 2, 4, byrow = TRUE)
  z <- sample(1:2, n, TRUE)
  X <- matrix(rbinom(n * 4, 1, p[z, ]), n, 4)
  colnames(X) <- paste0("q", seq_len(4))
  X
}

.ll01 <- function(X)
  fit_mixture(X, n_classes = 2, measurement = "binary", n_init = 5,
              random_state = 1)$metrics$ll

test_that("every two-valued coding gives the fit the 0/1 data would have", {
  X   <- .recode_fixture()
  ref <- .ll01(X)

  # 1/2, the coding every major package's example data ships with.
  expect_equal(suppressMessages(.ll01(as.data.frame(X + 1))), ref)
  # An arbitrary numeric pair: the lower value is the 0.
  expect_equal(suppressMessages(.ll01(as.data.frame(ifelse(X == 1, 5, 2)))), ref)
  # Logical, and a two-level factor.
  expect_equal(suppressMessages(.ll01(as.data.frame(X == 1))), ref)
  fac <- as.data.frame(lapply(as.data.frame(X), function(v)
    factor(ifelse(v == 1, "yes", "no"), levels = c("no", "yes"))))
  expect_equal(suppressMessages(.ll01(fac)), ref)
})

test_that("data already in 0/1 is left alone and says nothing", {
  X <- .recode_fixture()
  expect_message(fit <- fit_mixture(X, n_classes = 2, measurement = "binary",
                                    n_init = 3, random_state = 1), NA)
  expect_null(fit$binary_recode)
})

test_that("the recoding is reported so a probability is unambiguous", {
  X   <- .recode_fixture()
  fac <- as.data.frame(lapply(as.data.frame(X), function(v)
    factor(ifelse(v == 1, "yes", "no"), levels = c("no", "yes"))))

  expect_message(
    fit <- fit_mixture(fac, n_classes = 2, measurement = "binary",
                       n_init = 3, random_state = 1),
    "Recoded binary indicators")
  expect_equal(fit$binary_recode$q1$zero, "no")
  expect_equal(fit$binary_recode$q1$one,  "yes")
  # And the level shows up where the probabilities are read.
  out <- capture.output(suppressMessages(measurement_summary(fit)))
  expect_true(any(grepl("Probabilities are of.*q1 = yes", out)))
})

test_that("a binary block inside a mixed spec is recoded too", {
  # This is the case the old {0,1} guard skipped entirely: it was reachable only
  # through `is.character(measurement)`, so a list spec was never checked.
  X <- .recode_fixture()
  D <- as.data.frame(cbind(X[, 1:2] + 1, X[, 3:4] + 1))
  colnames(D) <- paste0("q", seq_len(4))

  expect_message(
    fit <- fit_mixture(D, n_classes = 2,
                       measurement = list(binary = 1:2, binary = 3:4),
                       n_init = 3, random_state = 1),
    "Recoded binary indicators")
  expect_named(fit$binary_recode, paste0("q", 1:4))
})

test_that("three or more values is refused, and points somewhere useful", {
  X  <- .recode_fixture()
  D3 <- as.data.frame(X)
  D3$q1 <- sample(1:3, nrow(X), TRUE)
  expect_error(
    suppressMessages(fit_mixture(D3, n_classes = 2, measurement = "binary",
                                 n_init = 2)),
    "categorical")
})

test_that("a 0-based categorical item is refused rather than silently wrong", {
  # Category 0 of item j indexed the last category of item j-1, so a 0-based
  # coding used to fit without complaint and return the wrong likelihood.
  set.seed(8)
  X <- matrix(sample(0:2, 300 * 3, TRUE), 300, 3)
  expect_error(fit_mixture(X, n_classes = 2, measurement = "categorical",
                           n_init = 2),
               "1-based")
  expect_no_error(suppressWarnings(
    fit_mixture(X + 1, n_classes = 2, measurement = "categorical", n_init = 2)))
})
