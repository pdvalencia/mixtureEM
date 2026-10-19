# `alongside` fits are shown beside a comparison and never ranked in it: the
# sweep's table, BF, cmP and best K must be exactly what they are without them.

count_data <- function() {
  set.seed(11)
  cl <- rbinom(150, 1, 0.4)
  sapply(1:4, function(j) rpois(150, ifelse(cl == 1, 4, 1)))
}

quiet <- function(expr) suppressWarnings(suppressMessages(expr))

test_that("a fit from another family is shown, labelled and not ranked", {
  X <- count_data()
  normal <- quiet(fit_mixture(X, n_classes = 2, measurement = "continuous",
                              n_init = 2, random_state = 1))
  out0 <- NULL
  utils::capture.output(out0 <- quiet(compare_mixtures(
    X, 1:2, measurement = "count", n_init = 2, random_state = 1)))
  txt <- utils::capture.output(out <- quiet(compare_mixtures(
    X, 1:2, measurement = "count", n_init = 2, random_state = 1,
    alongside = list(normal = normal))))

  expect_equal(out$fit_table, out0$fit_table)
  expect_identical(out$best_k, out0$best_k)
  expect_identical(attr(out, "family"), "Poisson")
  expect_identical(out$alongside$Model, "normal")
  expect_identical(out$alongside$Family, "normal")
  expect_equal(out$alongside$BIC, normal$metrics$bic)
  expect_true(any(grepl("Shown alongside, not ranked", txt)))
  expect_true(any(grepl("Not ranked against the table above", txt)))
  expect_length(grep("Best model", txt), 1L)
})

test_that("a fit from the same family is shown without the family note", {
  X <- count_data()
  other <- quiet(fit_mixture(X, n_classes = 3, measurement = "count",
                             n_init = 2, random_state = 1))
  txt <- utils::capture.output(out <- quiet(compare_mixtures(
    X, 1:2, measurement = "count", n_init = 2, random_state = 1,
    alongside = other)))
  expect_identical(out$alongside$Model, "alongside 1")
  expect_true(any(grepl("Shown alongside, not ranked", txt)))
  expect_false(any(grepl("Not ranked against the table above", txt)))
})

test_that("anything but a fitted model is refused before any fitting", {
  expect_error(compare_mixtures(count_data(), 1:2, measurement = "count",
                                alongside = list(1)),
               "must be a fitted model")
})

test_that("compare_longitudinal() labels a growth fit by its family", {
  set.seed(5)
  n  <- 120
  cl <- rbinom(n, 1, 0.5)
  Y  <- sapply(0:3, function(t) rpois(n, exp(0.3 + ifelse(cl == 1, 0.3, -0.1) * t)))
  gauss <- quiet(fit_lcga(Y, n_classes = 2, times = 4, family = "gaussian",
                          n_init = 2, random_state = 1))
  txt <- utils::capture.output(out <- quiet(compare_longitudinal(
    Y, k_range = 1:2, model = "lcga", times = 4, family = "poisson",
    n_init = 2, random_state = 1, alongside = list(gaussian = gauss))))
  expect_identical(attr(out, "family"), "Poisson")
  expect_identical(out$alongside$Family, "normal")
  expect_true(any(grepl("Not ranked against the table above", txt)))
})
