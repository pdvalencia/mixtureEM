# Summary functions return the printed numbers in tidy form, so users never
# need to reach into the fitted object's internals.

.sim_summary_data <- function(n = 200, seed = 11) {
  set.seed(seed)
  cl    <- rbinom(n, 1, 0.4)
  items <- sapply(1:4, function(j) rbinom(n, 1, ifelse(cl == 1, 0.85, 0.15)))
  colnames(items) <- paste0("item", 1:4)
  list(items = items,
       age   = rnorm(n) + cl,
       bmi   = rnorm(n, 25 + 2 * cl))
}

test_that("measurement_summary returns a long data frame of item parameters", {
  d    <- .sim_summary_data()
  fit  <- suppressMessages(fit_mixture(d$items, n_classes = 2,
                                       measurement = "binary",
                                       n_init = 3,
                                       random_state = 2))
  msdf <- expect_invisible(measurement_summary(fit))

  expect_s3_class(msdf, "data.frame")
  expect_named(msdf, c("block", "parameter", "item", "category", "class",
                       "estimate", "overall"))
  expect_identical(nrow(msdf), 8L)              # 4 items x 2 classes
  expect_setequal(unique(msdf$item), colnames(d$items))
  expect_true(all(msdf$parameter == "probability"))

  # The returned estimates are the fitted item-response probabilities.
  expect_equal(msdf$estimate[msdf$item == "item1"],
               as.vector(fit$mm$parameters$pis[, "item1"]),
               tolerance = 1e-12)

  # `overall` is the observed marginal, not a fitted quantity, and is constant
  # down the class rows so the frame stays joinable on `class`.
  expect_equal(msdf$overall[msdf$item == "item1"],
               rep(mean(d$items[, "item1"]), 2L), tolerance = 1e-12)
  expect_equal(msdf$overall[msdf$class == 1L],
               unname(colMeans(d$items)), tolerance = 1e-12)
})

test_that("scale = \"probability\" is unchanged by the argument's existence", {
  d   <- .sim_summary_data()
  fit <- suppressMessages(fit_mixture(d$items, n_classes = 2,
                                      measurement = "binary",
                                      n_init = 3, random_state = 2))
  out_default <- utils::capture.output(df_default <- measurement_summary(fit))
  out_explicit <- utils::capture.output(
    df_explicit <- measurement_summary(fit, scale = "probability"))
  expect_identical(out_default, out_explicit)
  expect_identical(df_default, df_explicit)
})

test_that("the effect-coded class deviations sum to zero within each item", {
  d   <- .sim_summary_data()
  fit <- suppressMessages(fit_mixture(d$items, n_classes = 2,
                                      measurement = "binary",
                                      n_init = 3, random_state = 2))
  out <- utils::capture.output(msdf <- measurement_summary(fit, scale = "effect"))
  for (it in unique(msdf$item))
    expect_lt(abs(sum(msdf$estimate[msdf$item == it])), 1e-8)

  # The logit scale is qlogis() of the probability-scale estimate.
  out_plain <- utils::capture.output(plain <- measurement_summary(fit))
  out_lg <- utils::capture.output(logit <- measurement_summary(fit, scale = "logit"))
  expect_equal(logit$estimate, qlogis(plain$estimate), tolerance = 1e-10)

  # Neither alternative scale prints an Overall column, since the observed
  # marginal has no meaningful transform there.
  expect_false(any(grepl("Overall", out, fixed = TRUE)))
  expect_false(any(grepl("Overall", out_lg, fixed = TRUE)))
})

test_that("scale = \"effect\" refuses a polytomous item", {
  set.seed(3)
  X <- cbind(poly = sample(1:3, 200, replace = TRUE),
            bin  = rbinom(200, 1, 0.5))
  fit <- fit_mixture(X, n_classes = 2,
                     measurement = list(categorical = "poly", binary = "bin"),
                     n_init = 2)
  expect_error(measurement_summary(fit, scale = "effect"), "polytomous")
  # The logit scale has no such restriction.
  out <- utils::capture.output(msdf <- measurement_summary(fit, scale = "logit"))
  expect_s3_class(msdf, "data.frame")
})

test_that("the observed marginal is weighted, and dropped when unavailable", {
  set.seed(11)
  X <- cbind(a = rnorm(150), b = rnorm(150) + 2)
  w <- runif(150, 0.5, 2)
  fit <- fit_mixture(X, n_classes = 2, measurement = "continuous", n_init = 2,
                     weights = w)
  msdf <- suppressWarnings(measurement_summary(fit))
  ws   <- fit$sample_weights

  # The case weights belong in the benchmark as much as in the estimates.
  expect_equal(msdf$overall[msdf$class == 1L],
               unname(colSums(ws * X) / sum(ws)), tolerance = 1e-12)
  expect_false(isTRUE(all.equal(msdf$overall[msdf$class == 1L],
                                unname(colMeans(X)))))

  # No stored indicators, no benchmark: the column is NA and the table says so
  # rather than printing a column of blanks.
  bare <- fit
  bare$data <- NULL
  out <- utils::capture.output(bare_df <- measurement_summary(bare))
  expect_true(all(is.na(bare_df$overall)))
  expect_false(any(grepl("| Overall", out, fixed = TRUE)))
  expect_true(any(grepl("Overall column", out, fixed = TRUE)))
  expect_true(any(grepl("is omitted above", out, fixed = TRUE)))
})

test_that("class_sizes reports proportions and counts", {
  d   <- .sim_summary_data()
  fit <- suppressMessages(fit_mixture(d$items, n_classes = 2,
                                      measurement = "binary",
                                      n_init = 3,
                                      random_state = 2))
  cs  <- class_sizes(fit)

  expect_s3_class(cs, "data.frame")
  expect_named(cs, c("class", "proportion", "n_expected", "n_modal"))
  expect_equal(cs$proportion, as.vector(fit$weights), tolerance = 1e-12)
  expect_equal(sum(cs$n_expected), nrow(d$items), tolerance = 1e-8)
  expect_equal(sum(cs$n_modal), nrow(d$items), tolerance = 1e-8)
})

test_that("summary returns coefficient and outcome tables invisibly", {
  d    <- .sim_summary_data()
  fit  <- suppressMessages(fit_mixture(d$items, n_classes = 2,
                                       measurement = "binary",
                                       n_init = 3,
                                       random_state = 2))

  # No structural model: NULL, with the notice.
  expect_null(suppressMessages(withVisible(summary(fit)))$value)

  # Covariate model: coefficients table matches coef().
  fitc <- suppressMessages(add_covariates(fit, d$age))
  sv   <- NULL
  invisible(capture.output(sv <- summary(fitc)))
  # (The omnibus element is absent here: with two classes and single-column
  # covariates each omnibus test would duplicate a z-test already shown.)
  expect_true(all(c("ref_class", "n_classes", "n_steps", "correction",
                    "coefficients") %in% names(sv)))
  expect_s3_class(sv$coefficients, "data.frame")
  ors <- coef(fitc)
  expect_equal(sv$coefficients$OR[sv$coefficients$term == "age"],
               unname(ors[2, "age"]), tolerance = 1e-10)

  # Continuous distal outcome: class means with omnibus test.
  fito <- suppressMessages(add_outcome(fit, d$bmi))
  sv2  <- NULL
  invisible(capture.output(sv2 <- summary(fito)))
  expect_identical(sv2$outcome$type, "continuous")
  expect_s3_class(sv2$outcome$means, "data.frame")
  expect_equal(sv2$outcome$means$mean,
               as.vector(fito$sm$parameters$means), tolerance = 1e-12)
  expect_s3_class(sv2$outcome$omnibus, "data.frame")
})
