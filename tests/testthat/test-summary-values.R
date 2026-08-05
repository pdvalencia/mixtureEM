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
  fit  <- suppressMessages(fit_mixture(d$items, n_classes = 2, n_init = 3,
                                       random_state = 2))
  msdf <- expect_invisible(measurement_summary(fit))

  expect_s3_class(msdf, "data.frame")
  expect_named(msdf, c("block", "parameter", "item", "category", "class",
                       "estimate"))
  expect_identical(nrow(msdf), 8L)              # 4 items x 2 classes
  expect_setequal(unique(msdf$item), colnames(d$items))
  expect_true(all(msdf$parameter == "probability"))

  # The returned estimates are the fitted item-response probabilities.
  expect_equal(msdf$estimate[msdf$item == "item1"],
               as.vector(fit$mm$parameters$pis[, "item1"]),
               tolerance = 1e-12)
})

test_that("class_sizes reports proportions and counts", {
  d   <- .sim_summary_data()
  fit <- suppressMessages(fit_mixture(d$items, n_classes = 2, n_init = 3,
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
  fit  <- suppressMessages(fit_mixture(d$items, n_classes = 2, n_init = 3,
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
