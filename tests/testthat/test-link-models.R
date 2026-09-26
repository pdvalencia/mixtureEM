# link_models(): the three-step between two different mixture models.

.link_sim <- function(n = 600, seed = 11) {
  set.seed(seed)
  s1 <- sample(1:3, n, replace = TRUE)
  s2 <- ifelse(runif(n) < c(0.85, 0.5, 0.15)[s1], 1L, 2L)
  X1 <- sapply(1:6, function(j) rbinom(n, 1, c(0.9, 0.5, 0.1)[s1]))
  X2 <- sapply(1:4, function(j) rnorm(n, c(0, 2.5)[s2]))
  list(s1 = s1, s2 = s2, X1 = X1, X2 = X2)
}

test_that("link_models() fits the saturated K1 x K2 table", {
  d <- .link_sim()
  a <- suppressMessages(fit_mixture(d$X1, n_classes = 3, measurement = "binary",
                                    n_init = 5, random_state = 1))
  b <- suppressMessages(fit_mixture(d$X2, n_classes = 2,
                                    measurement = "continuous",
                                    n_init = 5, random_state = 1))
  fit <- link_models(list(early = a, late = b), random_state = 1)

  expect_s3_class(fit, "linked_model")
  expect_equal(fit$n_parameters, 2L + 1L + 2L)
  expect_equal(dim(fit$transitions[[1]]), c(3L, 2L))
  expect_equal(unname(rowSums(fit$transitions[[1]])), rep(1, 3))
  expect_equal(sum(fit$initial), 1)
  expect_length(fit$se, 5L)
  expect_s3_class(logLik(fit), "logLik")
  expect_output(print(fit), "3 -> 2")

  # With no correction, step 3 believes the assignments, and the saturated
  # model reproduces the cross-table of the modal classes.
  naive <- link_models(list(a, b), correction = "none")
  tab <- prop.table(table(fit$modal[, 1], fit$modal[, 2]), 1)
  expect_equal(unname(naive$transitions[[1]]), unname(unclass(tab)),
               tolerance = 1e-5)

  # A covariate adds one slope per initial logit and one per destination.
  x <- rnorm(nrow(d$X1))
  cov <- link_models(list(a, b), predictors_initial = x,
                     predictors_transition = x, standard_errors = FALSE)
  expect_equal(cov$n_parameters, 5L + 2L + 1L)
})

test_that("link_models() refuses models fitted on different cases", {
  d <- .link_sim(n = 200)
  a <- suppressMessages(fit_mixture(d$X1, n_classes = 2, measurement = "binary",
                                    n_init = 2, random_state = 1))
  b <- suppressMessages(fit_mixture(d$X2[1:150, ], n_classes = 2,
                                    measurement = "continuous",
                                    n_init = 2, random_state = 1))
  expect_error(link_models(list(a, b)), "same cases")
  expect_error(link_models(a), "at least two")
})
