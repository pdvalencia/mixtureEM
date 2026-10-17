# measurement_estimates() and the error bars on plot(type = "bar").

.me_zero <- list(latent = 0, categorical = 0, poisson = 0, variances = 0)

# Two classes far apart, so every case is assigned with certainty and a class
# mean's standard error is the textbook sd / sqrt(n_k).
.me_lpa <- function(variances_equal = FALSE) {
  set.seed(11)
  X <- rbind(matrix(rnorm(300, 0, 1), 150), matrix(rnorm(150, 10, 2), 75))
  colnames(X) <- c("a", "b")
  fit_mixture(X, n_classes = 2, measurement = "continuous",
              variances_equal = variances_equal, bayes_constants = .me_zero,
              n_init = 3, random_state = 1, n_cores = 1)
}

test_that("a separated class mean's standard error is sd / sqrt(n_k)", {
  fit <- .me_lpa()
  me  <- measurement_estimates(fit, "hessian")
  mu  <- me[me$parameter == "mean", ]
  sdk <- sqrt(fit$mm$parameters$covariances)
  nk  <- fit$weights * nrow(fit$data)
  expect_equal(mu$se, as.vector(sdk / sqrt(nk)), tolerance = 1e-3)
  # A variance's is var * sqrt(2 / n_k) by the same argument.
  va <- me[me$parameter == "variance", ]
  expect_equal(va$se, as.vector(sdk^2 * sqrt(2 / nk)), tolerance = 1e-3)
  expect_named(me, c("block", "parameter", "item", "category", "class",
                     "estimate", "se", "z", "p_value", "boundary"))
  expect_equal(dim(attr(me, "vcov")), c(nrow(me), nrow(me)))
})

test_that("robust and hessian agree on a correctly specified model", {
  fit <- .me_lpa()
  h <- measurement_estimates(fit, "hessian")
  r <- measurement_estimates(fit, "robust")
  expect_equal(r$estimate, h$estimate)
  expect_true(all(abs(r$se / h$se - 1) < 0.3))
})

test_that("a variance held equal across classes has one standard error", {
  me <- measurement_estimates(.me_lpa(variances_equal = TRUE), "hessian")
  va <- me[me$parameter == "variance", ]
  for (it in unique(va$item))
    expect_equal(diff(va$se[va$item == it]), 0)
})

test_that("a probability at the boundary is flagged, held fixed, and named", {
  set.seed(3)
  X <- rbind(cbind(1, matrix(rbinom(600, 1, 0.8), 150)),
             cbind(rbinom(150, 1, 0.5), matrix(rbinom(600, 1, 0.2), 150)))
  fit <- fit_mixture(X, n_classes = 2, measurement = "binary",
                     bayes_constants = .me_zero, n_init = 3,
                     random_state = 1, n_cores = 1)
  expect_message(me <- measurement_estimates(fit, "hessian"), "boundary")
  b <- me[me$boundary, ]
  expect_true(nrow(b) >= 1L)
  expect_true(all(b$parameter == "probability"))
  # The rest keep finite standard errors.
  expect_true(all(is.finite(me$se[!me$boundary])))
})

test_that("growth fits and other objects are refused with a pointer", {
  expect_error(measurement_estimates(structure(list(), class = c("gmm", "mixture_model"))),
               "growth_estimates")
  expect_error(measurement_estimates(list()), "fit_mixture")
})

test_that("the bar heights carry their standard errors on the bar scale", {
  fit <- .me_lpa()
  H   <- .profile_bar_heights(fit, se = TRUE)
  me  <- measurement_estimates(fit)
  sj  <- apply(fit$data, 2, function(x) sqrt(mean((x - mean(x))^2)))
  expect_equal(attr(H, "se"),
               sweep(matrix(me$se[me$parameter == "mean"], 2), 2, sj, "/"),
               ignore_attr = TRUE)
  expect_null(attr(.profile_bar_heights(fit), "se"))

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_invisible(plot(fit, type = "bar", error_bars = TRUE))
  expect_error(plot(fit, type = "line", error_bars = TRUE), "bar")
})

test_that("the prior's curvature counts at default priors and vanishes without one", {
  set.seed(5)
  X <- rbind(matrix(rnorm(120, 0, 1), 40), matrix(rnorm(120, 1.5, 0.6), 40))
  plain <- function(fit) {
    v <- .reported_vcov(fit, .measurement_reported, "measurement model",
                        drop_walls = TRUE, prior = FALSE)
    sqrt(diag(v$hessian))
  }
  on  <- fit_mixture(X, n_classes = 2, measurement = "continuous",
                     variances_equal = FALSE, n_init = 3, random_state = 1,
                     n_cores = 1)
  off <- fit_mixture(X, n_classes = 2, measurement = "continuous",
                     variances_equal = FALSE, bayes_constants = .me_zero,
                     n_init = 3, random_state = 1, n_cores = 1)
  s_on <- measurement_estimates(on, "hessian")$se
  expect_gt(max(abs(s_on / plain(on) - 1)), 1e-3)
  expect_equal(measurement_estimates(off, "hessian")$se, plain(off),
               tolerance = 1e-10)
})
