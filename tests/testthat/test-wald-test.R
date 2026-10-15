# ==============================================================================
# wald_test(): linear constraints on named parameters
# ==============================================================================

.wt_gmm <- local({
  fit <- NULL
  function() {
    if (is.null(fit)) {
      set.seed(1)
      n <- 300
      cls <- rbinom(n, 1, 0.5)
      icept <- ifelse(cls == 1, 3, 1) + rnorm(n, 0, 1)
      slope <- ifelse(cls == 1, 1.0, 0.2) + rnorm(n, 0, 0.4)
      y <- outer(icept, rep(1, 4)) + outer(slope, 0:3) +
        matrix(rnorm(n * 4, 0, 0.7), n, 4)
      fit <<- suppressWarnings(fit_gmm(y, times = 4, n_classes = 2, n_init = 3,
                                       random_state = 1))
    }
    fit
  }
})

test_that("wald_test() is the quadratic form in growth_estimates()'s covariance", {
  fit <- .wt_gmm()
  ge  <- growth_estimates(fit, se = "hessian")
  V   <- attr(ge, "vcov")
  i1  <- which(ge$parameter == "growth_mean" & ge$item == "linear" & ge$class == 1)
  i2  <- which(ge$parameter == "growth_mean" & ge$item == "linear" & ge$class == 2)
  d   <- ge$estimate[i1] - ge$estimate[i2] - 0.5
  ref <- d^2 / (V[i1, i1] + V[i2, i2] - 2 * V[i1, i2])

  w <- wald_test(fit, list(c("growth_mean[linear, 1]" = 1,
                             "growth_mean[linear, 2]" = -1)),
                 rhs = 0.5, se = "hessian")
  expect_s3_class(w, "mixture_wald")
  expect_equal(w$Wald_Chi2, ref, tolerance = 1e-10)
  expect_equal(w$df, 1L)

  # A label pair is the same constraint; the matrix form too.
  w2 <- wald_test(fit, list(c("growth_mean[linear, 1]", "growth_mean[linear, 2]")),
                  rhs = 0.5, se = "hessian")
  R  <- matrix(c(1, -1), 1, dimnames = list(NULL, c("growth_mean[linear, 1]",
                                                    "growth_mean[linear, 2]")))
  w3 <- wald_test(fit, R, rhs = 0.5, se = "hessian")
  expect_equal(w2$Wald_Chi2, ref, tolerance = 1e-10)
  expect_equal(w3$Wald_Chi2, ref, tolerance = 1e-10)
  expect_output(print(w), "CONSTRAINTS")
})

test_that("wald_test() refuses unknown labels and redundant constraints", {
  fit <- .wt_gmm()
  expect_error(wald_test(fit, list(c("growth_mean[slope, 1]",
                                     "growth_mean[linear, 2]"))),
               "Unknown parameter label")
  a <- "growth_mean[intercept, 1]"; b <- "growth_mean[intercept, 2]"
  expect_error(wald_test(fit, list(c(a, b), c(b, a))), "redundant")
  expect_error(wald_test(fit, list(c(a, b)), rhs = c(0, 1)), "rhs")
})

test_that("on a covariate model, wald_test() reproduces analytical_wald_test()", {
  set.seed(3)
  n <- 400
  z <- rnorm(n)
  cls <- ifelse(runif(n) < plogis(0.8 * z), 1, 2)
  p <- rbind(c(.85, .8, .8, .2, .15), c(.15, .2, .2, .8, .85))
  X <- t(vapply(cls, function(k) rbinom(5, 1, p[k, ]), numeric(5)))
  Z <- matrix(z, ncol = 1, dimnames = list(NULL, "age"))
  fit <- suppressMessages(suppressWarnings(fit_mixture(X, Y = Z, n_components = 2,
                                      measurement = "binary",
                                      structural = "covariate", n_steps = 3,
                                      correction = "ML", n_init = 3,
                                      random_state = 1)))
  labs <- grep(":age$", rownames(vcov(fit)), value = TRUE)
  w <- wald_test(fit, lapply(labs, function(l) stats::setNames(1, l)))
  a <- analytical_wald_test(fit, "age")
  expect_equal(w$Wald_Chi2, a$Wald_Chi2, tolerance = 1e-3)
  expect_equal(w$df, a$df)
})

test_that("a one-step continuous outcome's means are labelled for wald_test()", {
  set.seed(2)
  n   <- 300
  cls <- rbinom(n, 1, 0.5)
  y   <- outer(ifelse(cls == 1, 3, 1) + rnorm(n, 0, 0.8), rep(1, 4)) +
    outer(ifelse(cls == 1, 1, 0.2), 0:3) + matrix(rnorm(n * 4, 0, 0.6), n, 4)
  dist <- rnorm(n, ifelse(cls == 1, 2, 0), 1)
  fit <- suppressWarnings(fit_gmm(y, times = 4, n_classes = 2, n_init = 3,
                                  random_state = 1, outcome = dist,
                                  outcome_type = "continuous", n_steps = 1))
  ge <- growth_estimates(fit)
  om <- ge[ge$parameter == "outcome_mean", ]
  ov <- ge[ge$parameter == "outcome_variance", ]
  expect_equal(nrow(om), 2L)
  expect_equal(om$estimate, as.vector(fit$sm$parameters$means), tolerance = 1e-10)
  # One shared variance, repeated across the classes.
  expect_equal(ov$estimate[1], ov$estimate[2])
  lab <- sprintf("outcome_mean[%s, %d]", om$item, om$class)
  V   <- attr(ge, "vcov")
  i   <- which(ge$parameter == "outcome_mean")
  ref <- diff(om$estimate)^2 / (V[i[1], i[1]] + V[i[2], i[2]] - 2 * V[i[1], i[2]])
  expect_equal(wald_test(fit, list(lab))$Wald_Chi2, ref, tolerance = 1e-10)
})
