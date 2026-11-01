# ==============================================================================
# derived_estimates() and joint_probabilities()
# ==============================================================================
#
# A weighted sum's standard error is exact (c'Vc), a smooth function's is the
# delta method's; a single linear restriction's Wald statistic is the square
# of its z. The conditional probabilities are the softmax of the coefficient
# sums they are built from, and sum to one over `y`.

sim_moderated <- function(n, seed) {
  set.seed(seed)
  w  <- rbinom(n, 1, 0.5)
  x <- sample(1:2, n, TRUE)
  y <- ifelse(runif(n) < plogis(-1 + 2 * (x == 1) + w * (x == 1)), 1L, 2L)
  A  <- sapply(1:4, function(j) rbinom(n, 1, ifelse(x == 1, 0.85, 0.15)))
  S  <- sapply(1:4, function(j) rbinom(n, 1, ifelse(y == 1, 0.85, 0.15)))
  list(A = A, S = S, z = data.frame(w = w))
}

test_that("linear and smooth derived parameters of a joint model", {
  d <- sim_moderated(800, seed = 41)
  f <- fit_joint_classes(d$A, d$S, predictors_y = d$z, y_effects = "by_x",
                         n_init = 3, random_state = 1, smoothing = 0,
                         bayes_constants = list(categorical = 0))
  pe <- .wald_parameters(f, "robust")
  expect_setequal(names(pe$estimate), c("Y#1 Intercept", "Y#1 ON X#1",
                                        "Y#1 ON w", "Y#1 ON w | X#1"))
  # The rows' own SEs come back on the diagonal.
  expect_equal(sqrt(diag(pe$vcov)), f$association$se, ignore_attr = TRUE,
               tolerance = 1e-10)

  cw <- c("Y#1 ON X#1" = 1, "Y#1 ON w | X#1" = 1)
  d1 <- derived_estimates(f, assoc_w1 = cw,
    logit = function(p) p[["Y#1 Intercept"]] + p[["Y#1 ON X#1"]])
  V  <- pe$vcov[names(cw), names(cw)]
  expect_equal(d1$estimate[1], sum(pe$estimate[names(cw)]))
  expect_equal(d1$se[1], sqrt(sum(V)), tolerance = 1e-10)
  # A linear function given as a function: the same SE, to the differencing.
  expect_equal(d1$se[2], sqrt(sum(pe$vcov[c("Y#1 Intercept", "Y#1 ON X#1"),
                                          c("Y#1 Intercept", "Y#1 ON X#1")])),
               tolerance = 1e-6)
  expect_equal(d1$upper - d1$estimate, qnorm(0.975) * d1$se)
  wt <- wald_test(f, list(cw))
  expect_equal(wt$Wald_Chi2, d1$z[1]^2, tolerance = 1e-8)

  # P(y | x, w): the softmax of the coefficient sums, rows summing to one.
  jp <- joint_probabilities(f, at = list(w = 1))
  expect_equal(nrow(jp), 4)
  expect_equal(tapply(jp$estimate, jp$x, sum), c("X#1" = 1, "X#2" = 1),
               ignore_attr = TRUE)
  e <- pe$estimate
  expect_equal(jp$estimate[jp$x == "X#1" & jp$y == "Y#1"],
               plogis(e[["Y#1 Intercept"]] + e[["Y#1 ON X#1"]] +
                      e[["Y#1 ON w"]] + e[["Y#1 ON w | X#1"]]),
               tolerance = 1e-10)
  jp0 <- joint_probabilities(f)
  expect_equal(jp0$estimate[jp0$x == "X#2" & jp0$y == "Y#1"],
               plogis(e[["Y#1 Intercept"]]), tolerance = 1e-10)
  # Two classes of y: the two probabilities of a row share one SE.
  expect_equal(jp$se[1], jp$se[2], tolerance = 1e-6)

  expect_error(derived_estimates(f, c("Y#1 ON X#1" = 1)), "named argument")
  expect_error(derived_estimates(f, x = c(nope = 1)), "Unknown parameter")
  expect_error(joint_probabilities(f, at = list(v = 1)), "covariates of `y`")
})
