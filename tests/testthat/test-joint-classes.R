# ==============================================================================
# Two latent class variables in one model: fit_joint_classes()
# ==============================================================================
#
# Structural checks. The fit is the padded two-occasion LTA, so it must keep
# that fit's likelihood exactly while counting only the parameters that touch
# data, and its coefficients must be the reference-class contrasts of the
# fitted table.

sim_joint <- function(n, K1, K2, Ja, Js, seed) {
  set.seed(seed)
  ca <- sample(seq_len(K1), n, TRUE)
  cs <- vapply(ca, function(a) sample(seq_len(K2), 1, prob = c(4, rep(1, K2 - 1))[
    (seq_len(K2) - a) %% K2 + 1]), integer(1))
  A <- sapply(seq_len(Ja), function(j) rbinom(n, 1, c(0.85, 0.15, 0.5)[ca]))
  S <- sapply(seq_len(Js), function(j) rbinom(n, 1, c(0.85, 0.15, 0.5)[cs]))
  list(A = A, S = S)
}

test_that("unequal classes and items: same likelihood as the padded fit, padding not counted", {
  d <- sim_joint(500, 3, 2, 5, 4, seed = 21)
  f <- fit_joint_classes(d$A, d$S, n_ca = 3, n_cs = 2, n_init = 5,
                         random_state = 2, smoothing = 0,
                         bayes_constants = list(categorical = 0))
  X <- cbind(d$A, d$S, NA)
  forb <- matrix(FALSE, 3, 3); forb[, 3] <- TRUE
  g <- fit_lta(X, n_statuses = 3, times = 2, measurement = "binary",
               measurement_invariance = "none", forbidden_transitions = forb,
               n_init = 5, random_state = 2, smoothing = 0,
               bayes_constants = list(categorical = 0), standard_errors = FALSE)
  expect_equal(f$metrics$ll, g$metrics$ll, tolerance = 1e-6)
  # (3 - 1) ca + 3 x (2 - 1) association + 3 x 5 + 2 x 4 items = 28; the
  # padded fit also counts the padded item (3) and the padded cs class (4).
  expect_equal(f$n_params, 28)
  expect_equal(g$n_params - f$n_params, 7)
  expect_equal(f$metrics$bic, -2 * f$metrics$ll + log(500) * 28)

  expect_equal(unname(rowSums(f$conditional)), rep(1, 3))
  expect_equal(sum(f$joint), 1, tolerance = 1e-8)
})

test_that("the coefficients are the reference-class contrasts of the fitted table", {
  d <- sim_joint(500, 3, 3, 4, 4, seed = 22)
  f <- fit_joint_classes(d$A, d$S, n_ca = 3, n_cs = 3, n_init = 5,
                         random_state = 3)
  a  <- f$association
  lt <- log(f$conditional)
  lo <- function(k, j) lt[k, j] - lt[k, 3]
  expect_equal(nrow(a), 6)
  expect_equal(a$estimate[a$cs == "CS#2" & a$term == "Intercept"], lo(3, 2),
               tolerance = 1e-8)
  expect_equal(a$estimate[a$cs == "CS#1" & a$term == "CA#2"],
               lo(2, 1) - lo(3, 1), tolerance = 1e-8)
  expect_true(all(is.finite(a$se) & a$se > 0))
  expect_output(print(f), "P\\(cs \\| ca\\)")
})

test_that("cs_effects names covariates and is checked", {
  d <- sim_joint(300, 2, 2, 4, 4, seed = 23)
  z <- data.frame(w = rnorm(300), x = rnorm(300))
  f <- fit_joint_classes(d$A, d$S, n_ca = 2, n_cs = 2, predictors_cs = z,
                         cs_effects = c(w = "by_ca"), n_init = 3,
                         random_state = 4)
  a <- f$association
  expect_setequal(a$term, c("Intercept", "CA#1", "w", "x"))
  expect_equal(sum(a$term == "w"), 2)          # base and the CA#1 difference
  expect_error(fit_joint_classes(d$A, d$S, predictors_cs = z,
                                 cs_effects = c(v = "by_ca")), "not in")
  expect_error(fit_joint_classes(d$A, d$S, n_ca = 2, n_cs = 3), "swap")
})
