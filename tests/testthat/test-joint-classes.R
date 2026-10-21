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

test_that("covariates on a smaller cs: the padded class has no category", {
  d <- sim_joint(400, 3, 2, 5, 4, seed = 24)
  z <- data.frame(w = rbinom(400, 1, 0.5))
  f <- fit_joint_classes(d$A, d$S, n_ca = 3, n_cs = 2, predictors_cs = z,
                         n_init = 3, random_state = 5)
  # (3 - 1) ca + 3 x 5 + 2 x 4 items + one cs logit on [1, 2 CA, w] = 29.
  expect_equal(f$n_params, 29)
  expect_equal(dim(f$tau_beta[[1]]), c(2L, 4L))
  lt <- .lta_log_tau(f)[[1]]
  expect_true(all(vapply(lt, function(M) all(exp(M[, 3]) < 1e-290), TRUE)))
  expect_equal(unname(rowSums(f$conditional)), rep(1, 3))
  a <- f$association
  expect_equal(a$term, c("Intercept", "CA#1", "CA#2", "w"))
  expect_true(all(is.finite(a$se) & a$se > 0))
  # a forbidden move that is not a whole destination still has no logit form
  forb <- matrix(FALSE, 3, 3); forb[1, 3] <- TRUE
  expect_error(fit_lta(cbind(d$A, d$S, NA), n_statuses = 3, times = 2,
                       forbidden_transitions = forb, predictors_transition = z,
                       n_init = 1), "whole destinations")
})

test_that("ordinal items with unequal item counts: the padded item stays finite", {
  # The padding is an all-missing item column; at zero priors its counts are
  # 0/0, which once turned every start's log-likelihood into NaN.
  d <- sim_joint(400, 3, 2, 5, 4, seed = 25)
  d$A[, 1:2] <- d$A[, 1:2] + rbinom(400 * 2, 1, 0.3)     # three categories
  f <- fit_joint_classes(d$A, d$S, n_ca = 3, n_cs = 2,
                         measurement_ca = "ordinal", measurement_cs = "ordinal",
                         n_init = 2, random_state = 6, smoothing = 0,
                         bayes_constants = list(categorical = 0),
                         standard_errors = FALSE)
  expect_true(is.finite(f$metrics$ll))
  expect_false(anyNA(f$mm$models[[2]]$parameters$pis))
})

test_that("two-step: step 1 is the two separate LCAs, step 2 keeps their items", {
  d <- sim_joint(500, 3, 2, 5, 4, seed = 25)
  z <- data.frame(w = rbinom(500, 1, 0.5))
  zero <- list(categorical = 0)
  f <- fit_joint_classes(d$A, d$S, n_ca = 3, n_cs = 2, n_steps = 2,
                         predictors_cs = z, cs_effects = c(w = "by_ca"),
                         n_init = 5, random_state = 6, smoothing = 0,
                         bayes_constants = zero)
  # With no association the joint likelihood factorises, so step 1's is the
  # sum of the two separate analyses'.
  la <- fit_mixture(d$A, n_classes = 3, measurement = "binary", n_init = 5,
                    random_state = 6, smoothing = 0, bayes_constants = zero)
  ls <- fit_mixture(d$S, n_classes = 2, measurement = "binary", n_init = 5,
                    random_state = 6, smoothing = 0, bayes_constants = zero)
  expect_equal(f$step1$metrics$ll, la$metrics$ll + ls$metrics$ll,
               tolerance = 1e-5)
  expect_identical(f$mm$models, f$step1$mm$models)
  expect_gte(f$metrics$ll, f$step1$metrics$ll)
  expect_true(isTRUE(f$se$twostep))
  idx2 <- .lta_par_split(.lta_par_layout(f))$structural
  expect_true(all(diag(f$se$vcov)[idx2] >= diag(f$se$twostep_V2) - 1e-12))
  expect_output(print(f), "Two-step")
  expect_error(fit_joint_classes(d$A, d$S, n_steps = 3), "1 or 2")
})

test_that("item slopes: step 1 is the two separate LCAs with their DIF, step 2 keeps them", {
  d <- sim_joint(600, 2, 2, 4, 4, seed = 27)
  set.seed(28)
  x <- rbinom(600, 1, 0.5)
  # Shift one item on each side for x = 1, whatever the class.
  d$A[, 2] <- ifelse(x == 1 & runif(600) < 0.3, 1, d$A[, 2])
  d$S[, 3] <- ifelse(x == 1 & runif(600) < 0.3, 1, d$S[, 3])
  colnames(d$A) <- paste0("a", 1:4); colnames(d$S) <- paste0("s", 1:4)
  px <- data.frame(x = x)
  zero <- list(categorical = 0)
  pi <- list(a2 = px, s3 = px)
  f <- fit_joint_classes(d$A, d$S, n_steps = 2, predictors_ca = px,
                         predictors_cs = px, predictors_items = pi,
                         n_init = 5, random_state = 3, smoothing = 0,
                         bayes_constants = zero)
  # Step 1 regresses each variable on x beside the item slopes, and without
  # the association the likelihood factorises into the two LCAs.
  lca <- function(Y, it) fit_mixture(Y + 1, n_classes = 2, measurement = "ordinal",
                                     cats = rep(2L, 4),
                                     predictors = ~ x, data = px,
                                     predictors_items = setNames(list(~ x), it),
                                     n_init = 5, random_state = 3, smoothing = 0,
                                     bayes_constants = zero)
  la <- lca(d$A, "a2"); ls <- lca(d$S, "s3")
  expect_equal(f$step1$metrics$ll, la$metrics$ll + ls$metrics$ll,
               tolerance = 1e-5)
  expect_identical(f$dif$beta, f$step1$dif$beta)
  expect_true(all(f$dif$beta[, c(2, 7), 1] != 0))
  expect_true(isTRUE(f$se$twostep))
  idx2 <- .lta_par_split(.lta_par_layout(f))$structural
  expect_true(all(diag(f$se$vcov)[idx2] >= diag(f$se$twostep_V2) - 1e-12))

  # One step: two slopes beyond the model without them.
  g1 <- fit_joint_classes(d$A, d$S, predictors_ca = px, predictors_cs = px,
                          predictors_items = pi, n_init = 5, random_state = 3,
                          smoothing = 0, bayes_constants = zero)
  g0 <- fit_joint_classes(d$A, d$S, predictors_ca = px, predictors_cs = px,
                          n_init = 5, random_state = 3, smoothing = 0,
                          bayes_constants = zero, standard_errors = FALSE)
  expect_equal(g1$n_params, g0$n_params + 2)
  expect_gte(g1$metrics$ll, f$metrics$ll - 1e-6)

  expect_error(fit_joint_classes(d$A, d$S, predictors_items = px), "named list")
  expect_error(fit_joint_classes(unname(d$A), d$S, predictors_items = pi),
               "column names")
})
