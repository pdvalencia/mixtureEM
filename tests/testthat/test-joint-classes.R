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
})

test_that("a smaller ca: without covariates, the same model as the swap", {
  d <- sim_joint(500, 2, 3, 4, 5, seed = 25)
  args <- list(n_init = 5, random_state = 2, smoothing = 0,
               bayes_constants = list(categorical = 0), standard_errors = FALSE)
  f <- do.call(fit_joint_classes, c(list(d$A, d$S, n_ca = 2, n_cs = 3), args))
  g <- do.call(fit_joint_classes, c(list(d$S, d$A, n_ca = 3, n_cs = 2), args))
  expect_lt(abs(f$metrics$ll - g$metrics$ll), 1e-4)
  # (2 - 1) ca + 2 x (3 - 1) association + 2 x 4 + 3 x 5 items = 28, and the
  # swap counts the same model the other way round.
  expect_equal(f$n_params, 28)
  expect_equal(g$n_params, 28)
  expect_equal(f$delta[3], 0)
  expect_equal(dim(f$joint), c(2L, 3L))
  expect_equal(sort(as.vector(f$joint)), sort(as.vector(g$joint)),
               tolerance = 1e-3)
})

test_that("a smaller ca with covariates on both: the padded class is no category", {
  d <- sim_joint(400, 2, 3, 4, 5, seed = 26)
  z <- data.frame(w = rbinom(400, 1, 0.5))
  f <- fit_joint_classes(d$A, d$S, n_ca = 2, n_cs = 3, predictors_ca = z,
                         predictors_cs = z, cs_effects = c(w = "by_ca"),
                         n_init = 3, random_state = 5)
  # ca logit on [1, w] = 2; items 2 x 4 + 3 x 5 = 23; two cs logits on
  # [1, CA#1, w, w x CA#1] = 8. 33 in all.
  expect_equal(f$n_params, 33)
  expect_equal(dim(f$delta_beta), c(2L, 2L))
  expect_equal(dim(f$tau_beta[[1]]), c(3L, 4L))
  expect_true(all(exp(.lta_log_delta(f))[, 3] < 1e-290))
  expect_equal(unname(rowSums(f$conditional)), rep(1, 2))
  a <- f$association
  expect_equal(nrow(a), 8)
  expect_true(all(is.finite(a$se) & a$se > 0))
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
  d$S[, 1:2] <- d$S[, 1:2] + rbinom(400 * 2, 1, 0.3)     # paired alike
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

test_that("a binary item paired with a three-category one: no empty category", {
  # Paired by position, the binary `cs` items would inherit the third
  # category of the `ca` items beside them; the stacked layout gives each
  # item its own.
  d <- sim_joint(400, 3, 2, 5, 4, seed = 30)
  d$A[, 1:2] <- d$A[, 1:2] + rbinom(400 * 2, 1, 0.3)
  f <- fit_joint_classes(d$A, d$S, n_ca = 3, n_cs = 2,
                         measurement_ca = "ordinal", measurement_cs = "binary",
                         n_init = 2, random_state = 7, smoothing = 0,
                         bayes_constants = list(categorical = 0),
                         standard_errors = FALSE)
  expect_true(f$joint_stacked)
  # (3 - 1) ca + 3 x (2 - 1) association + 3 x (2 x 2 + 3) + 2 x 4 items.
  expect_equal(f$n_params, 34)
})

test_that("no association: the two separate LCAs, tested on (K1 - 1)(K2 - 1) df", {
  d <- sim_joint(500, 3, 2, 5, 4, seed = 29)
  zero <- list(categorical = 0)
  f0 <- fit_joint_classes(d$A, d$S, n_ca = 3, n_cs = 2, association = "none",
                          n_init = 5, random_state = 8, smoothing = 0,
                          bayes_constants = zero)
  f1 <- fit_joint_classes(d$A, d$S, n_ca = 3, n_cs = 2, n_init = 5,
                          random_state = 8, smoothing = 0,
                          bayes_constants = zero)
  la <- fit_mixture(d$A, n_classes = 3, measurement = "binary", n_init = 5,
                    random_state = 8, smoothing = 0, bayes_constants = zero)
  ls <- fit_mixture(d$S, n_classes = 2, measurement = "binary", n_init = 5,
                    random_state = 8, smoothing = 0, bayes_constants = zero)
  expect_equal(f0$metrics$ll, la$metrics$ll + ls$metrics$ll, tolerance = 1e-5)
  # (3 - 1) ca + (2 - 1) cs + 3 x 5 + 2 x 4 items.
  expect_equal(f0$n_params, 26)
  lt <- lr_test(f0, f1)
  expect_equal(lt$df, 2)
  expect_equal(lt$statistic, 2 * (f1$metrics$ll - f0$metrics$ll))
  # The cs share is one logit, its SE the observed information's.
  expect_equal(f0$association$term, "Intercept")
  expect_equal(f0$association$estimate,
               log(f0$conditional[1, 1] / f0$conditional[1, 2]), tolerance = 1e-8)
  expect_true(is.finite(f0$association$se) && f0$association$se > 0)
  expect_equal(unname(f0$conditional[1, ]), unname(f0$conditional[3, ]))
  expect_output(print(f0), "No association")
  z <- data.frame(w = rbinom(500, 1, 0.5))
  expect_error(fit_joint_classes(d$A, d$S, n_ca = 3, n_cs = 2,
                                 predictors_cs = z, cs_effects = "by_ca",
                                 association = "none"),
               "needs the association")
})

test_that("no association, two-step: tested against the association on one step 1", {
  d <- sim_joint(500, 3, 2, 5, 4, seed = 31)
  z <- data.frame(w = rbinom(500, 1, 0.5))
  zero <- list(categorical = 0)
  fit <- function(...) fit_joint_classes(d$A, d$S, n_ca = 3, n_cs = 2,
                                         n_steps = 2, n_init = 5,
                                         random_state = 9, smoothing = 0,
                                         bayes_constants = zero, ...)
  f1 <- fit(predictors_cs = z)
  f0 <- fit(predictors_cs = z, association = "none")
  expect_identical(f0$mm$models, f1$mm$models)
  expect_true(isTRUE(f0$se$twostep))
  expect_setequal(f0$association$term, c("Intercept", "w"))
  expect_equal(lr_test(f0, f1)$df, 2)
  # Step 1 is the model with every structural parameter at zero, counted
  # without the padding: the association and w's slope.
  expect_equal(lr_test(f1$step1, f1)$df, 3)
  # Without covariates the no-association model is step 1 itself.
  g0 <- fit(association = "none")
  expect_equal(g0$metrics$ll, f1$step1$metrics$ll, tolerance = 1e-8)
  expect_equal(lr_test(g0, f1)$df, 3)
})

test_that("two-step: a step 1 handed over is reused, or continued from the two LCAs", {
  d <- sim_joint(500, 3, 2, 5, 4, seed = 33)
  z <- data.frame(w = rbinom(500, 1, 0.5))
  zero <- list(categorical = 0)
  fit <- function(...) fit_joint_classes(d$A, d$S, n_ca = 3, n_cs = 2,
                                         n_steps = 2, n_init = 3,
                                         random_state = 4, smoothing = 0,
                                         bayes_constants = zero, ...)
  f1 <- fit(predictors_cs = z)
  # Step 1 searched as the two LCAs is still the joint step-1 model at their
  # solution, and an EM run from it moves nothing.
  la <- fit_mixture(d$A, n_classes = 3, measurement = "binary", n_init = 3,
                    random_state = 4, smoothing = 0,
                    bayes_constants = list(categorical = 0, latent = 0))
  ls <- fit_mixture(d$S, n_classes = 2, measurement = "binary", n_init = 3,
                    random_state = 4, smoothing = 0,
                    bayes_constants = list(categorical = 0, latent = 0))
  expect_equal(f1$step1$metrics$ll, la$metrics$ll + ls$metrics$ll,
               tolerance = 1e-6)
  expect_lte(f1$step1$n_iter, 5)

  # A previous fit's step 1, reused: the same fit as one searching its own.
  f2  <- fit(predictors_cs = z, cs_effects = c(w = "by_ca"))
  f2r <- fit(predictors_cs = z, cs_effects = c(w = "by_ca"), step1 = f1)
  expect_identical(f2r$step1$mm$models, f2$step1$mm$models)
  expect_equal(f2r$metrics$ll, f2$metrics$ll, tolerance = 1e-10)
  expect_equal(f2r$n_params, f2$n_params)
  expect_equal(f2r$step1$n_params, f1$step1$n_params)
  expect_equal(lr_test(f1, f2r)$df, lr_test(f1, f2)$df)

  # The two LCAs, handed over.
  f3 <- fit(predictors_cs = z, step1 = list(ca = la, cs = ls))
  expect_equal(f3$metrics$ll, f1$metrics$ll, tolerance = 1e-6)
  # The redirected no-association fit is step 1 itself.
  g0 <- fit(association = "none", step1 = f1)
  expect_equal(g0$metrics$ll, f1$step1$metrics$ll, tolerance = 1e-6)

  expect_error(fit(predictors_cs = z, step1 = list(ca = ls, cs = la)), "n_ca")
  expect_error(fit_joint_classes(d$A[, -1], d$S, n_ca = 3, n_cs = 2,
                                 n_steps = 2, step1 = f1), "other items")
  expect_error(fit_joint_classes(d$A, d$S, n_ca = 3, n_cs = 2, step1 = f1),
               "n_steps = 2")
})

test_that("the default prior is spread over each variable's own classes, never the padding", {
  # Both directions: a smaller `cs` (padded destination) with fewer `cs`
  # items, and a smaller `ca` (padded origin) with fewer `ca` items. The
  # penalty's value and gradient must be the log-prior's, the fit must be
  # stationary on the penalised objective, and moving a padded cell -- a
  # padded class's item probabilities or a padded item column -- must leave
  # the log-prior where it was.
  cases <- list(list(K1 = 3, K2 = 2, Ja = 5, Js = 4),
                list(K1 = 2, K2 = 3, Ja = 4, Js = 5))
  for (cc in cases) {
    d   <- sim_joint(400, cc$K1, cc$K2, cc$Ja, cc$Js, seed = 31)
    fit <- fit_joint_classes(d$A, d$S, n_ca = cc$K1, n_cs = cc$K2, n_init = 3,
                             random_state = 1, standard_errors = FALSE)
    X   <- as.matrix(fit$data)
    w   <- fit$weights_vec
    lay <- .lta_par_layout(fit)
    par <- .lta_par_pack(fit, lay)
    lp  <- function(p) .lta_log_prior(.lta_par_unpack(p, fit, lay), X, 1)
    obj <- function(p) sum(w * .lta_ll_case(fit, X, p, lay)) + lp(p)
    fd  <- function(f) vapply(seq_along(par), function(j) {
      e <- numeric(length(par)); e[j] <- 1e-6
      (f(par + e) - f(par - e)) / 2e-6 }, 0)
    info <- paste(cc$K1, "x", cc$K2)
    expect_lt(max(abs(fd(obj))), 1e-3, label = info)
    pen <- .lta_penalty(fit, X, lay, 1)
    expect_equal(pen$value, lp(par), tolerance = 1e-10, info = info)
    expect_lt(max(abs(pen$gradient - fd(lp))), 1e-6, label = info)

    # The smaller variable also has the shorter item list here, so one
    # occasion holds both kinds of padding: its last status and last item.
    K  <- fit$n_statuses
    t_pad <- if (cc$K1 < K) 1L else 2L
    moved <- fit
    pis <- moved$mm$models[[t_pad]]$parameters$pis
    pis[K, ] <- 0.37
    pis[, max(cc$Ja, cc$Js)] <- 0.21
    moved$mm$models[[t_pad]]$parameters$pis <- pis
    expect_equal(.lta_log_prior(moved, X, 1), .lta_log_prior(fit, X, 1),
                 tolerance = 1e-12, info = info)
  }
})
