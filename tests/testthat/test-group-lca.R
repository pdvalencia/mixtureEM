# Multiple-group LCA (`group=`/`group_effects=` on fit_mixture()), Collins &
# Lanza (2010) sec. 5.7-5.12. Validated against synthetic data with a known
# per-group generating model: two groups with different item-response
# probabilities AND different class prevalences.

.make_group_data <- function(seed = 7, n = 800) {
  set.seed(seed)
  J <- 4L; K <- 2L
  grp  <- factor(sample(c("A", "B"), n, replace = TRUE))
  prev <- list(A = c(.7, .3), B = c(.3, .7))
  pis  <- list(
    A = list(c1 = rep(.1, J), c2 = rep(.9, J)),
    B = list(c1 = rep(.2, J), c2 = rep(.8, J))
  )
  X <- matrix(NA_real_, n, J)
  for (i in seq_len(n)) {
    g  <- as.character(grp[i])
    cl <- sample(1:K, 1, prob = prev[[g]])
    p  <- if (cl == 1) pis[[g]]$c1 else pis[[g]]$c2
    X[i, ] <- rbinom(J, 1, p)
  }
  colnames(X) <- paste0("item", 1:J)
  list(X = X, grp = grp, J = J, K = K)
}

test_that("group_effects='prevalence' matches predictors= (Collins & Lanza sec. 6.10.2)", {
  d <- .make_group_data()
  fit_group <- fit_mixture(d$X, n_classes = d$K, measurement = "binary",
                           group = d$grp, group_effects = "prevalence",
                           n_steps = 1, n_init = 5, random_state = 1)
  fit_pred  <- fit_mixture(d$X, n_classes = d$K, measurement = "binary",
                           predictors = d$grp,
                           n_steps = 1, n_init = 5, random_state = 1)
  expect_equal(fit_group$metrics$ll, fit_pred$metrics$ll, tolerance = 1e-6)
  expect_equal(fit_group$metrics$n_params, fit_pred$metrics$n_params)
})

test_that("n_params accounting matches hand-derived formulas for each group_effects mode", {
  d <- .make_group_data()
  J <- d$J; K <- d$K; G <- 2L

  fit_none <- fit_mixture(d$X, n_classes = K, measurement = "binary",
                          group = d$grp, group_effects = "none",
                          n_steps = 1, n_init = 3, random_state = 1)
  expect_equal(fit_none$metrics$n_params, K * J + (K - 1))

  fit_prev <- fit_mixture(d$X, n_classes = K, measurement = "binary",
                          group = d$grp, group_effects = "prevalence",
                          n_steps = 1, n_init = 3, random_state = 1)
  expect_equal(fit_prev$metrics$n_params, K * J + (K - 1) * G)

  fit_meas <- fit_mixture(d$X, n_classes = K, measurement = "binary",
                          group = d$grp, group_effects = "measurement",
                          n_steps = 1, n_init = 3, random_state = 1)
  expect_equal(fit_meas$metrics$n_params, K * J * G + (K - 1))

  fit_both <- fit_mixture(d$X, n_classes = K, measurement = "binary",
                          group = d$grp, group_effects = "both",
                          n_steps = 1, n_init = 3, random_state = 1)
  expect_equal(fit_both$metrics$n_params, K * J * G + (K - 1) * G)

  fit_partial <- fit_mixture(d$X, n_classes = K, measurement = "binary",
                             group = d$grp, group_effects = "both",
                             group_invariant_items = 1:2,
                             n_steps = 1, n_init = 3, random_state = 1)
  n_inv <- 2L
  expect_equal(fit_partial$metrics$n_params,
               K * (n_inv + (J - n_inv) * G) + (K - 1) * G)
})

test_that("longitudinal_lrt() rejects both false null hypotheses with correct df", {
  d <- .make_group_data()
  K <- d$K; J <- d$J; G <- 2L

  fit_both <- fit_mixture(d$X, n_classes = K, measurement = "binary",
                          group = d$grp, group_effects = "both",
                          n_steps = 1, n_init = 10, random_state = 1)
  fit_prev <- fit_mixture(d$X, n_classes = K, measurement = "binary",
                          group = d$grp, group_effects = "prevalence",
                          n_steps = 1, n_init = 10, random_state = 1)
  fit_none <- fit_mixture(d$X, n_classes = K, measurement = "binary",
                          group = d$grp, group_effects = "none",
                          n_steps = 1, n_init = 10, random_state = 1)

  # Measurement-invariance test (sec. 5.8): item-response probabilities
  # genuinely differ by group in the generating model, so this must reject.
  inv_test <- longitudinal_lrt(fit_prev, fit_both)
  expect_equal(inv_test$df, K * J * (G - 1))
  expect_lt(inv_test$p_value, 0.001)

  # Prevalence-equivalence test (sec. 5.11): prevalences genuinely differ by
  # group too, so this must also reject.
  prev_test <- longitudinal_lrt(fit_none, fit_prev)
  expect_equal(prev_test$df, K - 1)
  expect_lt(prev_test$p_value, 0.001)
})

test_that("group_effects='both' recovers each group's true class-membership probabilities", {
  d <- .make_group_data()
  fit <- fit_mixture(d$X, n_classes = d$K, measurement = "binary",
                     group = d$grp, group_effects = "both",
                     n_steps = 1, n_init = 10, random_state = 1)

  Z      <- cbind(Intercept = 1, fit$group_info$design)
  beta   <- fit$sm$parameters$beta
  logits <- Z %*% t(beta)
  probs  <- exp(logits) / rowSums(exp(logits))

  # Match each group's fitted classes to the truth by item-response pattern
  # (class labels are not guaranteed to align across group blocks; see the
  # group_effects docs), then compare recovered prevalences to the truth.
  is_high_endorsing <- function(block) {
    which.max(rowMeans(block$parameters$pis))
  }
  block_names <- names(fit$mm$models)
  truth <- list(A = c(low = .7, high = .3), B = c(low = .3, high = .7))

  for (i in seq_along(block_names)) {
    g <- substring(block_names[i], 2, 2)  # "G1" -> "1" is the factor level index
    level <- fit$group_info$levels[as.integer(g)]
    hi <- is_high_endorsing(fit$mm$models[[i]])
    row1 <- which(as.character(d$grp) == level)[1]
    p_hi <- unname(probs[row1, hi])
    expect_equal(p_hi, unname(truth[[level]]["high"]), tolerance = 0.1)
  }
})
