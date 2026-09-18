# ==============================================================================
# transition_patterns(), and the fit indices that arrived with it
# ==============================================================================
#
# The oracles here are internal consistency, not an external target: a joint
# pattern table has to reproduce the one- and two-occasion summaries the
# package already computes when it is collapsed back down to them, and the
# model-implied and posterior tables have to agree with each other on a fit
# that describes its data well. Both are exact identities, not tolerances on
# a benchmark.

.sim_tp <- function(seed, n, Tn, J, K, sep, stay) {
  set.seed(seed)
  S <- matrix(0L, n, Tn); S[, 1] <- sample(seq_len(K), n, TRUE)
  for (t in 2:Tn)
    S[, t] <- ifelse(stats::runif(n) < stay, S[, t - 1],
                     sample(seq_len(K), n, TRUE))
  p <- matrix(0, K, J)
  for (k in seq_len(K)) p[k, ] <- (1 - sep) / 2 + sep * (k - 1) / max(K - 1, 1)
  X <- matrix(0, n, Tn * J)
  for (t in seq_len(Tn)) for (j in seq_len(J))
    X[, (t - 1) * J + j] <- stats::rbinom(n, 1, p[S[, t], j])
  X
}

.tp_fit <- function(seed = 11, n = 400, Tn = 3, K = 2)
  fit_lta(.sim_tp(seed, n, Tn, 4, K, 0.8, 0.75), n_statuses = K, times = Tn,
          measurement = "binary", n_init = 5, random_state = seed,
          standard_errors = FALSE)


test_that("every type returns a complete, normalised pattern table", {
  fit <- .tp_fit()
  N   <- sum(fit$weights_vec)

  for (ty in c("model", "posterior", "modal")) {
    tab <- transition_patterns(fit, type = ty)
    expect_equal(nrow(tab), fit$n_statuses^fit$n_times)
    expect_equal(names(tab),
                 c(fit$longitudinal$time_labels, "count", "proportion"))
    expect_equal(sum(tab$proportion), 1)
    expect_equal(sum(tab$count), N)
    expect_true(all(tab$count >= 0))
  }
})

test_that("patterns are sorted by occasion 1, then 2, and so on", {
  tab <- transition_patterns(.tp_fit(), type = "model")
  labs <- setdiff(names(tab), c("count", "proportion"))
  expect_false(is.unsorted(do.call(order, as.list(tab[labs]))))
  # The sort is the identity permutation only if the rows are already ordered.
  expect_equal(do.call(order, as.list(tab[labs])), seq_len(nrow(tab)))
})

test_that("collapsing the joint table reproduces the one-occasion summaries", {
  fit <- .tp_fit()
  labs <- fit$longitudinal$time_labels

  # Each type's marginal at occasion t is that occasion's own prevalence, on
  # the matching convention: "model" against the model-implied prevalences,
  # "posterior" against the posterior ones.
  for (ty in c("model", "posterior")) {
    tab <- transition_patterns(fit, type = ty)
    P   <- status_prevalences(fit, type = ty)
    for (t in seq_along(labs)) {
      marg <- as.numeric(tapply(tab$proportion, tab[[labs[t]]], sum))
      expect_equal(marg, as.numeric(P[t, ]))
    }
  }
})

test_that("collapsing the model table over two occasions gives the transitions", {
  fit  <- .tp_fit()
  labs <- fit$longitudinal$time_labels
  tab  <- transition_patterns(fit, type = "model")
  K    <- fit$n_statuses

  # P(S_t = i, S_t+1 = j) / P(S_t = i) is the transition matrix, exactly.
  for (t in seq_len(fit$n_times - 1L)) {
    joint <- tapply(tab$proportion,
                    list(tab[[labs[t]]], tab[[labs[t + 1L]]]), sum)
    joint <- matrix(as.numeric(joint), K, K)
    expect_equal(joint / rowSums(joint),
                 unname(transition_matrix(fit, occasion = t)))
  }
})

test_that("the model-implied and posterior tables agree on a well-fitting model", {
  # Not an identity -- they are different quantities -- but on data generated
  # from a well-separated model of this shape they must be close, and a large
  # gap is the signature of a path posterior built wrong. On the LTA-FAQ
  # benchmark the two agree to five decimals; this is the simulated
  # counterpart of that check.
  fit <- .tp_fit(n = 1200)
  m <- transition_patterns(fit, type = "model")
  p <- transition_patterns(fit, type = "posterior")
  expect_lt(max(abs(m$proportion - p$proportion)), 0.01)
})

test_that("the modal table is the Viterbi decode, cross-tabulated", {
  fit  <- .tp_fit()
  labs <- fit$longitudinal$time_labels
  tab  <- transition_patterns(fit, type = "modal")
  path <- class_assignments(fit, "viterbi")

  key  <- do.call(paste, c(as.data.frame(path), sep = "_"))
  gkey <- do.call(paste, c(as.list(tab[labs]), sep = "_"))
  expect_equal(tab$count, as.numeric(table(factor(key, levels = gkey))))
})

test_that("a mixture over chains weights the classes and refuses the posterior", {
  X <- .sim_tp(3, 500, 3, 4, 2, 0.85, 0.8)
  fit <- fit_lta(X, n_statuses = 2, times = 3, measurement = "binary",
                 n_classes = 2, n_init = 5, random_state = 3,
                 standard_errors = FALSE)
  skip_if(is.null(fit$n_classes) || fit$n_classes < 2L)

  pooled <- transition_patterns(fit, type = "model")
  per <- lapply(seq_len(fit$n_classes),
                function(k) transition_patterns(fit, type = "model", class = k))
  mix <- Reduce(`+`, Map(function(p, w) w * p$proportion,
                         per, fit$class_weights))
  expect_equal(pooled$proportion, mix)

  expect_error(transition_patterns(fit, type = "posterior"),
               "not available for a mixture over chains")
})

test_that("the K^T cap is enforced before the table is enumerated", {
  # A stub, not a fit: the point is that the refusal happens before anything
  # of size K^T is allocated, so it must be reachable without one existing.
  # K = 5, T = 20 is 95 trillion patterns; if the guard is ever moved back
  # below the grid construction this test stops erroring and starts
  # exhausting memory.
  stub <- structure(
    list(n_statuses = 5L, n_classes = 1L,
         longitudinal = list(n_times = 20L, time_labels = paste0("T", 1:20))),
    class = "lta_model")

  for (ty in c("model", "posterior", "modal"))
    expect_error(transition_patterns(stub, type = ty), "refused above 10000")
})

test_that("transition_patterns() rejects a non-LTA object", {
  expect_error(transition_patterns(structure(list(), class = "mixture_model")),
               "must be a fitted latent transition model")
})


# ------------------------------------------------------------------------------
# CAIC, AIC3 and ICL
# ------------------------------------------------------------------------------

test_that("fit_lta() reports CAIC, AIC3 and ICL on the documented formulas", {
  m <- .tp_fit()$metrics
  expect_equal(m$caic, -2 * m$ll + (log(m$n_eff) + 1) * m$n_params)
  expect_equal(m$aic3, -2 * m$ll + 3 * m$n_params)
  # ICL is BIC penalised further by classification entropy, so it can only be
  # worse -- equal exactly when the classification is perfect.
  expect_gte(m$icl, m$bic)
})

test_that("entropy_by_occasion is one number per occasion, and not the headline", {
  fit <- .tp_fit()
  eo  <- fit$metrics$entropy_by_occasion
  expect_length(eo, fit$n_times)
  expect_true(all(eo >= 0 & eo <= 1))
  expect_true(fit$metrics$entropy >= 0 && fit$metrics$entropy <= 1)
  # Summing the per-occasion marginal entropies double-counts uncertainty the
  # joint path posterior does not carry, so the headline (joint) figure is the
  # higher relative entropy of the two conventions. Documented in ?fit_lta.
  expect_gt(fit$metrics$entropy, mean(eo) - 1e-8)
})

test_that("the comparison tables carry the new columns and still choose by BIC", {
  X <- .sim_tp(21, 300, 2, 4, 2, 0.85, 0.8)
  cmp <- suppressWarnings(suppressMessages(
    compare_longitudinal(X, k_range = 2:3, times = 2, measurement = "binary",
                         n_init = 3, random_state = 21,
                         standard_errors = FALSE)))
  expect_true(all(c("CAIC", "AIC3", "ICL") %in% names(cmp$fit_table)))
  expect_equal(cmp$best_k,
               cmp$fit_table$Classes[which.min(cmp$fit_table$BIC)])
  # CAIC penalises each parameter more heavily than BIC; AIC3 more than AIC.
  expect_true(all(cmp$fit_table$CAIC > cmp$fit_table$BIC))
  expect_true(all(cmp$fit_table$AIC3 > cmp$fit_table$AIC))
})
