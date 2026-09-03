# A random-intercept LTA (Muthen & Asparouhov 2022). `.gauss_hermite()` has
# its own file, test-quadrature.R; the six tests here are 2-6 of
# internal/ROADMAP.md ### 14.10.8 (item 1, the quadrature, is in that file).
# `.lta_refine_sim()`, `.ml`, and `.lta_ri_sim()` are in helper-lta-packing.R.

# --- 2. Parameter counts, from structure alone: no estimator needed --------

test_that("RI parameter counts match the paper's Table 5 shapes (K=2,R=2,T=4)", {
  # Stationary transitions, matching the roadmap's "stationary, structural = 3"
  # shape these targets are computed against -- without it T = 4's default
  # (a separate matrix per interval) adds parameters the target does not count.
  X <- .lta_refine_sim(n = 30, K = 2, Tn = 4, J = 2, seed = 1)
  fit_reg <- suppressWarnings(fit_lta(X, n_statuses = 2, times = 4,
    measurement = "binary", transition_invariance = "full",
    n_init = 1, max_iter = 1, random_state = 1, standard_errors = FALSE))
  fit_bin <- suppressWarnings(fit_lta(X, n_statuses = 2, times = 4,
    measurement = "binary", transition_invariance = "full",
    n_init = 1, max_iter = 1, random_state = 1, standard_errors = FALSE,
    random_intercept = "binary", n_ri = 2))
  fit_con <- suppressWarnings(fit_lta(X, n_statuses = 2, times = 4,
    measurement = "binary", transition_invariance = "full",
    n_init = 1, max_iter = 1, random_state = 1, standard_errors = FALSE,
    random_intercept = "continuous", n_quadrature = 5))
  expect_equal(fit_reg$n_params, 7L)
  expect_equal(fit_bin$n_params, 10L)
  expect_equal(fit_con$n_params, 9L)
})

test_that("RI parameter counts match the LTA-FAQ shape (K=4,R=5,T=2)", {
  X <- .lta_refine_sim(n = 30, K = 4, Tn = 2, J = 5, seed = 1)
  fit_reg <- suppressWarnings(fit_lta(X, n_statuses = 4, times = 2,
    measurement = "binary", n_init = 1, max_iter = 1, random_state = 1,
    standard_errors = FALSE))
  fit_con <- suppressWarnings(fit_lta(X, n_statuses = 4, times = 2,
    measurement = "binary", n_init = 1, max_iter = 1, random_state = 1,
    standard_errors = FALSE, random_intercept = "continuous",
    n_quadrature = 20))
  fit_bin <- suppressWarnings(fit_lta(X, n_statuses = 4, times = 2,
    measurement = "binary", n_init = 1, max_iter = 1, random_state = 1,
    standard_errors = FALSE, random_intercept = "binary", n_ri = 2))
  expect_equal(fit_reg$n_params, 35L)
  expect_equal(fit_con$n_params, 40L)
  expect_equal(fit_bin$n_params, 41L)
})

# --- 3. The refusals --------------------------------------------------------

test_that("`random_intercept` refuses partial/no invariance", {
  X <- .lta_refine_sim(n = 30, K = 2, Tn = 3, J = 4, seed = 1)
  expect_error(
    fit_lta(X, n_statuses = 2, times = 3, measurement = "binary",
           measurement_invariance = "none", random_intercept = "continuous"),
    "measurement_invariance")
})

test_that("`random_intercept` refuses non-binary measurement", {
  X <- .lta_gaussian_refine_sim()
  expect_error(
    fit_lta(X, n_statuses = 3, times = 2, measurement = "continuous",
           random_intercept = "continuous"),
    "binary indicators only")
})

# Mover-Stayer x RI-LTA (roadmap Part 14 Phase B, "a double outer loop over
# class and node, and no new mathematics"): the factor (`A`/`L`/`mass`) is
# shared across classes, only `delta`/`tau`/`class_weights` are class-
# specific, mirroring how the measurement model is already pooled across
# classes for a non-RI mover-stayer fit. This is a cheap smoke test only --
# n_init = 1, max_iter tiny, made-up data -- not a match against real-data
# reference targets, which stays a separate, more expensive follow-up.
test_that("mover-stayer combines with a random intercept (continuous)", {
  X <- .lta_refine_sim(n = 30, K = 2, Tn = 4, J = 3, seed = 1)
  fit0 <- suppressWarnings(fit_lta(X, n_statuses = 2, times = 4,
    measurement = "binary", mover_stayer = TRUE,
    n_init = 1, max_iter = 3, random_state = 1, standard_errors = FALSE))
  fit_ri <- suppressWarnings(fit_lta(X, n_statuses = 2, times = 4,
    measurement = "binary", mover_stayer = TRUE,
    random_intercept = "continuous", n_quadrature = 3,
    n_init = 1, max_iter = 3, random_state = 1, standard_errors = FALSE))
  expect_true(is.finite(fit_ri$loglik))
  # The only new parameters over the plain mover-stayer fit are the R
  # loadings (M = 1 for the continuous variant); everything else is either
  # shared (the factor's own A, folded into the already-counted pis) or
  # already scales with C in n_params (see .lta_n_parameters()).
  expect_equal(fit_ri$n_params, fit0$n_params + 3L)
  expect_equal(length(fit_ri$class_weights), 2L)
  expect_equal(dim(fit_ri$ri$L), c(3L, 1L))
})

test_that("mover-stayer combines with a random intercept (binary)", {
  X <- .lta_refine_sim(n = 30, K = 2, Tn = 4, J = 3, seed = 1)
  fit0 <- suppressWarnings(fit_lta(X, n_statuses = 2, times = 4,
    measurement = "binary", mover_stayer = TRUE,
    n_init = 1, max_iter = 3, random_state = 1, standard_errors = FALSE))
  fit_ri <- suppressWarnings(fit_lta(X, n_statuses = 2, times = 4,
    measurement = "binary", mover_stayer = TRUE,
    random_intercept = "binary", n_ri = 2,
    n_init = 1, max_iter = 3, random_state = 1, standard_errors = FALSE))
  expect_true(is.finite(fit_ri$loglik))
  # For n_ri = 2, the binary variant adds R loadings plus one free mass
  # (length(mass) - 1) over the plain mover-stayer fit.
  expect_equal(fit_ri$n_params, fit0$n_params + 3L + 1L)
  expect_equal(length(fit_ri$class_weights), 2L)
})

test_that("tie_initial_status ties the initial distribution across classes", {
  X <- .lta_refine_sim(n = 30, K = 2, Tn = 4, J = 3, seed = 1)
  fit0 <- suppressWarnings(fit_lta(X, n_statuses = 2, times = 4,
    measurement = "binary", mover_stayer = TRUE,
    n_init = 1, max_iter = 3, random_state = 1, standard_errors = FALSE))
  fit_tied <- suppressWarnings(fit_lta(X, n_statuses = 2, times = 4,
    measurement = "binary", mover_stayer = TRUE, tie_initial_status = TRUE,
    n_init = 1, max_iter = 3, random_state = 1, standard_errors = FALSE))
  expect_true(is.finite(fit_tied$loglik))
  # Tying collapses C free (K-1)-vectors into one shared one.
  expect_equal(fit_tied$n_params, fit0$n_params - 1L)
  expect_equal(fit_tied$delta_c[[1]], fit_tied$delta_c[[2]])
  expect_true(.lta_scores_full(fit_tied))
})

test_that("the refinement's gradient matches finite differences for a tied mover-stayer fit", {
  # A tied initial-status distribution collapses the C per-class delta blocks
  # in .lta_score_matrix()'s C > 1 branch into one shared block (the roadmap's
  # 14.13 derivation) -- this is the same check test-lta-refine.R runs for an
  # untied mixture over chains, on a model where that collapse is exercised.
  X <- .lta_refine_sim(n = 30, K = 2, Tn = 4, J = 3, seed = 1)
  fit <- suppressMessages(suppressWarnings(
    fit_lta(X, n_statuses = 2, times = 4, measurement = "binary",
            mover_stayer = TRUE, tie_initial_status = TRUE,
            smoothing = 0, bayes_constants = .ml, n_init = 2,
            random_state = 1, refine = FALSE, standard_errors = FALSE)))
  expect_true(.lta_scores_full(fit))

  layout <- .lta_par_layout(fit)
  par0   <- .lta_par_pack(fit, layout)
  w      <- fit$weights_vec

  expect_equal(length(par0), ncol(.lta_score_matrix(fit, X)$S))
  expect_equal(length(par0), fit$n_params)

  rt <- .lta_par_unpack(par0, fit, layout)
  expect_equal(rt$delta_c[[1]], fit$delta_c[[1]], tolerance = 1e-10)
  expect_equal(rt$delta_c[[2]], fit$delta_c[[2]], tolerance = 1e-10)
  expect_equal(rt$delta_c[[1]], rt$delta_c[[2]], tolerance = 1e-10)
  expect_equal(rt$tau_c, fit$tau_c, tolerance = 1e-10)
  expect_equal(rt$class_weights, fit$class_weights, tolerance = 1e-10)

  obj <- function(p) {
    st <- .lta_par_unpack(p, fit, layout)
    sum(w * .lta_score_matrix(st, X)$ll) +
      .lta_penalty(st, X, layout, 0)$value
  }
  ana <- colSums(sweep(.lta_score_matrix(fit, X)$S, 1, w, "*")) +
    .lta_penalty(fit, X, layout, 0)$gradient
  eps <- 1e-5
  fd <- vapply(seq_along(par0), function(i) {
    e <- numeric(length(par0)); e[i] <- eps
    (obj(par0 + e) - obj(par0 - e)) / (2 * eps)
  }, numeric(1))

  expect_lt(max(abs(ana - fd)) / max(1, max(abs(fd))), 1e-5)

  # The penalty check: at a tightly converged fit with the default priors on,
  # the penalised gradient (EM score + .lta_penalty()) should vanish. This is
  # the one that catches a wrong tied prior constant in .lta_penalty() -- with
  # alpha/(C*K) in place of alpha/K the delta component alone fails to vanish.
  fit2 <- suppressMessages(suppressWarnings(
    fit_lta(X, n_statuses = 2, times = 4, measurement = "binary",
            mover_stayer = TRUE, tie_initial_status = TRUE,
            n_init = 2, random_state = 1, tol = 1e-14, max_iter = 20000,
            refine = FALSE, standard_errors = FALSE)))
  layout2 <- .lta_par_layout(fit2)
  g <- colSums(sweep(.lta_score_matrix(fit2, X)$S, 1, fit2$weights_vec, "*")) +
    .lta_penalty(fit2, X, layout2, 1)$gradient
  expect_lt(max(abs(g)), 1e-2)
})

test_that("`random_intercept` still refuses covariate-driven classes", {
  X <- .lta_refine_sim(n = 30, K = 2, Tn = 3, J = 4, seed = 1)
  expect_error(
    fit_lta(X, n_statuses = 2, times = 3, measurement = "binary",
           n_classes = 2, predictors_initial = matrix(rnorm(30), 30, 1),
           random_intercept = "continuous"),
    "not yet available")
})

test_that("`random_intercept` refuses covariates", {
  sim <- .lta_cov_refine_sim()
  expect_error(
    fit_lta(sim$X, n_statuses = 2, times = 3, measurement = "binary",
           predictors_initial = sim$Z, random_intercept = "binary"),
    "covariates")
})

test_that("`n_quadrature`/`n_ri` are validated", {
  X <- .lta_refine_sim(n = 30, K = 2, Tn = 3, J = 4, seed = 1)
  expect_error(
    fit_lta(X, n_statuses = 2, times = 3, measurement = "binary",
           random_intercept = "continuous", n_quadrature = 0),
    "n_quadrature")
  expect_error(
    fit_lta(X, n_statuses = 2, times = 3, measurement = "binary",
           random_intercept = "binary", n_ri = 1),
    "n_ri")
})

# --- 4. Q = 1 reduces exactly to regular LTA --------------------------------

test_that("n_quadrature = 1 reproduces regular LTA exactly", {
  # A single deterministic start is enough here: this is an algebraic
  # identity (a zero node contributes nothing regardless of the loadings it
  # is multiplied by), not a statistical recovery claim, so it holds at
  # whatever point n_init = 1 happens to converge to, not only at the global
  # optimum -- confirmed at n = 100 to ~1e-9, well inside the asserted bound.
  X <- .lta_refine_sim(n = 100, K = 2, Tn = 3, J = 4, seed = 4)
  fit_reg <- suppressWarnings(fit_lta(X, n_statuses = 2, times = 3,
                     measurement = "binary",
                     smoothing = 0, bayes_constants = .ml, n_init = 1,
                     random_state = 3, tol = 1e-12, max_iter = 5000,
                     standard_errors = FALSE))
  fit_ri1 <- suppressWarnings(fit_lta(X, n_statuses = 2, times = 3,
                     measurement = "binary",
                     smoothing = 0, bayes_constants = .ml, n_init = 1,
                     random_state = 3, tol = 1e-12, max_iter = 5000,
                     standard_errors = FALSE,
                     random_intercept = "continuous", n_quadrature = 1))
  expect_lt(abs(fit_ri1$loglik - fit_reg$loglik), 1e-6)
  expect_lt(max(abs(fit_ri1$tau[[1]] - fit_reg$tau[[1]])), 1e-5)
})

# --- 5. Nesting check: RI fitted to regular-LTA data must find loadings ~0 -
#
# A spurious loading always improves the fit somewhat on finite data -- an
# extra free parameter fits sampling noise even where none of the true
# structure calls for it -- so this needs enough cases that the noise the
# loadings could chase is small next to the signal. Measured on this exact
# fixture (seed 4, n_init = 5, refine_from the converged regular fit as the
# RI search's own start -- 14.10.7's accelerator, not a relaxation of the
# check): loadings max 0.018/0.10/0.18/0.16 and an LL gap of 0.60 at n = 3000,
# against 0.65 max loading and a 2.5 LL gap at n = 250. That shrinkage with n
# and not with a better start is what distinguishes finite-sample noise from
# the bug this test exists to catch (14.10.10 failure mode 1: a wrong
# posterior would not shrink with n at all).

test_that("RI-LTA fitted to regular-LTA data recovers near-zero loadings", {
  skip_on_cran()
  X <- .lta_refine_sim(n = 3000, K = 2, Tn = 3, J = 4, seed = 4)
  fit_reg <- suppressWarnings(fit_lta(X, n_statuses = 2, times = 3,
                    measurement = "binary",
                    smoothing = 0, bayes_constants = .ml, n_init = 5,
                    random_state = 3, tol = 1e-11, max_iter = 3000,
                    standard_errors = FALSE))
  fit_ri <- suppressWarnings(fit_lta(X, n_statuses = 2, times = 3,
                    measurement = "binary",
                    smoothing = 0, bayes_constants = .ml,
                    refine_from = fit_reg, tol = 1e-11, max_iter = 3000,
                    standard_errors = FALSE,
                    random_intercept = "continuous", n_quadrature = 15))
  expect_lt(max(abs(fit_ri$ri$L)), 0.25)
  expect_lt(abs(fit_ri$loglik - fit_reg$loglik), 1)
})

# --- 6. Recovery check: one replication of the paper's own MC design -------

test_that("RI-LTA recovers loadings and transitions on its own generator", {
  skip_on_cran()
  sim <- .lta_ri_sim()
  fit <- suppressWarnings(fit_lta(sim$X, n_statuses = 2, times = 3,
                measurement = "binary",
                smoothing = 0, bayes_constants = .ml,
                random_intercept = "continuous", n_quadrature = 15,
                n_init = 1, random_state = 1, standard_errors = FALSE))
  expect_lt(mean(abs(fit$ri$L - 2)), 0.3)
  trans11 <- max(fit$tau[[1]][1, 1], fit$tau[[1]][2, 2])
  expect_lt(abs(trans11 - 0.622), 0.08)
})

# --- 7. The staged search and its coarse ranking grid ----------------------
#
# An RI fit stages its restarts (fit_lta(), `staged <- C > 1L || !is.null(ri)`)
# and ranks them on five quadrature nodes before promoting the survivors to
# `n_quadrature`. Both are search devices, so both are graded on recovery of
# the generator's truth rather than on a log-likelihood one of them optimises.
# Test 6 above already runs the ladder at `n_quadrature = 15`; what is left to
# assert is that the promotion actually happens, and that ranking a real pool
# of restarts on the coarse grid still lands in the right basin.

test_that("the coarse ranking grid never reaches the reported fit", {
  # The ladder is on by default at 5 nodes; set it explicitly anyway so this
  # test's intent (exercising the mechanism it guards) does not silently stop
  # doing so if the default ever changes again.
  old <- getOption("mixtureEM.ri_rank_nodes")
  on.exit(options(mixtureEM.ri_rank_nodes = old), add = TRUE)
  options(mixtureEM.ri_rank_nodes = 5L)
  X <- .lta_refine_sim(n = 200, K = 2, Tn = 3, J = 4, seed = 4)
  fit <- suppressWarnings(fit_lta(X, n_statuses = 2, times = 3,
                measurement = "binary", random_intercept = "continuous",
                n_quadrature = 20, n_init = 2, max_iter = 30,
                random_state = 1, standard_errors = FALSE))
  # 5 here rather than 20 is the whole failure this guards: it would mean the
  # survivors were run on, and the model reported, at the ranking accuracy.
  expect_equal(length(fit$ri$mass), 20L)
  expect_equal(nrow(fit$ri$Dnode), 20L)
})

test_that("a staged multi-start RI search recovers the generator's truth", {
  skip_on_cran()
  sim <- .lta_ri_sim()
  fit <- suppressWarnings(fit_lta(sim$X, n_statuses = 2, times = 3,
                measurement = "binary",
                smoothing = 0, bayes_constants = .ml,
                random_intercept = "continuous", n_quadrature = 15,
                n_init = 6, random_state = 1, standard_errors = FALSE))
  expect_lt(mean(abs(fit$ri$L - 2)), 0.3)
  trans11 <- max(fit$tau[[1]][1, 1], fit$tau[[1]][2, 2])
  expect_lt(abs(trans11 - 0.622), 0.08)
})

# --- 8. The status ordering has to carry the random intercept with it -------

# `order_by_size` relabels the statuses at the very end of fit_lta(). For a
# regular LTA the measurement model it has to carry is `pis`; for an RI fit
# `pis` is only the integral of `A` over the nodes, and `A` is the thing being
# estimated. Permuting one and not the other leaves the returned object
# describing two different models, with the log-likelihood it reports (computed
# before the relabelling) belonging to neither. That is invisible in every
# summary the fit prints, so it is asserted here directly: score the fit's own
# stored parameters and see whether the number comes back.
test_that("order_by_size relabels the random intercept's intercepts too", {
  X <- .lta_refine_sim(n = 150, K = 3, Tn = 2, J = 4, seed = 1)
  unsorted <- suppressWarnings(fit_lta(X, n_statuses = 3, times = 2,
                measurement = "binary", random_intercept = "continuous",
                n_quadrature = 5, n_init = 1, max_iter = 25, random_state = 1,
                standard_errors = FALSE, order_by_size = FALSE))
  # This fixture's Time-1 prevalences are 0.23/0.15/0.62, so the relabelling
  # below is a real permutation and not the identity.
  expect_false(identical(order(unsorted$delta, decreasing = TRUE), 1:3))

  sorted <- .sort_lta_statuses(unsorted)
  expect_equal(.lta_em(sorted, X, max_iter = 0L, alpha = 1)$loglik,
               unsorted$loglik, tolerance = 1e-10)
  # `pis` is derived from `A`, so the two must still agree afterwards.
  expect_equal(sorted$mm$models[[1]]$parameters$pis,
               .lta_ri_integrated_pis(sorted$ri, 3L, 4L),
               tolerance = 1e-12, ignore_attr = TRUE)
})

# --- 9. Standard errors on the loadings -------------------------------------

# The random intercept's loading is the headline estimate the source paper
# reports with a standard error; W10 adds it. The L-BFGS polish and the robust
# sandwich stay off for RI (`.lta_scores_full()`), so these only check the
# default (empirical-information) estimator.
test_that("standard_errors = TRUE returns finite loading SEs (continuous)", {
  X <- .lta_refine_sim(n = 150, K = 2, Tn = 4, J = 3, seed = 1)
  fit <- suppressWarnings(fit_lta(X, n_statuses = 2, times = 4,
    measurement = "binary", random_intercept = "continuous",
    n_quadrature = 10, n_init = 1, max_iter = 25, random_state = 1,
    standard_errors = TRUE))
  expect_true(.lta_scores_supported(fit))
  expect_false(.lta_scores_full(fit))
  expect_false(is.null(fit$se))
  expect_equal(dim(fit$se$loading_se), c(3L, 1L))
  expect_true(all(is.finite(fit$se$loading_se)))
  expect_true(all(fit$se$loading_se > 0))
  b <- fit$se$blocks[[which(vapply(fit$se$blocks, function(x) x$name, "") ==
                             "alpha[item 1]")]]
  alpha_se <- sqrt(diag(fit$se$vcov))[b$cols]
  expect_true(all(is.finite(alpha_se)) && all(alpha_se > 0))
})

test_that("standard_errors = TRUE returns finite loading SEs (binary)", {
  X <- .lta_refine_sim(n = 150, K = 2, Tn = 4, J = 3, seed = 1)
  fit <- suppressWarnings(fit_lta(X, n_statuses = 2, times = 4,
    measurement = "binary", random_intercept = "binary", n_ri = 2,
    n_init = 1, max_iter = 25, random_state = 1, standard_errors = TRUE))
  expect_false(is.null(fit$se))
  expect_equal(dim(fit$se$loading_se), c(3L, 1L))
  expect_true(all(is.finite(fit$se$loading_se)) && all(fit$se$loading_se > 0))
})

test_that('standard_errors = "robust" falls back silently for RI fits', {
  X <- .lta_refine_sim(n = 60, K = 2, Tn = 4, J = 2, seed = 1)
  fit <- suppressWarnings(fit_lta(X, n_statuses = 2, times = 4,
    measurement = "binary", random_intercept = "continuous", n_quadrature = 5,
    n_init = 1, max_iter = 10, random_state = 1, standard_errors = "robust"))
  expect_false(is.null(fit$se))
  expect_false(fit$se$robust)
})

# --- Mover-stayer x RI standard errors (roadmap ### 14.14) ------------------

test_that("standard errors exist for a mover-stayer RI fit", {
  X <- .lta_refine_sim(n = 60, K = 2, Tn = 4, J = 3, seed = 1)
  fit_con <- suppressWarnings(fit_lta(X, n_statuses = 2, times = 4,
    measurement = "binary", mover_stayer = TRUE,
    random_intercept = "continuous", n_quadrature = 5,
    n_init = 1, max_iter = 25, random_state = 1, standard_errors = TRUE))
  expect_true(.lta_scores_supported(fit_con))
  expect_false(.lta_scores_full(fit_con))
  expect_false(is.null(fit_con$se))
  expect_true(all(is.finite(fit_con$se$loading_se)))
  expect_true(all(fit_con$se$loading_se > 0))
  expect_true(all(is.finite(fit_con$se$prob_se[["class"]])))
  expect_equal(length(fit_con$se$prob_se[["class"]]), 2L)

  fit_bin <- suppressWarnings(fit_lta(X, n_statuses = 2, times = 4,
    measurement = "binary", mover_stayer = TRUE,
    random_intercept = "binary", n_ri = 2,
    n_init = 1, max_iter = 25, random_state = 1, standard_errors = TRUE))
  expect_false(is.null(fit_bin$se))
  expect_true(all(is.finite(fit_bin$se$loading_se)))
  expect_true(all(is.finite(fit_bin$se$prob_se[["ri_mass"]])))
  expect_true(all(fit_bin$se$prob_se[["ri_mass"]] > 0))
  expect_equal(length(fit_bin$se$prob_se[["ri_mass"]]), 2L)
})

test_that("ncol(S) == n_params for all four RI shapes", {
  X <- .lta_refine_sim(n = 60, K = 2, Tn = 4, J = 3, seed = 1)
  fit_con <- suppressWarnings(fit_lta(X, n_statuses = 2, times = 4,
    measurement = "binary", random_intercept = "continuous", n_quadrature = 5,
    n_init = 1, max_iter = 25, random_state = 1, standard_errors = FALSE))
  fit_bin <- suppressWarnings(fit_lta(X, n_statuses = 2, times = 4,
    measurement = "binary", random_intercept = "binary", n_ri = 2,
    n_init = 1, max_iter = 25, random_state = 1, standard_errors = FALSE))
  fit_con_ms <- suppressWarnings(fit_lta(X, n_statuses = 2, times = 4,
    measurement = "binary", mover_stayer = TRUE,
    random_intercept = "continuous", n_quadrature = 5,
    n_init = 1, max_iter = 25, random_state = 1, standard_errors = FALSE))
  fit_bin_ms <- suppressWarnings(fit_lta(X, n_statuses = 2, times = 4,
    measurement = "binary", mover_stayer = TRUE,
    random_intercept = "binary", n_ri = 2,
    n_init = 1, max_iter = 25, random_state = 1, standard_errors = FALSE))
  for (fit in list(fit_con, fit_bin, fit_con_ms, fit_bin_ms))
    expect_equal(ncol(.lta_score_matrix(fit, X)$S), fit$n_params)
})

test_that("the RI score blocks match finite differences (single-class and mover-stayer)", {
  check_ri_gradient <- function(fit, X) {
    sc    <- .lta_score_matrix(fit, X)
    S     <- sc$S
    w     <- fit$weights_vec
    cols  <- function(name)
      Find(function(b) identical(b$name, name), sc$blocks)$cols
    ll_at <- function(st) sum(st$weights_vec * .lta_score_matrix(st, X)$ll)
    eps   <- 1e-5
    fd    <- function(f) (ll_at(f(eps)) - ll_at(f(-eps))) / (2 * eps)
    bump_A <- function(k) function(h) { st <- fit; st$ri$A[k, 1] <- st$ri$A[k, 1] + h; st }
    bump_L <- function(m) function(h) { st <- fit; st$ri$L[1, m] <- st$ri$L[1, m] + h; st }
    bump_m <- function(q) function(h) {
      st <- fit; Q <- length(st$ri$mass)
      v <- log(st$ri$mass[-Q]) - log(st$ri$mass[Q]); v[q] <- v[q] + h
      p <- exp(c(v, 0) - max(c(v, 0))); st$ri$mass <- p / sum(p); st
    }
    for (k in seq_len(nrow(fit$ri$A)))
      expect_lt(abs(sum(w * S[, cols("alpha[item 1]")[k]]) - fd(bump_A(k))) /
                  max(1, abs(fd(bump_A(k)))), 1e-5)
    for (m in seq_len(ncol(fit$ri$L)))
      expect_lt(abs(sum(w * S[, cols("lambda[item 1]")[m]]) - fd(bump_L(m))) /
                  max(1, abs(fd(bump_L(m)))), 1e-5)
    if (identical(fit$ri$kind, "binary") && length(fit$ri$mass) > 1L)
      for (q in seq_len(length(fit$ri$mass) - 1L))
        expect_lt(abs(sum(w * S[, cols("ri_mass")[q]]) - fd(bump_m(q))) /
                    max(1, abs(fd(bump_m(q)))), 1e-5)
  }

  X <- .lta_refine_sim(n = 60, K = 2, Tn = 4, J = 3, seed = 1)
  fit_bin <- suppressWarnings(fit_lta(X, n_statuses = 2, times = 4,
    measurement = "binary", random_intercept = "binary", n_ri = 2,
    n_init = 1, max_iter = 25, random_state = 1, standard_errors = FALSE))
  check_ri_gradient(fit_bin, X)

  fit_bin_ms <- suppressWarnings(fit_lta(X, n_statuses = 2, times = 4,
    measurement = "binary", mover_stayer = TRUE,
    random_intercept = "binary", n_ri = 2,
    n_init = 1, max_iter = 25, random_state = 1, standard_errors = FALSE))
  check_ri_gradient(fit_bin_ms, X)
})
