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

test_that("`random_intercept` refuses continuous measurement", {
  X <- .lta_gaussian_refine_sim()
  expect_error(
    fit_lta(X, n_statuses = 3, times = 2, measurement = "continuous",
           random_intercept = "continuous"),
    "binary or ordinal indicators only")
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

test_that("`random_intercept` combines with covariates on the initial status and transitions", {
  # Roadmap ### 14.17: the guard that used to refuse this combination is
  # lifted; this is the structural replacement for the old refusal test.
  sim <- .lta_cov_refine_sim()
  fit <- suppressWarnings(fit_lta(sim$X, n_statuses = 2, times = 3,
    measurement = "binary", predictors_initial = sim$Z,
    predictors_transition = sim$Z, random_intercept = "binary", n_ri = 2,
    n_init = 1, max_iter = 5, random_state = 1, standard_errors = FALSE))
  expect_true(is.finite(fit$loglik))
  expect_false(is.null(fit$delta_beta))
  expect_false(is.null(fit$tau_beta))
  expect_false(is.null(fit$ri))
})

test_that("`random_intercept` combines with `group` (multiple-group RI-LTA)", {
  # `group` is implemented as covariates on delta/tau (R/lta.R), so this is
  # the same guard as the test above, exercised through the other interface.
  sim <- .lta_cov_refine_sim()
  g <- factor(ifelse(sim$Z$z > 0, "a", "b"))
  fit <- suppressWarnings(fit_lta(sim$X, n_statuses = 2, times = 3,
    measurement = "binary", group = g, group_effects = "both",
    random_intercept = "continuous", n_quadrature = 5,
    n_init = 1, max_iter = 5, random_state = 1, standard_errors = FALSE))
  expect_true(is.finite(fit$loglik))
  expect_false(is.null(fit$ri))
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
  # This is an algebraic identity (a zero node contributes nothing
  # regardless of the loadings it is multiplied by), not a statistical
  # recovery claim, but it is an identity between two DIFFERENT
  # parametrizations of the same surface, so a single restart of each can
  # still land in different basins of it. n_init = 1 confirmed this at
  # ~1e-9 back when the RI restart pool's second half was seeded from a
  # plain-LTA pre-fit, which happened to bias that one restart toward the
  # same basin fit_reg finds; with an unbiased second construction
  # (Part 47's W4) that bias is gone and n_init = 1 lands ~7.8 away. A
  # handful of restarts finds the shared optimum exactly again (confirmed
  # at n_init = 3, 5 and 10, all agreeing to ~1e-9), which is what this test
  # asks for.
  X <- .lta_refine_sim(n = 100, K = 2, Tn = 3, J = 4, seed = 4)
  fit_reg <- suppressWarnings(fit_lta(X, n_statuses = 2, times = 3,
                     measurement = "binary",
                     smoothing = 0, bayes_constants = .ml, n_init = 1,
                     random_state = 3, tol = 1e-12, max_iter = 5000,
                     standard_errors = FALSE))
  fit_ri1 <- suppressWarnings(fit_lta(X, n_statuses = 2, times = 3,
                     measurement = "binary",
                     smoothing = 0, bayes_constants = .ml, n_init = 3,
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

# How many stage-1 winners get promoted to the full grid. `fit$metrics$n_starts`
# is `length(final_lls)`, which on the staged path is populated only by the
# promoted-survivor loop (R/lta.R), so it doubles as an observable count of
# `n_survivors` without reaching into an internal.
test_that("the survivor count is proportional, floored at 3, and overridable", {
  X <- .lta_refine_sim(n = 120, K = 2, Tn = 2, J = 3, seed = 2)
  fit_default <- suppressWarnings(fit_lta(X, n_statuses = 2, times = 2,
                measurement = "binary", random_intercept = "continuous",
                n_quadrature = 5, n_init = 10, max_iter = 5,
                random_state = 1, standard_errors = FALSE))
  # 10% of 10 rounds up to 1, below the floor of 3.
  expect_equal(fit_default$metrics$n_starts, 3L)

  old <- getOption("mixtureEM.lta_survivors")
  on.exit(options(mixtureEM.lta_survivors = old), add = TRUE)
  options(mixtureEM.lta_survivors = 7L)
  fit_opt <- suppressWarnings(fit_lta(X, n_statuses = 2, times = 2,
                measurement = "binary", random_intercept = "continuous",
                n_quadrature = 5, n_init = 10, max_iter = 5,
                random_state = 1, standard_errors = FALSE))
  expect_equal(fit_opt$metrics$n_starts, 7L)
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
# reports with a standard error; W10 adds it. The L-BFGS polish now covers RI
# fits too (`.lta_scores_full()`, `### 14.15` W8); the robust sandwich also
# does (`### 14.15` W6), so the tests below check the default
# (empirical-information) estimator except where noted.
test_that("standard_errors = TRUE returns finite loading SEs (continuous)", {
  X <- .lta_refine_sim(n = 150, K = 2, Tn = 4, J = 3, seed = 1)
  fit <- suppressWarnings(fit_lta(X, n_statuses = 2, times = 4,
    measurement = "binary", random_intercept = "continuous",
    n_quadrature = 10, n_init = 1, max_iter = 25, random_state = 1,
    standard_errors = TRUE))
  expect_true(.lta_scores_supported(fit))
  expect_true(.lta_scores_full(fit))
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

test_that('standard_errors = "robust" now works for RI fits (roadmap ### 14.15)', {
  X <- .lta_refine_sim(n = 60, K = 2, Tn = 4, J = 2, seed = 1)
  fit <- suppressWarnings(fit_lta(X, n_statuses = 2, times = 4,
    measurement = "binary", random_intercept = "continuous", n_quadrature = 5,
    n_init = 1, max_iter = 10, random_state = 1, standard_errors = "robust"))
  expect_false(is.null(fit$se))
  expect_true(fit$se$robust)
})

# --- Mover-stayer x RI standard errors (roadmap ### 14.14) ------------------

test_that("standard errors exist for a mover-stayer RI fit", {
  X <- .lta_refine_sim(n = 60, K = 2, Tn = 4, J = 3, seed = 1)
  fit_con <- suppressWarnings(fit_lta(X, n_statuses = 2, times = 4,
    measurement = "binary", mover_stayer = TRUE,
    random_intercept = "continuous", n_quadrature = 5,
    n_init = 1, max_iter = 25, random_state = 1, standard_errors = TRUE))
  expect_true(.lta_scores_supported(fit_con))
  expect_true(.lta_scores_full(fit_con))
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

test_that("the RI penalty branch matches finite differences on the RI log-prior (binary and continuous)", {
  # ### 14.15 W8 step 1: .lta_penalty() grew alpha/lambda/ri_mass branches
  # mirroring .lta_ri_log_prior() term for term. This checks the whole
  # penalised gradient -- delta/tau/alpha/lambda/ri_mass together -- against
  # central differences on .lta_log_prior(), which for an RI fit dispatches to
  # .lta_ri_log_prior() (R/lta_ri.R:294-316).
  check_ri_penalty <- function(fit, X) {
    layout <- .lta_par_layout(fit)
    par0   <- .lta_par_pack(fit, layout)
    obj <- function(p) .lta_log_prior(.lta_par_unpack(p, fit, layout), X, 1)
    ana <- .lta_penalty(fit, X, layout, 1)$gradient
    eps <- 1e-6
    fd <- vapply(seq_along(par0), function(i) {
      e <- numeric(length(par0)); e[i] <- eps
      (obj(par0 + e) - obj(par0 - e)) / (2 * eps)
    }, numeric(1))
    expect_lt(max(abs(ana - fd)) / max(1, max(abs(fd))), 1e-6)
  }

  X <- .lta_refine_sim(n = 60, K = 2, Tn = 4, J = 3, seed = 1)
  # Binary RI: exercises the ri_mass block as well as alpha/lambda.
  fit_bin <- suppressWarnings(fit_lta(X, n_statuses = 2, times = 4,
    measurement = "binary", random_intercept = "binary", n_ri = 2,
    n_init = 1, max_iter = 25, random_state = 1, standard_errors = FALSE))
  check_ri_penalty(fit_bin, X)

  # Continuous RI: node weights are fixed Gauss-Hermite quadrature, so only
  # alpha/lambda get a block -- no ri_mass.
  fit_con <- suppressWarnings(fit_lta(X, n_statuses = 2, times = 4,
    measurement = "binary", random_intercept = "continuous", n_quadrature = 5,
    n_init = 1, max_iter = 25, random_state = 1, standard_errors = FALSE))
  check_ri_penalty(fit_con, X)
})

# --- 10. The packing trio (roadmap ### 14.15 W2/W3/W5) ----------------------
#
# .lta_par_layout()/.lta_par_pack()/.lta_par_unpack() now have alpha/lambda/
# ri_mass cases, and .lta_ll_case() dispatches an RI fit to .lta_ri_ll_case().
# These four assertions are what would have caught a mis-packed vector: it
# still round-trips and still produces plausible, wrong standard errors, which
# is how the 32.2 defect survived its own unit tests (roadmap ### 32.2).

test_that("the packed vector describes the same thing the score matrix does, on all four RI shapes", {
  X <- .lta_refine_sim(n = 60, K = 2, Tn = 4, J = 3, seed = 1)
  shapes <- list(
    con    = list(random_intercept = "continuous", n_quadrature = 5),
    bin    = list(random_intercept = "binary", n_ri = 2),
    con_ms = list(random_intercept = "continuous", n_quadrature = 5, mover_stayer = TRUE),
    bin_ms = list(random_intercept = "binary", n_ri = 2, mover_stayer = TRUE))

  for (shape in shapes) {
    fit <- suppressWarnings(do.call(fit_lta, c(list(
      X, n_statuses = 2, times = 4, measurement = "binary",
      n_init = 1, max_iter = 25, random_state = 1, standard_errors = FALSE),
      shape)))

    # Assertion 1: layout and score-matrix block widths agree.
    layout <- mixtureEM:::.lta_par_layout(fit)
    sc     <- mixtureEM:::.lta_score_matrix(fit, X)
    expect_equal(vapply(layout, function(b) b$len, integer(1)),
                 vapply(sc$blocks, function(b) length(b$cols), integer(1)))
    par <- mixtureEM:::.lta_par_pack(fit, layout)
    expect_equal(length(par), ncol(sc$S))
    expect_equal(ncol(sc$S), fit$n_params)

    # Assertion 2: pack -> unpack -> pack round-trips, and pis is restored.
    st2  <- mixtureEM:::.lta_par_unpack(par, fit, layout)
    par2 <- mixtureEM:::.lta_par_pack(st2, layout)
    expect_equal(par2, par, tolerance = 1e-12)
    expect_equal(st2$mm$models[[1]]$parameters$pis,
                 mixtureEM:::.lta_ri_integrated_pis(st2$ri, 2L, 3L),
                 tolerance = 1e-12, ignore_attr = TRUE)

    # Assertion 3: the likelihood the Hessian differentiates matches both the
    # E-step's own ll and the fitted loglik.
    ll_case <- mixtureEM:::.lta_ll_case(fit, X, par, layout)
    ll_e    <- mixtureEM:::.lta_ri_e_step(fit, X, fit$weights_vec)$ll
    expect_equal(ll_case, ll_e, tolerance = 1e-10)
    expect_equal(sum(fit$weights_vec * ll_case), fit$loglik, tolerance = 1e-8)
  }
})

# --- 11. RI x covariates / groups (roadmap ### 14.17) -----------------------

test_that("n_quadrature = 1 reduces RI-plus-covariates exactly to plain covariate LTA", {
  # A single node at 0 with mass 1 contributes nothing regardless of the
  # loadings, so the RI model started AT the plain covariate optimum is
  # already (up to reoptimisation noise) at a fixed point of its own EM --
  # refine_from() makes this close to an identity check rather than a claim
  # that two independent random starts converge to the same basin (they need
  # not, with covariates in play). Only `predictors_initial` is used: with
  # `predictors_transition` too, this fixture's second transition matrix
  # sits close to separation (coefficients drift past +-90 on the logit
  # scale), where two optimiser paths that agree on the likelihood to 0.006
  # can disagree on raw coefficients by a lot -- test-lta-refine.R's
  # covariate gradient check documents the same trap. The tolerances below
  # are optimiser-convergence slack, not algebraic-identity slack.
  sim <- .lta_cov_refine_sim()
  X <- sim$X; Z <- sim$Z
  fit_cov <- suppressWarnings(fit_lta(X, n_statuses = 2, times = 3,
    measurement = "binary", predictors_initial = Z,
    smoothing = 0, bayes_constants = .ml, n_init = 1, random_state = 3,
    tol = 1e-12, max_iter = 5000, standard_errors = FALSE))
  fit_ri1 <- suppressWarnings(fit_lta(X, n_statuses = 2, times = 3,
    measurement = "binary", predictors_initial = Z,
    smoothing = 0, bayes_constants = .ml, refine_from = fit_cov,
    tol = 1e-12, max_iter = 5000, standard_errors = FALSE,
    random_intercept = "continuous", n_quadrature = 1))
  expect_lt(abs(fit_ri1$loglik - fit_cov$loglik), 1e-3)
  expect_lt(max(abs(fit_ri1$delta_beta - fit_cov$delta_beta)), 1e-2)
})

test_that("standard errors and the L-BFGS polish switch themselves on for RI-plus-covariates (roadmap 14.17 W5)", {
  # .lta_scores_supported()/.lta_par_packable() add no covariate test, so
  # lifting the guard in W2 makes this combination packable/scores-full
  # without anyone deciding it should be -- confirm that on purpose here.
  sim <- .lta_cov_refine_sim()
  X <- sim$X; Z <- sim$Z
  fit <- suppressWarnings(fit_lta(X, n_statuses = 2, times = 3,
    measurement = "binary", predictors_initial = Z, predictors_transition = Z,
    random_intercept = "continuous", n_quadrature = 5,
    n_init = 1, max_iter = 25, random_state = 1, standard_errors = TRUE))
  expect_true(.lta_par_packable(fit))
  expect_true(.lta_scores_full(fit))
  expect_false(is.null(fit$se))
  bn  <- vapply(fit$se$blocks, function(x) x$name, "")
  se  <- sqrt(diag(fit$se$vcov))
  db  <- fit$se$blocks[[which(bn == "delta_beta")]]
  expect_true(all(is.finite(se[db$cols])) && all(se[db$cols] > 0))
  for (i in grep("^tau_beta", bn))
    expect_true(all(is.finite(se[fit$se$blocks[[i]]$cols])) &&
                  all(se[fit$se$blocks[[i]]$cols] > 0))

  fit_norefine <- suppressWarnings(fit_lta(X, n_statuses = 2, times = 3,
    measurement = "binary", predictors_initial = Z, predictors_transition = Z,
    random_intercept = "continuous", n_quadrature = 5,
    n_init = 1, max_iter = 25, random_state = 1, standard_errors = FALSE,
    refine = FALSE))
  fit_refine <- suppressWarnings(fit_lta(X, n_statuses = 2, times = 3,
    measurement = "binary", predictors_initial = Z, predictors_transition = Z,
    random_intercept = "continuous", n_quadrature = 5,
    n_init = 1, max_iter = 25, random_state = 1, standard_errors = FALSE,
    refine = TRUE))
  # The standing veto: the polish may only ever raise the log-likelihood.
  expect_gte(fit_refine$loglik, fit_norefine$loglik - 1e-8)
})

test_that("a zero true covariate effect is recovered under RI-LTA", {
  skip_on_cran()
  set.seed(5)
  sim <- .lta_cov_refine_sim()
  noise <- data.frame(z = stats::rnorm(nrow(sim$X)))
  fit <- suppressWarnings(fit_lta(sim$X, n_statuses = 2, times = 3,
    measurement = "binary", predictors_initial = noise,
    smoothing = 0, bayes_constants = .ml,
    random_intercept = "continuous", n_quadrature = 10,
    n_init = 1, random_state = 1, standard_errors = TRUE))
  expect_false(is.null(fit$se))
  b <- fit$se$blocks[[which(vapply(fit$se$blocks, function(x) x$name, "") ==
                             "delta_beta")]]
  est <- fit$delta_beta[1, 2]
  se  <- sqrt(diag(fit$se$vcov))[b$cols][2]
  expect_true(is.finite(se) && se > 0)
  expect_lt(abs(est) / se, 3.5)
})

test_that("the RI score blocks match finite differences for covariate delta/tau", {
  sim <- .lta_cov_refine_sim()
  X <- sim$X; Z <- sim$Z
  fit <- suppressMessages(suppressWarnings(
    fit_lta(X, n_statuses = 2, times = 3, measurement = "binary",
            predictors_initial = Z, predictors_transition = Z,
            random_intercept = "continuous", n_quadrature = 5,
            smoothing = 0, bayes_constants = .ml, n_init = 2,
            random_state = 1, refine = FALSE, standard_errors = FALSE)))
  expect_true(.lta_scores_full(fit))

  layout <- .lta_par_layout(fit)
  par0   <- .lta_par_pack(fit, layout)
  w      <- fit$weights_vec
  expect_equal(length(par0), ncol(.lta_score_matrix(fit, X)$S))
  expect_equal(length(par0), fit$n_params)

  rt <- .lta_par_unpack(par0, fit, layout)
  expect_equal(rt$delta_beta, fit$delta_beta, tolerance = 1e-10)
  expect_equal(rt$tau_beta, fit$tau_beta, tolerance = 1e-10)

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
})

test_that("the finite-difference Hessian is finite and negative definite (smoke, continuous RI)", {
  # refine = FALSE: this checks the packing/Hessian plumbing at whatever EM
  # converges to, not the L-BFGS polish (`### 14.15` W8 turned the polish on
  # by default for RI fits too) -- a small n = 60 smoke fixture is not the
  # place to require the polished optimum's plain-likelihood Hessian to be
  # cleanly negative definite.
  X <- .lta_refine_sim(n = 60, K = 2, Tn = 4, J = 3, seed = 1)
  fit <- suppressWarnings(fit_lta(X, n_statuses = 2, times = 4,
    measurement = "binary", random_intercept = "continuous", n_quadrature = 5,
    n_init = 1, max_iter = 25, random_state = 1, standard_errors = FALSE,
    refine = FALSE))
  layout <- mixtureEM:::.lta_par_layout(fit)
  par    <- mixtureEM:::.lta_par_pack(fit, layout)
  ll_fun <- function(p)
    sum(fit$weights_vec * mixtureEM:::.lta_ll_case(fit, X, p, layout))
  H <- mixtureEM:::.step1_fd_hessian(ll_fun, par)
  expect_true(all(is.finite(H)))
  expect_true(all(eigen(-H, symmetric = TRUE, only.values = TRUE)$values > 0))
})

# ------------------------------------------------------------------------------
# The wide restart search, options(mixtureEM.lta_ri_search = "wide")
# ------------------------------------------------------------------------------

test_that("one Newton step of the binomial M-step never lowers the aggregated log-likelihood", {
  set.seed(4)
  D <- cbind(diag(3)[rep(1:3, times = 4), ], rep(c(-2, -0.5, 0.5, 2), each = 3))
  y <- runif(12); w <- runif(12, 1, 20)
  wll <- function(b) { eta <- as.vector(D %*% b); sum(w * (y * eta - log1p(exp(eta)))) }
  for (s in 1:5) {
    start <- rnorm(4, sd = 4)                      # wild, as a wide start is
    out <- mixtureEM:::.wglm_newton_step(D, y, w, start)
    expect_true(all(is.finite(out$coefficients)))
    expect_gte(wll(out$coefficients), wll(start) - 1e-10)
  }
  # An empty design keeps its start.
  out <- mixtureEM:::.wglm_newton_step(D, y, w * 0, c(0, 0, 0, 0))
  expect_equal(out$coefficients, c(0, 0, 0, 0))
})

test_that("the wide search runs, hands survivors back to the full M-step and leaves the default path alone", {
  X <- .lta_refine_sim(n = 120, K = 2, Tn = 3, J = 3, seed = 2)
  call_fit <- function() suppressWarnings(fit_lta(X, n_statuses = 2, times = 3,
    measurement = "binary", random_intercept = "continuous", n_quadrature = 6,
    n_init = 4, max_iter = 60, random_state = 5, standard_errors = FALSE,
    refine = FALSE))
  old <- options(mixtureEM.lta_ri_search = "current"); on.exit(options(old), add = TRUE)
  f_cur1 <- call_fit()
  options(mixtureEM.lta_ri_search = "wide")
  f_wide <- call_fit()
  expect_true(is.finite(f_wide$loglik))
  expect_true(all(is.finite(f_wide$ri$L)))
  expect_null(f_wide$ri$gem)                  # cleared at promotion
  expect_equal(f_wide$n_params, f_cur1$n_params)
  # The option is read at fit time, so the default path is untouched by it.
  options(mixtureEM.lta_ri_search = "current")
  f_cur2 <- call_fit()
  expect_identical(f_cur2$loglik, f_cur1$loglik)
  # The wide construction itself: every start finite, loadings wide, and the
  # even/odd schemes distinct.
  st <- mixtureEM:::.lta_random_start(f_cur1, X)
  set.seed(1); s_odd  <- mixtureEM:::.lta_ri_wide_start(st, X, 1L)
  set.seed(1); s_even <- mixtureEM:::.lta_ri_wide_start(st, X, 2L)
  expect_true(all(is.finite(s_odd$ri$A)) && all(is.finite(s_even$ri$A)))
  expect_false(isTRUE(all.equal(s_odd$ri$A, s_even$ri$A)))
  expect_true(max(abs(s_odd$ri$L)) > 0.8)
})
