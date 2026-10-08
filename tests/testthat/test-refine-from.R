# ==============================================================================
# refine_from - continuing a fitted model's own EM
# ==============================================================================
#
# The claim being tested is an identity, not an approximation: a run stopped by
# a loose tolerance and then continued under a tight one climbs the same hill as
# a run held to the tight one from the start, because it *is* that run, resumed.
# Anything else would mean the two paths left the same basin, which EM cannot do
# from the same starting values.

.refine_data <- function() {
  set.seed(20)
  K <- 2L
  s1 <- sample.int(K, 500, replace = TRUE, prob = c(0.6, 0.4))
  s2 <- ifelse(s1 == 1L, sample.int(K, 500, TRUE, prob = c(0.8, 0.2)),
                          sample.int(K, 500, TRUE, prob = c(0.3, 0.7)))
  emit <- function(s) {
    p <- ifelse(s == 1L, 0.2, 0.8)
    vapply(seq_len(3), function(j) stats::rbinom(length(s), 1, p),
           numeric(length(s)))
  }
  cbind(emit(s1), emit(s2))
}

test_that("a refined fit reaches the same optimum as a from-scratch tight fit", {
  X <- .refine_data()
  args <- list(indicators = X, n_statuses = 2, times = 2,
               measurement = "binary", n_init = 3, random_state = 5,
               standard_errors = FALSE)

  loose <- suppressWarnings(do.call(fit_lta, args))
  tight <- suppressWarnings(do.call(fit_lta, c(args, list(tol = 1e-13,
                                                          max_iter = 20000))))
  refined <- suppressWarnings(
    fit_lta(X, n_statuses = 2, times = 2, measurement = "binary",
            refine_from = loose, tol = 1e-13, max_iter = 20000,
            random_state = 5, standard_errors = FALSE))

  expect_equal(refined$loglik, tight$loglik, tolerance = 1e-6)
  # Not `- 1e-10`. What the fit maximises is the PENALISED log-posterior, and
  # with the default priors in force that is not the same function as the plain
  # log-likelihood asserted here: measured on this fit, the refined solution is
  # 4.0e-7 *above* the loose one on the penalised objective while sitting
  # 1.7e-5 below it on the raw log-likelihood. Both are behaving correctly and
  # trading one against the other is what the prior is for. The tight rule only
  # held before the L-BFGS refinement existed, when both fits stopped at the
  # same EM fixed point and the two objectives could not separate.
  expect_gte(refined$loglik, loose$loglik - 1e-4)
})

test_that("a refined fit records that it ran no pool", {
  X <- .refine_data()
  loose <- suppressWarnings(
    fit_lta(X, n_statuses = 2, times = 2, measurement = "binary",
            n_init = 3, random_state = 5, standard_errors = FALSE))
  refined <- suppressWarnings(
    fit_lta(X, n_statuses = 2, times = 2, measurement = "binary",
            refine_from = loose, tol = 1e-13, standard_errors = FALSE))

  expect_identical(refined$metrics$n_requested, 1L)
  expect_true(isTRUE(refined$refined_from))
  expect_false(isTRUE(loose$refined_from))
})

test_that("fit_mixture() refines its own fit and refuses the wrong donor", {
  set.seed(21)
  z <- sample(2, 400, replace = TRUE)
  X <- vapply(seq_len(5), function(j)
    stats::rbinom(400, 1, ifelse(z == 1L, 0.2, 0.8)), numeric(400))

  a <- suppressWarnings(fit_mixture(X, measurement = "binary", n_classes = 2,
                                    n_init = 5, random_state = 2))
  b <- suppressWarnings(fit_mixture(X, measurement = "binary", n_classes = 2,
                                    refine_from = a, random_state = 2))

  # `a` already sits at the maximum (the Newton-type finish, R/qn_finish.R),
  # so refining it can move the reported log-likelihood by rounding only.
  expect_gte(b$metrics$ll, a$metrics$ll - 1e-8)
  expect_equal(b$metrics$ll, a$metrics$ll, tolerance = 1e-6)
  expect_identical(b$metrics$n_requested, 1L)
  expect_true(isTRUE(b$refined_from))

  expect_error(fit_mixture(X, measurement = "binary", n_classes = 3,
                           refine_from = a), "classes")
  expect_error(fit_mixture(X, measurement = "binary", n_classes = 2,
                           refine_from = a, n_init = 5), "nothing to size")
  expect_error(fit_mixture(X, measurement = "binary", n_classes = 2,
                           refine_from = list()), "fitted by fit_mixture")
})

test_that("refine_from refuses a donor of the wrong shape, and refuses n_init", {
  X <- .refine_data()
  loose <- suppressWarnings(
    fit_lta(X, n_statuses = 2, times = 2, measurement = "binary",
            n_init = 3, random_state = 5, standard_errors = FALSE))

  expect_error(
    fit_lta(X, n_statuses = 3, times = 2, measurement = "binary",
            refine_from = loose, standard_errors = FALSE),
    "statuses")
  expect_error(
    fit_lta(X, n_statuses = 2, times = 3, measurement = "binary",
            refine_from = loose, standard_errors = FALSE),
    "occasions")
  expect_error(
    fit_lta(X, n_statuses = 2, times = 2, measurement = "binary",
            refine_from = list(), standard_errors = FALSE),
    "fitted by fit_lta")
  expect_error(
    fit_lta(X, n_statuses = 2, times = 2, measurement = "binary",
            refine_from = loose, n_init = 5, standard_errors = FALSE),
    "nothing to size")
})

# A random-intercept donor hands over its intercepts (or ordinal thresholds) and
# loadings, not only the integrated item table they summarise. With no EM
# iteration the refit is the donor, so it must report the donor's likelihood;
# before the fix the loadings came from a random draw and this was ~150 units
# off on simulated data.
test_that("refine_from continues a random-intercept donor exactly", {
  sim <- .lta_ri_sim(n = 300, Tn = 3, J = 4, seed = 7)
  donor <- suppressWarnings(
    fit_lta(sim$X, n_statuses = 2, times = 3, measurement = "binary",
            random_intercept = "continuous", n_quadrature = 7, n_init = 2,
            n_cores = 1, random_state = 1, standard_errors = FALSE))
  again <- suppressWarnings(
    fit_lta(sim$X, n_statuses = 2, times = 3, measurement = "binary",
            random_intercept = "continuous", n_quadrature = 7,
            refine_from = donor, max_iter = 0, refine = FALSE,
            random_state = 99, standard_errors = FALSE))
  expect_equal(again$loglik, donor$loglik, tolerance = 1e-8)
  expect_equal(again$ri$L, donor$ri$L)
  expect_equal(again$ri$A, donor$ri$A)

  Xo <- .lta_ordinal_sim(n = 150, Tn = 3, K = 2, cats = c(3L, 3L, 2L), seed = 1)
  donor_o <- suppressWarnings(
    fit_lta(Xo, n_statuses = 2, times = 3, measurement = "ordinal",
            random_intercept = "continuous", n_quadrature = 5, n_init = 2,
            n_cores = 1, random_state = 1, max_iter = 200,
            standard_errors = FALSE))
  again_o <- suppressWarnings(
    fit_lta(Xo, n_statuses = 2, times = 3, measurement = "ordinal",
            random_intercept = "continuous", n_quadrature = 5,
            refine_from = donor_o, max_iter = 0, refine = FALSE,
            random_state = 99, standard_errors = FALSE))
  expect_equal(again_o$loglik, donor_o$loglik, tolerance = 1e-8)
  expect_equal(again_o$ri$theta, donor_o$ri$theta)
})

test_that("refine_from refuses a random intercept of the other kind", {
  sim <- .lta_ri_sim(n = 200, Tn = 3, J = 4, seed = 7)
  donor <- suppressWarnings(
    fit_lta(sim$X, n_statuses = 2, times = 3, measurement = "binary",
            random_intercept = "binary", n_ri = 2, n_init = 1, max_iter = 50,
            n_cores = 1, random_state = 1, standard_errors = FALSE))
  expect_error(
    fit_lta(sim$X, n_statuses = 2, times = 3, measurement = "binary",
            random_intercept = "continuous", n_quadrature = 5,
            refine_from = donor, standard_errors = FALSE),
    "binary random intercept")
})

# A donor fitted without predictors, continued with them: the new regressions
# start at the coefficients that reproduce the donor's own probabilities, so
# with no EM iteration the refit is still the donor. Started from zero
# coefficients, the first M-step on a donor with an empty transition cell fell
# below the donor and the run settled on a worse solution than the model it
# nests.
test_that("refine_from starts added predictors at the donor's probabilities", {
  sim <- .lta_cov_refine_sim()
  ml <- list(smoothing = 0, bayes_constants = list(categorical = 0))
  donor <- fit_lta(sim$X, n_statuses = 2, times = 3, measurement = "binary",
                   n_init = 3, n_cores = 1, random_state = 1,
                   standard_errors = FALSE, smoothing = ml$smoothing,
                   bayes_constants = ml$bayes_constants)
  for (effects in c("common", "by_origin")) {
    again <- suppressWarnings(fit_lta(sim$X, n_statuses = 2, times = 3,
                     measurement = "binary",
                     predictors_initial = sim$Z, predictors_transition = sim$Z,
                     transition_effects = effects, refine_from = donor,
                     max_iter = 0, refine = FALSE, standard_errors = FALSE,
                     smoothing = ml$smoothing,
                     bayes_constants = ml$bayes_constants))
    expect_equal(again$loglik, donor$loglik, tolerance = 1e-8)
    # The regressions exist from the start and reproduce the donor's tables:
    # every slope zero, the averaged probabilities the donor's own.
    expect_false(is.null(again$tau_beta))
    expect_false(is.null(again$delta_beta))
    expect_equal(again$tau, donor$tau, tolerance = 1e-8)
    expect_equal(again$delta, donor$delta, tolerance = 1e-8)
  }
  full <- fit_lta(sim$X, n_statuses = 2, times = 3, measurement = "binary",
                  predictors_initial = sim$Z, predictors_transition = sim$Z,
                  refine_from = donor, standard_errors = FALSE,
                  smoothing = ml$smoothing,
                  bayes_constants = ml$bayes_constants)
  expect_gte(full$loglik, donor$loglik)
})
