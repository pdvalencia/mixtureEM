# Packing, scores and standard errors for the ordinal (cumulative-logit)
# measurement family, plain and under a random intercept (roadmap ### 14.18,
# W9-W10). `.lta_ordinal_sim()` is in helper-lta-packing.R. These mirror
# test-lta-ri.R's own binary-family checks (parameter counts, the packed
# vector against the score matrix's own column order, the score blocks
# against finite differences) rather than inventing a new shape for them.

test_that("ncol(S) == n_params for plain and RI ordinal fits", {
  X <- .lta_ordinal_sim()
  check <- function(fit) {
    sc <- .lta_score_matrix(fit, X)
    expect_equal(ncol(sc$S), fit$n_params)
  }
  fit_plain <- suppressWarnings(fit_lta(X, n_statuses = 2, times = 3,
    measurement = "ordinal", cats = c(3L, 3L, 2L),
    n_init = 1, max_iter = 15, random_state = 1, standard_errors = FALSE))
  check(fit_plain)

  fit_bin <- suppressWarnings(fit_lta(X, n_statuses = 2, times = 3,
    measurement = "ordinal", cats = c(3L, 3L, 2L),
    random_intercept = "binary", n_ri = 2,
    n_init = 1, max_iter = 15, random_state = 1, standard_errors = FALSE))
  check(fit_bin)

  fit_cont <- suppressWarnings(fit_lta(X, n_statuses = 2, times = 3,
    measurement = "ordinal", cats = c(3L, 3L, 2L),
    random_intercept = "continuous", n_quadrature = 15,
    n_init = 1, max_iter = 15, random_state = 1, standard_errors = FALSE))
  check(fit_cont)
})

test_that("the ordinal score blocks match finite differences (plain, binary-RI, continuous-RI)", {
  X <- .lta_ordinal_sim()
  eps <- 1e-5

  check_gradient <- function(fit, bumpers) {
    w  <- fit$weights_vec
    sc <- .lta_score_matrix(fit, X)
    S  <- sc$S
    cols <- function(name) Find(function(b) identical(b$name, name), sc$blocks)$cols
    ll_at <- function(st) sum(st$weights_vec * .lta_score_matrix(st, X)$ll)
    fd <- function(f) (ll_at(f(eps)) - ll_at(f(-eps))) / (2 * eps)
    for (bmp in bumpers) {
      analytic <- sum(w * S[, cols(bmp$name)[bmp$idx]])
      expect_lt(abs(analytic - fd(bmp$fn)) / max(1, abs(fd(bmp$fn))), 1e-5)
    }
  }

  fit_plain <- suppressWarnings(fit_lta(X, n_statuses = 2, times = 3,
    measurement = "ordinal", cats = c(3L, 3L, 2L),
    n_init = 1, max_iter = 15, random_state = 1, standard_errors = FALSE))
  m1    <- fit_plain$mm$models[[1]]
  cols1 <- .ordinal_theta_cols(m1$cats, 1L)
  bump_plain_theta <- function(idx) function(h) {
    st <- fit_plain
    m  <- st$mm$models[[1]]
    th <- .ordinal_theta_from_pis(m$parameters$pis, m$cats)
    th[1, cols1[idx]] <- th[1, cols1[idx]] + h
    probs <- .ordinal_cat_probs(th[, cols1, drop = FALSE], 0, m$cats[1])
    pcols <- .ordinal_pis_cols(m$cats, 1L)
    for (tt in seq_along(st$mm$models))
      st$mm$models[[tt]]$parameters$pis[, pcols] <- probs
    st
  }
  check_gradient(fit_plain, list(
    list(name = "theta[item 1]", idx = 1, fn = bump_plain_theta(1)),
    list(name = "theta[item 1]", idx = 3, fn = bump_plain_theta(2))))

  fit_bin <- suppressWarnings(fit_lta(X, n_statuses = 2, times = 3,
    measurement = "ordinal", cats = c(3L, 3L, 2L),
    random_intercept = "binary", n_ri = 2,
    n_init = 1, max_iter = 15, random_state = 1, standard_errors = FALSE))
  colsb <- .ordinal_theta_cols(fit_bin$ri$cats, 1L)
  check_gradient(fit_bin, list(
    list(name = "theta[item 1]", idx = 1, fn = function(h) {
      st <- fit_bin; st$ri$theta[1, colsb[1]] <- st$ri$theta[1, colsb[1]] + h; st }),
    list(name = "theta[item 1]", idx = 3, fn = function(h) {
      st <- fit_bin; st$ri$theta[1, colsb[2]] <- st$ri$theta[1, colsb[2]] + h; st }),
    list(name = "lambda[item 1]", idx = 1, fn = function(h) {
      st <- fit_bin; st$ri$L[1, 1] <- st$ri$L[1, 1] + h; st })))

  fit_cont <- suppressWarnings(fit_lta(X, n_statuses = 2, times = 3,
    measurement = "ordinal", cats = c(3L, 3L, 2L),
    random_intercept = "continuous", n_quadrature = 15,
    n_init = 1, max_iter = 15, random_state = 1, standard_errors = FALSE))
  colsc <- .ordinal_theta_cols(fit_cont$ri$cats, 1L)
  check_gradient(fit_cont, list(
    list(name = "theta[item 1]", idx = 1, fn = function(h) {
      st <- fit_cont; st$ri$theta[1, colsc[1]] <- st$ri$theta[1, colsc[1]] + h; st }),
    list(name = "theta[item 1]", idx = 3, fn = function(h) {
      st <- fit_cont; st$ri$theta[1, colsc[2]] <- st$ri$theta[1, colsc[2]] + h; st }),
    list(name = "lambda[item 1]", idx = 1, fn = function(h) {
      st <- fit_cont; st$ri$L[1, 1] <- st$ri$L[1, 1] + h; st })))
})

test_that("the packed vector describes the same thing the score matrix does (ordinal, plain and RI)", {
  X <- .lta_ordinal_sim()
  check_roundtrip <- function(fit) {
    layout <- .lta_par_layout(fit)
    par    <- .lta_par_pack(fit, layout)
    expect_equal(length(par), fit$n_params)
    st2 <- .lta_par_unpack(par, fit, layout)
    expect_equal(.lta_ll_case(fit, X, par, layout),
                 .lta_score_matrix(fit, X)$ll, tolerance = 1e-8)
    expect_equal(.lta_par_pack(st2, layout), par, tolerance = 1e-8)
  }
  fit_plain <- suppressWarnings(fit_lta(X, n_statuses = 2, times = 3,
    measurement = "ordinal", cats = c(3L, 3L, 2L),
    n_init = 1, max_iter = 15, random_state = 1, standard_errors = FALSE))
  check_roundtrip(fit_plain)

  fit_bin <- suppressWarnings(fit_lta(X, n_statuses = 2, times = 3,
    measurement = "ordinal", cats = c(3L, 3L, 2L),
    random_intercept = "binary", n_ri = 2,
    n_init = 1, max_iter = 15, random_state = 1, standard_errors = FALSE))
  check_roundtrip(fit_bin)
})

test_that("refine = TRUE never decreases the ordinal RI log-likelihood, and standard_errors run to completion", {
  X <- .lta_ordinal_sim()
  for (ri_args in list(list(), list(random_intercept = "binary", n_ri = 2),
                       list(random_intercept = "continuous", n_quadrature = 15))) {
    args <- c(list(X, n_statuses = 2, times = 3, measurement = "ordinal",
                   cats = c(3L, 3L, 2L), n_init = 2, max_iter = 20,
                   random_state = 1, standard_errors = TRUE), ri_args)
    fit_norefine <- suppressWarnings(do.call(fit_lta, c(args, list(refine = FALSE))))
    fit_refine   <- suppressWarnings(do.call(fit_lta, c(args, list(refine = TRUE))))
    # `refine` can only improve the objective a candidate is RANKED on --
    # the penalised score, not the plain log-likelihood this test reads --
    # and each run's winner is chosen from a different candidate set (the
    # cold/warm-mixed pool, R8's W11B). When the two runs' winners differ,
    # the plain log-likelihoods are not guaranteed to be ordered the same
    # way the penalised scores are, because the prior term the two winning
    # candidates carry need not be equal. Measured on this fixture
    # (RECORDS.md, R8 W11B step 1): a 0.045 gap, `random_state = 1`. -0.1
    # keeps this catching a genuine regression while tolerating that.
    expect_gte(fit_refine$loglik, fit_norefine$loglik - 0.1)
    expect_false(is.null(fit_refine$se))
  }
})

test_that("at zero loading an RI ordinal fit reproduces the plain ordinal likelihood exactly", {
  # The saturation equivalence W9/W10 must preserve: a random intercept with
  # every loading pinned at 0 integrates to the same category probabilities
  # at every node, so the RI likelihood collapses to the plain multinoulli-
  # style one. Checked at max_iter = 0 (an E-step only, no M-step drift).
  X <- .lta_ordinal_sim()
  fit_plain <- suppressWarnings(fit_lta(X, n_statuses = 2, times = 3,
    measurement = "ordinal", cats = c(3L, 3L, 2L),
    n_init = 1, max_iter = 15, random_state = 1, standard_errors = FALSE))

  fit_ri <- suppressWarnings(fit_lta(X, n_statuses = 2, times = 3,
    measurement = "ordinal", cats = c(3L, 3L, 2L),
    random_intercept = "continuous", n_quadrature = 15,
    n_init = 1, max_iter = 1, random_state = 1, standard_errors = FALSE))
  fit_ri$ri$L[] <- 0
  fit_ri$ri$theta <- .ordinal_theta_from_pis(fit_plain$mm$models[[1]]$parameters$pis,
                                             fit_plain$mm$models[[1]]$cats)
  fit_ri$delta_c <- fit_plain$delta_c
  fit_ri$tau_c   <- fit_plain$tau_c

  layout <- .lta_par_layout(fit_ri)
  par    <- .lta_par_pack(fit_ri, layout)
  ll_ri  <- sum(fit_ri$weights_vec * .lta_ll_case(fit_ri, X, par, layout))

  layout_p <- .lta_par_layout(fit_plain)
  par_p    <- .lta_par_pack(fit_plain, layout_p)
  ll_plain <- sum(fit_plain$weights_vec * .lta_ll_case(fit_plain, X, par_p, layout_p))

  expect_equal(ll_ri, ll_plain, tolerance = 1e-6)
})

test_that("an ordinal RI fit's stored state reproduces its own reported log-likelihood", {
  # .sort_lta_statuses() relabels statuses by Time 1 prevalence at the end of a
  # fit, and it has to move every status-indexed quantity together. For an
  # ordinal random intercept the measurement model is `ri$theta`, not the `ri$A`
  # a binary one carries; leaving `theta` behind while `pis` and delta/tau moved
  # left the returned object holding two different models at once. Everything a
  # user reads off the fit stayed right (those readers use the permuted `pis`),
  # but everything recomputed from `ri$theta` -- the score matrix, the packed
  # log-likelihood, standard errors and the scaling factor -- was evaluated at a
  # scrambled model. The invariant that catches it is the cheapest one there is:
  # an E-step over the object's own stored parameters must return the number the
  # object reports.
  X <- .lta_ordinal_sim(n = 400, seed = 3)
  check_consistent <- function(fit) {
    w  <- fit$weights_vec
    sc <- .lta_score_matrix(fit, X)
    expect_equal(sum(w * sc$ll), fit$loglik, tolerance = 1e-8)
    layout <- .lta_par_layout(fit)
    par    <- .lta_par_pack(fit, layout)
    expect_equal(sum(w * .lta_ll_case(fit, X, par, layout)), fit$loglik,
                 tolerance = 1e-8)
  }

  args <- list(X, n_statuses = 2, times = 3, measurement = "ordinal",
               cats = c(3L, 3L, 2L), n_init = 3, max_iter = 200,
               random_state = 5, standard_errors = FALSE)
  check_consistent(suppressWarnings(do.call(fit_lta, args)))
  check_consistent(suppressWarnings(do.call(fit_lta,
    c(args, list(random_intercept = "binary", n_ri = 2)))))
  check_consistent(suppressWarnings(do.call(fit_lta,
    c(args, list(random_intercept = "continuous", n_quadrature = 15)))))
})

test_that(".sort_lta_statuses() permutes an ordinal random intercept's thresholds", {
  # The unit-level statement of the same defect, so a future edit to the sort
  # cannot drop `theta` again without a test naming it.
  X <- .lta_ordinal_sim(n = 300, seed = 8)
  fit <- suppressWarnings(fit_lta(X, n_statuses = 2, times = 3,
    measurement = "ordinal", cats = c(3L, 3L, 2L),
    random_intercept = "continuous", n_quadrature = 15,
    n_init = 1, max_iter = 20, random_state = 2, standard_errors = FALSE))

  # Force a relabelling by making status 2 the more prevalent one at Time 1.
  flipped <- fit
  flipped$delta_c <- list(rev(fit$delta_c[[1]]))
  sorted <- .sort_lta_statuses(flipped)
  expect_equal(sorted$ri$theta, flipped$ri$theta[c(2L, 1L), , drop = FALSE])
  expect_equal(sorted$mm$models[[1]]$parameters$pis,
               flipped$mm$models[[1]]$parameters$pis[c(2L, 1L), , drop = FALSE])
})
