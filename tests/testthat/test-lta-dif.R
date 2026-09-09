# Direct covariate effects on the indicators -- `predictors_items`, the
# package's test of measurement invariance with respect to a covariate
# (part-r9-item-dif.md, W10). One proportional-odds slope per latent status per
# item per covariate, shared across occasions.
#
# The fixture is `.lta_dif_sim()` in helper-lta-packing.R: three ordinal items
# with ragged category counts, two statuses, three occasions, one binary
# covariate, and a true slope of 0.8 on item 1 in status 1 only.
#
# Every fit in this file is memoised at file level and capped at a small
# `max_iter`, because the structural identities below are properties of the
# code at whatever point EM has reached, not of a converged maximum. The one
# exception is the recovery check at the end, which needs a real optimum. The
# same budgeting took test-lta-transition-slopes.R from 12.3 minutes to 2.

.dif_dat <- local({
  d <- NULL
  function() {
    if (is.null(d)) d <<- .lta_dif_sim(n = 600, seed = 3)
    d
  }
})

.dif_fit <- local({
  cache <- list()
  function(key, ...) {
    if (is.null(cache[[key]])) {
      d <- .dif_dat()
      cache[[key]] <<- suppressWarnings(suppressMessages(
        fit_lta(d$X, n_statuses = 2, times = 3, measurement = "ordinal",
                n_init = 2, n_cores = 1, random_state = 11, max_iter = 150,
                ...)))
    }
    cache[[key]]
  }
})

# The four model shapes the packing has to describe.
.fit_dif      <- function() .dif_fit("dif", predictors_items = .dif_dat()$Z,
                                     standard_errors = TRUE)
.fit_dif_ri   <- function() .dif_fit("dif_ri", predictors_items = .dif_dat()$Z,
                                     random_intercept = "continuous",
                                     n_quadrature = 5, standard_errors = FALSE)
.fit_plain    <- function() .dif_fit("plain", standard_errors = FALSE)
.fit_ri       <- function() .dif_fit("ri", random_intercept = "continuous",
                                     n_quadrature = 5, standard_errors = FALSE)

.ll_of <- function(st, X) sum(st$weights_vec * .lta_ri_ll_case(st, X))

# ------------------------------------------------------------------------------
# 1. At d = 0 the model IS the model without the covariate
# ------------------------------------------------------------------------------

test_that("zero DIF slopes reproduce the no-DIF likelihood exactly", {
  fit <- .fit_dif()
  X   <- fit$data

  zeroed <- fit
  zeroed$dif$beta[] <- 0
  stripped <- fit
  stripped$dif <- NULL

  expect_lt(abs(.ll_of(zeroed, X) - .ll_of(stripped, X)), 1e-10)
})

# ------------------------------------------------------------------------------
# 2. The degenerate factor is regular LTA, in both the likelihood and the count
# ------------------------------------------------------------------------------
#
# This is what proves the parameter count of the reference program's
# 76-parameter model before any expensive fit is run: a `predictors_items` fit
# borrows the random-intercept machinery at one node with the loading fixed at
# zero, and that borrowing must cost exactly nothing.

test_that("the degenerate one-node factor reduces exactly to a plain ordinal fit", {
  fit  <- .fit_plain()
  X    <- fit$data
  cats <- fit$mm$models[[1]]$cats

  degen <- fit
  degen$ri <- list(
    kind = "continuous", Dnode = matrix(0, 1L, 1L), mass = 1,
    A = NULL, L = matrix(0, fit$n_items, 1L),
    theta = .ordinal_theta_from_pis(fit$mm$models[[1]]$parameters$pis, cats),
    cats = cats, loading_free = FALSE)

  expect_lt(abs(.ll_of(degen, X) - fit$loglik), 1e-8)
  expect_identical(.lta_n_parameters(degen), fit$n_params)
  # ... and with the loading treated as free it would cost three more, which is
  # exactly the 79-versus-76 difference the `loading_free` flag exists to make.
  free <- degen
  free$ri$loading_free <- TRUE
  expect_identical(.lta_n_parameters(free), fit$n_params + 3L)
})

# ------------------------------------------------------------------------------
# 3. The two likelihood paths agree, with DIF switched on
# ------------------------------------------------------------------------------
#
# .lta_ri_e_step() and .lta_ri_ll_case() build the emission tables separately;
# the finite-difference Hessian differentiates only the second, so a
# disagreement would make every standard error wrong while every reported
# log-likelihood stayed right.

test_that("the E-step and the case-likelihood paths agree under DIF", {
  for (fit in list(.fit_dif(), .fit_dif_ri())) {
    X <- fit$data
    expect_lt(abs(.lta_score_matrix(fit, X)$ll[1] -
                    .lta_ri_ll_case(fit, X)[1]), 1e-10)
    expect_lt(abs(sum(fit$weights_vec * .lta_score_matrix(fit, X)$ll) -
                    .ll_of(fit, X)), 1e-10)
  }
})

# ------------------------------------------------------------------------------
# 4. The score blocks against finite differences
# ------------------------------------------------------------------------------

test_that("the DIF score block matches finite differences", {
  eps <- 1e-5
  for (fit in list(.fit_dif(), .fit_dif_ri())) {
    X      <- fit$data
    w      <- fit$weights_vec
    layout <- .lta_par_layout(fit)
    par    <- .lta_par_pack(fit, layout)
    sc     <- .lta_score_matrix(fit, X)
    expect_equal(ncol(sc$S), length(par))

    analytic <- colSums(sweep(sc$S, 1, w, "*"))
    at <- function(v) sum(w * .lta_ll_case(fit, X, v, layout))

    dif_cols <- unlist(lapply(sc$blocks, function(b)
      if (grepl("^dif\\[", b$name)) b$cols else NULL))
    expect_gt(length(dif_cols), 0L)

    # Every DIF column, and a spread of the rest of the vector.
    check <- union(dif_cols, unique(round(seq(1, length(par), length.out = 8))))
    for (i in check) {
      up <- dn <- par; up[i] <- up[i] + eps; dn[i] <- dn[i] - eps
      fd <- (at(up) - at(dn)) / (2 * eps)
      expect_lt(abs(analytic[i] - fd) / max(1, abs(fd)), 1e-4)
    }
  }
})

# ------------------------------------------------------------------------------
# 5. The prior in the M-step and the prior in the penalty are one objective
# ------------------------------------------------------------------------------
#
# The measurement prior is applied at DIF = 0, on grid rows with a zero
# covariate design. If the M-step's prior rows and .lta_penalty()'s branches
# were not the same objective, EM's stopping rule and the L-BFGS polish would
# be climbing different hills, and this is the only check that would say so.

test_that("EM's fixed point is a stationary point of the penalised objective", {
  d   <- .dif_dat()
  fit <- suppressWarnings(suppressMessages(
    fit_lta(d$X, n_statuses = 2, times = 3, measurement = "ordinal",
            predictors_items = d$Z, n_init = 1, n_cores = 1,
            random_state = 11, max_iter = 600, tol = 1e-10,
            smoothing = 0.5, bayes_constants = list(categorical = 1),
            refine = FALSE, standard_errors = FALSE)))
  X      <- fit$data
  w      <- fit$weights_vec
  layout <- .lta_par_layout(fit)
  sc     <- .lta_score_matrix(fit, X)
  pen    <- .lta_penalty(fit, X, layout, 0.5)
  grad   <- colSums(sweep(sc$S, 1, w, "*")) + pen$gradient
  expect_lt(max(abs(grad)), 0.05)
})

# ------------------------------------------------------------------------------
# 6. Relabelling: dif$beta is status-indexed and must travel with the labels
# ------------------------------------------------------------------------------
#
# `ri$theta` was left unpermuted once and the fit then carried two different
# models at once -- everything a user read stayed right and everything
# recomputed was evaluated at a scrambled model. `dif$beta` is the second
# status-indexed quantity outside `pis`.

test_that(".sort_lta_statuses() permutes dif$beta with the statuses", {
  fit <- .fit_dif()

  st <- fit
  st$delta_c[[1]] <- c(0.2, 0.8)        # forces the ordering to swap
  before <- st$dif$beta
  after  <- .sort_lta_statuses(st)$dif$beta
  expect_equal(after[1, , , drop = FALSE], before[2, , , drop = FALSE])
  expect_equal(after[2, , , drop = FALSE], before[1, , , drop = FALSE])
})

test_that("a relabelled fit still reproduces its own reported log-likelihood", {
  fit <- .fit_dif()
  expect_lt(abs(.ll_of(fit, fit$data) - fit$loglik), 1e-8)
})

# ------------------------------------------------------------------------------
# 6b. The invariance test itself
# ------------------------------------------------------------------------------
#
# Fitting the model with and without `predictors_items` and comparing the two
# with lr_test() IS the test of measurement invariance, and the documentation
# says so, so the pair has to actually nest in the way lr_test() reads.

test_that("the no-DIF and DIF fits nest, so lr_test() tests the invariance", {
  a <- .fit_plain()
  b <- .fit_dif()
  res <- suppressWarnings(lr_test(a, b))
  expect_equal(res$df, b$n_params - a$n_params)
  expect_equal(res$df, 2 * 3 * 1)          # statuses x items x covariates
  expect_true(is.finite(res$statistic))
})

# ------------------------------------------------------------------------------
# 7. Monotonicity
# ------------------------------------------------------------------------------

test_that("EM does not go downhill with DIF on", {
  d <- .dif_dat()
  objective <- vapply(c(2L, 5L, 10L, 20L, 35L, 50L), function(m) {
    f <- suppressWarnings(suppressMessages(
      fit_lta(d$X, n_statuses = 2, times = 3, measurement = "ordinal",
              predictors_items = d$Z, n_init = 1, n_cores = 1,
              random_state = 11, max_iter = m, refine = FALSE,
              order_by_size = FALSE, standard_errors = FALSE)))
    f$loglik + .lta_log_prior(f, f$data, 1.0)
  }, numeric(1))
  expect_true(all(diff(objective) > -1e-8))
})

# ------------------------------------------------------------------------------
# 8. Packing invariants, on all four shapes
# ------------------------------------------------------------------------------

test_that("the packed vector describes every DIF and RI shape", {
  shapes <- list(`DIF only` = .fit_dif(), `DIF + RI` = .fit_dif_ri(),
                 `RI only` = .fit_ri(), `neither` = .fit_plain())
  for (nm in names(shapes)) {
    fit    <- shapes[[nm]]
    X      <- fit$data
    layout <- .lta_par_layout(fit)
    par    <- .lta_par_pack(fit, layout)

    # (a) the vector is the model's free parameters, no more and no fewer
    expect_equal(length(par), fit$n_params, info = nm)
    # (b) the score matrix describes the same vector
    expect_equal(ncol(.lta_score_matrix(fit, X)$S), length(par), info = nm)
    # (c) unpack(pack(state)) is the same model
    round_trip <- .lta_par_unpack(par, fit, layout)
    expect_lt(abs(sum(fit$weights_vec *
                        .lta_score_matrix(round_trip, X)$ll) - fit$loglik), 1e-8)
  }
})

# ------------------------------------------------------------------------------
# 9. Recovery of a known slope
# ------------------------------------------------------------------------------

test_that("a known DIF slope is recovered", {
  skip_on_cran()
  d <- .lta_dif_sim(n = 3000, seed = 17)
  fit <- suppressWarnings(suppressMessages(
    fit_lta(d$X, n_statuses = 2, times = 3, measurement = "ordinal",
            predictors_items = d$Z, n_init = 4, n_cores = 1,
            random_state = 2, standard_errors = FALSE)))
  # Status labels are arbitrary, so the truth is "one status of item 1 carries
  # 0.8 and everything else is near zero", not "row 1 does".
  b <- fit$dif$beta[, , 1]
  expect_lt(min(abs(b[, 1] - 0.8)), 0.15)
  expect_lt(max(abs(b[, 2:3])), 0.3)
})

# ------------------------------------------------------------------------------
# 10. The pattern cap is a refusal, not a warning
# ------------------------------------------------------------------------------

test_that("too many covariate patterns are refused at the argument", {
  n <- 200
  ok  <- matrix(rep_len(seq_len(16), n), ncol = 1L, dimnames = list(NULL, "z"))
  bad <- matrix(rep_len(seq_len(17), n), ncol = 1L, dimnames = list(NULL, "z"))

  expect_silent(.lta_dif_init(ok, n, 2L, 3L))
  expect_error(.lta_dif_init(bad, n, 2L, 3L),
               "17 distinct covariate patterns, above the limit of 16")
  # Raising the option deliberately is the documented escape hatch.
  withr::with_options(list(mixtureEM.dif_max_patterns = 20L),
                      expect_silent(.lta_dif_init(bad, n, 2L, 3L)))
})
