# `predictors_items` as a named list: DIF on the items named, uniform (one
# slope shared by every status) unless the item is in
# `predictors_items_by_status`, and inside the three-step's step 1, where the
# status at every occasion is regressed on the same covariates. The external
# grade is internal/validation-tests/test-part52-lta-dif-validation.R.
#
# Fixture: .lta_dif_sim() (helper-lta-packing.R), three ordinal items, two
# statuses, three occasions, one binary covariate `x`.

.dii_dat <- local({
  d <- NULL
  function() {
    if (is.null(d)) d <<- .lta_dif_sim(n = 500, seed = 5)
    d
  }
})

.dii_fit <- local({
  cache <- list()
  function(key, ...) {
    if (is.null(cache[[key]])) {
      d <- .dii_dat()
      cache[[key]] <<- suppressWarnings(suppressMessages(
        fit_lta(d$X, n_statuses = 2, times = 3, measurement = "ordinal",
                n_init = 2, n_cores = 1, random_state = 11, max_iter = 150,
                ...)))
    }
    cache[[key]]
  }
})

# Item 1 uniform, item 3 by status, item 2 without DIF.
.dii_list <- function() {
  z  <- as.data.frame(.dii_dat()$Z)
  nm <- .dii_fit("plain")$longitudinal$item_names
  stats::setNames(list(z, z), nm[c(1, 3)])
}
.dii_masked <- function() .dii_fit("masked", predictors_items = .dii_list(),
                                   predictors_items_by_status =
                                     .dii_fit("plain")$longitudinal$item_names[3],
                                   standard_errors = TRUE)

test_that("the list form frees exactly the slopes it names", {
  plain  <- .dii_fit("plain", standard_errors = FALSE)
  masked <- .dii_masked()
  K <- 2L
  # one uniform slope on item 1, K on item 3, none on item 2
  expect_equal(masked$n_params, plain$n_params + 1L + K)
  b <- masked$dif$beta
  expect_equal(b[1, 1, 1], b[2, 1, 1])
  expect_true(all(b[, 2, ] == 0))
  expect_false(isTRUE(all.equal(b[1, 3, 1], b[2, 3, 1])))
})

test_that("the masked vector packs, unpacks and scores consistently", {
  fit    <- .dii_masked()
  X      <- fit$data
  w      <- fit$weights_vec
  layout <- .lta_par_layout(fit)
  par    <- .lta_par_pack(fit, layout)
  expect_equal(length(par), fit$n_params)
  back <- .lta_par_unpack(par, fit, layout)
  expect_equal(back$dif$beta, fit$dif$beta)

  sc <- .lta_score_matrix(fit, X)
  expect_equal(ncol(sc$S), length(par))
  analytic <- colSums(sweep(sc$S, 1, w, "*"))
  at <- function(v) sum(w * .lta_ll_case(fit, X, v, layout))
  dif_cols <- unlist(lapply(sc$blocks, function(b)
    if (grepl("^dif\\[", b$name)) b$cols else NULL))
  expect_length(dif_cols, 3L)
  eps <- 1e-5
  for (i in dif_cols) {
    up <- dn <- par; up[i] <- up[i] + eps; dn[i] <- dn[i] - eps
    fd <- (at(up) - at(dn)) / (2 * eps)
    expect_lt(abs(analytic[i] - fd) / max(1, abs(fd)), 1e-4)
  }
})

test_that("the list form refuses what it cannot read", {
  d  <- .dii_dat()
  z  <- as.data.frame(d$Z)
  nm <- .dii_fit("plain")$longitudinal$item_names
  f  <- function(...) fit_lta(d$X, n_statuses = 2, times = 3,
                              measurement = "ordinal", n_init = 1,
                              n_cores = 1, max_iter = 2, ...)
  expect_error(f(predictors_items = list(nope = z)), "not among the items")
  expect_error(f(predictors_items = stats::setNames(list(d$Z[, 1]), nm[1])),
               "column names")
  expect_error(f(predictors_items = stats::setNames(list(z), nm[1]),
                 predictors_items_by_status = nm[2]), "no entry")
  expect_error(f(predictors_items = z, predictors_items_by_status = nm[1]),
               "named list")
})

test_that("a covariate on the transition-free model is origin-free", {
  d  <- .dii_dat()
  z  <- as.data.frame(d$Z)
  fit <- suppressWarnings(fit_lta(d$X, n_statuses = 2, times = 3,
    measurement = "ordinal", n_init = 1, n_cores = 1, random_state = 2,
    max_iter = 50, predictors_initial = z, predictors_transition = z,
    .transition_free = TRUE, standard_errors = FALSE))
  lt <- .lta_log_tau(fit)
  for (t in seq_along(lt)) expect_equal(lt[[t]][[1]], lt[[t]][[2]])
  # one intercept and one slope per non-reference destination, per occasion
  expect_equal(fit$tau_n_params, (3L - 1L) * (2L - 1L) * 2L)
})

test_that("the three-step accepts DIF and keeps it in step 1 only", {
  d  <- .dii_dat()
  nm <- .dii_fit("plain")$longitudinal$item_names
  fit <- suppressWarnings(fit_lta(d$X, n_statuses = 2, times = 3,
    measurement = "ordinal", n_init = 2, n_cores = 1, random_state = 3,
    max_iter = 150, n_steps = 3, assignment = "modal",
    predictors_items = stats::setNames(list(as.data.frame(d$Z)), nm[1])))
  expect_equal(fit$n_steps, 3L)
  expect_false(is.null(fit$step1$dif))
  expect_false(is.null(fit$step1$delta_beta))
  expect_false(is.null(fit$step1$tau_beta))
  expect_null(fit$dif)
  expect_null(fit$delta_beta)
})
