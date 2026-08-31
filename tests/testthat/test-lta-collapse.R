# fit_lta() runs its search on the response-pattern table where it can. The
# claim is that this is a saving and not a change: the likelihood is a weighted
# sum over cases either way, and the starts never see the data. These tests are
# what stands behind that claim.

.coll_sim <- function(seed = 3, n = 1200, K = 3, Tn = 2, J = 4) {
  set.seed(seed)
  tau <- matrix(c(.8, .15, .05,
                  .10, .80, .10,
                  .05, .15, .80), 3, 3, byrow = TRUE)
  rho <- matrix(c(.85, .75, .20, .25,
                  .75, .25, .80, .20,
                  .20, .25, .25, .85), 3, 4, byrow = TRUE)
  s <- sample(K, n, TRUE)
  b <- vector("list", Tn)
  for (t in seq_len(Tn)) {
    if (t > 1L) s <- vapply(s, function(k) sample(K, 1L, prob = tau[k, ]),
                            integer(1))
    b[[t]] <- matrix(stats::rbinom(n * J, 1L, rho[s, ]), n, J)
  }
  X <- do.call(cbind, b)
  colnames(X) <- paste0("t", rep(seq_len(Tn), each = J),
                        "_i", rep(seq_len(J), Tn))
  X
}

# The skeleton a real fit builds, obtained the cheap way, so the eligibility
# predicate can be asked about a state the driver would actually have made.
.coll_state <- function(X, ...) suppressMessages(suppressWarnings(
  fit_lta(X, n_statuses = 3, times = 2, measurement = "binary",
          n_init = 1, max_iter = 1L, refine = FALSE,
          standard_errors = FALSE, order_by_size = FALSE, ...)))

test_that(".pattern_index() groups rows exactly, missing cells included", {
  X <- rbind(c(1, 0, 1),
             c(0, 0, 1),
             c(1, 0, 1),      # duplicate of row 1
             c(1, NA, 1),     # NA is its own value, not a wildcard
             c(1, NA, 1))
  pat <- .pattern_index(X)

  expect_identical(pat$rep_row, c(TRUE, TRUE, FALSE, TRUE, FALSE))
  expect_identical(pat$idx, c(1L, 2L, 1L, 3L, 3L))
  # idx indexes the retained rows in first-appearance order.
  expect_equal(X[pat$rep_row, , drop = FALSE][pat$idx, ], X)
})

test_that(".pattern_index() agrees with the radix key it replaced", {
  # The radix key is exact below 2^53 and was the package's definition of "the
  # same response pattern" until this refactor. Where both are valid they must
  # give the same grouping, or fit_mixture()'s collapsed path changed meaning.
  set.seed(21)
  for (mv in c(1, 3)) {
    X <- matrix(sample(c(0:mv, NA), 400 * 6, TRUE), 400, 6)
    key <- numeric(nrow(X))
    for (j in seq_len(ncol(X)))
      key <- key * (mv + 2) + ifelse(is.na(X[, j]), 0, X[, j] + 1)
    expect_identical(.pattern_index(X)$idx, match(key, key[!duplicated(key)]))
  }
})

test_that(".pattern_index() survives a table too wide for any radix", {
  # Eight items over five occasions is forty columns; the radix key overflowed
  # 2^53 at about a dozen, which is why fit_lta() could not use it.
  set.seed(22)
  X <- matrix(stats::rbinom(600 * 40, 1L, 0.5), 600, 40)
  X <- rbind(X, X[1:50, ])
  pat <- .pattern_index(X)
  expect_equal(X[pat$rep_row, , drop = FALSE][pat$idx, ], X)
  expect_lt(sum(pat$rep_row), nrow(X))
})

test_that(".lta_collapse() accepts a repetitive binary panel and declines the rest", {
  X  <- .coll_sim()
  st <- .coll_state(X)

  coll <- .lta_collapse(st, X)
  expect_false(is.null(coll))
  expect_lt(nrow(coll$X), nrow(X) / 2)
  expect_equal(sum(coll$w), nrow(X))
  expect_equal(coll$X[.pattern_index(X)$idx, ], X)

  # Per-case designs cannot be collapsed: two people with the same answers and
  # different covariates are not one row.
  with_z <- st; with_z$Z_delta <- matrix(1, nrow(X), 1)
  expect_null(.lta_collapse(with_z, X))
  with_t <- st; with_t$Z_tau <- matrix(1, nrow(X), 1)
  expect_null(.lta_collapse(with_t, X))
  with_g <- st; with_g$group_info <- list(design = matrix(1, nrow(X), 1))
  expect_null(.lta_collapse(with_g, X))
  with_d <- st; with_d$has_survey_design <- TRUE
  expect_null(.lta_collapse(with_d, X))

  # Nothing to gain when the rows are nearly all distinct.
  expect_null(.lta_collapse(st, matrix(stats::rnorm(nrow(X) * 12), nrow(X))))
})

test_that("the collapsed search returns the fit the pattern table returns", {
  X   <- .coll_sim()
  pat <- .pattern_index(X)
  Xc  <- X[pat$rep_row, , drop = FALSE]
  wc  <- as.vector(rowsum(rep(1, nrow(X)), pat$idx))

  args <- list(n_statuses = 3, times = 2, measurement = "binary",
               n_init = 4, random_state = 5, standard_errors = FALSE)

  # Collapsed internally.
  full <- suppressMessages(suppressWarnings(
    do.call(fit_lta, c(list(X), args))))
  # Handed the same table by the user instead. Its rows are all distinct, so
  # .lta_collapse() declines it and this arm genuinely does not collapse --
  # which is what makes it a control rather than the same computation twice.
  expect_null(.lta_collapse(.coll_state(Xc, weights = wc,
                                        weight_type = "frequency"), Xc))
  hand <- suppressMessages(suppressWarnings(
    do.call(fit_lta, c(list(Xc), args,
                       list(weights = wc, weight_type = "frequency")))))

  # Bit-for-bit, not "close": the two arms see identical starting values and
  # identical weighted sufficient statistics.
  expect_identical(full$loglik,  hand$loglik)
  expect_identical(full$delta_c, hand$delta_c)
  expect_identical(full$tau_c,   hand$tau_c)
  expect_identical(full$mm$models, hand$mm$models)
  expect_identical(full$n_iter,  hand$n_iter)
  expect_identical(full$metrics$bic, hand$metrics$bic)
  expect_identical(full$metrics$entropy, hand$metrics$entropy)
})

test_that("everything per case comes back on the full sample", {
  X <- .coll_sim()
  fit <- suppressMessages(suppressWarnings(
    fit_lta(X, n_statuses = 3, times = 2, measurement = "binary",
            n_init = 3, random_state = 5)))

  # The search ran on a few hundred patterns; the fit describes 1,200 people.
  expect_identical(nrow(fit$data), nrow(X))
  expect_length(fit$ll_case, nrow(X))
  expect_identical(nrow(fit$gamma[[1]]), nrow(X))
  expect_equal(sum(fit$weights_vec), nrow(X))
  expect_equal(fit$loglik, sum(fit$weights_vec * fit$ll_case))
  # Standard errors are built from case-level scores, so they must have been
  # computed after the expansion.
  expect_false(is.null(fit$se))
  expect_true(all(is.finite(unlist(fit$se$delta))))
})
