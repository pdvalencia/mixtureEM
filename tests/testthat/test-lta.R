# ==============================================================================
# Latent transition analysis - the single-chain forward-backward engine
# ==============================================================================
#
# Structural checks that need no reference program: internal consistency of
# the forward-backward recursions, weighting, and restriction arithmetic.

# ------------------------------------------------------------------------------
# Structural checks, which need no reference program
# ------------------------------------------------------------------------------

test_that("frequency weights agree with the expanded case-level data", {
  set.seed(11)
  pat <- as.matrix(expand.grid(rep(list(0:1), 3)))
  w   <- c(40, 12, 9, 15, 11, 7, 6, 30)
  expanded <- pat[rep(seq_len(nrow(pat)), w), , drop = FALSE]

  a <- fit_lta(pat, n_statuses = 2, times = 3, measurement = "binary",
               weights = w, weight_type = "frequency", n_init = 10,
               random_state = 3, standard_errors = FALSE)
  b <- fit_lta(expanded, n_statuses = 2, times = 3, measurement = "binary",
               n_init = 10, random_state = 3, standard_errors = FALSE)

  expect_equal(a$loglik, b$loglik, tolerance = 1e-6)
  expect_equal(a$n_params, b$n_params)
})

test_that("the forward-backward posteriors are proper distributions", {
  set.seed(12)
  X <- matrix(rbinom(400 * 4, 1, 0.5), ncol = 4)
  fit <- fit_lta(X, n_statuses = 2, times = 4, measurement = "binary",
                 n_init = 5, random_state = 4, standard_errors = FALSE)

  for (g in fit$gamma) expect_equal(unname(rowSums(g)), rep(1, nrow(X)),
                                    tolerance = 1e-8)
  # Each pairwise table sums to the effective sample size, being a joint
  # distribution over consecutive statuses summed over cases.
  for (xi in fit$xi) expect_equal(sum(xi), nrow(X), tolerance = 1e-6)
  # And its row sums are the occasion's status counts.
  for (t in seq_along(fit$xi))
    expect_equal(unname(rowSums(fit$xi[[t]])),
                 unname(colSums(fit$gamma[[t]])), tolerance = 1e-6)
})

test_that("forbidden transitions hold the ruled-out cells at exactly zero", {
  set.seed(13)
  X <- matrix(rbinom(300 * 4, 1, 0.5), ncol = 4)
  # No return from status 2 to status 1: an absorbing second stage.
  fit <- fit_lta(X, n_statuses = 2, times = 4, measurement = "binary",
                 forbidden_transitions = matrix(c(0, 0, 1, 0), 2, 2,
                                                byrow = TRUE),
                 n_init = 5, random_state = 5, standard_errors = FALSE)

  for (m in fit$tau) expect_equal(m[2, 1], 0)
  for (m in fit$tau) expect_equal(unname(rowSums(m)), c(1, 1), tolerance = 1e-10)
  # The forbidden cell is not charged a parameter: 1 (delta) + 3 x 1 (row 1 of
  # each transition matrix, row 2 being fixed) + 2 (the one item's response
  # probabilities, held equal across occasions by the default invariance).
  expect_equal(fit$n_params, 1L + 3L + 2L)
})

# ------------------------------------------------------------------------------
# Two-level indicators that are not coded 0/1
# ------------------------------------------------------------------------------

test_that("1/2-coded indicators fit the same model as their 0/1 recoding", {
  set.seed(21)
  n   <- 300
  st  <- sample(1:2, n, TRUE)
  p   <- rbind(c(0.8, 0.75, 0.85), c(0.2, 0.25, 0.15))
  X01 <- cbind(matrix(rbinom(n * 3, 1, p[st, ]), n, 3),
               matrix(rbinom(n * 3, 1, p[st, ]), n, 3))
  colnames(X01) <- paste0("i", rep(1:3, 2), "@T", rep(1:2, each = 3))

  a <- fit_lta(X01, n_statuses = 2, times = 2, measurement = "binary",
               n_init = 5, random_state = 7, standard_errors = FALSE)
  expect_message(
    b <- fit_lta(X01 + 1, n_statuses = 2, times = 2, measurement = "binary",
                 n_init = 5, random_state = 7, standard_errors = FALSE),
    "Recoded binary indicators")

  expect_equal(a$loglik, b$loglik, tolerance = 1e-8)
})

test_that("the binary recode pools the occasions rather than deciding per column", {
  set.seed(22)
  n  <- 200
  X  <- matrix(rbinom(n * 4, 1, 0.5), n, 4) + 1
  colnames(X) <- paste0("i", rep(1:2, 2), "@T", rep(1:2, each = 2))
  # Item 1 happens to observe only the upper level at the second occasion. A
  # per-column recode would map that 2 to 0, contradicting the first occasion.
  X[, 3] <- 2

  expect_silent(
    fit <- suppressMessages(
      fit_lta(X, n_statuses = 2, times = 2, measurement = "binary",
              n_init = 3, random_state = 8, standard_errors = FALSE)))
  expect_true(is.finite(fit$loglik))

  rec <- suppressMessages(.recode_binary_blocks(X, n_items = 2, n_blocks = 2)$X)
  expect_equal(unname(rec[, 3]), rep(1, n))
  expect_true(all(rec %in% c(0, 1)))
})

# ------------------------------------------------------------------------------
# The two priors are separate: `smoothing` for the chain, `bayes_constants` for
# the measurement model
# ------------------------------------------------------------------------------
#
# Four binary items, the fourth of them constant within status, so unsmoothed ML
# drives its response probabilities to exactly 0 and 1 and the measurement prior
# is the only thing holding them off the boundary.
.lta_prior_fixture <- function() {
  set.seed(11)
  n <- 120; J <- 4
  s1 <- sample(1:2, n, TRUE, c(.6, .4))
  s2 <- ifelse(runif(n) < .8, s1, 3 - s1)
  mk <- function(s) {
    m <- matrix(0, n, J)
    for (j in 1:J) {
      p <- ifelse(s == 1, c(.95, .95, .95, 1.0)[j], c(.05, .05, .05, 0.0)[j])
      m[, j] <- rbinom(n, 1, p)
    }
    m
  }
  X <- cbind(mk(s1), mk(s2))
  colnames(X) <- c(paste0("i", 1:J, "_t1"), paste0("i", 1:J, "_t2"))
  X
}

.lta_rho <- function(f)
  sort(as.numeric(unlist(lapply(f$mm$models, function(m) m$parameters$pis))))

test_that("bayes_constants reaches the measurement model", {
  X  <- .lta_prior_fixture()
  go <- function(...) fit_lta(X, n_statuses = 2, times = 2, n_init = 3,
                              random_state = 1, standard_errors = FALSE, ...)

  a <- .lta_rho(go())
  b <- .lta_rho(go(bayes_constants = list(categorical = 5)))

  # `smoothing` used to be passed into the measurement model's M-step, where it
  # won the fallback to `bayes_constants$categorical` and made the documented
  # argument inert: these two were identical to the last digit.
  expect_false(isTRUE(all.equal(a, b)))
  expect_lt(max(b), max(a))
})

test_that("smoothing does not supply the measurement model's prior", {
  X <- .lta_prior_fixture()
  f <- fit_lta(X, n_statuses = 2, times = 2, n_init = 3, random_state = 1,
               smoothing = 0, standard_errors = FALSE)
  r <- .lta_rho(f)

  # `smoothing = 0` asks for unsmoothed transitions. It used to strip the
  # measurement prior with them, sending rho to 1 and to 3.6e-10 - the boundary
  # estimates that make interval estimates for those parameters meaningless.
  # The default `bayes_constants$categorical` now holds them off it.
  expect_lt(max(r), 1 - 1e-6)
  expect_gt(min(r), 1e-4)
  # The transitions themselves are unsmoothed, which is what was asked for.
  expect_lt(min(vapply(f$tau, min, numeric(1))),
            min(vapply(fit_lta(X, n_statuses = 2, times = 2, n_init = 3,
                               random_state = 1, standard_errors = FALSE)$tau,
                       min, numeric(1))))
})
