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
