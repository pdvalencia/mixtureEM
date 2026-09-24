# ==============================================================================
# The step-two classification-error matrix
# ==============================================================================
#
# `.classification_error()` is the only genuinely new arithmetic in the
# bias-adjusted three-step estimator: the table that says how often a case
# belonging to class `c` gets filed under class `r`. Everything downstream
# treats that table as known, so an error here is invisible and propagates.
#
# The main check is a published worked example rather than a simulation.
# Nylund-Gibson, Grimm, Quirk and Furlong (2014), Appendix B, print the inputs
# (average posterior probabilities by assigned class, and the class counts),
# the resulting classification probabilities, and the fixed logits Appendix C
# then feeds to the third step. Reproducing those numbers grades the formula
# and -- more importantly -- the orientation of the table, which is the part
# that fails silently.

# ------------------------------------------------------------------------------
# The published example
# ------------------------------------------------------------------------------

# Average latent class probabilities for most likely class membership (row) by
# latent class (column), Appendix B. The paper's note says cells that were zero
# became .0001; its own worked denominator for column 1 uses .0001 where the
# printed table shows .001, so the [2, 1] cell is written that way here.
.appendix_b_A <- matrix(
  c(0.842, 0.001, 0.039, 0.064, 0.054,
    1e-4,  0.934, 0.061, 0.005, 0.001,
    0.029, 0.068, 0.822, 0.081, 0.001,
    0.051, 0.006, 0.106, 0.817, 0.020,
    0.067, 0.001, 0.001, 0.020, 0.913),
  nrow = 5, byrow = TRUE)

.appendix_b_N <- c(280, 734, 520, 367, 271)

# One case per person, each carrying the average posterior vector of the class
# it was filed under. Every row is diagonal-dominant, so modal assignment puts
# it back in that row and the weighted sums collapse to the paper's own
# arithmetic.
.appendix_b_resp <- function() {
  .appendix_b_A[rep(seq_len(5), times = .appendix_b_N), , drop = FALSE]
}

# Columns 1 and 4 of the published q table. Columns 2, 3 and 5 are deliberately
# absent: they cannot be reproduced from the paper's own printed inputs --
# back-solving column 2's denominator from its published q needs an input cell
# of 5.2e-6 where the table prints .001 -- so asserting on them would be
# asserting on a typesetting error. Do not "fix" this by adding them.
.appendix_b_q1 <- c(0.819216, 0.000255, 0.052400, 0.065038, 0.063092)
.appendix_b_q4 <- c(0.048568, 0.009947, 0.114156, 0.81264,  0.01469)

# The same two columns as fixed logits, from Appendix C's third-step input,
# which prints them to nine digits. Four rows, the fifth assigned category
# being the reference.
.appendix_c_logit1 <- c(2.563758177, -5.510887505, -0.185686799, 0.030376041)
.appendix_c_logit4 <- c(1.195821592, -0.389904153,  2.050426872, 4.013149848)


test_that("the published classification probabilities reproduce", {
  ce <- .classification_error(.appendix_b_resp(), assignment = "modal")

  expect_equal(dim(ce$D), c(5L, 5L))
  expect_lt(max(abs(ce$D[, 1] - .appendix_b_q1)), 1e-6)
  expect_lt(max(abs(ce$D[, 4] - .appendix_b_q4)), 1e-6)
})

test_that("the published fixed logits reproduce", {
  ce <- .classification_error(.appendix_b_resp(), assignment = "modal")

  expect_lt(max(abs(ce$logits[1:4, 1] - .appendix_c_logit1)), 1e-6)
  expect_lt(max(abs(ce$logits[1:4, 4] - .appendix_c_logit4)), 1e-6)

  # The last assigned category is the reference, so its own logit is zero.
  expect_equal(unname(ce$logits[5, ]), rep(0, 5))
})

# ------------------------------------------------------------------------------
# Orientation
# ------------------------------------------------------------------------------

test_that("columns are the true class and sum to one", {
  ce <- .classification_error(.appendix_b_resp(), assignment = "modal")
  expect_equal(unname(colSums(ce$D)), rep(1, 5))

  # Not symmetric, which is why the next test exists at all.
  expect_gt(max(abs(ce$D - t(ce$D))), 0.01)
})

test_that("D is the transpose of the cross-sectional correction's table", {
  # fit_ml() builds the same numbers the other way up, rows the true class
  # (R/corrections.R). Under proportional assignment both expressions use the
  # posteriors on both sides, so the two must agree exactly after a transpose.
  set.seed(11)
  n    <- 200
  raw  <- matrix(stats::runif(n * 3), n, 3)^3
  resp <- raw / rowSums(raw)
  w    <- stats::runif(n, 0.5, 2)

  ce <- .classification_error(resp, assignment = "proportional", weights = w)

  C_prop     <- t(resp * w) %*% resp
  Nk         <- colSums(resp * w)
  C_row_norm <- sweep(C_prop, 1, Nk, "/")

  expect_lt(max(abs(ce$D - t(C_row_norm))), 1e-12)
})

# ------------------------------------------------------------------------------
# The zero-cell floor
# ------------------------------------------------------------------------------

# Three cases, three classes, arranged so that the case assigned to class 3
# carries no posterior mass in class 1 at all. The [3, 1] cell is then exactly
# zero, and it is the reference row, so every logit in that column depends on
# what the floor is.
.zero_cell_resp <- function() {
  matrix(c(0.9, 0.1, 0.0,
           0.1, 0.9, 0.0,
           0.0, 0.2, 0.8), nrow = 3, byrow = TRUE)
}

test_that("a zero cell is floored and leaves the logits finite", {
  ce <- .classification_error(.zero_cell_resp(), assignment = "modal",
                              zero_floor = 1e-6)

  expect_lt(abs(ce$D[3, 1] - 1e-6), 1e-8)
  expect_true(all(is.finite(ce$logits)))
  expect_equal(unname(colSums(ce$D)), rep(1, 3))
})

test_that("the floor is the only thing that moves a floored column's logits", {
  fine   <- .classification_error(.zero_cell_resp(), assignment = "modal",
                                  zero_floor = 1e-6)
  coarse <- .classification_error(.zero_cell_resp(), assignment = "modal",
                                  zero_floor = 1e-4)

  # Column 1's reference cell is the floored one, so raising the floor by a
  # factor of 100 lowers every logit in the column by exactly log(100) -- the
  # probabilities move by about 1e-4 and the parameter table by 4.6. That is
  # why the value is an argument and not a constant.
  expect_lt(max(abs((fine$logits[1:2, 1] - coarse$logits[1:2, 1]) - log(100))),
            1e-6)

  # A column with no floored cell does not notice.
  expect_lt(max(abs(fine$D[, 2] - coarse$D[, 2])), 1e-6)

  expect_error(.classification_error(.zero_cell_resp(), zero_floor = 0))
})

# ------------------------------------------------------------------------------
# Weights and the assignment rule
# ------------------------------------------------------------------------------

test_that("integer weights match replicating the rows", {
  resp <- .zero_cell_resp()
  w    <- c(2, 1, 3)

  weighted   <- .classification_error(resp, assignment = "modal", weights = w)
  replicated <- .classification_error(resp[rep(1:3, times = w), ],
                                      assignment = "modal")

  expect_lt(max(abs(weighted$D - replicated$D)), 1e-12)
})

test_that("the proportional table is the softer of the two", {
  set.seed(24)
  n    <- 400
  cl   <- sample(1:3, n, replace = TRUE)
  raw  <- matrix(stats::runif(n * 3, 0.05, 0.35), n, 3)
  raw[cbind(seq_len(n), cl)] <- stats::runif(n, 1, 3)
  resp <- raw / rowSums(raw)

  modal <- .classification_error(resp, assignment = "modal")
  prop  <- .classification_error(resp, assignment = "proportional")

  # Modal assignment throws away the uncertainty in each case's class, so its
  # table sits closer to the identity than the proportional one does.
  expect_true(all(diag(modal$D) > diag(prop$D)))
  expect_equal(unname(colSums(prop$D)), rep(1, 3))
})

test_that("a fixed modal assignment is used as given", {
  # The three-step variance differentiates the table with the assignments held
  # where step 2 put them; a posterior nudged across a boundary must not move
  # its case to another row.
  resp <- rbind(c(.6, .3, .1), c(.2, .7, .1), c(.1, .2, .7), c(.45, .44, .11))
  a    <- max.col(resp)
  expect_equal(.classification_error(resp, "modal", assigned = a),
               .classification_error(resp, "modal"))
  nudged <- resp
  nudged[4, ] <- c(.44, .45, .11)
  held <- .classification_error(nudged, "modal", assigned = a)$D
  free <- .classification_error(nudged, "modal")$D
  expect_false(isTRUE(all.equal(held, free)))
  # Row 1 of `held` still holds case 4's mass: column 2 gains it from case 4.
  expect_gt(held[1, 2], free[1, 2])
})
