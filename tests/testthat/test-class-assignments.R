# class_assignments(): the accessor for the per-case classification.
#
# The contract is that it agrees with the hand-written form it replaces, and
# with get_modal_resp(), which is what the modal three-step correction is a
# table of. Those two must never disagree about a tie.

test_that("modal assignment matches max.col, and 'both' has K + 2 columns", {
  set.seed(1)
  X   <- matrix(rbinom(500, 1, 0.5), nrow = 100)
  fit <- suppressWarnings(fit_mixture(X, n_classes = 2,
                                      measurement = "binary",
                                      n_init = 5,
                                      random_state = 1))

  expect_identical(class_assignments(fit),
                   max.col(fit$log_resp, ties.method = "first"))
  expect_identical(class_assignments(fit, "modal"), class_assignments(fit))

  post <- class_assignments(fit, "posterior")
  expect_equal(dim(post), c(nrow(X), fit$n_components))
  expect_equal(unname(rowSums(post)), rep(1, nrow(X)), tolerance = 1e-10)

  both <- class_assignments(fit, "both")
  expect_s3_class(both, "data.frame")
  expect_identical(ncol(both), as.integer(fit$n_components) + 2L)
  expect_identical(both$class, class_assignments(fit))
  # The probability column is the assigned class's own posterior, so it is the
  # row maximum by construction.
  expect_equal(both$probability, apply(post, 1, max), tolerance = 1e-12)
})

test_that("an LTA assigns a status per occasion", {
  set.seed(22)
  X   <- matrix(rbinom(150 * 3, 1, 0.5), ncol = 3)
  fit <- suppressWarnings(
    fit_lta(X, n_statuses = 2, times = 3, measurement = "binary", n_init = 2,
            random_state = 2, standard_errors = FALSE))

  a <- class_assignments(fit)
  expect_equal(dim(a), c(150L, 3L))
  expect_identical(a[, 2], class_assignments(fit, occasion = 2))

  # One occasion behaves exactly as the cross-sectional case does.
  both <- class_assignments(fit, "both", occasion = 1)
  expect_identical(ncol(both), as.integer(fit$n_statuses) + 2L)
  expect_error(class_assignments(fit, "both"), "single `occasion`")
})
