# ==============================================================================
# The multi-start report: counts, the warning, and the guards on it
# ==============================================================================
#
# The predicate and the warning are exercised directly wherever a real fit would
# only be a slow way of building the same list. What a fit is needed for is the
# plumbing: that the counts reach `metrics` at all, and that the comparison
# functions stay quiet.

# ------------------------------------------------------------------------------
# The predicate
# ------------------------------------------------------------------------------

test_that("a lone replication counts only when enough starts were requested", {
  # The threshold is on what was asked for, not on what converged: a staged
  # search runs three survivors however many restarts it was given.
  expect_true(.is_unreplicated(list(n_replicated = 1L, n_requested = 20L,
                                    n_starts = 3L)))
  expect_false(.is_unreplicated(list(n_replicated = 1L, n_requested = 5L,
                                     n_starts = 5L)))
  expect_false(.is_unreplicated(list(n_replicated = 6L, n_requested = 20L,
                                     n_starts = 20L)))
  # No counts at all is unknown, not clean.
  expect_true(is.na(.is_unreplicated(list())))
})

# ------------------------------------------------------------------------------
# The warning and the two guards
# ------------------------------------------------------------------------------

test_that("an unreplicated maximum warns, names 100, and reports both counts", {
  fit <- list(metrics = list(n_replicated = 1L, n_requested = 50L,
                             n_starts = 3L))
  w <- tryCatch(.check_replication(fit), warning = function(w) w)
  expect_s3_class(w, "mixtureEM_replication")
  expect_match(conditionMessage(w), "n_init = 100", fixed = TRUE)
  expect_match(conditionMessage(w), "1 of 3 starts that ran to convergence")
  expect_match(conditionMessage(w), "out of 50 requested")
})

test_that("the advice scales with how many starts were asked for", {
  # Telling a user who ran n_init = 200 to refit with 100 is worse than saying
  # nothing, so above the threshold the message changes what it recommends.
  low <- list(metrics = list(n_replicated = 1L, n_requested = 20L,
                             n_starts = 20L))
  w <- tryCatch(.check_replication(low), warning = function(w) w)
  expect_match(conditionMessage(w), "n_init = 100", fixed = TRUE)

  high <- list(metrics = list(n_replicated = 1L, n_requested = 200L,
                              n_starts = 200L))
  w2 <- tryCatch(.check_replication(high), warning = function(w) w)
  expect_false(grepl("n_init = 100", conditionMessage(w2), fixed = TRUE))
  expect_match(conditionMessage(w2), "about the specification")

  # The printed note reads off the same helper.
  expect_match(paste(capture.output(.print_replication_note(low)),
                     collapse = " "),
               "n_init = 100", fixed = TRUE)
  expect_false(grepl("n_init = 100",
                     paste(capture.output(.print_replication_note(high)),
                           collapse = " "), fixed = TRUE))
})

test_that("the specification reading waits for starts that were run out", {
  # The strong reading is a claim about restarts that reached convergence. A
  # staged search refines a fraction of what it was asked for, so a rule on the
  # requested count alone would assert a hundred-start test that never ran.
  expect_match(.replication_advice(50L, 50L), "n_init = 100", fixed = TRUE)

  staged <- .replication_advice(200L, 40L)
  expect_match(staged, "40 restarts run out to convergence", fixed = TRUE)
  expect_match(staged, "200 requested", fixed = TRUE)
  expect_match(staged, "raise n_init further", fixed = TRUE)
  expect_false(grepl("about the specification", staged, fixed = TRUE))

  full <- .replication_advice(200L, 200L)
  expect_match(full, "about the specification", fixed = TRUE)
  # The claim is hedged rather than asserted: more starts can still be the fix.
  expect_match(full, "More starts can still help", fixed = TRUE)

  # A caller with only the requested count is treated as an unstaged search.
  expect_equal(.replication_advice(200L), full)

  # Singular where there is one survivor, and the counts come off the fit.
  fit <- list(metrics = list(n_replicated = 1L, n_requested = 150L,
                             n_starts = 1L))
  w <- tryCatch(.check_replication(fit), warning = function(w) w)
  expect_match(conditionMessage(w), "1 restart run out to convergence",
               fixed = TRUE)
})

test_that("a fit below the threshold stays silent", {
  fit <- list(metrics = list(n_replicated = 1L, n_requested = 5L,
                             n_starts = 5L))
  expect_no_warning(.check_replication(fit))
})

test_that("a collapsed variance or a growth boundary suppresses it", {
  # Both of those warnings say raising n_init can make matters worse, and this
  # one says raise it. They must never both fire on one fit.
  base <- list(metrics = list(n_replicated = 1L, n_requested = 50L,
                              n_starts = 3L))

  degenerate <- base
  degenerate$degenerate <- list(classes = 1L)
  expect_no_warning(.check_replication(degenerate))

  boundary <- base
  boundary$growth <- list(boundary = list("residual variance at occasion 5"))
  expect_no_warning(.check_replication(boundary))
})

# ------------------------------------------------------------------------------
# Non-convergence
# ------------------------------------------------------------------------------

test_that("the non-convergence warning names the doubled iteration budget", {
  w <- tryCatch(.warn_non_convergence(1000L), warning = function(w) w)
  expect_match(conditionMessage(w), "max_iter = 2000", fixed = TRUE)
})

test_that("fit_lta() warns when it stops at the iteration cap", {
  set.seed(21)
  X <- matrix(rbinom(120 * 3, 1, 0.5), ncol = 3)
  expect_warning(
    fit_lta(X, n_statuses = 2, times = 3, measurement = "binary", n_init = 1,
            max_iter = 2, random_state = 1, standard_errors = FALSE),
    "max_iter = 4", fixed = TRUE)
})

# ------------------------------------------------------------------------------
# The counts reach the fit
# ------------------------------------------------------------------------------

test_that("fit_lta() records how many restarts found its solution", {
  set.seed(22)
  X <- matrix(rbinom(150 * 3, 1, 0.5), ncol = 3)
  fit <- fit_lta(X, n_statuses = 2, times = 3, measurement = "binary",
                 n_init = 4, random_state = 2, standard_errors = FALSE)
  expect_equal(fit$metrics$n_requested, 4L)
  expect_true(fit$metrics$n_starts >= 1L)
  expect_true(fit$metrics$n_replicated <= fit$metrics$n_starts)
})

test_that("fit_mixture() carries the requested count alongside the converged one", {
  set.seed(23)
  X <- matrix(rbinom(120 * 4, 1, 0.5), ncol = 4)
  fit <- fit_mixture(X, n_classes = 2, measurement = "binary", n_init = 3)
  expect_equal(fit$metrics$n_requested, 3)
  expect_equal(fit$metrics$n_starts, 3)
})

# ------------------------------------------------------------------------------
# The comparison tables
# ------------------------------------------------------------------------------

test_that("compare_longitudinal() reports replication in the table, not by warning", {
  set.seed(24)
  X <- matrix(rbinom(150 * 3, 1, 0.5), ncol = 3)
  expect_no_warning(
    res <- compare_longitudinal(X, k_range = 2L, model = "lta", times = 3,
                                measurement = "binary", n_init = 2,
                                standard_errors = FALSE, verbose = FALSE))
  expect_true("Unreplicated" %in% names(res$fit_table))
  expect_type(res$fit_table$Unreplicated, "logical")
})

test_that("the BLRT warns when 100 draws cannot resolve the decision", {
  # Within one Monte Carlo step of .05 at 100 draws, and outside it.
  w <- tryCatch(.blrt_check_granularity(0.0495, 100L), warning = function(w) w)
  expect_s3_class(w, "warning")
  expect_match(conditionMessage(w), "n_reps = 999", fixed = TRUE)
  expect_no_warning(.blrt_check_granularity(0.4, 100L))
  # At 999 draws there is no larger number to recommend, so it goes quiet.
  expect_no_warning(.blrt_check_granularity(0.0495, 999L))
})

test_that("the BLRT warns on a negative bootstrap statistic", {
  w <- tryCatch(.blrt_check_negative(3L, 100L, 2L, 3L), warning = function(w) w)
  expect_match(conditionMessage(w), "n_init_boot = 50", fixed = TRUE)
  expect_no_warning(.blrt_check_negative(0L, 100L, 2L, 3L))
})

test_that("compare_mixtures() gains the same column and still prints", {
  set.seed(25)
  X <- matrix(rbinom(150 * 4, 1, 0.5), ncol = 4)
  out <- utils::capture.output(
    res <- compare_mixtures(X, k_range = 1:2, measurement = "binary",
                            n_init = 2))
  expect_true("Unreplicated" %in% names(res$fit_table))
  # A logical column used to break the rounding on the way to the console.
  expect_true(any(grepl("Unreplicated", out)))
})
