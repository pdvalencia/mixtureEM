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

test_that("the printed note survives a fit where nothing converged", {
  # n_replicated goes NA when every restart in the pool scored NaN (R/lta.R),
  # while n_starts stays a finite count -- the shape that used to crash
  # print.lta_model() with "missing value where TRUE/FALSE needed" instead of
  # reporting anything.
  fit <- list(metrics = list(n_starts = 3L, n_replicated = NA_real_))
  expect_no_error(out <- capture.output(.print_replication_note(fit)))
  expect_match(paste(out, collapse = " "), "none of 3 starts ran to convergence")
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
  # EM alone (`refine = FALSE`). A refined fit is taken on from the cap by the
  # Newton-type finish (R/qn_finish.R), and its flag then says whether that
  # ended at a stationary point.
  expect_warning(
    fit_lta(X, n_statuses = 2, times = 3, measurement = "binary", n_init = 1,
            max_iter = 2, random_state = 1, standard_errors = FALSE,
            refine = FALSE),
    "max_iter = 4", fixed = TRUE)
  fin <- suppressWarnings(
    fit_lta(X, n_statuses = 2, times = 3, measurement = "binary", n_init = 1,
            max_iter = 2, random_state = 1, standard_errors = FALSE))
  if (!is.null(fin$qn_finish))
    expect_identical(fin$converged, fin$qn_finish$converged)
})

# ------------------------------------------------------------------------------
# What counts as the same solution
# ------------------------------------------------------------------------------

# EM leaves restarts on one peak short of its maximum by different amounts, so
# their scores can differ by a tenth of a unit while every case's
# log-likelihood barely moves; a different peak moves the cases by far more.
test_that("restarts are counted on their per-case log-likelihoods", {
  set.seed(31)
  n <- 400
  w <- rep(1, n)
  peak <- -runif(n, 0.5, 3)
  same_a <- peak - 0.1 / n                 # 0.1 below, spread evenly
  same_b <- peak - abs(rnorm(n, 0, 0.004)) # a little noise per case
  other  <- peak + rnorm(n, 0, 0.1)        # another peak, near in total
  expect_identical(.qn_n_replicated(list(same_a, same_b, other), peak, w), 2L)
  # A missing or misshapen signature is no count at all, not a match.
  expect_null(.qn_n_replicated(list(same_a, NULL), peak, w))
  expect_null(.qn_n_replicated(list(same_a, peak[-1]), peak, w))

  # When the finish has to travel far from EM's end point -- a growth model
  # whose maximum lies outside the space EM can reach -- the restarts match the
  # end point the winner was finished from, and that counts as the same peak.
  far    <- peak - abs(rnorm(n, 0, 0.06))    # the finish moved every case
  near_o <- far - abs(rnorm(n, 0, 0.004))    # a restart on the same peak
  expect_identical(.qn_n_replicated(list(far, near_o, other), peak, w), 0L)
  expect_identical(.qn_n_replicated(list(far, near_o, other), peak, w,
                                    origin = far), 2L)
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

test_that("summary() repeats the replication line print() shows", {
  # Hipp & Bauer's diagnostic belongs wherever a reader might start, so both
  # summary() methods carry the same sentence, read off the same counts.
  best_line <- function(out) grep("Best solution", out, value = TRUE)

  set.seed(23)
  X <- matrix(rbinom(120 * 4, 1, 0.5), ncol = 4)
  fit <- fit_mixture(X, n_classes = 2, measurement = "binary", n_init = 3)
  p <- best_line(capture.output(print(fit)))
  expect_length(p, 1L)
  expect_identical(best_line(capture.output(summary(fit))), p)

  set.seed(22)
  X <- matrix(rbinom(150 * 3, 1, 0.5), ncol = 3)
  lta <- fit_lta(X, n_statuses = 2, times = 3, measurement = "binary",
                 n_init = 4, random_state = 2, standard_errors = FALSE)
  p <- best_line(capture.output(print(lta)))
  expect_length(p, 1L)
  expect_identical(best_line(capture.output(summary(lta))), p)
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
  expect_true(any(grepl("Model Selection Summary", out)))
  # Printed as a mark on the log-likelihood with a footnote, not a column.
  res$fit_table$Unreplicated <- c(FALSE, TRUE)
  out <- utils::capture.output(print(res))
  expect_true(any(grepl("^ +2 +-[0-9.]+\\*", out)))
  expect_true(any(grepl("Unreplicated maximum at K = 2", out, fixed = TRUE)))
})

# ------------------------------------------------------------------------------
# The spike check: the score with the single best-fitted case removed
# ------------------------------------------------------------------------------

test_that(".k_deleted_ll() removes the best-fitted case, in part when weighted", {
  d <- c(-1.2, -0.3, -2.5, -0.9)
  expect_equal(.k_deleted_ll(d, rep(1, 4)), sum(d) - max(d))
  # A pattern of frequency f is f identical cases: deleting one removes one
  # unit of its weight, which is the fractional boundary term.
  w <- c(2, 3, 1, 4)
  expect_equal(.k_deleted_ll(d, w), .k_deleted_ll(rep(d, w), rep(1, sum(w))))
  expect_equal(.k_deleted_ll(d, c(0.4, 0.5, 1, 1)),
               sum(c(0.4, 0.5, 1, 1) * d) - (0.5 * -0.3 + 0.5 * -0.9))
  # Any score built on the log-likelihood can be passed as the total.
  expect_equal(.k_deleted_ll(d, rep(1, 4), total = -10), -10 - max(d))
  expect_true(is.na(.k_deleted_ll(d[1], 1)))
  expect_equal(.k_deleted_ll(c(d, NaN, -Inf), c(rep(1, 4), 1, 1)),
               sum(d) - max(d))
})

test_that("every fit carries the spike check, and a clean fit stays silent", {
  set.seed(31)
  X <- matrix(rnorm(300 * 3, rep(c(-2, 0, 2), each = 100)), ncol = 3)
  fit <- fit_mixture(X, n_classes = 3, measurement = "continuous", n_init = 5,
                     random_state = 1)
  expect_true(is.finite(fit$metrics$kdel_ll))
  expect_false(fit$metrics$kdel_disagrees)
  expect_false(any(grepl("Spike check", capture.output(print(fit)))))

  # The growth models reach it through fit_em() as fit_mixture() does.
  set.seed(32)
  cls <- rbinom(200, 1, 0.5)
  Y   <- outer(ifelse(cls == 1, 1, -1), rep(1, 4)) +
         outer(ifelse(cls == 1, 0.8, -0.2), 0:3) + matrix(rnorm(800, sd = 0.5), 200)
  g <- suppressWarnings(fit_gmm(Y, n_classes = 2, times = 4, n_init = 3,
                                random_state = 1))
  expect_true(is.finite(g$metrics$kdel_ll))
  l <- fit_lcga(Y, times = 4, n_classes = 2, family = "gaussian", n_init = 3,
                random_state = 1)
  expect_true(is.finite(l$metrics$kdel_ll))

  set.seed(22)
  B <- matrix(rbinom(150 * 3, 1, 0.5), ncol = 3)
  lta <- fit_lta(B, n_statuses = 2, times = 3, measurement = "binary",
                 n_init = 4, random_state = 2, standard_errors = FALSE)
  expect_true(is.finite(lta$metrics$kdel_ll))
  expect_true(lta$metrics$kdel_best >= lta$metrics$kdel_ll)
})

test_that("the spike check fires when one case decides the ranking", {
  # One class too many on 40 cases, fitted without the variance prior: the
  # winner has a class shrunk onto a point inside the data, and another start
  # is better once each loses its best-fitted case.
  set.seed(303)
  z <- rbinom(40, 1, 0.4)
  d <- data.frame(y = ifelse(z == 1, rnorm(40, 6.25, sqrt(0.27)),
                             rnorm(40, 4.33, sqrt(0.14))))
  fit <- suppressWarnings(fit_mixture(d, 3, measurement = "continuous",
           variances_equal = FALSE, bayes_constants = list(variances = 0, latent = 0),
           n_init = 100, random_state = 3))
  # The search first, so a failure below is the flag's and not a moved maximum.
  expect_equal(fit$metrics$ll, -42.9960, tolerance = 1e-4)
  expect_true(fit$metrics$kdel_disagrees)
  expect_output(print(fit), "Spike check")
})
