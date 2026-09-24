# ==============================================================================
# A three-step LTA fit's log-likelihood is step 3's: the likelihood of the
# assigned statuses under fixed classification errors, not of the items. It is
# on a different scale from a one- or two-step fit's, so every surface that sets
# likelihoods side by side either refuses the pair or reads step 1 instead.
# ==============================================================================

.tsc_cache <- new.env(parent = emptyenv())

.tsc_sim <- function() {
  if (!is.null(.tsc_cache$X)) return(.tsc_cache$X)
  set.seed(11)
  K <- 2; Tn <- 3; J <- 4; n <- 300
  rho  <- rbind(c(.85, .80, .90, .75), c(.15, .20, .10, .25))
  prev <- rbind(c(.7, .3), c(.5, .5), c(.3, .7))
  X <- matrix(NA_integer_, n, J * Tn)
  for (t in seq_len(Tn)) {
    s <- sample.int(K, n, TRUE, prob = prev[t, ])
    for (j in seq_len(J))
      X[, (t - 1L) * J + j] <- rbinom(n, 1, rho[s, j])
  }
  .tsc_cache$X <- X
  X
}

.tsc_fit <- function(n_steps, ...) {
  key <- paste(n_steps, ..., sep = "_")
  if (!is.null(.tsc_cache[[key]])) return(.tsc_cache[[key]])
  .tsc_cache[[key]] <- fit_lta(.tsc_sim(), n_statuses = 2, times = 3,
                               measurement = "binary", n_init = 2,
                               random_state = 1, n_cores = 1,
                               standard_errors = FALSE, n_steps = n_steps, ...)
}

test_that("lr_test() refuses a three-step fit beside a one-step one", {
  one   <- .tsc_fit(1)
  three <- .tsc_fit(3, assignment = "modal")
  expect_error(lr_test(three, one), "One model is a three-step fit")
  expect_error(lr_test(one, three), "One model is a three-step fit")
})

test_that("lr_test() refuses two three-step fits on different step-3 data", {
  modal <- .tsc_fit(3, assignment = "modal")
  prop  <- .tsc_fit(3, assignment = "proportional")
  expect_error(lr_test(modal, prop), "different step-1 results")
})

test_that("print() reports step 1's criteria for a three-step fit", {
  three <- .tsc_fit(3, assignment = "modal")
  out <- paste(capture.output(print(three)), collapse = "\n")
  expect_match(out, "Log-Likelihood (Step 1)", fixed = TRUE)
  expect_match(out, sprintf("%.2f", three$step1$metrics$ll), fixed = TRUE)
  expect_false(grepl(sprintf("%.2f", three$metrics$ll), out, fixed = TRUE))

  one <- .tsc_fit(1)
  expect_false(grepl("(Step 1)", paste(capture.output(print(one)),
                                       collapse = "\n"), fixed = TRUE))
})

test_that("compare_longitudinal(n_steps = 3) tabulates step 1 on every row", {
  cmp <- compare_longitudinal(.tsc_sim(), k_range = 2:3, model = "lta",
                              times = 3, verbose = FALSE,
                              measurement = "binary", n_init = 2,
                              random_state = 1, n_cores = 1,
                              standard_errors = FALSE, n_steps = 3,
                              assignment = "modal")
  expect_equal(cmp$fit_table$LL,
               unname(vapply(cmp$models, function(m) m$step1$metrics$ll,
                             numeric(1))))
  expect_equal(cmp$fit_table$Params,
               unname(vapply(cmp$models, function(m)
                 as.numeric(m$step1$metrics$n_params), numeric(1))))
})
