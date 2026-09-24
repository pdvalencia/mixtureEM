# A distal outcome off the last occasion of the three-step latent transition
# model (fit_lta(n_steps = 3, distal = )). A smoke test on simulated data: the
# outcome's status-specific parameters are recovered, counted and given
# standard errors, the pairwise contrasts run, and the unsupported
# combinations are refused.

.lta_distal_sim <- function() {
  set.seed(11)
  n <- 600; J <- 5
  rho <- rbind(rep(.85, J), rep(.15, J))
  s1 <- sample.int(2, n, TRUE, prob = c(.6, .4))
  tr <- rbind(c(.8, .2), c(.3, .7))
  s2 <- vapply(s1, function(k) sample.int(2, 1, prob = tr[k, ]), integer(1))
  X <- cbind(matrix(rbinom(n * J, 1, rho[s1, ]), n),
             matrix(rbinom(n * J, 1, rho[s2, ]), n))
  distal <- data.frame(score = rnorm(n, c(10, 14)[s2], c(2, 3)[s2]),
                       event = rbinom(n, 1, c(.2, .7)[s2]))
  list(X = X, distal = distal, s2 = s2)
}

test_that("a distal outcome on the last occasion is recovered with SEs", {
  sim <- .lta_distal_sim()
  fit <- suppressWarnings(fit_lta(
    sim$X, n_statuses = 2, times = 2, measurement = "binary",
    n_init = 4, random_state = 1, n_cores = 1, n_steps = 3,
    assignment = "modal", distal = sim$distal))

  # Status labels follow step 1; read them off the recovered means.
  hi <- which.max(fit$mm$distal$mean[, "score"])
  lo <- 3L - hi
  expect_equal(fit$mm$distal$mean[c(lo, hi), "score"], c(10, 14),
               tolerance = 0.05)
  expect_equal(sqrt(fit$mm$distal$var[c(lo, hi), "score"]), c(2, 3),
               tolerance = 0.1)
  expect_equal(fit$mm$distal$prob[c(lo, hi), "event"], c(.2, .7),
               tolerance = 0.2)

  # 1 initial + 2 transition logits + 2 x (mean, variance) + 2 probabilities.
  expect_equal(fit$n_params, 1L + 2L + 4L + 2L)
  expect_equal(nrow(fit$distal), 6L)
  expect_true(all(is.finite(fit$distal$se) & fit$distal$se > 0))

  cs <- outcome_contrasts(fit, outcome = "score")
  expect_equal(nrow(cs), 1L)
  expect_equal(abs(cs$estimate),
               abs(diff(fit$mm$distal$mean[, "score"])), tolerance = 1e-10)
  expect_true(cs$se > 0)
  ce <- outcome_contrasts(fit, outcome = "event")
  expect_true(all(c("OR", "OR_lower", "OR_upper") %in% names(ce)))
})

test_that("distal is refused outside modal three-step", {
  sim <- .lta_distal_sim()
  expect_error(fit_lta(sim$X, n_statuses = 2, times = 2, distal = sim$distal),
               "n_steps = 3")
  expect_error(fit_lta(sim$X, n_statuses = 2, times = 2, n_steps = 3,
                       distal = sim$distal), "modal")
})
