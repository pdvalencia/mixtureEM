# blrt() under sampling weights and a survey design: the i.i.d. replicates are
# rescaled by the null model's design effect instead of being refused.

.blrt_survey_sim <- function(n = 300) {
  cls <- 1L + (runif(n) > 0.55)
  p   <- rbind(c(.85, .80, .80, .75, .80), c(.20, .25, .15, .25, .20))
  X   <- t(vapply(cls, function(k) rbinom(5, 1, p[k, ]), numeric(5)))
  list(X = X, w = ifelse(cls == 2L, 2.5, 1) * exp(rnorm(n, 0, 0.4)),
       st = rep(1:2, each = n / 2), cl = rep(seq_len(n / 10), each = 10))
}

test_that("an unweighted BLRT is not rescaled", {
  set.seed(31)
  d <- .blrt_survey_sim()
  res <- suppressWarnings(blrt(d$X, 2, 3, measurement = "binary", n_reps = 3,
                               n_init_base = 3, n_init_boot = 1,
                               verbose = FALSE))
  expect_false(res$scaled)
  expect_true(is.na(res$scaling_factor))
  # All-one sampling weights are no design at all.
  set.seed(31)
  d <- .blrt_survey_sim()
  ones <- suppressWarnings(blrt(d$X, 2, 3, measurement = "binary", n_reps = 3,
                                n_init_base = 3, n_init_boot = 1,
                                verbose = FALSE, weights = rep(1, 300)))
  expect_false(ones$scaled)
  expect_equal(ones$null_dist, res$null_dist)
})

test_that("a weighted, clustered BLRT rescales its draws by the null fit's design effect", {
  set.seed(32)
  d <- .blrt_survey_sim()
  set.seed(5)
  res <- suppressWarnings(blrt(d$X, 2, 3, measurement = "binary", n_reps = 3,
                               n_init_base = 3, n_init_boot = 1,
                               verbose = FALSE, weights = d$w,
                               strata = d$st, cluster = d$cl))
  expect_true(res$scaled)

  # The factor is the k_small fit's own tr(A^-1 B) / p: the same first fit,
  # from the same random-number state.
  set.seed(5)
  f0 <- fit_mixture_internal(d$X, n_components = 2, measurement = "binary",
                             n_init = 3, n_cores = 1L, weights = d$w,
                             strata = d$st, cluster = d$cl)
  expect_equal(f0$metrics$ll, res$ll_small)
  c0 <- .scaling_pieces(.nested_fit_info(f0))$c
  expect_equal(res$scaling_factor, c0)
  expect_gt(c0, 1)
  expect_equal(res$p_value,
               (sum(res$null_dist >= res$obs_diff) + 1) / (res$n_reps + 1))
  expect_output(print(res), "Design effect")
})
