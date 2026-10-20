# A per-covariate `transition_effects` (R/lta_covariates.R, .lta_tau_cross()):
# the "common" transition regression with the covariates named "by_origin"
# crossed with the origin dummies. The external grade is RECORDS "C4"; these
# are the identities that need no outside number, and the refusals.

.mix_sim <- function() {
  s <- .lta_cov_refine_sim(n = 300, K = 2, Tn = 2, seed = 11)
  set.seed(5)
  s$Z$g <- factor(sample(c("a", "b", "c"), 300, TRUE))
  s
}

test_that("names resolve to design columns, factors by their variable name", {
  Z <- cbind(Intercept = 1, z = 0, gb = 0, gc = 0)
  expect_null(.lta_tau_cross(NULL, Z))
  expect_null(.lta_tau_cross(c(z = "common"), Z))
  expect_identical(.lta_tau_cross(c(g = "by_origin"), Z),
                   c(z = FALSE, gb = TRUE, gc = TRUE))
  expect_identical(.lta_tau_cross(c(z = "common", .default = "by_origin"), Z),
                   c(z = FALSE, gb = TRUE, gc = TRUE))
  expect_error(.lta_tau_cross(c(w = "by_origin"), Z), "not among")
  expect_error(.lta_tau_cross(c(z = "both"), Z), "named vector")
  expect_error(.lta_tau_cross(c(z = "by_origin"), NULL), "predictors_transition")
})

test_that("one covariate named by_origin is the by-origin model; all common is common", {
  s   <- .mix_sim()
  arg <- list(s$X, n_statuses = 2, times = 2, measurement = "binary",
              predictors_transition = s$Z["z"], smoothing = 0,
              bayes_constants = list(categorical = 0), n_init = 3,
              random_state = 1, n_cores = 1, max_iter = 2000,
              standard_errors = FALSE)
  f_by  <- do.call(fit_lta, c(arg, list(transition_effects = "by_origin")))
  f_mix <- do.call(fit_lta, c(arg, list(transition_effects = c(z = "by_origin"))))
  expect_lt(abs(f_by$loglik - f_mix$loglik), 1e-4)
  expect_identical(f_by$metrics$n_params, f_mix$metrics$n_params)
  expect_true("from:1 x z" %in% colnames(.lta_tau_design(f_mix, 1L)))

  f_c  <- do.call(fit_lta, c(arg, list(transition_effects = "common")))
  f_c2 <- do.call(fit_lta, c(arg, list(transition_effects = c(z = "common"))))
  expect_identical(f_c$loglik, f_c2$loglik)
})

test_that("a by-origin covariate beside a common one counts its own slopes", {
  s <- .mix_sim()
  f <- fit_lta(s$X, n_statuses = 2, times = 2, measurement = "binary",
               predictors_transition = s$Z, transition_effects = c(g = "by_origin"),
               n_init = 2, random_state = 1, n_cores = 1, max_iter = 100,
               standard_errors = FALSE)
  cn <- colnames(.lta_tau_design(f, 1L))
  expect_identical(cn, c("Intercept", "from:1", "z", "g.b", "g.c",
                         "from:1 x g.b", "from:1 x g.c"))
  expect_output(lta_covariate_summary(f), "from:1 x g.b")
  expect_error(fit_lta(s$X, n_statuses = 2, times = 2, measurement = "binary",
                       predictors_transition = s$Z["z"], group = s$Z$g,
                       transition_effects = c(z = "by_origin"), n_init = 1),
               "cannot be combined")
})
