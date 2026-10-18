# The prior `smoothing` puts on the initial-status and transition regressions:
# smoothing / (K U) pseudo-cases per status and covariate pattern for the
# initial status, smoothing / (K^2 U) per origin, destination and pattern for
# each estimated transition matrix (U distinct covariate patterns).

.lcp_sim <- function(n = 300, Tn = 2, seed = 11, no_return = FALSE,
                     p_item = 0.9) {
  set.seed(seed)
  x <- rbinom(n, 1, 0.5)
  S <- matrix(0L, n, Tn)
  S[, 1] <- 1L + rbinom(n, 1, plogis(0.3 - 0.8 * x))
  for (t in 2:Tn) {
    move12 <- plogis(-1.6 + 0.6 * x)
    move21 <- if (no_return) ifelse(x == 1, 0, plogis(-1.9)) else
      plogis(-1.9 - 0.5 * x)
    u <- runif(n)
    S[, t] <- ifelse(S[, t - 1] == 1L, ifelse(u < move12, 2L, 1L),
                     ifelse(u < move21, 1L, 2L))
  }
  J <- 5
  X <- matrix(0L, n, J * Tn)
  for (t in seq_len(Tn)) for (j in seq_len(J))
    X[, (t - 1) * J + j] <- rbinom(n, 1,
                                   ifelse(S[, t] == 1L, p_item, 1 - p_item))
  list(X = X, x = matrix(x, ncol = 1, dimnames = list(NULL, "x")), Tn = Tn)
}

.lcp_fit <- function(s, effects = "by_origin", inv = "none", smoothing = 1) {
  fit_lta(s$X, n_statuses = 2, times = s$Tn, measurement = "binary",
          measurement_invariance = "full", n_init = 4, random_state = 1,
          n_cores = 1, predictors_initial = s$x, predictors_transition = s$x,
          transition_effects = effects, transition_invariance = inv,
          smoothing = smoothing)
}

# The kernel written from its definition, on one row per covariate pattern,
# sharing no code with the pieces the package prices it from.
.lcp_kernel <- function(fit, alpha = 1) {
  K  <- fit$n_statuses
  u  <- !duplicated(fit$Z_tau)
  U  <- sum(u)
  ld <- .lta_log_delta(fit)
  lt <- .lta_log_tau(fit)
  # One matrix is estimated under "full"; under "slopes" every occasion has
  # its own intercepts, so each occasion's matrix carries its own mass.
  slopes <- isTRUE(fit$tau_occasion_free_intercepts) && fit$n_times > 2
  mats <- if (isTRUE(fit$tau_homogeneous) && !slopes) 1L else seq_along(lt)
  val <- alpha / (K * U) * sum(ld[u, ])
  for (t in mats) for (k in seq_len(K))
    val <- val + alpha / (K^2 * U) * sum(lt[[t]][[k]][u, ])
  val
}

test_that("the covariate prior is the per-pattern kernel in every design", {
  s2 <- .lcp_sim(Tn = 2)
  s3 <- .lcp_sim(Tn = 3, seed = 12)
  cases <- list(list(s2, "by_origin", "none"), list(s2, "common", "none"),
                list(s3, "common", "none"), list(s3, "common", "full"),
                list(s3, "common", "slopes"), list(s3, "by_origin", "full"))
  for (cs in cases) {
    fit <- .lcp_fit(cs[[1]], cs[[2]], cs[[3]])
    expect_equal(.lta_log_prior_dt(fit, 1), .lcp_kernel(fit), tolerance = 1e-10,
                 info = paste(cs[[2]], cs[[3]], cs[[1]]$Tn))
  }
})

test_that("the fit is stationary on the penalised objective, and the penalty's gradient is exact", {
  s   <- .lcp_sim(Tn = 3, seed = 12)
  fit <- .lcp_fit(s, "common", "none")
  X   <- as.matrix(fit$data)
  w   <- fit$weights_vec
  lay <- .lta_par_layout(fit)
  par <- .lta_par_pack(fit, lay)
  lp  <- function(p) .lta_log_prior(.lta_par_unpack(p, fit, lay), X, 1)
  obj <- function(p) sum(w * .lta_ll_case(fit, X, p, lay)) + lp(p)
  fd  <- function(f) vapply(seq_along(par), function(j) {
    e <- numeric(length(par)); e[j] <- 1e-6
    (f(par + e) - f(par - e)) / 2e-6 }, 0)
  expect_lt(max(abs(fd(obj))), 1e-4)
  pen <- .lta_penalty(fit, X, lay, 1)
  expect_equal(pen$value, lp(par), tolerance = 1e-10)
  expect_lt(max(abs(pen$gradient - fd(lp))), 1e-6)
})

test_that("a pattern that never makes a transition keeps a finite coefficient", {
  # Nobody with x = 1 returns from status 2 to status 1, and the items leave
  # little doubt about who is where. Without a prior the
  # by-origin slope runs towards minus infinity; the prior holds it.
  s  <- .lcp_sim(n = 400, Tn = 2, seed = 21, no_return = TRUE, p_item = 0.99)
  on  <- suppressWarnings(.lcp_fit(s, smoothing = 1))
  off <- suppressWarnings(.lcp_fit(s, smoothing = 0))
  slope <- function(fit) {
    B <- fit$tau_beta[[1]]
    max(abs(vapply(B, function(b) b[1, 2], 0)))
  }
  expect_lt(slope(on), 8)
  expect_gt(slope(off), slope(on) + 2)
})

test_that("a one-step fit prints the observed-information standard errors", {
  # The M-step's complete-data Hessian treats the statuses as observed, so
  # its standard errors are too small; the printed table must read se$vcov.
  fit <- .lcp_fit(.lcp_sim(), effects = "common")
  out <- capture.output(lta_covariate_summary(fit, digits = 6))
  blk <- Find(function(b) identical(b$name, "delta_beta"), fit$se$blocks)
  se_v <- sqrt(diag(fit$se$vcov)[blk$cols])
  se_h <- sqrt(diag(pinv(-fit$delta_hessian)))
  row  <- grep("^ *Status 1 +x ", out, value = TRUE)
  expect_length(row, 1L)
  printed <- as.numeric(strsplit(trimws(row), " +")[[1]][5])
  expect_equal(printed, round(se_v[2], 6))
  expect_gt(se_v[2], se_h[2])
})
