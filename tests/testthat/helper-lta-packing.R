# Shared fixtures for the files that exercise the LTA's unconstrained parameter
# packing -- test-lta-refine.R (the L-BFGS refinement) and test-lta-scaling.R
# (the MLR scaling factor and the robust standard errors). Both drive the same
# packing through the same small synthetic fit, so the simulator lives here
# rather than in either file: testthat gives each test file its own environment
# and does not share definitions between them, but it sources every `helper-`
# file into the environment they all inherit from.

# A three-status, three-occasion, four-item binary LTA with a sticky transition
# matrix. Small enough that a finite-difference Hessian over its parameter
# vector costs seconds.
.lta_refine_sim <- function(n = 250, K = 3, Tn = 3, J = 4, seed = 4) {
  set.seed(seed)
  rho <- matrix(stats::runif(K * J, 0.15, 0.85), K, J)
  s <- sample(K, n, TRUE)
  out <- list()
  for (t in seq_len(Tn)) {
    if (t > 1) s <- ifelse(stats::runif(n) < 0.7, s, sample(K, n, TRUE))
    out[[t]] <- matrix(stats::rbinom(n * J, 1, rho[s, ]), n, J)
  }
  X <- do.call(cbind, out)
  colnames(X) <- paste0("t", rep(seq_len(Tn), each = J), "_i",
                        rep(seq_len(J), Tn))
  X
}

# Priors off, so the fitted objective is the plain log-likelihood. Every
# identity these two files assert -- the gradient, the fixed point, the
# case-level likelihood -- is stated on the unpenalised objective.
.ml <- list(categorical = 0, latent = 0, variances = 0)

# A three-status, two-occasion, three-item continuous LTA with well-separated
# status means, so the Gaussian variance blocks stay well away from the
# degenerate near-zero-variance region.
.lta_gaussian_refine_sim <- function(n = 300, K = 3, Tn = 2, J = 3, seed = 21) {
  set.seed(seed)
  mu <- matrix(seq(-3, 3, length.out = K), K, J) +
    matrix(stats::runif(K * J, -0.3, 0.3), K, J)
  sd0 <- matrix(stats::runif(K * J, 0.6, 1.2), K, J)
  s <- sample(K, n, TRUE)
  out <- list()
  for (t in seq_len(Tn)) {
    if (t > 1) s <- ifelse(stats::runif(n) < 0.75, s, sample(K, n, TRUE))
    out[[t]] <- matrix(stats::rnorm(n * J, mu[s, ], sd0[s, ]), n, J)
  }
  X <- do.call(cbind, out)
  colnames(X) <- paste0("t", rep(seq_len(Tn), each = J), "_i",
                        rep(seq_len(J), Tn))
  X
}

# A covariate-driven LTA where the covariate actually predicts delta and tau,
# with a moderate true effect. Unlike bolting an unrelated random covariate
# onto .lta_refine_sim()'s covariate-free data, this keeps the MLE well away
# from separation: delta_beta/tau_beta carry no prior of their own (by
# design -- see .lta_penalty()/.lta_log_prior()), so an underidentified or
# noise-only covariate drives the fitted coefficients to the +/-25 box the
# L-BFGS refinement clips to, or past the point where a softmax underflows to
# its 1e-300 floor -- both of which make a finite-difference gradient check
# fail for reasons that are about the fixture, not the score formula.
.lta_cov_refine_sim <- function(n = 300, K = 2, Tn = 3, J = 3, seed = 11) {
  set.seed(seed)
  z <- stats::rnorm(n)
  delta1 <- stats::plogis(0.3 + 0.6 * z)
  s <- rbinom(n, 1, delta1) + 1L
  tau_stay <- function(k, z) stats::plogis(0.4 + (if (k == 1) 0.5 else -0.5) * z)
  S <- matrix(0L, n, Tn); S[, 1] <- s
  for (t in 2:Tn) {
    p_stay <- ifelse(S[, t - 1] == 1, tau_stay(1, z), tau_stay(2, z))
    u <- stats::runif(n)
    S[, t] <- ifelse(u < p_stay, S[, t - 1], 3L - S[, t - 1])
  }
  rho <- matrix(stats::runif(K * J, 0.25, 0.75), K, J)
  X <- array(0L, dim = c(n, Tn, J))
  for (t in seq_len(Tn)) for (j in seq_len(J))
    X[, t, j] <- stats::rbinom(n, 1, rho[S[, t], j])
  Xwide <- do.call(cbind, lapply(seq_len(Tn), function(t) X[, t, ]))
  colnames(Xwide) <- paste0("t", rep(seq_len(Tn), each = J), "_i",
                            rep(seq_len(J), Tn))
  list(X = Xwide, Z = data.frame(z = z))
}

# Generated from RI-LTA itself, at the published article's own Monte Carlo
# settings (Muthen & Asparouhov 2022, Table 3/4 design): 5 binary items, a
# continuous factor loading of 2 on every item, item logits +1/-1 by status,
# delta = (.5, .5), tau row 1 = (.622, .378), row 2 = (.5, .5).
.lta_ri_sim <- function(n = 2000, Tn = 3, J = 5, seed = 99) {
  set.seed(seed)
  K <- 2L
  delta <- c(0.5, 0.5)
  tau   <- matrix(c(0.622, 0.378, 0.500, 0.500), K, K, byrow = TRUE)
  A <- matrix(c(1, -1), K, J)
  L <- rep(2, J)
  z <- stats::rnorm(n)

  S <- matrix(0L, n, Tn)
  S[, 1] <- sample.int(K, n, replace = TRUE, prob = delta)
  for (t in 2:Tn)
    for (k in seq_len(K)) {
      i <- S[, t - 1] == k
      if (any(i)) S[i, t] <- sample.int(K, sum(i), replace = TRUE, prob = tau[k, ])
    }

  X <- matrix(NA_real_, n, Tn * J)
  for (t in seq_len(Tn)) for (j in seq_len(J)) {
    eta <- A[S[, t], j] + L[j] * z
    X[, (t - 1) * J + j] <- stats::rbinom(n, 1, stats::plogis(eta))
  }
  colnames(X) <- paste0("t", rep(seq_len(Tn), each = J), "_i", rep(seq_len(J), Tn))
  list(X = X, z = z, S = S, delta = delta, tau = tau)
}

# A small ordinal (cumulative-logit) LTA, ragged category counts by item
# (### 14.18, W9's own test fixture) -- two 3-category items and one binary
# item, matching the article's own Dating-data shape without its size or
# extreme thresholds, so a finite-difference check over it costs seconds.
# The same ordinal shape, generated from the cumulative-logit model itself so
# that a known direct covariate effect is in the data: item j in status k has
# free thresholds and gains `x_i * dif[k, j]` on its linear predictor, shared
# across occasions -- exactly what `predictors_items` estimates. One binary
# covariate, so there are two covariate patterns and the emission cost is two
# tables per node instead of one.
.lta_dif_sim <- function(n = 600, Tn = 3, K = 2L, cats = c(3L, 3L, 2L),
                         dif = NULL, seed = 3) {
  set.seed(seed)
  J <- length(cats)
  if (is.null(dif)) {
    dif <- matrix(0, K, J)
    dif[1, 1] <- 0.8            # item 1 behaves differently in status 1 only
  }
  # Thresholds on the P(U >= s) scale, decreasing across s within an item, and
  # separated by status so the statuses are recoverable.
  th <- lapply(seq_len(J), function(j) {
    base <- seq(0.8, -0.8, length.out = cats[j] - 1L)
    # matrix(), not t(vapply()): a two-category item has one threshold, and
    # vapply() would hand back a vector that t() turns into a 1 x K matrix.
    matrix(rep(base, each = K) + rep(ifelse(seq_len(K) == 1L, 0.9, -0.9),
                                     times = cats[j] - 1L),
           K, cats[j] - 1L)
  })
  x <- stats::rbinom(n, 1L, 0.5)
  delta <- rep(1 / K, K)
  tau   <- matrix(0.2 / (K - 1L), K, K); diag(tau) <- 0.8

  S <- matrix(0L, n, Tn)
  S[, 1] <- sample.int(K, n, replace = TRUE, prob = delta)
  for (t in 2:Tn) for (k in seq_len(K)) {
    i <- S[, t - 1] == k
    if (any(i)) S[i, t] <- sample.int(K, sum(i), replace = TRUE, prob = tau[k, ])
  }

  X <- matrix(NA_integer_, n, Tn * J)
  for (t in seq_len(Tn)) for (j in seq_len(J)) {
    Sj <- cats[j]
    for (i in seq_len(n)) {
      k <- S[i, t]
      F_ <- c(1, stats::plogis(th[[j]][k, ] + x[i] * dif[k, j]), 0)
      p  <- F_[seq_len(Sj)] - F_[seq_len(Sj) + 1L]
      X[i, (t - 1L) * J + j] <- sample.int(Sj, 1L, prob = p)
    }
  }
  colnames(X) <- paste0("t", rep(seq_len(Tn), each = J), "_i",
                        rep(seq_len(J), Tn))
  list(X = X, Z = matrix(x, ncol = 1L, dimnames = list(NULL, "x")),
       dif = dif, S = S)
}

.lta_ordinal_sim <- function(n = 150, Tn = 3, K = 2, cats = c(3L, 3L, 2L),
                             seed = 1) {
  set.seed(seed)
  J <- length(cats)
  delta <- c(0.5, 0.5)
  tau <- matrix(c(0.8, 0.2, 0.2, 0.8), K, K, byrow = TRUE)
  X <- matrix(NA_integer_, n, Tn * J)
  for (i in seq_len(n)) {
    s <- sample(K, 1, prob = delta)
    for (t in seq_len(Tn)) {
      if (t > 1) s <- sample(K, 1, prob = tau[s, ])
      for (j in seq_len(J)) {
        Sj <- cats[j]
        probs <- if (s == 1) rev(seq_len(Sj)) else seq_len(Sj)
        probs <- probs / sum(probs)
        X[i, (t - 1) * J + j] <- sample(Sj, 1, prob = probs)
      }
    }
  }
  colnames(X) <- paste0("t", rep(seq_len(Tn), each = J), "_i", rep(seq_len(J), Tn))
  X
}
