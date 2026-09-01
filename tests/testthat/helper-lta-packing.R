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
