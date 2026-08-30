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
.ml <- list(categorical = 0, latent = 0)
