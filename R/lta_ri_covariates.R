# ==============================================================================
# Latent transition analysis - covariates on the continuous random intercept
# ==============================================================================
#
# f_i | x_i ~ N(x_i'beta, 1): the article's own `f ON x`. The residual variance
# stays fixed at 1 and there is no intercept -- a constant shift of the factor
# is absorbed exactly by the free per-class thresholds, so a location shift of
# `f` is not identified and a constant column is refused below.
#
# The fixed Gauss-Hermite grid is never shifted per case (that would make the
# shared emission tables person-specific and the fit unrunnable); instead the
# node prior is reweighted by the exact importance identity
# phi(z - mu) / phi(z) = exp(mu*z - mu^2/2), so
#
#   log mass[i, q] = log w_q + mu_i * z_q - mu_i^2 / 2.
#
# The prior part of the expected complete-data log-likelihood is then
# sum_i w_i (mu_i * fbar_i - mu_i^2 / 2), with fbar_i the posterior factor
# mean, so its maximiser is ordinary weighted least squares -- no optimiser.

# ------------------------------------------------------------------------------
# Design matrix
# ------------------------------------------------------------------------------

# No intercept column: the residual mean of the factor is fixed at 0, so an
# intercept is not identified. This is the one deliberate difference from
# .lta_design() (R/lta_covariates.R), which always prepends
# cbind(Intercept = 1, Z).
.lta_ri_design <- function(predictors, n, label) {
  if (is.null(predictors)) return(NULL)
  Z <- prepare_covariates(.as_named_covariates(predictors, NULL, label))
  Z <- complete_covariates(as.matrix(Z))
  if (nrow(Z) != n)
    stop(sprintf("`%s` must have one row per case.", label), call. = FALSE)
  if (ncol(Z) == 0L)
    stop(sprintf("`%s` must have at least one column.", label), call. = FALSE)
  const <- apply(Z, 2L, function(col) length(unique(col[!is.na(col)])) <= 1L)
  if (any(const))
    stop(paste0(
      "`", label, "` must not contain a constant column: the random ",
      "intercept's residual mean is fixed at zero, so there is no intercept ",
      "to estimate and a constant column is not identified."), call. = FALSE)
  Z
}

# ------------------------------------------------------------------------------
# The person-specific node prior
# ------------------------------------------------------------------------------

# n x Q matrix of log node masses. Without covariates every row is the same
# fixed Gauss-Hermite weight vector; with them, the exact importance-
# reweighting identity for f_i ~ N(x_i'beta, 1) on the residual grid.
.lta_ri_log_mass <- function(state, n) {
  base <- log(pmax(state$ri$mass, 1e-300))
  M <- matrix(base, n, length(base), byrow = TRUE)
  if (is.null(state$Z_ri) || is.null(state$ri_beta)) return(M)
  mu <- as.vector(state$Z_ri %*% state$ri_beta)
  M + outer(mu, state$ri$Dnode[, 1L]) - mu^2 / 2
}

# ------------------------------------------------------------------------------
# The M-step: closed-form weighted least squares
# ------------------------------------------------------------------------------

# Exact ECM for beta. The prior part of the expected complete-data
# log-likelihood is sum_i w_i (mu_i fbar_i - mu_i^2 / 2) with fbar_i the
# posterior factor mean, so the maximiser is ordinary weighted least squares.
# No optimiser, no .fit_mnl(): this is a linear model in closed form.
.lta_ri_mstep_beta <- function(state, E) {
  if (is.null(state$Z_ri)) return(state)
  Z    <- state$Z_ri
  w    <- state$weights_vec
  fbar <- as.vector(E$ri$pq %*% state$ri$Dnode[, 1L])
  Zw   <- Z * w
  b <- tryCatch(solve(crossprod(Z, Zw), crossprod(Zw, fbar)),
                error = function(e) state$ri_beta)
  state$ri_beta <- matrix(as.vector(b), ncol = 1L)
  state
}
