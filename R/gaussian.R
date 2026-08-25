# ==============================================================================
# S3 Gaussian Models (Continuous Data)
# ==============================================================================

#' Constructor for Gaussian models
#'
#' @description
#' Sets up the initial state and class structure for continuous emission models
#' (like Gaussian with diagonal or unit variance).
#'
#' @param n_components Integer. The number of latent classes/components to estimate.
#' @param type Character. The specific variance structure, e.g., "gaussian_diag" or "gaussian_unit".
#' @param variances_equal Logical. Hold each item's variance equal across the
#'   classes, so the classes differ in location only (the homoscedastic latent
#'   profile model, and the default parameterisation of several commercial
#'   programs). Ignored by the unit-variance types, which have no variances to
#'   estimate. The estimated variance is still stored once per class — a
#'   \code{K x J} matrix with identical rows — so every reader of the parameters
#'   (likelihood, plotting, alignment, boundary checks) is unchanged; the
#'   constraint lives in the M-step and in \code{n_parameters()}.
#'
#' @return A list object containing the model state.
#' @export
gaussian_model <- function(n_components, type = "gaussian_unit",
                           variances_equal = FALSE) {
  state <- list(
    n_components    = n_components,
    variances_equal = isTRUE(variances_equal),
    parameters      = list()
  )
  class(state) <- c(type, "emission")
  return(state)
}

# ------------------------------------------------------------------------------
# 1. Gaussian Unit (Fixed Variance = 1) S3 Methods
# ------------------------------------------------------------------------------

#' @exportS3Method
init_params.gaussian_unit <- function(model_state, X, resp, random_state = NULL, ...) {
  if (!is.null(random_state)) set.seed(random_state)
  idx <- sample.int(nrow(X), model_state$n_components)
  model_state$parameters$means <- X[idx, , drop = FALSE]
  return(model_state)
}

#' @exportS3Method
m_step.gaussian_unit <- function(model_state, X, resp, weights = NULL, ...) {
  if (!is.null(weights)) resp <- sweep(resp, 1, weights, "*")

  means <- t(resp) %*% X
  sum_resp <- colSums(resp)
  model_state$parameters$means <- sweep(means, 1, sum_resp, "/")
  return(model_state)
}

#' @exportS3Method
log_likelihood.gaussian_unit <- function(model_state, X, ...) {
  mu  <- model_state$parameters$means
  out <- X %*% t(mu) +
         rep(-0.5 * rowSums(mu * mu) - 0.5 * ncol(X) * log(2 * pi), each = nrow(X)) -
         0.5 * rowSums(X * X)
  dimnames(out) <- NULL
  return(out)
}

#' @exportS3Method
n_parameters.gaussian_unit <- function(model_state, ...) {
  return(length(model_state$parameters$means))
}

# ------------------------------------------------------------------------------
# 2. Gaussian Unit NaN (Fixed Variance = 1, Missing Data FIML) S3 Methods
# ------------------------------------------------------------------------------

#' @exportS3Method
init_params.gaussian_unit_nan <- function(model_state, X, resp, random_state = NULL, ...) {
  if (!is.null(random_state)) set.seed(random_state)
  idx   <- sample.int(nrow(X), model_state$n_components)
  means <- X[idx, , drop = FALSE]
  # Seed rows may themselves contain NA; falling back to the column mean over
  # observed cases keeps every component's starting mean well defined and stops
  # an NA seed from masking out otherwise-valid likelihood contributions.
  if (anyNA(means)) {
    col_means <- colMeans(X, na.rm = TRUE)
    na_cells  <- which(is.na(means), arr.ind = TRUE)
    means[na_cells] <- col_means[na_cells[, "col"]]
  }
  model_state$parameters$means <- means
  return(model_state)
}

#' @exportS3Method
m_step.gaussian_unit_nan <- function(model_state, X, resp, weights = NULL, ...) {
  if (!is.null(weights)) resp <- sweep(resp, 1, weights, "*")

  # dimnames carried explicitly: this M-step allocates rather than deriving the
  # result from a matrix product, so without them the item names are lost and
  # summaries fall back to "Item_1, Item_2, ...".
  means <- matrix(0, nrow = model_state$n_components, ncol = ncol(X),
                  dimnames = list(NULL, colnames(X)))
  for (j in seq_len(ncol(X))) {
    valid <- !is.na(X[, j])
    if (any(valid)) {
      resp_valid <- resp[valid, , drop = FALSE]
      means[, j] <- t(resp_valid) %*% X[valid, j] / colSums(resp_valid)
    }
  }
  model_state$parameters$means <- means
  return(model_state)
}

#' @exportS3Method
log_likelihood.gaussian_unit_nan <- function(model_state, X, ...) {
  mu  <- model_state$parameters$means
  X0  <- X
  X0[is.na(X0)] <- 0
  M   <- (!is.na(X)) * 1
  out <- X0 %*% t(mu) + M %*% t(-0.5 * mu * mu - 0.5 * log(2 * pi)) -
         0.5 * rowSums(X0 * X0)
  dimnames(out) <- NULL
  return(out)
}

#' @exportS3Method n_parameters gaussian_unit_nan
n_parameters.gaussian_unit_nan <- n_parameters.gaussian_unit

# ------------------------------------------------------------------------------
# 3. Gaussian Diag (Estimated feature-specific variance) S3 Methods
# ------------------------------------------------------------------------------

#' @exportS3Method
init_params.gaussian_diag <- function(model_state, X, resp, random_state = NULL, ...) {
  model_state <- init_params.gaussian_unit(model_state, X, resp, random_state)
  model_state$parameters$covariances <- matrix(1, nrow = model_state$n_components, ncol = ncol(X))
  return(model_state)
}

# Weighted marginal variance of each column, over that column's observed cells.
#
# This is the centre of the prior on the class variances: the variance item j
# would have if the classes told us nothing about it. Centring on the observed
# marginal rather than on a fixed constant is what makes the prior scale-free —
# it says "a class variance is not far below the item's overall spread", which
# is a statement about the item, not about the units it happens to be measured
# in.
#
# The marginal is *weighted*, matching the weights the M-step itself uses. An
# unweighted marginal would centre the prior somewhere the sufficient statistics
# never go, so the same data supplied as a frequency-weighted table and as
# expanded case-level rows would not give the same answer.
#
# A column with no usable spread (all cells missing, zero total weight, or a
# constant) falls back to 1. Such a column carries no information about any
# class anyway; the fallback exists so the prior stays finite rather than to
# say anything substantive.
# Turn class-wise sums of squares into the class variances, applying the
# homoscedastic constraint if the emission carries one.
#
# `ss` and `n` are K x J (or K-vectors recycled over the columns): the weighted
# sum of squares about each class's own mean, and the weight behind it.
# `prior_obs` is the alpha/K pseudo-observations at the marginal `s2` that the
# free update already used.
#
# Under `variances_equal` the J parameters are shared by the K classes, so the
# constrained optimum sums both the data and the prior over the classes: each
# class still contributes its own alpha/K pseudo-observations to the objective,
# and they now all bear on one parameter, which is K * (alpha/K) = alpha in
# total. Writing the result back into all K rows keeps the stored shape
# rectangular, which is what lets every other method ignore the constraint.
.pool_variances_over_classes <- function(ss, n, prior_obs, s2, variances_equal) {
  K <- nrow(ss)
  n <- if (is.matrix(n)) n else matrix(n, nrow = K, ncol = ncol(ss))
  if (!isTRUE(variances_equal)) {
    out <- (ss + prior_obs * rep(s2, each = K)) / (n + prior_obs)
  } else {
    pooled <- (colSums(ss) + K * prior_obs * s2) / (colSums(n) + K * prior_obs)
    out <- matrix(pooled, nrow = K, ncol = ncol(ss), byrow = TRUE)
  }
  out
}

.marginal_var <- function(X, weights = NULL) {
  w <- if (is.null(weights)) rep(1, nrow(X)) else weights
  vapply(seq_len(ncol(X)), function(j) {
    ok <- !is.na(X[, j])
    if (!any(ok)) return(1)
    wj <- w[ok]
    sw <- sum(wj)
    if (!is.finite(sw) || sw <= 0) return(1)
    xj <- X[ok, j]
    mu <- sum(wj * xj) / sw
    v  <- sum(wj * (xj - mu)^2) / sw
    if (!is.finite(v) || v <= 0) 1 else v
  }, numeric(1))
}

#' @exportS3Method
m_step.gaussian_diag <- function(model_state, X, resp, weights = NULL, ...) {
  if (!is.null(weights)) resp <- sweep(resp, 1, weights, "*")

  # 1. Update Means
  means <- t(resp) %*% X
  sum_resp <- colSums(resp)
  means <- sweep(means, 1, sum_resp, "/")
  model_state$parameters$means <- means

  # 2. Update Variances
  #
  # The posterior mode under a truncated inverse-Wishart prior centred on the
  # item's observed marginal variance, expressed as alpha/K pseudo-observations
  # per class:
  #
  #     sigma^2_kj = (SS_kj + (alpha/K) s^2_j) / (n_k + alpha/K)
  #
  # This replaces the additive `+ 1e-6` floor that used to sit here, and the
  # difference is not cosmetic. A floor is a constant in the units of the data,
  # so it is invisible on an item scored 1-5 and enormous on one scored in
  # thousands; worse, it is applied *after* the maximisation, so nothing stops
  # the objective from wanting a smaller value, and refine_lbfgs() — which
  # re-optimises log(sd) with no floor at all — would walk straight through it
  # and climb an unbounded likelihood. The prior is part of the objective
  # instead, so both stages agree about what is being maximised and the
  # likelihood the polish climbs actually has a maximum.
  #
  # The prior is the *truncated* inverse-Wishart: the usual
  # -((K+1)/2) log|Sigma| normalising term is deliberately omitted. That term is
  # what would keep pulling variances upward even at alpha = 0; without it,
  # alpha = 0 recovers the unpenalised ML update exactly, which is what makes
  # `bayes_constants = list(variances = 0)` a usable escape hatch rather than an
  # approximate one.
  alpha     <- .bayes_alpha(model_state, "variances")
  prior_obs <- alpha / model_state$n_components
  s2        <- .marginal_var(X, weights)

  rx <- t(resp) %*% X                       # this is what `means` was built from
  ss <- t(resp) %*% (X * X) - 2 * means * rx + means^2 * sum_resp
  covariances <- .pool_variances_over_classes(ss, sum_resp, prior_obs, s2,
                                              model_state$variances_equal)
  # Numerical guard only, and reachable only at alpha = 0, where the user has
  # asked for exactly this behaviour. It stops an exactly-zero variance from
  # producing a non-finite density; it is not a statistical regularisation and
  # is far too small to act as one.
  model_state$parameters$covariances <- pmax(covariances, 1e-12)
  return(model_state)
}

#' @exportS3Method
log_likelihood.gaussian_diag <- function(model_state, X, ...) {
  mu  <- model_state$parameters$means
  v   <- model_state$parameters$covariances
  Cd  <- -0.5 * (log(2 * pi * v) + mu * mu / v)
  out <- -0.5 * ((X * X) %*% t(1 / v)) + X %*% t(mu / v) +
         rep(rowSums(Cd), each = nrow(X))
  dimnames(out) <- NULL
  return(out)
}

#' @exportS3Method
n_parameters.gaussian_diag <- function(model_state, ...) {
  n_var <- length(model_state$parameters$covariances)
  # Stored per class even when shared by them; only the free ones are counted.
  if (isTRUE(model_state$variances_equal))
    n_var <- n_var / model_state$n_components
  return(length(model_state$parameters$means) + n_var)
}

# ------------------------------------------------------------------------------
# 4. Gaussian Diag NaN (Missing Data FIML) S3 Methods
# ------------------------------------------------------------------------------

#' @exportS3Method
init_params.gaussian_diag_nan <- function(model_state, X, resp, random_state = NULL, ...) {
  model_state <- init_params.gaussian_unit_nan(model_state, X, resp, random_state)
  model_state$parameters$covariances <- matrix(1, nrow = model_state$n_components, ncol = ncol(X))
  return(model_state)
}

#' @exportS3Method n_parameters gaussian_diag_nan
n_parameters.gaussian_diag_nan <- n_parameters.gaussian_diag

#' @exportS3Method
m_step.gaussian_diag_nan <- function(model_state, X, resp, weights = NULL, ...) {
  if (!is.null(weights)) resp <- sweep(resp, 1, weights, "*")

  # Same prior as the complete-data M-step above; see the comment there for the
  # form and for why alpha = 0 recovers ML. Under FIML both the prior centre and
  # the class counts it competes against are computed *per column* over that
  # column's observed cells, so an item observed on a tenth of the sample gets a
  # prior calibrated to the tenth that saw it rather than to the whole.
  alpha     <- .bayes_alpha(model_state, "variances")
  prior_obs <- alpha / model_state$n_components
  s2        <- .marginal_var(X, weights)

  K <- model_state$n_components
  means <- matrix(0, nrow = K, ncol = ncol(X),
                  dimnames = list(NULL, colnames(X)))
  ss    <- matrix(0, nrow = K, ncol = ncol(X))
  nk    <- matrix(0, nrow = K, ncol = ncol(X))
  empty <- logical(ncol(X))

  for (j in seq_len(ncol(X))) {
    valid <- !is.na(X[, j])
    if (any(valid)) {
      resp_valid <- resp[valid, , drop=FALSE]
      sum_resp <- colSums(resp_valid)
      nk[, j] <- sum_resp

      # Means
      means[, j] <- t(resp_valid) %*% X[valid, j] / sum_resp

      # Sums of squares about each class's own mean
      for (c in seq_len(K))
        ss[c, j] <- sum(resp_valid[, c] * (X[valid, j] - means[c, j])^2)
    } else {
      # No observed cell for this item: the data say nothing, so the prior is
      # the whole of the posterior and the update returns its centre.
      empty[j] <- TRUE
    }
  }
  covariances <- .pool_variances_over_classes(ss, nk, prior_obs, s2,
                                              model_state$variances_equal)
  if (any(empty)) covariances[, empty] <- rep(s2[empty], each = K)
  model_state$parameters$means <- means
  model_state$parameters$covariances <- pmax(covariances, 1e-12)
  return(model_state)
}

#' @exportS3Method
log_likelihood.gaussian_diag_nan <- function(model_state, X, ...) {
  mu  <- model_state$parameters$means
  v   <- model_state$parameters$covariances
  Cd  <- -0.5 * (log(2 * pi * v) + mu * mu / v)
  X0  <- X
  X0[is.na(X0)] <- 0
  M   <- (!is.na(X)) * 1
  out <- -0.5 * ((X0 * X0) %*% t(1 / v)) + X0 %*% t(mu / v) + M %*% t(Cd)
  dimnames(out) <- NULL
  return(out)
}
