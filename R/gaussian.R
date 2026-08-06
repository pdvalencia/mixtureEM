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
#'
#' @return A list object containing the model state.
#' @export
gaussian_model <- function(n_components, type = "gaussian_unit") {
  state <- list(
    n_components = n_components,
    parameters = list()
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
  n <- nrow(X)
  log_eps <- matrix(0, nrow = n, ncol = model_state$n_components)
  for (c in seq_len(model_state$n_components)) {
    mean_c <- matrix(model_state$parameters$means[c, ], nrow = n, ncol = ncol(X), byrow = TRUE)
    ll_matrix <- dnorm(X, mean = mean_c, sd = 1, log = TRUE)
    log_eps[, c] <- rowSums(ll_matrix)
  }
  return(log_eps)
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

  means <- matrix(0, nrow = model_state$n_components, ncol = ncol(X))
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
  n <- nrow(X)
  log_eps <- matrix(0, nrow = n, ncol = model_state$n_components)
  for (c in seq_len(model_state$n_components)) {
    mean_c <- matrix(model_state$parameters$means[c, ], nrow = n, ncol = ncol(X), byrow = TRUE)
    ll_matrix <- dnorm(X, mean = mean_c, sd = 1, log = TRUE)
    ll_matrix[is.na(ll_matrix)] <- 0 # FIML: missing cells drop out of the sum
    log_eps[, c] <- rowSums(ll_matrix)
  }
  return(log_eps)
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

  covariances <- matrix(0, nrow = model_state$n_components, ncol = ncol(X))
  for (c in seq_len(model_state$n_components)) {
    diff_sq <- sweep(X, 2, means[c, ], "-")^2
    covariances[c, ] <- (colSums(resp[, c] * diff_sq) + prior_obs * s2) /
      (sum_resp[c] + prior_obs)
  }
  # Numerical guard only, and reachable only at alpha = 0, where the user has
  # asked for exactly this behaviour. It stops an exactly-zero variance from
  # producing a non-finite density; it is not a statistical regularisation and
  # is far too small to act as one.
  model_state$parameters$covariances <- pmax(covariances, 1e-12)
  return(model_state)
}

#' @exportS3Method
log_likelihood.gaussian_diag <- function(model_state, X, ...) {
  n <- nrow(X)
  log_eps <- matrix(0, nrow = n, ncol = model_state$n_components)
  for (c in seq_len(model_state$n_components)) {
    mean_c <- matrix(model_state$parameters$means[c, ], nrow = n, ncol = ncol(X), byrow = TRUE)
    sd_c <- matrix(sqrt(model_state$parameters$covariances[c, ]), nrow = n, ncol = ncol(X), byrow = TRUE)

    ll_matrix <- dnorm(X, mean = mean_c, sd = sd_c, log = TRUE)
    log_eps[, c] <- rowSums(ll_matrix)
  }
  return(log_eps)
}

#' @exportS3Method
n_parameters.gaussian_diag <- function(model_state, ...) {
  return(length(model_state$parameters$means) + length(model_state$parameters$covariances))
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

  means <- matrix(0, nrow = model_state$n_components, ncol = ncol(X))
  covariances <- matrix(0, nrow = model_state$n_components, ncol = ncol(X))

  for (j in seq_len(ncol(X))) {
    valid <- !is.na(X[, j])
    if (any(valid)) {
      resp_valid <- resp[valid, , drop=FALSE]
      sum_resp <- colSums(resp_valid)

      # Means
      means[, j] <- t(resp_valid) %*% X[valid, j] / sum_resp

      # Variances
      for (c in seq_len(model_state$n_components)) {
        diff_sq <- (X[valid, j] - means[c, j])^2
        covariances[c, j] <- (sum(resp_valid[, c] * diff_sq) + prior_obs * s2[j]) /
          (sum_resp[c] + prior_obs)
      }
    } else {
      # No observed cell for this item: the data say nothing, so the prior is
      # the whole of the posterior and the update returns its centre.
      covariances[, j] <- s2[j]
    }
  }
  model_state$parameters$means <- means
  model_state$parameters$covariances <- pmax(covariances, 1e-12)
  return(model_state)
}

#' @exportS3Method
log_likelihood.gaussian_diag_nan <- function(model_state, X, ...) {
  n <- nrow(X)
  log_eps <- matrix(0, nrow = n, ncol = model_state$n_components)
  for (c in seq_len(model_state$n_components)) {
    mean_c <- matrix(model_state$parameters$means[c, ], nrow = n, ncol = ncol(X), byrow = TRUE)
    sd_c <- matrix(sqrt(model_state$parameters$covariances[c, ]), nrow = n, ncol = ncol(X), byrow = TRUE)

    ll_matrix <- dnorm(X, mean = mean_c, sd = sd_c, log = TRUE)
    ll_matrix[is.na(ll_matrix)] <- 0 # FIML Masking
    log_eps[, c] <- rowSums(ll_matrix)
  }
  return(log_eps)
}
