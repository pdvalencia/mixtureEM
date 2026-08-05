# ==============================================================================
# S3 Poisson Models (Count Data)
# ==============================================================================

#' Constructor for Poisson (count) models
#'
#' @description
#' Sets up the initial state and class structure for count emission models, in
#' which each item is Poisson-distributed within class with its own rate.
#'
#' @param n_components Integer. The number of latent classes/components to estimate.
#' @param type Character. The specific distribution type, "poisson" or the
#'   missing-data variant "poisson_nan".
#' @param ... Additional arguments passed to the method.
#'
#' @return A list object of class \code{c(type, "emission")} containing the model state.
#' @export
poisson_model <- function(n_components, type = "poisson", ...) {
  state <- list(
    n_components = n_components,
    parameters = list()
  )
  class(state) <- c(type, "emission")
  return(state)
}

# ------------------------------------------------------------------------------
# 1. Poisson (Count) S3 Methods
# ------------------------------------------------------------------------------

# Starting rates are the item marginal means multiplied by a random factor in
# [0.5, 1.5]. Counts have no natural scale — an item may average 0.3 events or
# 300 — so a scale-free draw like the Bernoulli model's runif(0.25, 0.75) would
# start every class in the wrong place for all but one item. Perturbing the
# observed marginal keeps the starting values on the data's own scale while
# still separating the classes enough for EM to move them apart.
#' @exportS3Method
init_params.poisson <- function(model_state, X, resp, random_state = NULL, ...) {
  if (!is.null(random_state)) set.seed(random_state)
  K <- model_state$n_components
  marginal <- colMeans(X, na.rm = TRUE)
  # An item that is all-missing or entirely zero would otherwise seed a rate of
  # zero, which is an absorbing state: the M-step can never move it away.
  marginal[!is.finite(marginal) | marginal <= 0] <- 1e-2
  model_state$parameters$rates <-
    matrix(runif(K * ncol(X), 0.5, 1.5), nrow = K) *
    matrix(marginal, nrow = K, ncol = ncol(X), byrow = TRUE)
  return(model_state)
}

# The weighted-mean update is the Poisson M-step, with the same conservative
# Gamma prior the Bernoulli model uses in Dirichlet form: prior_obs = alpha / K
# pseudo-observations centred on the item's marginal mean. Without it a class
# that empties out during EM drives a rate to exactly zero, and dpois(x, 0) is
# -Inf for any x > 0 — an absorbing state that no later iteration can leave.
#' @exportS3Method
m_step.poisson <- function(model_state, X, resp, weights = NULL, alpha = 1.0, ...) {
  if (!is.null(weights)) {
    resp <- sweep(resp, 1, weights, "*")
    marginal <- colSums(sweep(X, 1, weights, "*"), na.rm = TRUE) / sum(weights)
  } else {
    marginal <- colMeans(X, na.rm = TRUE)
  }
  marginal[!is.finite(marginal)] <- 0

  K         <- model_state$n_components
  prior_obs <- alpha / K

  rates    <- t(resp) %*% X
  rates    <- sweep(rates, 2, prior_obs * marginal, "+")
  sum_resp <- colSums(resp) + prior_obs
  rates    <- sweep(rates, 1, sum_resp, "/")

  model_state$parameters$rates <- pmax(rates, 1e-10)
  return(model_state)
}

#' @exportS3Method
log_likelihood.poisson <- function(model_state, X, ...) {
  n <- nrow(X)
  log_eps <- matrix(0, nrow = n, ncol = model_state$n_components)
  for (c in seq_len(model_state$n_components)) {
    rate_c <- matrix(model_state$parameters$rates[c, ],
                     nrow = n, ncol = ncol(X), byrow = TRUE)
    log_eps[, c] <- rowSums(dpois(X, lambda = rate_c, log = TRUE))
  }
  return(log_eps)
}

#' @exportS3Method
n_parameters.poisson <- function(model_state, ...) {
  return(length(model_state$parameters$rates))
}

# ------------------------------------------------------------------------------
# 2. Poisson NaN (Count with Missing Data, FIML) S3 Methods
# ------------------------------------------------------------------------------

#' @exportS3Method init_params poisson_nan
init_params.poisson_nan <- init_params.poisson
#' @exportS3Method n_parameters poisson_nan
n_parameters.poisson_nan <- n_parameters.poisson

#' @exportS3Method
m_step.poisson_nan <- function(model_state, X, resp, weights = NULL, alpha = 1.0, ...) {
  if (!is.null(weights)) resp <- sweep(resp, 1, weights, "*")

  K         <- model_state$n_components
  prior_obs <- alpha / K
  rates     <- matrix(0, nrow = K, ncol = ncol(X))

  # Each item is updated over the cases that observed it, so a missing cell
  # informs neither the numerator nor the effective count for that item.
  for (j in seq_len(ncol(X))) {
    valid <- !is.na(X[, j])
    if (any(valid)) {
      resp_valid <- resp[valid, , drop = FALSE]
      x_j        <- X[valid, j]

      marg_j <- if (!is.null(weights))
        sum(x_j * weights[valid]) / sum(weights[valid])
      else mean(x_j)

      rates[, j] <- (as.vector(t(resp_valid) %*% x_j) + prior_obs * marg_j) /
        (colSums(resp_valid) + prior_obs)
    } else {
      rates[, j] <- 1e-2   # item never observed: hold at the init_params floor
    }
  }
  model_state$parameters$rates <- pmax(rates, 1e-10)
  return(model_state)
}

#' @exportS3Method
log_likelihood.poisson_nan <- function(model_state, X, ...) {
  n <- nrow(X)
  log_eps <- matrix(0, nrow = n, ncol = model_state$n_components)
  for (c in seq_len(model_state$n_components)) {
    rate_c <- matrix(model_state$parameters$rates[c, ],
                     nrow = n, ncol = ncol(X), byrow = TRUE)
    ll_matrix <- dpois(X, lambda = rate_c, log = TRUE)
    ll_matrix[is.na(ll_matrix)] <- 0 # FIML: missing cells drop out of the sum
    log_eps[, c] <- rowSums(ll_matrix)
  }
  return(log_eps)
}
