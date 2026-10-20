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
m_step.poisson <- function(model_state, X, resp, weights = NULL, alpha = NULL,
                           prior_scale = 1, ...) {
  alpha <- alpha %||% .bayes_alpha(model_state, "poisson")
  if (!is.null(weights)) {
    resp <- sweep(resp, 1, weights, "*")
    marginal <- colSums(sweep(X, 1, weights, "*"), na.rm = TRUE) / sum(weights)
  } else {
    marginal <- colMeans(X, na.rm = TRUE)
  }
  marginal[!is.finite(marginal)] <- 0

  K         <- model_state$n_components
  prior_obs <- prior_scale * alpha / K

  rates    <- t(resp) %*% X
  rates    <- sweep(rates, 2, prior_obs * marginal, "+")
  sum_resp <- colSums(resp) + prior_obs
  rates    <- sweep(rates, 1, sum_resp, "/")

  model_state$parameters$rates <- pmax(rates, 1e-10)
  return(model_state)
}

#' @exportS3Method
log_likelihood.poisson <- function(model_state, X, ...) {
  lam <- model_state$parameters$rates
  out <- X %*% t(log(lam)) - rep(rowSums(lam), each = nrow(X)) -
         rowSums(lgamma(X + 1))
  dimnames(out) <- NULL
  return(out)
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
m_step.poisson_nan <- function(model_state, X, resp, weights = NULL, alpha = NULL,
                               prior_scale = 1, ...) {
  alpha <- alpha %||% .bayes_alpha(model_state, "poisson")
  if (!is.null(weights)) resp <- sweep(resp, 1, weights, "*")

  K         <- model_state$n_components
  prior_obs <- prior_scale * alpha / K
  rates     <- matrix(0, nrow = K, ncol = ncol(X),
                      dimnames = list(NULL, colnames(X)))

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
  lam <- model_state$parameters$rates
  X0  <- X
  X0[is.na(X0)] <- 0
  M   <- (!is.na(X)) * 1
  out <- X0 %*% t(log(lam)) - M %*% t(lam) - rowSums(lgamma(X0 + 1) * M)
  dimnames(out) <- NULL
  return(out)
}

# ------------------------------------------------------------------------------
# 3. Zero-inflated Poisson (Count with Excess Zeros) S3 Methods
# ------------------------------------------------------------------------------
#
# Each item, within each class, has a rate `rates` for the cases at risk and a
# probability `inflation` that a case is not at risk at all and records a
# structural zero. One implementation serves both `zip` and `zip_nan`: an
# unobserved cell drops out of every sum, which is a no-op on complete data.

#' @exportS3Method
init_params.zip <- function(model_state, X, resp, random_state = NULL, ...) {
  model_state <- init_params.poisson(model_state, X, resp, random_state)
  # Seed every class at the item's excess zeros -- the observed share of
  # zeros beyond what a Poisson at the item mean would produce -- kept off
  # 0 and 1 and spread by the same random factor as the rates.
  K  <- model_state$n_components
  m  <- colMeans(X, na.rm = TRUE)
  m[!is.finite(m)] <- 0
  p0  <- exp(-m)
  exc <- (colMeans(X == 0, na.rm = TRUE) - p0) / pmax(1 - p0, 1e-10)
  exc[!is.finite(exc)] <- 0
  exc <- pmin(pmax(exc, 0.05), 0.6)
  model_state$parameters$inflation <-
    matrix(runif(K * ncol(X), 0.5, 1.5), nrow = K) *
    matrix(exc, nrow = K, ncol = ncol(X), byrow = TRUE)
  model_state
}

# Inner E-step and closed-form M-step. z is the chance that an observed zero is
# structural, given the class and the current parameters; the rate is a
# weighted mean over the share 1 - z of each record that is at risk, under the
# Poisson model's gamma prior centred on the item mean, and the inflation is
# the posterior-weighted mean of z.
#' @exportS3Method
m_step.zip <- function(model_state, X, resp, weights = NULL, alpha = NULL,
                       prior_scale = 1, ...) {
  alpha <- alpha %||% .bayes_alpha(model_state, "poisson")
  if (!is.null(weights)) resp <- sweep(resp, 1, weights, "*")
  K   <- model_state$n_components
  obs <- !is.na(X)
  X0  <- replace(X, !obs, 0)
  zer <- obs & X0 == 0
  marg <- colSums(X0 * (if (is.null(weights)) 1 else weights)) /
    pmax(colSums(obs * (if (is.null(weights)) 1 else weights)), 1e-12)
  prior_obs <- prior_scale * alpha / K

  lam <- model_state$parameters$rates
  pz  <- model_state$parameters$inflation
  for (k in seq_len(K)) {
    L  <- matrix(lam[k, ], nrow(X), ncol(X), byrow = TRUE)
    P  <- matrix(pz[k, ],  nrow(X), ncol(X), byrow = TRUE)
    Z  <- ifelse(zer, P / (P + (1 - P) * exp(-L)), 0)
    rk <- resp[, k]
    at_risk  <- colSums(rk * obs * (1 - Z))
    lam[k, ] <- (colSums(rk * X0 * (1 - Z)) + prior_obs * marg) /
      (at_risk + prior_obs)
    pz[k, ]  <- colSums(rk * Z) / pmax(colSums(rk * obs), 1e-300)
  }
  model_state$parameters$rates     <- pmax(lam, 1e-10)
  model_state$parameters$inflation <- pmin(pmax(pz, 0), 1 - 1e-10)
  model_state
}

#' @exportS3Method
log_likelihood.zip <- function(model_state, X, ...) {
  lam <- model_state$parameters$rates
  pz  <- model_state$parameters$inflation
  obs <- !is.na(X)
  X0  <- replace(X, !obs, 0)
  zer <- X0 == 0
  out <- vapply(seq_len(model_state$n_components), function(k) {
    L  <- matrix(lam[k, ], nrow(X), ncol(X), byrow = TRUE)
    P  <- matrix(pz[k, ],  nrow(X), ncol(X), byrow = TRUE)
    ll <- ifelse(zer, log(P + (1 - P) * exp(-L)),
                 log1p(-P) + X0 * log(L) - L - lgamma(X0 + 1))
    rowSums(ll * obs)
  }, numeric(nrow(X)))
  matrix(out, nrow = nrow(X))
}

#' @exportS3Method
n_parameters.zip <- function(model_state, ...)
  length(model_state$parameters$rates) + length(model_state$parameters$inflation)

#' @exportS3Method init_params zip_nan
init_params.zip_nan <- init_params.zip
#' @exportS3Method m_step zip_nan
m_step.zip_nan <- m_step.zip
#' @exportS3Method log_likelihood zip_nan
log_likelihood.zip_nan <- log_likelihood.zip
#' @exportS3Method n_parameters zip_nan
n_parameters.zip_nan <- n_parameters.zip

# ------------------------------------------------------------------------------
# 4. Negative binomial (Overdispersed Count) S3 Methods
# ------------------------------------------------------------------------------
#
# NB2: each item, within each class, has a mean `rates`; the item's
# `overdispersion` alpha, variance = mean + alpha mean^2, is shared by the
# classes (stored K x items with equal rows, so class sorting needs no special
# case). One implementation serves `negbin` and `negbin_nan`.

#' @exportS3Method
init_params.negbin <- function(model_state, X, resp, random_state = NULL, ...) {
  model_state <- init_params.poisson(model_state, X, resp, random_state)
  # The item's moment estimate (var - mean) / mean^2, kept off the Poisson
  # boundary, from which the one-dimensional search could still leave.
  m  <- pmax(colMeans(X, na.rm = TRUE), 1e-2)
  v  <- apply(X, 2, stats::var, na.rm = TRUE)
  a0 <- (v - m) / m^2
  a0[!is.finite(a0)] <- 0.05
  a0 <- pmin(pmax(a0, 0.05), 10)
  model_state$parameters$overdispersion <-
    matrix(a0, nrow = model_state$n_components, ncol = ncol(X), byrow = TRUE)
  model_state
}

# ECM. The mean of an item in a class, at fixed alpha and under the Poisson
# model's gamma prior (prior_obs pseudo-cases at the item mean m), solves
#   (Sx - R mu) / (mu (1 + alpha mu)) + prior_obs (m - mu) / mu = 0,
# a quadratic with one positive root; with no prior it is the weighted mean,
# whatever alpha. Then each item's alpha, pooled over the classes.
#' @exportS3Method
m_step.negbin <- function(model_state, X, resp, weights = NULL, alpha = NULL,
                          prior_scale = 1, ...) {
  alpha <- alpha %||% .bayes_alpha(model_state, "poisson")
  if (!is.null(weights)) resp <- sweep(resp, 1, weights, "*")
  K   <- model_state$n_components
  obs <- !is.na(X)
  X0  <- replace(X, !obs, 0)
  marg <- colSums(X0 * (if (is.null(weights)) 1 else weights)) /
    pmax(colSums(obs * (if (is.null(weights)) 1 else weights)), 1e-12)
  p  <- prior_scale * alpha / K
  od <- model_state$parameters$overdispersion

  Sx <- t(resp) %*% X0                       # K x J
  R  <- t(resp) %*% obs
  M  <- matrix(marg, K, ncol(X), byrow = TRUE)
  A  <- p * od
  B  <- R + p - p * od * M
  C  <- Sx + p * M
  lam <- 2 * C / (B + sqrt(B^2 + 4 * A * C))
  lam[!is.finite(lam)] <- 1e-10
  lam <- pmax(lam, 1e-10)

  for (j in seq_len(ncol(X))) {
    o <- obs[, j]
    od[, j] <- .nb_alpha_fit(rep(X0[o, j], K), rep(lam[, j], each = sum(o)),
                             as.vector(resp[o, , drop = FALSE]))
  }
  model_state$parameters$rates          <- lam
  model_state$parameters$overdispersion <- od
  model_state
}

#' @exportS3Method
log_likelihood.negbin <- function(model_state, X, ...) {
  lam <- model_state$parameters$rates
  od  <- model_state$parameters$overdispersion
  obs <- !is.na(X)
  X0  <- replace(X, !obs, 0)
  out <- vapply(seq_len(model_state$n_components), function(k) {
    ll <- .nb_log_dens(X0, matrix(lam[k, ], nrow(X), ncol(X), byrow = TRUE),
                       matrix(od[k, ], nrow(X), ncol(X), byrow = TRUE))
    rowSums(ll * obs)
  }, numeric(nrow(X)))
  matrix(out, nrow = nrow(X))
}

#' @exportS3Method
n_parameters.negbin <- function(model_state, ...)
  length(model_state$parameters$rates) + ncol(model_state$parameters$overdispersion)

#' @exportS3Method init_params negbin_nan
init_params.negbin_nan <- init_params.negbin
#' @exportS3Method m_step negbin_nan
m_step.negbin_nan <- m_step.negbin
#' @exportS3Method log_likelihood negbin_nan
log_likelihood.negbin_nan <- log_likelihood.negbin
#' @exportS3Method n_parameters negbin_nan
n_parameters.negbin_nan <- n_parameters.negbin
