# ==============================================================================
# S3 Distal Continuous Outcome (Weighted ANOVA) - ROBUST SEs
# ==============================================================================

# `variances`: "class_specific" (one residual variance per class) or "equal"
# (one shared by every class). The means are weighted averages either way; only
# the variance step differs.
distal_continuous_model <- function(n_components, variances = "class_specific",
                                    outcome_variances = NULL, ...) {
  # fit_mixture() passes its choice under this name, because the same
  # arguments also reach the measurement model.
  if (!is.null(outcome_variances)) variances <- outcome_variances
  state <- list(n_components = n_components, variances = variances,
                parameters = list())
  class(state) <- c("distal_continuous", "emission")
  return(state)
}

#' @exportS3Method
init_params.distal_continuous <- function(model_state, X, resp, ...) {
  Y <- as.numeric(X[, 1])
  K <- model_state$n_components
  model_state$parameters$means <- matrix(mean(Y, na.rm=TRUE), nrow = K, ncol = 1)
  model_state$parameters$covariances <- matrix(var(Y, na.rm=TRUE), nrow = K, ncol = 1)
  return(model_state)
}

#' @exportS3Method
m_step.distal_continuous <- function(model_state, X, resp, weights = NULL, ...) {
  Y <- as.numeric(X[, 1])
  valid <- !is.na(Y)
  Y_v <- Y[valid]
  resp_v <- resp[valid, , drop=FALSE]

  if (!is.null(weights)) resp_v <- sweep(resp_v, 1, weights[valid], "*")

  K <- model_state$n_components
  means <- numeric(K)
  vars <- numeric(K)
  ses <- numeric(K)

  for (k in 1:K) {
    W_k <- resp_v[, k]
    Nk <- sum(W_k)
    Nk_abs <- sum(abs(W_k))

    if (Nk > 1e-5) {
      # Point Estimate
      mu_k <- sum(W_k * Y_v) / Nk
      resids <- Y_v - mu_k

      # Log-Likelihood Variance (Absolute weights)
      var_k <- sum(abs(W_k) * resids^2) / Nk_abs

      # Robust Sandwich SE (Meat / Bread^2)
      meat <- sum((W_k * resids)^2)
      se_k <- sqrt(meat) / abs(Nk)

    } else {
      mu_k <- 0
      var_k <- 1e-5
      se_k <- 1e-5
    }

    means[k] <- mu_k
    vars[k] <- max(var_k, 1e-5)
    ses[k] <- se_k
  }

  # One shared variance: the class sums of squares pooled over the classes,
  #   sigma^2 = sum_k sum_i w_ik (y_i - mu_k)^2 / sum_k sum_i w_ik,
  # with the same absolute weights as the class-specific form.
  if (identical(model_state$variances, "equal")) {
    A <- abs(resp_v)
    den <- sum(A)
    pooled <- if (den > 1e-5) sum(A * outer(Y_v, means, "-")^2) / den else 1e-5
    vars[] <- max(pooled, 1e-5)
  }

  model_state$parameters$means <- matrix(means, nrow = K, ncol = 1)
  model_state$parameters$covariances <- matrix(vars, nrow = K, ncol = 1)
  model_state$parameters$ses <- matrix(ses, nrow = K, ncol = 1)

  return(model_state)
}

#' @exportS3Method
log_likelihood.distal_continuous <- function(model_state, X, ...) {
  Y <- as.numeric(X[, 1])
  K <- model_state$n_components
  N <- length(Y)
  ll <- matrix(0, nrow = N, ncol = K)
  valid <- !is.na(Y)

  for (k in 1:K) {
    if (any(valid)) {
      ll[valid, k] <- dnorm(Y[valid],
                            mean = model_state$parameters$means[k, 1],
                            sd = sqrt(model_state$parameters$covariances[k, 1]),
                            log = TRUE)
    }
  }
  return(ll)
}

#' @exportS3Method
n_parameters.distal_continuous <- function(model_state, ...) {
  K <- model_state$n_components
  K + if (identical(model_state$variances, "equal")) 1L else K
}
