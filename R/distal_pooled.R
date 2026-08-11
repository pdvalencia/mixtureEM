# ==============================================================================
# S3 Distal Pooled Regression Model (Pooled Slopes, Class-Varying Intercepts)
# ==============================================================================
# distal_one_hot() and distal_forward() are defined in utils.R

distal_pooled_model <- function(n_components, tol = 1e-4, max_iter = 500,
                                method = "newton-raphson") {
  state <- list(n_components = n_components, tol = tol, max_iter = max_iter,
                method = method, parameters = list())
  class(state) <- c("distal_pooled", "emission")
  return(state)
}

# ------------------------------------------------------------------------------
# Internal helper: validate and normalise the categorical outcome column.
# Shared logic with distal_regression, defined here locally so distal_pooled.R
# is self-contained without a cross-file dependency on a private function.
# ------------------------------------------------------------------------------
.validate_pooled_Y <- function(Y_col, context = "distal_pooled") {
  Y <- as.numeric(Y_col)

  # Shift 0-indexed to 1-indexed
  if (!all(is.na(Y)) && min(Y, na.rm = TRUE) == 0) Y <- Y + 1

  # Guard against degenerate constant outcome
  M <- max(Y, na.rm = TRUE)
  if (M <= 1) {
    warning(sprintf(
      paste0("%s: outcome Y has only one unique category (M = %d). ",
             "A regression model cannot be estimated. ",
             "Returning the model state unchanged."),
      context, M
    ))
    return(NULL)
  }

  return(Y)
}

#' @exportS3Method
init_params.distal_pooled <- function(model_state, X, resp,
                                      random_state = NULL, ...) {
  if (!is.null(random_state)) set.seed(random_state)

  Y <- .validate_pooled_Y(X[, 1], "distal_pooled init_params")
  if (is.null(Y)) {
    model_state$max_val <- 1L
    model_state$parameters$beta_pooled <- matrix(0, nrow = 0, ncol = 0)
    return(model_state)
  }

  Z     <- X[, -1, drop = FALSE]
  K     <- model_state$n_components
  D_cov <- ncol(Z)
  M     <- max(Y, na.rm = TRUE)
  L     <- K + D_cov

  model_state$max_val <- M
  model_state$parameters$beta_pooled <-
    matrix(rnorm((M - 1) * L, 0, 1.0), nrow = M - 1, ncol = L)
  return(model_state)
}

#' @exportS3Method
m_step.distal_pooled <- function(model_state, X, resp, weights = NULL,
                                 alpha = NULL, ...) {
  Y <- .validate_pooled_Y(X[, 1], "distal_pooled m_step")
  if (is.null(Y)) return(model_state)  # constant outcome - skip

  Z <- complete_covariates(X[, -1, drop = FALSE])

  M     <- model_state$max_val
  K     <- model_state$n_components
  D_cov <- ncol(Z)
  L     <- K + D_cov
  N     <- nrow(X)

  w_case <- if (is.null(weights)) rep(1, N) else weights
  if (!is.null(weights)) resp <- sweep(resp, 1, weights, "*")

  U <- matrix(0, nrow = N * K, ncol = L)
  for (k in 1:K) {
    idx <- ((k - 1) * N + 1):(k * N)
    U[idx, k] <- 1
    if (D_cov > 0) U[idx, (K + 1):L] <- Z
  }

  W_flat       <- as.vector(resp)
  Y_oh_flat    <- do.call(rbind, replicate(K, distal_one_hot(Y, M),
                                           simplify = FALSE))
  valid_Y_flat <- rep(!is.na(Y), K)

  U_valid    <- U[valid_Y_flat, , drop = FALSE]
  W_valid    <- W_flat[valid_Y_flat]
  Y_oh_valid <- Y_oh_flat[valid_Y_flat, , drop = FALSE]
  n_data     <- nrow(U_valid)

  # Dirichlet prior on the class-specific response probabilities, in the same
  # alpha/K pseudo-observations form the categorical measurement M-step uses, so
  # that bayes_constants$categorical means the same thing on both engines.
  # Without it a class in which no member gives a particular response drives its
  # intercept to minus infinity instead of to a large finite logit.
  #
  # The prior enters as K extra observations, one per class, each carrying
  # alpha/K weight and the outcome's observed marginal distribution in place of
  # a one-hot response: the existing gradient and Hessian sums then pick it up
  # unchanged. That equivalence holds because the design row of a pseudo-
  # observation is the class indicator itself, which is only true when there are
  # no covariates -- with covariates there is no non-arbitrary place in
  # covariate space to put them, so the prior is confined to the no-covariate
  # (distal_categorical) engine.
  if (inherits(model_state, "distal_categorical") && D_cov == 0L) {
    alpha     <- alpha %||% .bayes_alpha(model_state, "categorical")
    prior_obs <- alpha / K
    obs_Y     <- !is.na(Y)
    if (prior_obs > 0 && any(obs_Y)) {
      q <- colSums(distal_one_hot(Y[obs_Y], M) * w_case[obs_Y]) / sum(w_case[obs_Y])
      U_prior <- matrix(0, nrow = K, ncol = L)
      U_prior[cbind(seq_len(K), seq_len(K))] <- 1
      U_valid    <- rbind(U_valid, U_prior)
      W_valid    <- c(W_valid, rep(prior_obs, K))
      Y_oh_valid <- rbind(Y_oh_valid, matrix(q, nrow = K, ncol = M, byrow = TRUE))
    }
  }

  beta_mat <- model_state$parameters$beta_pooled

  for (iter in seq_len(model_state$max_iter)) {
    P <- distal_forward(U_valid, beta_mat)

    G <- matrix(0, nrow = M - 1, ncol = L)
    for (m in 2:M) {
      diff       <- W_valid * (Y_oh_valid[, m] - P[, m])
      G[m - 1, ] <- t(U_valid) %*% diff
    }
    G_flat <- as.vector(t(G))
    if (max(abs(G_flat)) < model_state$tol) break

    H_nr <- matrix(0, nrow = (M - 1) * L, ncol = (M - 1) * L)
    for (m in 2:M) {
      for (m_prime in 2:M) {
        kronecker <- if (m == m_prime) 1 else 0
        r     <- W_valid * P[, m] * (kronecker - P[, m_prime])
        sub_H <- -t(U_valid) %*% sweep(U_valid, 1, r, "*")
        idx_m       <- ((m - 2) * L + 1):((m - 1) * L)
        idx_m_prime <- ((m_prime - 2) * L + 1):((m_prime - 1) * L)
        H_nr[idx_m, idx_m_prime] <- sub_H
      }
    }
    diag(H_nr) <- diag(H_nr) - 1e-4

    update_step <- pinv(-H_nr) %*% G_flat
    max_step    <- 1.0
    if (max(abs(update_step)) > max_step)
      update_step <- update_step * (max_step / max(abs(update_step)))

    beta_flat <- as.vector(t(beta_mat)) + update_step
    beta_mat  <- t(matrix(beta_flat, nrow = L, ncol = M - 1))
  }

  # ----------------------------------------------------------------
  # Always recompute H at the FINAL converged beta_mat (ensures H is
  # computed correctly even if the loop breaks early): if the gradient
  # was already small at iteration 1, the loop breaks before H_nr
  # is ever assigned, and the previous call's H would be silently returned.
  # The ridge penalty is omitted here — we want the true Fisher
  # information matrix at the MLE, not the NR stabilisation term.
  # The prior's pseudo-observations are dropped for the same reason: they are
  # part of the objective that was maximised, not part of the data whose
  # curvature the standard errors describe.
  # ----------------------------------------------------------------
  U_final <- U_valid[seq_len(n_data), , drop = FALSE]
  W_final <- W_valid[seq_len(n_data)]
  P_final <- distal_forward(U_final, beta_mat)
  H_final <- matrix(0, nrow = (M - 1) * L, ncol = (M - 1) * L)
  for (m in 2:M) {
    for (m_prime in 2:M) {
      kronecker <- if (m == m_prime) 1 else 0
      r     <- W_final * P_final[, m] * (kronecker - P_final[, m_prime])
      sub_H <- -t(U_final) %*% sweep(U_final, 1, r, "*")
      idx_m       <- ((m - 2) * L + 1):((m - 1) * L)
      idx_m_prime <- ((m_prime - 2) * L + 1):((m_prime - 1) * L)
      H_final[idx_m, idx_m_prime] <- sub_H
    }
  }

  model_state$parameters$beta_pooled <- beta_mat
  model_state$parameters$hessian     <- H_final
  return(model_state)
}

#' @exportS3Method
log_likelihood.distal_pooled <- function(model_state, X, ...) {
  Y <- .validate_pooled_Y(X[, 1], "distal_pooled log_likelihood")
  if (is.null(Y)) {
    # Constant outcome: contribute zero log-likelihood (uninformative)
    return(matrix(0, nrow = nrow(X), ncol = model_state$n_components))
  }

  Z <- complete_covariates(X[, -1, drop = FALSE])

  M     <- model_state$max_val
  K     <- model_state$n_components
  D_cov <- ncol(Z)
  L     <- K + D_cov
  N     <- nrow(X)

  U <- matrix(0, nrow = N * K, ncol = L)
  for (k in 1:K) {
    idx <- ((k - 1) * N + 1):(k * N)
    U[idx, k] <- 1
    if (D_cov > 0) U[idx, (K + 1):L] <- Z
  }

  P_stacked <- pmax(pmin(
    distal_forward(U, model_state$parameters$beta_pooled), 1 - 1e-15), 1e-15)

  prob_obs_stacked                   <- numeric(N * K)
  valid_Y_stacked                    <- rep(!is.na(Y), K)
  Y_stacked                          <- rep(Y, K)
  prob_obs_stacked[valid_Y_stacked]  <-
    P_stacked[cbind(which(valid_Y_stacked), Y_stacked[valid_Y_stacked])]
  prob_obs_stacked[!valid_Y_stacked] <- 1

  return(log(matrix(prob_obs_stacked, nrow = N, ncol = K)))
}

#' @exportS3Method
n_parameters.distal_pooled <- function(model_state, ...) {
  return(length(model_state$parameters$beta_pooled))
}
