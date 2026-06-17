# ==============================================================================
# S3 Covariate Model - FINAL (Firth Estimation + Unpenalized Inference)
# ==============================================================================

covariate_model <- function(n_components, tol = 1e-6, max_iter = 500, intercept = TRUE) {
  state <- list(n_components = n_components, tol = tol, max_iter = max_iter,
                intercept = intercept, parameters = list())
  class(state) <- c("covariate", "emission")
  return(state)
}

#' @exportS3Method
m_step.covariate <- function(model_state, X, resp, weights = NULL, ...) {
  # Missing covariates are completed under the class-invariant Gaussian marginal
  # (endogenous-constrained-x; Sterba, 2014) so every case is retained rather
  # than dropped or filled with an unconditional column mean. With complete
  # covariates this is a no-op.
  X <- complete_covariates(as.matrix(X))
  X_mat <- if(model_state$intercept) cbind(1, X) else X
  K <- model_state$n_components
  D <- ncol(X_mat)
  w_vec <- if(!is.null(weights)) weights else rep(1, nrow(X_mat))

  # ============================================================================
  # 1. ESTIMATION PASS (Data Augmentation)
  # ============================================================================
  # We add a tiny "ghost" observation to every class at the mean of X.
  # This makes Complete Separation mathematically impossible during optimization.
  pseudo_X <- matrix(colMeans(X_mat), nrow = K, ncol = D, byrow = TRUE)
  pseudo_resp <- diag(K)
  pseudo_w <- rep(0.01, K)

  X_aug <- rbind(X_mat, pseudo_X)
  resp_aug <- rbind(resp, pseudo_resp)
  w_aug <- c(w_vec, pseudo_w)

  nll_func <- function(pars) {
    B <- rbind(matrix(pars, K-1, D, byrow=TRUE), 0)
    logits <- X_aug %*% t(B)
    logits <- pmax(pmin(logits, 50), -50)
    max_l <- apply(logits, 1, max)
    prob <- exp(logits - max_l) / rowSums(exp(logits - max_l))
    -sum(w_aug * rowSums(resp_aug * log(prob + 1e-15)))
  }

  # Analytical gradient of nll_func wrt pars.
  # For multinomial logistic regression:
  #   ∂nll/∂B[k,d] = Σ_i w_i (prob[i,k] - resp[i,k]) * X[i,d]  for k < K
  # Vectorised: grad[k, :] = t(X_aug) %*% (w_aug * (prob[,k] - resp_aug[,k]))
  nll_grad <- function(pars) {
    B <- rbind(matrix(pars, K-1, D, byrow=TRUE), 0)
    logits <- X_aug %*% t(B)
    logits <- pmax(pmin(logits, 50), -50)
    max_l <- apply(logits, 1, max)
    prob <- exp(logits - max_l) / rowSums(exp(logits - max_l))
    # residual: (prob - resp_aug) weighted by w_aug, for free classes k = 1..K-1
    resid <- sweep(prob[, seq_len(K-1), drop=FALSE] -
                     resp_aug[, seq_len(K-1), drop=FALSE], 1, w_aug, "*")
    as.vector(t(X_aug) %*% resid)   # D × (K-1), flattened row-major
  }

  fit <- optim(par = rep(0, (K-1)*D), fn = nll_func, gr = nll_grad,
               method = "BFGS")
  beta_final <- rbind(matrix(fit$par, K-1, D, byrow=TRUE), 0)

  # ============================================================================
  # 2. INFERENCE PASS (Original Data Only)
  # ============================================================================
  # We throw away the pseudo-data and calculate the Hessian strictly on the
  # original X_mat and w_vec to ensure Standard Errors are not artificially shrunk.
  logits <- X_mat %*% t(beta_final)
  logits <- pmax(pmin(logits, 50), -50)
  max_l <- apply(logits, 1, max)
  prob <- exp(logits - max_l) / rowSums(exp(logits - max_l))

  H <- matrix(0, (K-1)*D, (K-1)*D)

  for(k in seq_len(K-1)) {
    for(j in seq_len(K-1)) {
      W <- prob[,k] * ((if(k==j) 1 else 0) - prob[,j]) * w_vec
      H_kj <- -t(X_mat) %*% sweep(X_mat, 1, W, "*")

      H[((k-1)*D+1):(k*D), ((j-1)*D+1):(j*D)] <- H_kj
    }
  }

  # Pad with zeros for the anchor class (identifiability constraint)
  H_full <- matrix(0, K*D, K*D)
  if ((K-1)*D > 0)
    H_full[1:((K-1)*D), 1:((K-1)*D)] <- H
  diag(H_full)[((K-1)*D+1):(K*D)] <- -1e8

  model_state$parameters$beta <- beta_final
  model_state$parameters$hessian <- H_full

  # ============================================================================
  # 3. SURVEY-ROBUST COVARIANCE (optional)
  # ============================================================================
  # When a complex survey design is attached to the sub-model, the naive
  # information-based variance is replaced by the linearization sandwich
  #   V = (-H)^{-1} B (-H)^{-1}
  # where B aggregates the multinomial-logistic scores to the PSU level within
  # strata. Computed only when the design vectors are present and row-aligned
  # with the data used in this M-step; otherwise downstream code falls back to
  # the Hessian-based variance.
  has_design   <- isTRUE(model_state$has_survey_design) &&
    !is.null(model_state$strata) && !is.null(model_state$cluster) &&
    length(model_state$strata)  == nrow(X_mat) &&
    length(model_state$cluster) == nrow(X_mat)

  if (has_design && (K-1)*D > 0) {
    # Individual score vectors for the free classes k = 1..K-1. The score for
    # the multinomial logit is X_i * (resp_ik - prob_ik), weighted by w_i, the
    # same weighting used to form the Hessian above.
    score_list <- vector("list", K - 1)
    for (k in seq_len(K - 1)) {
      resid_k         <- (resp[, k] - prob[, k]) * w_vec
      score_list[[k]] <- sweep(X_mat, 1, resid_k, "*")
    }
    score_mat <- do.call(cbind, score_list)   # N x ((K-1)*D)

    meat     <- compute_survey_B(score_mat, model_state$strata, model_state$cluster)
    B_inv    <- pinv(-H)
    V_robust <- B_inv %*% meat %*% B_inv

    # Pad with zeros for the anchor class to match the K*D layout of H_full.
    V_full <- matrix(0, K*D, K*D)
    V_full[1:((K-1)*D), 1:((K-1)*D)] <- V_robust
    model_state$parameters$V_robust <- V_full
  }

  return(model_state)
}

#' @exportS3Method
init_params.covariate <- function(model_state, X, resp, ...) {
  D <- ncol(X) + as.integer(model_state$intercept)
  model_state$parameters$beta <- matrix(0, model_state$n_components, D)
  return(model_state)
}

#' @exportS3Method
log_likelihood.covariate <- function(model_state, X, ...) {
  # 1. Prepare Matrix (completing any missing covariates under the shared,
  #    class-invariant Gaussian marginal; a no-op when covariates are complete).
  X <- complete_covariates(as.matrix(X))
  X_mat <- if(model_state$intercept) cbind(1, X) else X

  # 2. Compute Raw Logits
  logits <- X_mat %*% t(model_state$parameters$beta)

  # 3. Use the Log-Sum-Exp trick for the denominator
  # This uses your existing helper in utils.R
  log_denominator <- logsumexp(logits, MARGIN = 1)

  # 4. Compute Log-Probabilities: Log(Prob) = Logits - Log(Sum_Exp_Logits)
  # sweep subtracts the log_denominator from each row of logits
  # sweep() can silently drop the dim attribute in ALTREP system,
  # causing non-conformable-arrays when log_prob is added to the mm log-likelihood
  # in e_step. matrix() re-attaches explicit dimensions unconditionally.
  log_prob <- matrix(
    sweep(logits, 1, log_denominator, "-"),
    nrow = nrow(logits), ncol = ncol(logits)
  )

  return(log_prob)
}

#' @exportS3Method
n_parameters.covariate <- function(model_state, ...) {
  return((nrow(model_state$parameters$beta) - 1) * ncol(model_state$parameters$beta))
}
