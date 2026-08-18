# ==============================================================================
# S3 Covariate Model - FINAL (Firth Estimation + Unpenalized Inference)
# ==============================================================================

covariate_model <- function(n_components, tol = 1e-6, max_iter = 500, intercept = TRUE) {
  state <- list(n_components = n_components, tol = tol, max_iter = max_iter,
                intercept = intercept, parameters = list())
  class(state) <- c("covariate", "emission")
  return(state)
}

# ==============================================================================
# Weighted multinomial logistic regression (shared fitter)
# ==============================================================================
#
# Fits P(class = k | z) by maximising  Σ_i w_i Σ_k resp_ik log p_k(z_i)  with the
# last class anchored at zero. Used for the covariate structural model and, in
# latent transition models, for the regressions on the initial latent status and
# on each row of the transition matrix.
#
# `Z` must already contain any intercept column. `resp` is n x K and is treated
# as rows summing to one; callers whose responsibilities carry their own mass
# (e.g. the pairwise posteriors ξ of a latent Markov model, whose rows sum to
# γ rather than to 1) should normalise the rows and fold that mass into
# `weights`. The two are algebraically equivalent, and folding the mass into the
# weights keeps this fitter free of special cases.
#
# Estimation uses data augmentation, in one of two forms.
#
# With `alpha = 0` the augmentation is one "ghost" observation per class placed
# at the column means with weight 0.01, which makes complete separation
# impossible during optimisation. The Hessian is then recomputed on the original
# data only so that standard errors are not artificially shrunk by the
# pseudo-data.
#
# With `alpha > 0` the ghost is replaced by the Dirichlet prior on the class
# probabilities, written as fractional pseudo-data so this fitter needs no
# special case: one row per class per unique covariate pattern, each with weight
# alpha / (K * U0). That adds alpha / K cases to each class, spread evenly over
# the covariate patterns, which makes the class sizes slightly more equal and the
# covariate effects slightly smaller. Unlike the ghost, these rows *do* enter the
# Hessian: they are part of the objective being maximised, so they belong in its
# curvature, and the standard errors shrink accordingly.
#
# Returns the K x D coefficient matrix (last row zero), the ((K-1)D)^2 Hessian on
# the free parameters, and the fitted n x K probabilities.
#
# `start` warm-starts the optimiser from a previous coefficient matrix. Inside an
# EM loop successive M-steps move the coefficients very little, so starting from
# the last iteration's values cuts the work sharply without changing the optimum.
.fit_mnl <- function(Z, resp, weights = NULL, augment = TRUE, start = NULL,
                     alpha = 0) {
  K <- ncol(resp)
  D <- ncol(Z)
  w_vec <- if (!is.null(weights)) weights else rep(1, nrow(Z))

  probs_for <- function(B, M) softmax_rows(M %*% t(B))

  if (K < 2L || D < 1L) {
    return(list(beta = matrix(0, K, D), hessian = matrix(0, 0, 0),
                prob = matrix(1 / max(K, 1L), nrow(Z), K)))
  }

  # Does the prior enter the curvature? Only when it is the thing doing the
  # augmenting; the ghost never does.
  prior_rows <- isTRUE(augment) && alpha > 0

  if (prior_rows) {
    Zu      <- unique(Z)
    U0      <- nrow(Zu)
    w_prior <- alpha / (K * U0)
    Z_pri    <- Zu[rep(seq_len(U0), times = K), , drop = FALSE]
    resp_pri <- diag(K)[rep(seq_len(K), each = U0), , drop = FALSE]
    Z_aug    <- rbind(Z, Z_pri)
    resp_aug <- rbind(resp, resp_pri)
    w_aug    <- c(w_vec, rep(w_prior, K * U0))
  } else if (isTRUE(augment)) {
    Z_aug    <- rbind(Z, matrix(colMeans(Z), nrow = K, ncol = D, byrow = TRUE))
    resp_aug <- rbind(resp, diag(K))
    w_aug    <- c(w_vec, rep(0.01, K))
  } else {
    Z_aug <- Z; resp_aug <- resp; w_aug <- w_vec
  }

  nll_func <- function(pars) {
    B    <- rbind(matrix(pars, K - 1, D, byrow = TRUE), 0)
    prob <- probs_for(B, Z_aug)
    -sum(w_aug * rowSums(resp_aug * log(prob + 1e-15)))
  }

  # Analytical gradient of nll_func wrt pars.
  # For multinomial logistic regression:
  #   ∂nll/∂B[k,d] = Σ_i w_i (prob[i,k] - resp[i,k]) * Z[i,d]  for k < K
  # Vectorised: grad[k, :] = t(Z_aug) %*% (w_aug * (prob[,k] - resp_aug[,k]))
  nll_grad <- function(pars) {
    B     <- rbind(matrix(pars, K - 1, D, byrow = TRUE), 0)
    prob  <- probs_for(B, Z_aug)
    resid <- sweep(prob[, seq_len(K - 1), drop = FALSE] -
                     resp_aug[, seq_len(K - 1), drop = FALSE], 1, w_aug, "*")
    as.vector(t(Z_aug) %*% resid)   # D × (K-1), flattened row-major
  }

  par0 <- if (!is.null(start) && identical(dim(start), c(K, D)))
    as.vector(t(start[seq_len(K - 1), , drop = FALSE])) else rep(0, (K - 1) * D)

  fit  <- optim(par = par0, fn = nll_func, gr = nll_grad, method = "BFGS")
  beta <- rbind(matrix(fit$par, K - 1, D, byrow = TRUE), 0)

  prob <- probs_for(beta, Z)

  # Rows the Hessian is taken over: data alone under the ghost, data plus the
  # prior rows when the prior is what augmented the fit.
  #
  # Note that this Hessian therefore includes the prior rows while
  # .step3_hessian() (R/step3_variance.R) excludes them, so the same nominal
  # estimator is reachable by two routes whose answers differ by O(1/n). The
  # step-3 bread normally wins, because R/corrections.R overwrites
  # sm$parameters$hessian with the inverse of the step-3 sandwich; it is only on
  # the fallback paths — BCH, n_steps = 1, covariate-plus-distal — that nothing
  # overwrites it and confint()/the analytical Wald read the matrix built here.
  Z_h    <- if (prior_rows) Z_aug else Z
  w_h    <- if (prior_rows) w_aug else w_vec
  prob_h <- if (prior_rows) probs_for(beta, Z_aug) else prob

  H <- matrix(0, (K - 1) * D, (K - 1) * D)
  for (k in seq_len(K - 1)) {
    for (j in seq_len(K - 1)) {
      W    <- prob_h[, k] * ((if (k == j) 1 else 0) - prob_h[, j]) * w_h
      H_kj <- -t(Z_h) %*% sweep(Z_h, 1, W, "*")
      H[((k - 1) * D + 1):(k * D), ((j - 1) * D + 1):(j * D)] <- H_kj
    }
  }

  list(beta = beta, hessian = H, prob = prob)
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

  # `bayes_constants$latent` is the Dirichlet prior on the class probabilities.
  # It applies wherever those probabilities are estimated, and a covariate model
  # is where they are estimated as a regression rather than as K-1 free weights;
  # reading it here is what stops the prior from silently lapsing the moment
  # covariates enter.
  alpha      <- .bayes_alpha(model_state, "latent")
  mnl        <- .fit_mnl(X_mat, resp, weights = w_vec, augment = TRUE,
                         alpha = alpha)
  beta_final <- mnl$beta
  prob       <- mnl$prob
  H          <- mnl$hessian

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
    model_state$parameters$V_method <- "Survey-robust (linearization)"
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
