# ==============================================================================
# S3 Distal Continuous Pooled Model (Pooled + Class-Specific Slopes, Class-Varying
# Intercepts)
# ==============================================================================

# Expand the N x D_cov covariate matrix Z into the N*K x L design used to
# estimate intercepts, pooled slopes and class-specific ("moderated") slopes
# jointly. `mod` is an integer vector of column indices into Z (1-indexed)
# naming the covariates that get their own slope per class; every other
# covariate gets one slope shared across classes.
#
#     cols 1 .. K                    class intercepts
#     cols K+1 .. K+D_pool           pooled slopes, in covariate-column order
#     then, per moderated covariate, a block of K columns (class 1 .. K)
#
# Shared by m_step.distal_continuous_pooled() and the BCH covariance
# computation in R/corrections.R, which must reproduce this exact layout or
# the two would silently drift apart.
.distal_U <- function(Z, K, N, mod = integer(0)) {
  D_cov  <- ncol(Z)
  pooled <- setdiff(seq_len(D_cov), mod)
  D_pool <- length(pooled)
  L      <- K + D_pool + K * length(mod)

  U <- matrix(0, nrow = N * K, ncol = L)
  for (k in seq_len(K)) {
    idx <- ((k - 1L) * N + 1L):(k * N)
    U[idx, k] <- 1
    if (D_pool > 0) U[idx, (K + 1L):(K + D_pool)] <- Z[, pooled, drop = FALSE]
    for (j in seq_along(mod)) {
      col       <- K + D_pool + (j - 1L) * K + k
      U[idx, col] <- Z[, mod[j]]
    }
  }
  U
}

# Full-length permutation of a theta/ses/cov_theta index set under
# `sort_model_classes()`'s class reordering: intercepts and each moderated
# covariate's class block are permuted by `new_order`; pooled-slope columns
# are left in place. `L` is the total column count (ncol(beta_pooled) or
# nrow(cov_theta), which are the same value).
.distal_pooled_reorder_idx <- function(sm, new_order, L) {
  K      <- sm$n_components
  mod    <- sm$moderated %||% integer(0)
  D_mod  <- length(mod)
  D_pool <- L - K - K * D_mod

  idx <- new_order
  if (D_pool > 0) idx <- c(idx, K + seq_len(D_pool))
  for (j in seq_len(D_mod)) {
    block_start <- K + D_pool + (j - 1L) * K
    idx <- c(idx, block_start + new_order)
  }
  idx
}

# `variances`: "equal" (one residual variance shared by every class) or
# "class_specific" (one per class). The class-specific form has no closed-form
# M-step, because the pooled slopes are then a weighted least-squares fit whose
# weights are the class precisions; see m_step below.
distal_continuous_pooled_model <- function(n_components, moderated = integer(0),
                                           variances = "equal", ...) {
  state <- list(n_components = n_components, moderated = as.integer(moderated),
                variances = variances, parameters = list())
  class(state) <- c("distal_continuous_pooled", "emission")
  return(state)
}

#' @exportS3Method
init_params.distal_continuous_pooled <- function(model_state, X, resp, ...) {
  Y <- as.numeric(X[, 1])
  Z <- complete_covariates(X[, -1, drop = FALSE])

  K      <- model_state$n_components
  D_cov  <- ncol(Z)
  mod    <- model_state$moderated
  D_pool <- D_cov - length(mod)
  L      <- K + D_pool + K * length(mod)

  model_state$parameters$beta_pooled <- matrix(0, nrow = 1, ncol = L)
  # Initialize the intercepts to the global mean
  model_state$parameters$beta_pooled[1, 1:K] <- mean(Y, na.rm = TRUE)

  model_state$parameters$covariances <- matrix(var(Y, na.rm = TRUE), nrow = K, ncol = 1)
  model_state$parameters$ses <- matrix(0, nrow = 1, ncol = L)
  return(model_state)
}

#' @exportS3Method
m_step.distal_continuous_pooled <- function(model_state, X, resp, weights = NULL, ...) {
  Y <- as.numeric(X[, 1])
  Z <- complete_covariates(X[, -1, drop = FALSE])
  valid <- !is.na(Y)

  Y_v    <- Y[valid]
  Z_v    <- Z[valid, , drop = FALSE]
  resp_v <- resp[valid, , drop = FALSE]

  if (!is.null(weights)) resp_v <- sweep(resp_v, 1, weights[valid], "*")

  N_v <- nrow(Z_v)
  K   <- model_state$n_components
  mod <- model_state$moderated

  # 1. Expand design matrix for simultaneous intercept/slope estimation
  U <- .distal_U(Z_v, K, N_v, mod)

  W_flat <- as.vector(resp_v)
  Y_flat <- rep(Y_v, K)

  # Class-specific variances make the slopes a weighted least-squares fit with
  # each record also weighted by its class's precision 1/sigma2_k. One
  # conditional-maximisation cycle -- coefficients at the current variances,
  # then the variances at the new coefficients -- raises the expected
  # complete-data log-likelihood, which is all a generalised EM step needs
  # (Meng & Rubin, 1993). Under equal variances the precision is a common
  # factor and cancels, so it is left out.
  cls_spec  <- identical(model_state$variances, "class_specific")
  prec_flat <- if (cls_spec)
    rep(1 / as.vector(model_state$parameters$covariances), each = N_v) else 1
  W_fit <- W_flat * prec_flat

  # 2. Estimate intercepts, pooled slopes and class-specific slopes  (bread of
  #    the sandwich)
  UWU <- t(U) %*% sweep(U, 1, W_fit, "*")
  UWY <- t(U) %*% (W_fit * Y_flat)

  diag(UWU) <- diag(UWU) + 1e-6        # ridge penalty for stability
  B_inv <- pinv(UWU)
  theta <- B_inv %*% UWY

  # 3. Pooled residual variance
  #
  #    Use SIGNED BCH weights, not absolute weights.
  #
  #    sigma^2 = sum_k sum_i [ w_ik * (y_i - yhat_ik)^2 ]
  #              / sum_k sum_i w_ik
  #
  #    This calculates the pooled error variance under the
  #    homoskedastic model, giving a single scalar shared across all classes.
  #    Using abs() inflates sigma^2 and, consequently, all SEs.
  preds  <- U %*% theta
  resids <- as.vector(Y_flat - preds)

  N_total <- sum(W_flat)
  sigma2  <- if (abs(N_total) > 1e-5)
    sum(W_flat * resids^2) / N_total
  else
    1e-5
  sigma2 <- max(sigma2, 1e-5)

  # Store a K x 1 matrix of the pooled variance (same value for every class)
  # to keep the log_likelihood method working unchanged.
  vars <- matrix(sigma2, nrow = K, ncol = 1)
  if (cls_spec) {
    cls  <- rep(seq_len(K), each = N_v)
    Nk   <- as.vector(tapply(W_flat, cls, sum))
    SSk  <- as.vector(tapply(W_flat * resids^2, cls, sum))
    vars <- matrix(pmax(ifelse(abs(Nk) > 1e-5, SSk / Nk, 1e-5), 1e-5),
                   nrow = K, ncol = 1)
  }

  # 4. Standard errors: a sandwich clustered on the case.
  #
  #    The expanded data set carries K weighted records per person, so the
  #    records are not independent and a record-level sandwich would count
  #    each person K times over. The answer is to cluster the meat on the
  #    person -- sum each person's K record scores into one score before
  #    squaring -- not to abandon the sandwich, which is what this code did
  #    until the BCH standard errors were traced through by hand.
  #
  #    Dropping the meat is not conservative. B = U'WU is built from the
  #    *signed* BCH weights, and the whole point of inverting the
  #    classification table is that undoing classification error costs
  #    information; sigma^2 B^{-1} prices that loss as though this were an
  #    ordinary least-squares fit, and so understates every standard error.
  #    The clustered sandwich keeps it:
  #
  #      s_i = sum_k w_ik (y_i - u_ik' theta) u_ik
  #      V   = B^{-1} ( sum_i s_i s_i' ) B^{-1}
  #
  #    Note: We report SEs for the absolute intercepts rather than pairwise
  #    contrasts. Both are correct; they differ only in parameterisation. The
  #    contrast SE is recoverable as sqrt(V[k,k] + V[1,1] - 2*V[k,1]).
  scores  <- U * (W_fit * resids)
  # rowsum() over the case index folds each person's K records into one score.
  s_case  <- rowsum(scores, rep.int(seq_len(N_v), K), reorder = FALSE)
  meat    <- crossprod(s_case)
  cov_theta <- B_inv %*% meat %*% B_inv
  ses <- sqrt(pmax(diag(cov_theta), 1e-8))

  model_state$parameters$beta_pooled <- matrix(as.vector(theta), nrow = 1)
  model_state$parameters$covariances <- vars
  model_state$parameters$ses         <- matrix(ses, nrow = 1)
  # The omnibus Wald tests and the pairwise class contrasts read cov_theta.
  # It is built here, from the same U, W and mask the estimates came from, so
  # that the covariance and the coefficients it belongs to cannot drift apart.
  model_state$parameters$cov_theta   <- cov_theta

  return(model_state)
}

#' @exportS3Method
log_likelihood.distal_continuous_pooled <- function(model_state, X, ...) {
  Y <- as.numeric(X[, 1])
  Z <- complete_covariates(X[, -1, drop = FALSE])

  K      <- model_state$n_components
  mod    <- model_state$moderated
  D_cov  <- ncol(Z)
  pooled <- setdiff(seq_len(D_cov), mod)
  D_pool <- length(pooled)
  N      <- length(Y)
  valid  <- !is.na(Y)

  ll    <- matrix(0, nrow = N, ncol = K)
  theta <- as.vector(model_state$parameters$beta_pooled)

  for (k in 1:K) {
    if (any(valid)) {
      intercept_k <- theta[k]
      preds <- rep(intercept_k, sum(valid))
      if (D_pool > 0)
        preds <- preds + Z[valid, pooled, drop = FALSE] %*%
          theta[(K + 1L):(K + D_pool)]
      # That case's own moderated columns: class k's slope on each
      # class-specific covariate.
      for (j in seq_along(mod)) {
        col   <- K + D_pool + (j - 1L) * K + k
        preds <- preds + Z[valid, mod[j]] * theta[col]
      }

      ll[valid, k] <- dnorm(Y[valid],
                            mean = as.vector(preds),
                            sd   = sqrt(model_state$parameters$covariances[k, 1]),
                            log  = TRUE)
    }
  }
  return(ll)
}

#' @exportS3Method
n_parameters.distal_continuous_pooled <- function(model_state, ...) {
  # L free coefficients (intercepts + pooled slopes + class-specific slopes)
  # + 1 pooled variance, or K of them when they vary by class
  n_var <- if (identical(model_state$variances, "class_specific"))
    model_state$n_components else 1L
  return(ncol(model_state$parameters$beta_pooled) + n_var)
}
