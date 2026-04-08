# ==============================================================================
# S3 3-Step Corrections (BCH and ML)
# ==============================================================================

get_modal_resp <- function(resp) {
  n <- nrow(resp)
  K <- ncol(resp)
  modal <- matrix(0, nrow = n, ncol = K)
  assigned_classes <- max.col(resp, ties.method = "first")
  for (i in 1:n) {
    modal[i, assigned_classes[i]] <- 1
  }
  return(modal)
}

# Apply the BCH Correction
#
# Implementation follows Vermunt (2010) and LatentGOLD's proportional BCH:
#
#   Step 1: Build the classification error matrix C from proportional assignment.
#           C[k,j] = sum_i resp[i,k] * resp[i,j] / sum_i resp[i,j]
#           Columns index the true class and sum to 1.
#
#   Step 2: BCH weight matrix D = t(inv(C))
#
#   Step 3: Apply to proportional (posterior) weights:
#           bch_resp = resp %*% D
#
#   Negative weights are retained; they are a mathematically necessary
#   feature of the BCH correction for poorly-separated classes (Bakk et al.,
#   2013) and must not be clipped to zero.
#
fit_bch <- function(model_state, X, Y) {

  if (inherits(model_state$sm, c("covariate", "distal_regression", "distal_pooled"))) {
    warning(paste(
      "BCH correction is not recommended for covariates or categorical distal outcomes.",
      "Vermunt (2010) and Bakk et al. (2013) show ML correction is preferred in this case.",
      "Consider correction = 'ML' instead."
    ))
  }

  weights <- model_state$sample_weights
  e_res   <- e_step(model_state, X, NULL)
  resp    <- exp(e_res$log_resp)          # n x K posterior probabilities

  # Classification error matrix C (proportional base).
  # C[j, k] = P(assigned approx j | true = k), built from soft assignment.
  # Columns index the true class and sum to 1.
  C <- t(resp) %*% (resp * weights)
  C <- sweep(C, 2, colSums(resp * weights), "/")

  # BCH weight matrix: D = t(C^{-1})
  D <- t(pinv(C))

  # Apply to proportional weights. Negative weights are retained intentionally.
  bch_resp <- resp %*% D

  model_state$sm <- init_params(model_state$sm, Y, bch_resp)
  model_state$sm <- m_step(model_state$sm, Y, bch_resp)

  # Compute and store the full sandwich variance-covariance matrix of the class
  # means. This only applies to distal_continuous (which stores $means);
  # categorical models do not store $means and are excluded by the guard below.
  # Raw BCH weights are used here (before sample-weight multiplication), as
  # required by the Wald test for equality of means.
  if (inherits(model_state$sm, "distal_continuous") &&
      !is.null(model_state$sm$parameters$means)) {
    Y_vec    <- as.numeric(Y[, 1])
    K        <- model_state$n_components
    mu       <- as.vector(model_state$sm$parameters$means)
    Nk       <- colSums(bch_resp)
    Sigma_mu <- matrix(0, K, K)
    for (j in 1:K) {
      for (k in j:K) {
        cov_jk <- sum(bch_resp[, j] * bch_resp[, k] *
                        (Y_vec - mu[j]) * (Y_vec - mu[k])) /
          (Nk[j] * Nk[k])
        Sigma_mu[j, k] <- cov_jk
        Sigma_mu[k, j] <- cov_jk
      }
    }
    model_state$sm$parameters$Sigma_mu <- Sigma_mu
  }

  # For distal_continuous_pooled: store the full model-based covariance matrix
  # of theta = [intercepts, slopes] for use in the omnibus Wald test.
  if (inherits(model_state$sm, "distal_continuous_pooled") &&
      !is.null(model_state$sm$parameters$beta_pooled)) {
    sigma2_p <- model_state$sm$parameters$covariances[1, 1]
    L_p      <- length(as.vector(model_state$sm$parameters$beta_pooled))
    K_p      <- model_state$sm$n_components
    D_cov    <- L_p - K_p
    Y_p      <- as.numeric(Y[, 1])
    Z_p      <- if (D_cov > 0) as.matrix(Y[, -1, drop = FALSE]) else
      matrix(0, nrow = length(Y_p), ncol = 0)
    N_p      <- length(Y_p)
    U_p      <- matrix(0, N_p * K_p, L_p)
    for (k in seq_len(K_p)) {
      idx_p         <- ((k - 1L) * N_p + 1L):(k * N_p)
      U_p[idx_p, k] <- 1
      if (D_cov > 0) U_p[idx_p, (K_p + 1L):L_p] <- Z_p
    }
    W_p    <- as.vector(bch_resp)
    UWU_p  <- t(U_p) %*% sweep(U_p, 1, W_p, "*")
    diag(UWU_p) <- diag(UWU_p) + 1e-6
    B_inv_p <- pinv(UWU_p)
    model_state$sm$parameters$cov_theta <- sigma2_p * B_inv_p
  }

  return(model_state)
}

# Apply the Maximum Likelihood (ML) Correction (Vermunt, 2010; Bolck et al., 2004)
#
# The core identity (Vermunt 2010, eq. 13) is:
#   P(a_i | z_i) = sum_k P(x=k | z_i) * P(a_i | x=k)
#
# P(a_i | x=k) is fixed from the proportional classification table built from
# step-1 posteriors:
#   C_prop[j,k]     = sum_i w_i * resp1[i,j] * resp1[i,k]   (K x K, symmetric)
#   C_row_norm[j,k] = C_prop[j,k] / sum_k C_prop[j,k]        = P(a=k | x=j)
#
# For individual i under proportional (soft) assignment:
#   P(a_i | x=j) = sum_k resp1[i,k] * C_row_norm[j,k]
#
# The EM for the expanded dataset (K records per person, weight resp1[i,k]) gives:
#   Z_mat[i,k] = sum_j P(x=j | z_i) * C_row_norm[j,k]       (normaliser)
#   W[i,j]     = P(x=j | z_i) * (R %*% t(C_row_norm))[i,j]  (posterior weights)
#   where R[i,k] = resp1[i,k] / Z_mat[i,k]
#
# W sums to 1 per row and is passed to m_step as classification weights.
# Convergence LL = sum_i w_i * sum_k resp1[i,k] * log Z_mat[i,k].
fit_ml <- function(model_state, X, Y, max_iter = 1000, abs_tol = 1e-10, rel_tol = 1e-10) {

  if (inherits(model_state$sm, c("distal_continuous", "distal_continuous_regression"))) {
    warning(paste(
      "ML correction is not recommended for continuous distal outcomes.",
      "Bakk & Vermunt (2016) show BCH is preferred in this case,",
      "as ML requires strong distributional assumptions.",
      "Consider correction = 'BCH' instead."
    ))
  }

  keep    <- complete.cases(Y)
  X_clean <- X[keep, , drop = FALSE]
  Y_clean <- Y[keep, , drop = FALSE]
  w_clean <- model_state$sample_weights[keep]

  if (nrow(X_clean) == 0) stop("All rows have missing covariates.")

  K <- model_state$n_components

  # Step 1: compute and freeze posteriors from the measurement model.
  e_res_step1 <- e_step(model_state, X_clean, NULL)
  resp_step1  <- exp(e_res_step1$log_resp)        # n_clean x K

  # Proportional classification error matrix.
  C_prop     <- t(resp_step1 * w_clean) %*% resp_step1   # K x K (symmetric)
  Nk         <- colSums(resp_step1 * w_clean)
  C_row_norm <- sweep(C_prop, 1, Nk, "/")                 # K x K, row-normalised

  # Initialise structural model.
  model_state$sm <- init_params(model_state$sm, Y_clean, resp_step1)
  model_state$sm <- m_step(model_state$sm, Y_clean, resp_step1, weights = w_clean)

  prev_ll <- -Inf

  for (iter in seq_len(max_iter)) {

    log_sm <- log_likelihood(model_state$sm, Y_clean)

    # E-step.
    #
    # For discrete outcomes (distal_pooled, distal_regression), the correct
    # ML EM formulation following Vermunt (2010, eq. 14) uses raw outcome
    # probabilities P(o_i|x=j) = exp(log_sm)[i,j] without row-normalisation:
    #   Z_mat[i,k] = sum_j pi_j * P(o_i|x=j) * C_row_norm[j,k]
    #   W_eff[i,j] = pi_j * P(o_i|x=j) * sum_k [resp1[i,k]/Z[i,k]] * C_norm[j,k]
    #   LL         = sum_i sum_k resp1[i,k] * log Z_mat[i,k]
    #
    # For continuous outcomes (distal_continuous*), log_sm returns
    # log-likelihoods on the correct scale and row-normalisation is harmless
    # because the Gaussian density cancels in the ratio.

    if (inherits(model_state$sm, c("distal_pooled", "distal_regression"))) {
      po_given_x <- exp(log_sm)                      # n_clean x K
      pi_k_clean <- colSums(resp_step1 * w_clean) /
        sum(w_clean)                                  # K (weighted class props)
      Z_mat      <- sweep(po_given_x, 2, pi_k_clean, "*") %*% C_row_norm
      current_ll <- sum(w_clean * rowSums(
        resp_step1 * log(pmax(Z_mat, 1e-300))))
      RC <- resp_step1 / pmax(Z_mat, 1e-300)
      W  <- sweep(po_given_x, 2, pi_k_clean, "*") *
        (RC %*% t(C_row_norm))
    } else {
      lsh     <- apply(log_sm, 1, max)
      sm_prob <- exp(log_sm - lsh)
      sm_prob <- sm_prob / rowSums(sm_prob)
      Z_mat   <- sm_prob %*% C_row_norm
      current_ll <- sum(w_clean * rowSums(
        resp_step1 * log(pmax(Z_mat, 1e-300))))
      R <- resp_step1 / pmax(Z_mat, 1e-300)
      W <- sm_prob * (R %*% t(C_row_norm))
    }

    if (iter > 1) {
      change <- current_ll - prev_ll
      denom  <- max(abs(prev_ll), 1e-9)
      if (is.na(change) || abs(change) < abs_tol || abs(change / denom) < rel_tol) break
    }
    prev_ll <- current_ll

    model_state$sm <- m_step(model_state$sm, Y_clean, W, weights = w_clean)
  }

  # Robust sandwich standard errors for discrete structural models.
  #
  # For distal_pooled and distal_regression under ML step-3, the Q-function
  # Hessian stored by m_step reflects only the expected complete-data curvature,
  # not the marginal LL curvature, which leads to underestimation of variance.
  #
  # The robust sandwich estimator (Bakk, Oberski & Vermunt, 2014) is used:
  #   V_robust = B^{-1} M B^{-1}
  # where:
  #   B = -H_marg  (numerical Hessian of the marginal LL at convergence)
  #   M = sum_i s_i s_i^T  (outer product of person-level marginal scores)
  #
  # The marginal LL is: L(theta) = sum_i sum_k resp1[i,k] * log Z_mat[i,k].
  # The Hessian slot is replaced with -V_robust^{-1} so that downstream
  # inference code, which calls pinv(-hessian) to obtain the variance matrix,
  # recovers V_robust.
  if (inherits(model_state$sm, c("distal_pooled", "distal_regression"))) {
    theta0 <- as.vector(model_state$sm$parameters$beta_pooled)
    if (is.null(theta0))
      theta0 <- as.vector(model_state$sm$parameters$betas)
    n_theta <- length(theta0)
    eps_nd  <- 1e-4

    # Capture the true dimensions of beta_pooled before entering the closures.
    # beta_pooled is an (M-1) x (K + D_cov) matrix; restoring it from a flat
    # vector requires the correct row count (M-1), which equals 1 only for
    # binary outcomes and is larger for polytomous ones.
    beta_pooled_dim <- if (!is.null(model_state$sm$parameters$beta_pooled))
      dim(model_state$sm$parameters$beta_pooled)
    else NULL

    # Dimensions of the betas array for distal_regression: c(K, M-1, D).
    betas_dim <- if (inherits(model_state$sm, "distal_regression"))
      dim(model_state$sm$parameters$betas)
    else NULL

    # Marginal log-likelihood as a function of the structural parameters.
    marg_ll_fn <- function(theta) {
      sm_tmp <- model_state$sm
      if (!is.null(beta_pooled_dim)) {
        sm_tmp$parameters$beta_pooled <-
          matrix(theta, nrow = beta_pooled_dim[1], ncol = beta_pooled_dim[2])
      } else if (!is.null(betas_dim)) {
        sm_tmp$parameters$betas <- array(theta, dim = betas_dim)
      }
      log_s <- log_likelihood(sm_tmp, Y_clean)
      po_x  <- exp(log_s)
      Z_m   <- sweep(po_x, 2, pi_k_clean, "*") %*% C_row_norm
      sum(w_clean * rowSums(resp_step1 * log(pmax(Z_m, 1e-300))))
    }

    # Numerical Hessian of the marginal LL (central differences, O(eps^2)).
    H_marg <- matrix(0, n_theta, n_theta)
    for (j in seq_len(n_theta)) {
      for (k in j:n_theta) {
        ej <- rep(0, n_theta); ej[j] <- eps_nd
        ek <- rep(0, n_theta); ek[k] <- eps_nd
        H_marg[j, k] <- H_marg[k, j] <-
          (marg_ll_fn(theta0 + ej + ek) - marg_ll_fn(theta0 + ej - ek) -
             marg_ll_fn(theta0 - ej + ek) + marg_ll_fn(theta0 - ej - ek)) /
          (4 * eps_nd^2)
      }
    }

    # Person-level marginal scores (numerical gradient).
    score_mat <- matrix(0, nrow(Y_clean), n_theta)
    for (j in seq_len(n_theta)) {
      ej   <- rep(0, n_theta); ej[j] <- eps_nd
      sm_p <- model_state$sm
      sm_m <- model_state$sm
      if (!is.null(beta_pooled_dim)) {
        sm_p$parameters$beta_pooled <-
          matrix(theta0 + ej, nrow = beta_pooled_dim[1], ncol = beta_pooled_dim[2])
        sm_m$parameters$beta_pooled <-
          matrix(theta0 - ej, nrow = beta_pooled_dim[1], ncol = beta_pooled_dim[2])
      } else if (!is.null(betas_dim)) {
        sm_p$parameters$betas <- array(theta0 + ej, dim = betas_dim)
        sm_m$parameters$betas <- array(theta0 - ej, dim = betas_dim)
      }
      log_sp  <- log_likelihood(sm_p, Y_clean)
      log_sm_m <- log_likelihood(sm_m, Y_clean)
      Zp <- sweep(exp(log_sp),  2, pi_k_clean, "*") %*% C_row_norm
      Zm <- sweep(exp(log_sm_m), 2, pi_k_clean, "*") %*% C_row_norm
      # Person i score for parameter j:
      # s_i[j] = sum_k resp1[i,k] * (log Zp[i,k] - log Zm[i,k]) / (2 * eps)
      score_mat[, j] <- rowSums(
        resp_step1 * (log(pmax(Zp, 1e-300)) - log(pmax(Zm, 1e-300)))
      ) / (2 * eps_nd)
    }

    # Meat = sum_i s_i s_i^T
    meat <- t(score_mat) %*% score_mat

    # Robust sandwich variance: V = B^{-1} M B^{-1} where B = -H_marg.
    B_inv    <- pinv(-H_marg)
    V_robust <- B_inv %*% meat %*% B_inv

    # Store -V_robust^{-1} so that pinv(-hessian) returns V_robust downstream.
    model_state$sm$parameters$hessian <- -pinv(V_robust)
  }

  e_res_full <- e_step(model_state, X, NULL)
  resp1_full <- exp(e_res_full$log_resp)

  p_a_gvn_x_full <- resp1_full %*% t(C_row_norm)

  sm_logp_full <- matrix(0, nrow = nrow(X), ncol = K)
  if (any(keep)) {
    log_sm_f <- log_likelihood(model_state$sm, Y_clean)
    lsh_f    <- apply(log_sm_f, 1, max)
    sp_f     <- exp(log_sm_f - lsh_f)
    sp_f     <- sp_f / rowSums(sp_f)
    sm_logp_full[keep, ] <- log(pmax(sp_f, 1e-300))
  }

  log_comb_full <- log(pmax(p_a_gvn_x_full, 1e-300)) + sm_logp_full
  max_lcf       <- apply(log_comb_full, 1, max)
  log_norm_full <- max_lcf + log(rowSums(exp(sweep(log_comb_full, 1, max_lcf, "-"))))

  model_state$log_resp    <- sweep(log_comb_full, 1, log_norm_full, "-")
  model_state$lower_bound <- log_norm_full

  return(model_state)
}
