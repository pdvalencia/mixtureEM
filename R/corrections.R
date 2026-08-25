# ==============================================================================
# S3 3-Step Corrections (BCH and ML)
# ==============================================================================

get_modal_resp <- function(resp) {
  n <- nrow(resp)
  K <- ncol(resp)
  modal <- matrix(0, nrow = n, ncol = K)
  assigned_classes <- max.col(resp, ties.method = "first")
  # One-hot the modal class by index rather than a row loop. Also correct for
  # an empty `resp`, where `1:n` would have run over rows 1 and 0.
  modal[cbind(seq_len(n), assigned_classes)] <- 1
  return(modal)
}

# Apply the BCH Correction
#
# Implementation follows the proportional BCH correction of Vermunt (2010)
# and Bakk, Tekle & Vermunt (2013):
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
fit_bch <- function(model_state, X, Y, assignment = "proportional") {

  if (inherits(model_state$sm, c("covariate", "distal_regression", "distal_pooled"))) {
    warning(paste(
      "BCH correction is not recommended for covariates or categorical distal outcomes.",
      "Vermunt (2010) and Bakk et al. (2013) show ML correction is preferred in this case,",
      "Consider correction = 'ML' instead."
    ))
  }

  # Complete missing covariates under the class-invariant Gaussian marginal
  # (endogenous-constrained-x; Sterba, 2014) so all cases are retained; a
  # missing distal outcome stays NA and is masked by the structural likelihood.
  Y <- complete_structural_covariates(model_state$sm, Y)

  weights <- model_state$sample_weights
  e_res   <- e_step(model_state, X, NULL)
  resp    <- exp(e_res$log_resp)          # n x K posterior probabilities

  # The assignment rule decides what the classification-error matrix is a table
  # of. Under proportional assignment a case contributes its posterior
  # probability to every class; under modal it contributes 1 to its most likely
  # class and 0 elsewhere (Bolck, Croon & Hagenaars, 2004; Vermunt, 2010). Both
  # are consistent; they differ in finite samples because the modal table sits
  # closer to the identity, so inverting it applies a smaller correction.
  A <- if (assignment == "modal") get_modal_resp(resp) else resp

  # Classification error matrix C.
  # C[j, k] = P(assigned approx j | true = k). Columns index the true class and
  # sum to 1. The true class is never observed, so its column weights are always
  # the posteriors; only the assigned variable A changes with the rule. Writing
  # it as A crossed with `resp` rather than A crossed with itself matters: a
  # modal A has a single 1 per row, so t(A) %*% A is diagonal and would make C
  # the identity, turning the correction into a no-op. Under proportional
  # assignment A *is* `resp` and the two expressions are the same matrix.
  C <- t(A) %*% (resp * weights)
  C <- sweep(C, 2, colSums(resp * weights), "/")

  # BCH weight matrix: D = t(C^{-1})
  D <- t(pinv(C))

  # Apply to the assignment weights. Negative weights are retained intentionally.
  bch_resp <- A %*% D

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

    # Restrict to cases with an observed outcome. The m_step above masks the
    # same rows, so the means this covariance matrix belongs to are estimated
    # from these cases only; Nk must be recomputed on the mask rather than over
    # all rows, or the denominators would not match the weights in the sums.
    obs     <- is.finite(Y_vec)
    bch_obs <- bch_resp[obs, , drop = FALSE]
    Y_obs   <- Y_vec[obs]
    Nk      <- colSums(bch_obs)

    if (sum(obs) < 2L || any(abs(Nk) < 1e-8)) {
      # Leave Sigma_mu unset: .wald_omnibus_means() falls back to the SEs when
      # it is NULL, whereas a matrix of zeros would be taken at face value.
      model_state$sm$parameters$Sigma_mu <- NULL
    } else {
      Sigma_mu <- matrix(0, K, K)
      for (j in 1:K) {
        for (k in j:K) {
          cov_jk <- sum(bch_obs[, j] * bch_obs[, k] *
                          (Y_obs - mu[j]) * (Y_obs - mu[k])) /
            (Nk[j] * Nk[k])
          Sigma_mu[j, k] <- cov_jk
          Sigma_mu[k, j] <- cov_jk
        }
      }
      model_state$sm$parameters$Sigma_mu <- Sigma_mu
    }
  }

  # distal_continuous_pooled's cov_theta used to be rebuilt here, and that copy
  # of the design was built over every row rather than over the rows with an
  # observed outcome that m_step actually fits -- so on data with a partly
  # missing outcome the covariance carried far more information than the
  # estimates did. m_step() now builds it alongside the coefficients, from the
  # same U, W and mask, and there is no second copy to drift.

  return(model_state)
}

# Apply the Maximum Likelihood (ML) Correction (Vermunt, 2010; Bolck et al., 2004)
#
# The core identity (Vermunt 2010, eq. 13) is:
#   P(a_i | z_i) = sum_k P(x=k | z_i) * P(a_i | x=k)
#
# P(a_i | x=k) is fixed from the proportional classification table built from
# step-1 posteriors:
#   C_prop[j,k]     = sum_i w_i * resp1[i,j] * A1[i,k]      (K x K)
#   C_row_norm[j,k] = C_prop[j,k] / sum_k C_prop[j,k]        = P(a=k | x=j)
#
# A1 is the assigned-class variable: resp1 itself under proportional assignment
# (where C_prop is symmetric), or the modal indicator under modal.
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
fit_ml <- function(model_state, X, Y, max_iter = 1000, abs_tol = 1e-10,
                   rel_tol = 1e-10, se = "corrected",
                   assignment = "proportional") {

  if (inherits(model_state$sm, c("distal_continuous", "distal_continuous_regression"))) {
    warning(paste(
      "ML correction is not recommended for continuous distal outcomes.",
      "Bakk & Vermunt (2016) show BCH is preferred in this case,",
      "as ML requires strong distributional assumptions.",
      "Consider correction = 'BCH' instead."
    ))
  }

  # Retain every case. Missing covariates are completed under the class-invariant
  # Gaussian marginal (endogenous-constrained-x; Sterba, 2014) instead of being
  # listwise deleted, which sacrifices efficiency and biases estimates under MAR
  # (Little, 1992; Little & Zhang, 2011). A missing distal outcome is left as NA
  # and masked by the structural model's own FIML likelihood.
  Y       <- complete_structural_covariates(model_state$sm, Y)
  keep    <- rep(TRUE, nrow(Y))
  X_clean <- X[keep, , drop = FALSE]
  Y_clean <- Y[keep, , drop = FALSE]
  w_clean <- model_state$sample_weights[keep]

  # Keep the sub-model's survey design row-aligned with the retained cases so
  # that variance code running inside m_step (e.g. the covariate M-step) sees
  # strata and cluster vectors that match Y_clean.
  if (!is.null(model_state$sm$strata)) {
    model_state$sm$strata  <- model_state$sm$strata[keep]
    model_state$sm$cluster <- model_state$sm$cluster[keep]
  }

  if (nrow(X_clean) == 0) stop("All rows have missing covariates.")

  K <- model_state$n_components

  # Step 1: compute and freeze posteriors from the measurement model.
  e_res_step1 <- e_step(model_state, X_clean, NULL)
  resp_step1  <- exp(e_res_step1$log_resp)        # n_clean x K

  # The assigned-class variable whose classification error the correction
  # inverts. Under proportional assignment it is the posterior itself; under
  # modal it is the indicator of the most likely class (Bolck, Croon &
  # Hagenaars, 2004; Vermunt, 2010). Everything downstream of the frozen
  # posteriors reads A1.
  A1 <- if (assignment == "modal") get_modal_resp(resp_step1) else resp_step1

  # Classification error matrix. Rows index the true class, which is never
  # observed and so is always weighted by the posteriors; columns index the
  # assigned variable A1, which is what the rule changes. As in fit_bch(), A1
  # must be crossed with `resp_step1` and not with itself -- a modal A1 has one
  # 1 per row, so t(A1) %*% A1 is diagonal and C_row_norm would come out as the
  # identity, leaving the correction with nothing to invert. Under proportional
  # assignment A1 is `resp_step1` and this is the same matrix as before.
  C_prop     <- t(resp_step1 * w_clean) %*% A1            # K x K
  Nk         <- colSums(resp_step1 * w_clean)
  C_row_norm <- sweep(C_prop, 1, Nk, "/")                 # K x K, row-normalised

  # Initialise structural model.
  model_state$sm <- init_params(model_state$sm, Y_clean, A1)
  model_state$sm <- m_step(model_state$sm, Y_clean, A1, weights = w_clean)

  # Weighted class proportions; computed once for discrete models and reused
  # both in the EM loop and in the variance estimation block below.
  if (inherits(model_state$sm, c("distal_pooled", "distal_regression"))) {
    pi_k_clean <- colSums(A1 * w_clean) / sum(w_clean)
  }

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
      Z_mat      <- sweep(po_given_x, 2, pi_k_clean, "*") %*% C_row_norm
      current_ll <- sum(w_clean * rowSums(
        A1 * log(pmax(Z_mat, 1e-300))))
      RC <- A1 / pmax(Z_mat, 1e-300)
      W  <- sweep(po_given_x, 2, pi_k_clean, "*") *
        (RC %*% t(C_row_norm))
    } else {
      lsh     <- .row_max(log_sm)
      sm_prob <- exp(log_sm - lsh)
      sm_prob <- sm_prob / rowSums(sm_prob)
      Z_mat   <- sm_prob %*% C_row_norm
      current_ll <- sum(w_clean * rowSums(
        A1 * log(pmax(Z_mat, 1e-300))))
      R <- A1 / pmax(Z_mat, 1e-300)
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

  # ============================================================================
  # Variance estimation for structural models under ML step-3.
  #
  # The Q-function Hessian stored by m_step reflects only the expected
  # complete-data curvature, not the marginal LL curvature, which leads to
  # underestimation of variance.  The robust sandwich estimator is used instead
  # (Bakk, Oberski & Vermunt, 2014):
  #
  #   V = B^{-1} M B^{-1}
  #
  # where B = -H_marg is the Hessian of the marginal log-likelihood
  #   L(theta) = sum_i sum_k resp1[i,k] * log Z_mat[i,k]
  # and M = sum_i s_i s_i^T is the outer product of person-level scores.
  #
  # The covariate case is delegated to R/step3_variance.R, which computes B and
  # the scores in closed form rather than by numerical differentiation and adds
  # the second variance component the two distal branches below still omit: the
  # uncertainty carried over from step 1, which is estimated and not known.
  #
  # For distal_pooled the parameters form one joint vector and the corrected
  # variance is stored in $hessian (singular), which is where downstream
  # inference reads from.
  #
  # For distal_regression the parameters are class-specific; the marginal LL
  # Hessian is block-diagonal across classes.  The per-class sandwich is
  # computed independently for each class k and stored in $hessians[[k]]
  # (plural list), which is where downstream inference reads from.
  # ============================================================================

  # --- covariate (class-membership regression) --------------------------------
  #
  # Downstream inference reads $V_robust in preference to $hessian, so attaching
  # it here is enough.
  if (inherits(model_state$sm, "covariate")) {
    model_state <- .attach_step3_covariate_vcov(
      model_state, X_clean, Y_clean, A1, C_row_norm, w_clean, se = se)
  }

  # A covariate block inside a nested structural model shares its step-three
  # likelihood with the distal block, so the single-block formulas above do not
  # apply. The naive Hessian stands, but it is labelled as such rather than
  # reported as though it were the corrected variance.
  if (inherits(model_state$sm, "nested") &&
      !is.null(model_state$sm$models[["predictor"]]) &&
      inherits(model_state$sm$models$predictor, "covariate")) {
    model_state$sm$models$predictor$parameters$V_method <-
      "Q-function Hessian (uncorrected; covariate combined with a distal outcome)"
  }

  # --- distal_pooled (joint parameter vector, one shared Hessian) -------------

  if (inherits(model_state$sm, "distal_pooled")) {
    theta0          <- as.vector(model_state$sm$parameters$beta_pooled)
    n_theta         <- length(theta0)
    eps_nd          <- 1e-4
    beta_pooled_dim <- dim(model_state$sm$parameters$beta_pooled)

    marg_ll_fn <- function(theta) {
      sm_tmp <- model_state$sm
      sm_tmp$parameters$beta_pooled <-
        matrix(theta, nrow = beta_pooled_dim[1], ncol = beta_pooled_dim[2])
      po_x <- exp(log_likelihood(sm_tmp, Y_clean))
      Z_m  <- sweep(po_x, 2, pi_k_clean, "*") %*% C_row_norm
      sum(w_clean * rowSums(A1 * log(pmax(Z_m, 1e-300))))
    }

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

    score_mat <- matrix(0, nrow(Y_clean), n_theta)
    for (j in seq_len(n_theta)) {
      ej   <- rep(0, n_theta); ej[j] <- eps_nd
      sm_p <- model_state$sm; sm_m <- model_state$sm
      sm_p$parameters$beta_pooled <-
        matrix(theta0 + ej, nrow = beta_pooled_dim[1], ncol = beta_pooled_dim[2])
      sm_m$parameters$beta_pooled <-
        matrix(theta0 - ej, nrow = beta_pooled_dim[1], ncol = beta_pooled_dim[2])
      Zp <- sweep(exp(log_likelihood(sm_p, Y_clean)), 2, pi_k_clean, "*") %*% C_row_norm
      Zm <- sweep(exp(log_likelihood(sm_m, Y_clean)), 2, pi_k_clean, "*") %*% C_row_norm
      score_mat[, j] <- rowSums(
        A1 * (log(pmax(Zp, 1e-300)) - log(pmax(Zm, 1e-300)))
      ) / (2 * eps_nd)
    }

    meat     <- compute_survey_B(score_mat,
                                 model_state$strata[keep],
                                 model_state$cluster[keep])
    B_inv    <- pinv(-H_marg)
    V_robust <- B_inv %*% meat %*% B_inv

    model_state$sm$parameters$hessian <- -pinv(V_robust)
  }

  # --- distal_regression (class-specific parameters, per-class Hessians) ------
  #
  # The marginal LL is block-diagonal in the class-specific parameter blocks
  # because P(o_i | x=k) depends only on betas[k,,].  Each class k's sandwich
  # is therefore computed independently from the others.

  if (inherits(model_state$sm, "distal_regression")) {
    betas_dim   <- dim(model_state$sm$parameters$betas)   # c(K, M-1, D)
    K_dr        <- betas_dim[1]
    n_per_class <- betas_dim[2] * betas_dim[3]             # (M-1) * D
    eps_nd      <- 1e-4

    for (kk in seq_len(K_dr)) {
      theta_k <- as.vector(model_state$sm$parameters$betas[kk, , ])

      marg_ll_k <- function(theta) {
        sm_tmp <- model_state$sm
        sm_tmp$parameters$betas[kk, , ] <-
          array(theta, dim = betas_dim[2:3])
        po_x <- exp(log_likelihood(sm_tmp, Y_clean))
        Z_m  <- sweep(po_x, 2, pi_k_clean, "*") %*% C_row_norm
        sum(w_clean * rowSums(A1 * log(pmax(Z_m, 1e-300))))
      }

      H_k <- matrix(0, n_per_class, n_per_class)
      for (j in seq_len(n_per_class)) {
        for (l in j:n_per_class) {
          ej <- rep(0, n_per_class); ej[j] <- eps_nd
          el <- rep(0, n_per_class); el[l] <- eps_nd
          H_k[j, l] <- H_k[l, j] <-
            (marg_ll_k(theta_k + ej + el) - marg_ll_k(theta_k + ej - el) -
               marg_ll_k(theta_k - ej + el) + marg_ll_k(theta_k - ej - el)) /
            (4 * eps_nd^2)
        }
      }

      score_k <- matrix(0, nrow(Y_clean), n_per_class)
      for (j in seq_len(n_per_class)) {
        ej   <- rep(0, n_per_class); ej[j] <- eps_nd
        sm_p <- model_state$sm; sm_m <- model_state$sm
        sm_p$parameters$betas[kk, , ] <- array(theta_k + ej, dim = betas_dim[2:3])
        sm_m$parameters$betas[kk, , ] <- array(theta_k - ej, dim = betas_dim[2:3])
        Zp <- sweep(exp(log_likelihood(sm_p, Y_clean)), 2, pi_k_clean, "*") %*% C_row_norm
        Zm <- sweep(exp(log_likelihood(sm_m, Y_clean)), 2, pi_k_clean, "*") %*% C_row_norm
        score_k[, j] <- rowSums(
          A1 * (log(pmax(Zp, 1e-300)) - log(pmax(Zm, 1e-300)))
        ) / (2 * eps_nd)
      }

      meat_k  <- compute_survey_B(score_k,
                                  model_state$strata[keep],
                                  model_state$cluster[keep])
      B_inv_k <- pinv(-H_k)
      V_k     <- B_inv_k %*% meat_k %*% B_inv_k

      model_state$sm$parameters$hessians[[kk]] <- -pinv(V_k)
    }
  }

  e_res_full <- e_step(model_state, X, NULL)
  resp1_full <- exp(e_res_full$log_resp)

  p_a_gvn_x_full <- resp1_full %*% t(C_row_norm)

  sm_logp_full <- matrix(0, nrow = nrow(X), ncol = K)
  if (any(keep)) {
    log_sm_f <- log_likelihood(model_state$sm, Y_clean)
    lsh_f    <- .row_max(log_sm_f)
    sp_f     <- exp(log_sm_f - lsh_f)
    sp_f     <- sp_f / rowSums(sp_f)
    sm_logp_full[keep, ] <- log(pmax(sp_f, 1e-300))
  }

  log_comb_full <- log(pmax(p_a_gvn_x_full, 1e-300)) + sm_logp_full
  max_lcf       <- .row_max(log_comb_full)
  log_norm_full <- max_lcf + log(rowSums(exp(sweep(log_comb_full, 1, max_lcf, "-"))))

  model_state$log_resp    <- sweep(log_comb_full, 1, log_norm_full, "-")
  model_state$lower_bound <- log_norm_full

  return(model_state)
}
