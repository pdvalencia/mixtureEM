# ==============================================================================
# S3 User Wrappers and Pipeline Tools (LCA/LPA Mixture Engine)
# ==============================================================================

sort_model_classes <- function(model_state) {
  K <- model_state$n_components
  if (K <= 1) return(model_state)

  new_order <- order(model_state$weights, decreasing = TRUE)
  model_state$weights <- model_state$weights[new_order]

  # --- Sort flat measurement model parameters ---
  if (!is.null(model_state$mm$parameters[["pis"]])) {
    model_state$mm$parameters$pis <-
      model_state$mm$parameters$pis[new_order, , drop = FALSE]
  }
  # Gaussian LPA: also sort means and covariances
  if (!is.null(model_state$mm$parameters[["means"]])) {
    model_state$mm$parameters$means <-
      model_state$mm$parameters$means[new_order, , drop = FALSE]
  }
  if (!is.null(model_state$mm$parameters[["covariances"]])) {
    model_state$mm$parameters$covariances <-
      model_state$mm$parameters$covariances[new_order, , drop = FALSE]
  }

  # --- Sort nested measurement model sub-model parameters ---
  # The flat-parameter block above only touches model_state$mm$parameters, which
  # is empty for nested models.  Sub-model parameters live one level deeper at
  # model_state$mm$models[[name]]$parameters and must be sorted independently.
  if (inherits(model_state$mm, "nested")) {
    for (name in names(model_state$mm$models)) {
      sub <- model_state$mm$models[[name]]
      if (!is.null(sub$parameters[["pis"]]))
        sub$parameters$pis <- sub$parameters$pis[new_order, , drop = FALSE]
      if (!is.null(sub$parameters[["means"]]))
        sub$parameters$means <- sub$parameters$means[new_order, , drop = FALSE]
      if (!is.null(sub$parameters[["covariances"]]))
        sub$parameters$covariances <-
          sub$parameters$covariances[new_order, , drop = FALSE]
      model_state$mm$models[[name]] <- sub
    }
  }

  # --- Sort structural model parameters ---
  if (!is.null(model_state$sm)) {

    sort_sm_params <- function(sm) {
      if (!is.null(sm$parameters[["beta"]])) {
        sm$parameters$beta <- sm$parameters$beta[new_order, , drop = FALSE]
        if (!is.null(sm$parameters[["hessian"]])) {
          H       <- sm$parameters$hessian
          D       <- ncol(sm$parameters$beta)
          idx_map <- as.vector(sapply(new_order, function(k) ((k-1)*D + 1):(k*D)))
          sm$parameters$hessian <- H[idx_map, idx_map, drop = FALSE]
        }
        # The survey-robust covariance is blocked by class in the same layout
        # as the Hessian, so it is permuted with the identical index map.
        if (!is.null(sm$parameters[["V_robust"]])) {
          Vr      <- sm$parameters$V_robust
          D       <- ncol(sm$parameters$beta)
          idx_map <- as.vector(sapply(new_order, function(k) ((k-1)*D + 1):(k*D)))
          sm$parameters$V_robust <- Vr[idx_map, idx_map, drop = FALSE]
        }
      }
      if (!is.null(sm$parameters[["pis"]])) {
        sm$parameters$pis <- sm$parameters$pis[new_order, , drop = FALSE]
      }
      if (!is.null(sm$parameters[["means"]])) {
        sm$parameters$means <- sm$parameters$means[new_order, , drop = FALSE]
      }
      if (!is.null(sm$parameters[["covariances"]])) {
        sm$parameters$covariances <-
          sm$parameters$covariances[new_order, , drop = FALSE]
      }
      if (!is.null(sm$parameters[["ses"]])) {
        if (inherits(sm, "distal_continuous_pooled")) {
          sm$parameters$ses[1, 1:K] <- sm$parameters$ses[1, new_order]
        } else {
          sm$parameters$ses <- sm$parameters$ses[new_order, , drop = FALSE]
        }
      }
      if (!is.null(sm$parameters[["Sigma_mu"]])) {
        sm$parameters$Sigma_mu <-
          sm$parameters$Sigma_mu[new_order, new_order, drop = FALSE]
      }
      if (!is.null(sm$parameters[["cov_theta"]])) {
        # cov_theta is L x L where the first K rows/cols are intercepts.
        # Reorder only the intercept block; slope block stays unchanged.
        K_ct  <- sm$n_components
        L_ct  <- nrow(sm$parameters$cov_theta)
        D_ct  <- L_ct - K_ct
        idx   <- c(new_order, if (D_ct > 0) (K_ct + seq_len(D_ct)) else integer(0))
        sm$parameters$cov_theta <-
          sm$parameters$cov_theta[idx, idx, drop = FALSE]
      }
      # Reorder the robust-sandwich hessian for distal_pooled / distal_regression.
      # This hessian is K x K (or larger) and is indexed by class in the
      # pre-sort order.  After sort_model_classes reorders beta_pooled, the
      # hessian must be permuted to match so that Sigma = pinv(-H) is aligned.
      if (!is.null(sm$parameters[["hessian"]]) &&
          inherits(sm, c("distal_pooled", "distal_regression"))) {
        H_sm  <- sm$parameters$hessian
        K_sm  <- sm$n_components
        L_sm  <- nrow(H_sm)
        D_sm  <- L_sm - K_sm
        idx_h <- c(new_order,
                   if (D_sm > 0) (K_sm + seq_len(D_sm)) else integer(0))
        sm$parameters$hessian <- H_sm[idx_h, idx_h, drop = FALSE]
      }
      if (!is.null(sm$parameters[["betas"]])) {
        if (length(dim(sm$parameters$betas)) == 3) {
          sm$parameters$betas <- sm$parameters$betas[new_order, , , drop = FALSE]
        } else {
          sm$parameters$betas <- sm$parameters$betas[new_order, , drop = FALSE]
        }
      }
      if (!is.null(sm$parameters[["beta_pooled"]])) {
        if (inherits(sm, "distal_continuous_pooled")) {
          sm$parameters$beta_pooled[1, 1:K] <- sm$parameters$beta_pooled[1, new_order]
        } else {
          # Guard: when beta_pooled is a degenerate 0x0 placeholder (produced by
          # a constant-outcome distal_pooled model), there is nothing to reorder.
          if (nrow(sm$parameters$beta_pooled) > 0 &&
              ncol(sm$parameters$beta_pooled) >= K) {
            sm$parameters$beta_pooled[, 1:K] <-
              sm$parameters$beta_pooled[, new_order, drop = FALSE]
          }
        }
      }
      return(sm)
    }

    if (inherits(model_state$sm, "nested")) {
      for (name in names(model_state$sm$models))
        model_state$sm$models[[name]] <- sort_sm_params(model_state$sm$models[[name]])
    } else {
      model_state$sm <- sort_sm_params(model_state$sm)
    }
  }

  if (!is.null(model_state$log_resp))
    model_state$log_resp <- model_state$log_resp[, new_order, drop = FALSE]

  return(model_state)
}

#' Print Measurement Model Parameters
#'
#' @description
#' Prints a formatted table of the fitted measurement model parameters:
#' item-response probabilities for categorical models, or means for Gaussian
#' models. Results are broken down by latent class. Handles both flat and
#' nested (mixed) measurement models.
#'
#' @param object A fitted \code{mixture_model} object returned by
#'   \code{\link{fit_mixture}}.
#'
#' @return Invisibly returns \code{NULL}. Called for its printed side-effect.
#'
#' @examples
#' set.seed(1)
#' X <- matrix(rbinom(500, 1, 0.5), nrow = 100)
#' fit <- fit_mixture(X, n_components = 2, measurement = "binary")
#' measurement_summary(fit)
#'
#' @export
measurement_summary <- function(object) {
  K <- object$n_components
  cat("=========================================================\n")
  cat("             MEASUREMENT MODEL PARAMETERS                \n")
  cat("=========================================================\n")

  print_item_matrix <- function(mat, title, sub_model = NULL) {
    cat(sprintf("\n%s\n", title))
    item_names <- colnames(mat)

    if (is.null(item_names)) {
      if (!is.null(sub_model) && !is.null(sub_model$max_val)) {
        M       <- sub_model$max_val
        n_items <- ncol(mat) / M
        base    <- sub_model$item_names
        if (is.null(base) || length(base) != n_items)
          base <- paste0("Poly_Item_", seq_len(n_items))
        item_names <- paste0(rep(base, each = M),
                             " (Cat ", rep(seq_len(M), times = n_items), ")")
      } else {
        item_names <- paste0("Item_", 1:ncol(mat))
      }
    }

    label_w <- max(20L, max(nchar(item_names)))
    cat(sprintf("%-*s", label_w, "Indicator"))
    for (k in 1:K) cat(sprintf(" | Class %d", k))
    cat("\n")
    cat(paste0(rep("-", label_w + K * 10), collapse = ""), "\n")

    for (j in 1:ncol(mat)) {
      cat(sprintf("%-*s", label_w, item_names[j]))
      for (k in 1:K) cat(sprintf(" | %7.3f", mat[k, j]))
      cat("\n")
    }
  }

  mm <- object$mm
  if (inherits(mm, "nested")) {
    for (name in names(mm$models)) {
      sub_mm <- mm$models[[name]]
      if (!is.null(sub_mm$parameters$pis))
        print_item_matrix(sub_mm$parameters$pis,
                          paste("Categorical Probabilities:", toupper(name)), sub_mm)
      if (!is.null(sub_mm$parameters$means))
        print_item_matrix(sub_mm$parameters$means,
                          paste("Continuous Means:", toupper(name)), sub_mm)
    }
  } else {
    if (!is.null(mm$parameters$pis))
      print_item_matrix(mm$parameters$pis, "CATEGORICAL PROBABILITIES", mm)
    if (!is.null(mm$parameters$means))
      print_item_matrix(mm$parameters$means, "CONTINUOUS MEANS", mm)
  }
  if (!is.null(object$missing_data) && isTRUE(object$missing_data$any_missing)) {
    md <- object$missing_data
    cat(sprintf("\nMissing data: %d of %d cells (%.1f%%) across %d item%s, handled via %s.\n",
                md$n_missing, md$n_cells, 100 * md$prop_missing,
                md$n_items_affected, if (md$n_items_affected == 1L) "" else "s",
                md$handled_by))
  }
  cat("=========================================================\n")
}

#' Print Classification Diagnostics
#'
#' @description
#' Computes and prints the Average Posterior Probability (AvePP) matrix.
#' Each row corresponds to observations modally assigned to a given class;
#' each column shows the mean posterior probability for that class. High
#' values on the diagonal (and low values off it) indicate well-separated,
#' clearly-assigned classes.
#'
#' @param object A fitted \code{mixture_model} object returned by
#'   \code{\link{fit_mixture}}.
#'
#' @return A K x K numeric matrix of average posterior probabilities,
#'   returned invisibly. The matrix is also printed to the console.
#'
#' @examples
#' set.seed(1)
#' X <- matrix(rbinom(500, 1, 0.5), nrow = 100)
#' fit <- fit_mixture(X, n_components = 2, measurement = "binary")
#' classification_diagnostics(fit)
#'
#' @export
classification_diagnostics <- function(object) {
  resp       <- exp(object$log_resp)
  pred_class <- max.col(resp)
  K          <- object$n_components

  ave_pp <- matrix(0, nrow = K, ncol = K)
  for (k in 1:K) {
    idx <- which(pred_class == k)
    if (length(idx) > 0)
      ave_pp[k, ] <- colMeans(resp[idx, , drop = FALSE])
  }
  rownames(ave_pp) <- paste("Assigned Class", 1:K)
  colnames(ave_pp) <- paste("Prob C", 1:K)

  cat("=========================================================\n")
  cat("          AVERAGE POSTERIOR PROBABILITIES (AvePP)        \n")
  cat("=========================================================\n")
  cat("Rows: Modal Assignment | Columns: Mean Probability\n\n")
  print(round(ave_pp, 3))
  cat("=========================================================\n")
  invisible(ave_pp)
}

# ==============================================================================
# Internal helpers for summary.mixture_model
# ==============================================================================

# Format p-values to publication conventions.
# Values below .001 are shown as "< .001"; NaN/NA appear as a dash.
.fmt_pval <- function(p) {
  if (is.na(p) || is.nan(p)) return("       -")
  if (p < 0.001)              return("  < .001")
  sprintf("   %5.3f", p)
}

# Omnibus Wald test for equality of K means (continuous distal outcomes).
#
# When Sigma_mu (the full K x K sandwich variance-covariance matrix of the
# means) is available, uses the full-covariance formulation:
#   W = c^T V^{-1} c,  where c = R * mu,  V = R * Sigma_mu * R^T
#   R = contrast matrix [class k vs class 1, k = 2..K]  (df = K-1)
#
# This matches LatentGOLD's robust Wald statistic (Bakk, Oberski & Vermunt,
# 2014) and accounts for the cross-class covariance induced by BCH weights.
#
# Falls back to the diagonal (precision-weighted) approximation when
# Sigma_mu is not stored (e.g. for non-BCH structural models).
#
# Returns a list with stat, df, and p.
.wald_omnibus_means <- function(means, ses, K, Sigma_mu = NULL) {
  if (K <= 1L) return(list(stat = NA_real_, df = 0L, p = NA_real_))
  df <- K - 1L

  if (!is.null(Sigma_mu) && all(is.finite(Sigma_mu))) {
    # Full sandwich Wald: contrast matrix R (K-1 x K), class 2..K vs class 1
    R     <- cbind(-1, diag(K - 1L))
    V     <- R %*% Sigma_mu %*% t(R)
    theta <- R %*% means
    W     <- tryCatch(
      as.numeric(t(theta) %*% solve(V) %*% theta),
      error = function(e) NA_real_
    )
  } else {
    # Fallback: diagonal precision-weighted approximation
    prec   <- 1 / pmax(ses^2, 1e-15)
    mu_bar <- sum(prec * means) / sum(prec)
    W      <- sum(prec * (means - mu_bar)^2)
  }

  p <- if (is.na(W)) NA_real_ else pchisq(W, df = df, lower.tail = FALSE)
  list(stat = W, df = df, p = p)
}

# Omnibus Wald test for equality of class effects in distal_pooled.
# Tests H0: beta[m, k] = beta[m, ref] for all k != ref and all m.
# df = (M - 1) * (K - 1).
.wald_omnibus_pooled <- function(beta_mat, Hessian, K, D_cov, ref_class) {
  M_minus_1 <- nrow(beta_mat)
  L         <- K + D_cov
  if (M_minus_1 == 0L || K <= 1L)
    return(list(stat = NA_real_, df = 0L, p = NA_real_))

  Sigma   <- pinv(-Hessian)
  non_ref <- setdiff(seq_len(K), ref_class)
  n_ctr   <- M_minus_1 * (K - 1L)
  R       <- matrix(0, nrow = n_ctr, ncol = M_minus_1 * L)

  row_i <- 1L
  for (m in seq_len(M_minus_1)) {
    for (k in non_ref) {
      R[row_i, (m - 1L) * L + k]         <-  1
      R[row_i, (m - 1L) * L + ref_class] <- -1
      row_i <- row_i + 1L
    }
  }

  # beta_mat is (M-1) x L; vectorise row-major to match Hessian block ordering
  beta_vec <- as.vector(t(beta_mat))
  r_vec    <- R %*% beta_vec
  V        <- R %*% Sigma %*% t(R)
  W        <- tryCatch(as.numeric(t(r_vec) %*% pinv(V) %*% r_vec),
                       error = function(e) NA_real_)
  p        <- if (is.na(W)) NA_real_ else pchisq(W, df = n_ctr, lower.tail = FALSE)
  list(stat = W, df = n_ctr, p = p)
}

# Predicted outcome probabilities for one class in a distal_pooled model,
# evaluated at covariates = 0 (i.e., the class intercept only).
.pred_probs_pooled <- function(beta_mat, K, D_cov, k) {
  L        <- K + D_cov
  U_k      <- matrix(0, nrow = 1, ncol = L)
  U_k[1, k] <- 1                                   # one-hot class indicator
  as.vector(distal_forward(U_k, beta_mat))
}

# Predicted outcome probabilities for one class in a distal_regression model,
# evaluated at covariates = 0 (i.e., using the intercept column only).
.pred_probs_reg <- function(betas_k, D) {
  U       <- matrix(0, nrow = 1, ncol = D)
  U[1, 1] <- 1                                     # intercept; covariates at 0
  as.vector(distal_forward(U, betas_k))
}

#' Summarise a Fitted Mixture Model
#'
#' @description
#' Prints a detailed summary of the structural model parameters. Depending on
#' which structural model was fitted, this includes covariate regression
#' coefficients (as odds ratios with 95\% confidence intervals and p-values),
#' distal outcome means, or class-specific regression effects. If no structural
#' model is present, a notice is printed directing the user to
#' \code{\link{measurement_summary}}.
#'
#' @param object A fitted \code{mixture_model} object returned by
#'   \code{\link{fit_mixture}}.
#' @param ref_class Integer. The reference latent class for pairwise contrasts.
#'   Defaults to the first class (\code{1}).
#' @param ... Currently unused. Present for S3 method compatibility.
#'
#' @return Invisibly returns \code{NULL}. Called for its printed side-effect.
#'
#' @examples
#' set.seed(1)
#' X <- matrix(rbinom(500, 1, 0.5), nrow = 100)
#' fit <- fit_mixture(X, n_components = 2, measurement = "binary")
#' summary(fit)
#'
#' @export
summary.mixture_model <- function(object, ref_class = NULL, ...) {
  K <- object$n_components
  if (is.null(ref_class)) ref_class <- 1

  # Input validation: ref_class must be a valid class index.
  # Without this guard the function starts printing output, then crashes
  # mid-way with a cryptic "subscript out of bounds" error.
  if (!is.numeric(ref_class) || length(ref_class) != 1 ||
      ref_class < 1 || ref_class > K)
    stop(sprintf(
      "ref_class must be an integer between 1 and %d. Got: %s",
      K, ref_class
    ))

  if (is.null(object$sm)) {
    cat("Notice: No structural model found. Use measurement_summary() for item parameters.\n")
    return(invisible())
  }

  cat("=========================================================\n")
  cat("             STRUCTURAL MODEL SUMMARY                    \n")
  cat("=========================================================\n")

  # A. Covariate model
  sm_sub <- NULL
  if (inherits(object$sm, "covariate")) sm_sub <- object$sm
  if (inherits(object$sm, "nested") && "predictor" %in% names(object$sm$models))
    sm_sub <- object$sm$models$predictor

  if (!is.null(sm_sub) && !is.null(sm_sub$parameters$hessian)) {
    cat("\nCATEGORICAL LATENT VARIABLE REGRESSION (CLASS PREDICTORS)\n")
    cat(sprintf("Reference Class: %d\n", ref_class))
    cat("---------------------------------------------------------\n")
    cat("                     OR       [95% CI]        P-Value\n")

    betas     <- sm_sub$parameters$beta
    D         <- ncol(betas)
    var_names <- if (!is.null(colnames(betas))) colnames(betas) else paste0("V", 1:D)
    Sigma     <- if (!is.null(sm_sub$parameters$V_robust))
      sm_sub$parameters$V_robust else pinv(-sm_sub$parameters$hessian)

    for (c in setdiff(1:K, ref_class)) {
      cat(sprintf("\nClass %d ON\n", c))
      for (v in 1:D) {
        est     <- betas[c, v] - betas[ref_class, v]
        idx_c   <- (c - 1) * D + v
        idx_ref <- (ref_class - 1) * D + v
        var_diff <- Sigma[idx_c, idx_c] + Sigma[idx_ref, idx_ref] -
          2 * Sigma[idx_c, idx_ref]
        se    <- sqrt(max(0, var_diff))
        z_val <- est / se
        p_val <- 2 * (1 - pnorm(abs(z_val)))
        cat(sprintf("  %-15s %7.3f  [%6.3f, %6.3f]  %s\n",
                    var_names[v], exp(est),
                    exp(est - 1.96 * se), exp(est + 1.96 * se),
                    .fmt_pval(p_val)))
      }
    }
  }

  # B0. Categorical distal outcome with no covariate (distal_categorical).
  #     This section is checked before B because distal_categorical inherits
  #     from distal_pooled; without the explicit class check below, section B
  #     would match first and display the wrong header.
  cat_sub <- NULL
  if (inherits(object$sm, "distal_categorical"))
    cat_sub <- object$sm
  if (inherits(object$sm, "nested") && "distal" %in% names(object$sm$models) &&
      inherits(object$sm$models$distal, "distal_categorical"))
    cat_sub <- object$sm$models$distal

  if (!is.null(cat_sub) && !is.null(cat_sub$parameters$beta_pooled)) {
    cat("\nCATEGORICAL DISTAL OUTCOME (CLASS PROBABILITIES)\n")
    cat("---------------------------------------------------------\n")

    beta_mat  <- cat_sub$parameters$beta_pooled
    M_minus_1 <- nrow(beta_mat)
    K_distal  <- K
    # distal_categorical has no covariate columns; D_cov is always 0.
    D_cov <- ncol(beta_mat) - K_distal
    M     <- M_minus_1 + 1L

    if (M_minus_1 == 0L) {
      cat("  (Constant outcome - no parameters to display)\n")
    } else {
      Sigma <- pinv(-cat_sub$parameters$hessian)

      # --- Omnibus test ---
      omni <- .wald_omnibus_pooled(beta_mat, cat_sub$parameters$hessian,
                                   K_distal, D_cov, ref_class)
      if (!is.na(omni$stat)) {
        cat(sprintf(
          "\nOmnibus test (class differences): Wald chi^2(%d) = %.2f, p%s\n",
          omni$df, omni$stat, .fmt_pval(omni$p)))
      }

      # --- Predicted probabilities ---
      cat("\nPredicted Probabilities:\n")
      cat(sprintf("  %-12s", ""))
      for (m in seq_len(M)) cat(sprintf(" Cat %-4d", m))
      cat("\n")
      for (k in seq_len(K_distal)) {
        probs <- .pred_probs_pooled(beta_mat, K_distal, D_cov = 0L, k)
        cat(sprintf("  Class %-6d", k))
        for (m in seq_len(M)) cat(sprintf("  %6.3f ", probs[m]))
        cat("\n")
      }

      # --- Pairwise OR table ---
      cat(sprintf("\nPairwise Odds Ratios (Reference: Class %d)\n", ref_class))
      cat("                     OR       [95% CI]        P-Value\n")
      for (m in seq_len(M_minus_1)) {
        cat(sprintf("\nOutcome Category %d (vs Category 1) ON\n", m + 1L))
        cat("  Latent Class:\n")
        for (c in setdiff(seq_len(K_distal), ref_class)) {
          est      <- beta_mat[m, c] - beta_mat[m, ref_class]
          idx_c    <- (m - 1L) * K_distal + c
          idx_ref  <- (m - 1L) * K_distal + ref_class
          var_diff <- Sigma[idx_c, idx_c] + Sigma[idx_ref, idx_ref] -
            2 * Sigma[idx_c, idx_ref]
          se    <- sqrt(max(0, var_diff))
          z_val <- if (se > 0) est / se else NA_real_
          p_val <- if (!is.na(z_val)) 2 * (1 - pnorm(abs(z_val))) else NA_real_
          cat(sprintf("    Class %d        %7.3f  [%6.3f, %6.3f]  %s\n",
                      c, exp(est), exp(est - 1.96 * se), exp(est + 1.96 * se),
                      .fmt_pval(p_val)))
        }
      }
    }
  }

  # B. Categorical distal outcome with a pooled covariate slope (distal_pooled).
  #    Excludes distal_categorical, which is handled in section B0 above.
  pooled_sub <- NULL
  if (inherits(object$sm, "distal_pooled") &&
      !inherits(object$sm, "distal_categorical")) pooled_sub <- object$sm
  if (inherits(object$sm, "nested") && "distal" %in% names(object$sm$models) &&
      inherits(object$sm$models$distal, "distal_pooled") &&
      !inherits(object$sm$models$distal, "distal_categorical"))
    pooled_sub <- object$sm$models$distal

  if (!is.null(pooled_sub) && !is.null(pooled_sub$parameters$beta_pooled)) {
    cat("\nCATEGORICAL DISTAL OUTCOME (POOLED SLOPES)\n")
    cat("---------------------------------------------------------\n")

    beta_mat  <- pooled_sub$parameters$beta_pooled
    M_minus_1 <- nrow(beta_mat)
    K_distal  <- K
    D_cov     <- ncol(beta_mat) - K_distal
    M         <- M_minus_1 + 1L
    Sigma     <- pinv(-pooled_sub$parameters$hessian)
    var_names <- if (D_cov > 0) paste0("Z", seq_len(D_cov)) else character(0)

    if (M_minus_1 > 0L) {

      # --- Omnibus test ---
      omni <- .wald_omnibus_pooled(beta_mat, pooled_sub$parameters$hessian,
                                   K_distal, D_cov, ref_class)
      if (!is.na(omni$stat)) {
        cat(sprintf(
          "\nOmnibus test (class differences): Wald chi^2(%d) = %.2f, p%s\n",
          omni$df, omni$stat, .fmt_pval(omni$p)))
      }

      # --- Predicted probabilities (primary display) ---
      cov_note <- if (D_cov > 0) " (covariates held at zero)" else ""
      cat(sprintf("\nPredicted Probabilities%s:\n", cov_note))
      cat(sprintf("  %-12s", ""))
      for (m in seq_len(M)) cat(sprintf(" Cat %-4d", m))
      cat("\n")
      for (k in seq_len(K_distal)) {
        probs <- .pred_probs_pooled(beta_mat, K_distal, D_cov, k)
        cat(sprintf("  Class %-6d", k))
        for (m in seq_len(M)) cat(sprintf("  %6.3f ", probs[m]))
        cat("\n")
      }

      # --- Pairwise OR table (secondary) ---
      cat(sprintf("\nPairwise Odds Ratios (Reference: Class %d)\n", ref_class))
      cat("                     OR       [95% CI]        P-Value\n")
      for (m in seq_len(M_minus_1)) {
        cat(sprintf("\nOutcome Category %d (vs Category 1) ON\n", m + 1L))
        cat("  Latent Class:\n")
        for (c in setdiff(seq_len(K_distal), ref_class)) {
          est      <- beta_mat[m, c] - beta_mat[m, ref_class]
          idx_c    <- (m - 1L) * (K_distal + D_cov) + c
          idx_ref  <- (m - 1L) * (K_distal + D_cov) + ref_class
          var_diff <- Sigma[idx_c, idx_c] + Sigma[idx_ref, idx_ref] -
            2 * Sigma[idx_c, idx_ref]
          se    <- sqrt(max(0, var_diff))
          z_val <- if (se > 0) est / se else NA_real_
          p_val <- if (!is.na(z_val)) 2 * (1 - pnorm(abs(z_val))) else NA_real_
          cat(sprintf("    Class %d        %7.3f  [%6.3f, %6.3f]  %s\n",
                      c, exp(est), exp(est - 1.96 * se), exp(est + 1.96 * se),
                      .fmt_pval(p_val)))
        }
        if (D_cov > 0) {
          cat("  Covariates (Pooled Slope):\n")
          for (v in seq_len(D_cov)) {
            est   <- beta_mat[m, K_distal + v]
            idx   <- (m - 1L) * (K_distal + D_cov) + K_distal + v
            se    <- sqrt(max(0, Sigma[idx, idx]))
            z_val <- if (se > 0) est / se else NA_real_
            p_val <- if (!is.na(z_val)) 2 * (1 - pnorm(abs(z_val))) else NA_real_
            cat(sprintf("    %-13s %7.3f  [%6.3f, %6.3f]  %s\n",
                        var_names[v], exp(est),
                        exp(est - 1.96 * se), exp(est + 1.96 * se),
                        .fmt_pval(p_val)))
          }
        }
      }
    }
  }

  # C. Distal regression (moderated)
  distal_sub <- NULL
  if (inherits(object$sm, "distal_regression")) distal_sub <- object$sm
  if (inherits(object$sm, "nested") && "distal" %in% names(object$sm$models) &&
      inherits(object$sm$models$distal, "distal_regression"))
    distal_sub <- object$sm$models$distal

  if (!is.null(distal_sub) && !is.null(distal_sub$parameters$betas)) {
    cat("\nCATEGORICAL DISTAL OUTCOME (CLASS-SPECIFIC SLOPES)\n")
    cat("---------------------------------------------------------\n")

    distal_betas <- distal_sub$parameters$betas
    K_distal     <- dim(distal_betas)[1]
    M_minus_1    <- dim(distal_betas)[2]
    D_distal     <- dim(distal_betas)[3]
    M            <- M_minus_1 + 1L

    if (M_minus_1 == 0L) {
      cat("  (Constant outcome - no parameters to display)\n")
    } else {
      var_names <- c("Intercept", paste0("Z", seq_len(D_distal - 1L)))
      cov_note  <- if (D_distal > 1L) " (covariates held at zero)" else ""

      # --- Predicted probabilities (primary display) ---
      cat(sprintf("\nPredicted Probabilities%s:\n", cov_note))
      cat(sprintf("  %-12s", ""))
      for (m in seq_len(M)) cat(sprintf(" Cat %-4d", m))
      cat("\n")
      for (k in seq_len(K_distal)) {
        betas_k <- matrix(distal_betas[k, , ], nrow = M_minus_1, ncol = D_distal)
        probs   <- .pred_probs_reg(betas_k, D_distal)
        cat(sprintf("  Class %-6d", k))
        for (m in seq_len(M)) cat(sprintf("  %6.3f ", probs[m]))
        cat("\n")
      }

      # --- Class-specific OR tables (secondary) ---
      cat("\nClass-Specific Estimates\n")
      cat("                     OR       [95% CI]        P-Value\n")
      for (k in seq_len(K_distal)) {
        cat(sprintf("\nClass %d:\n", k))
        Sigma <- if (!is.null(distal_sub$parameters$hessians) &&
                     length(distal_sub$parameters$hessians) >= k)
          pinv(-distal_sub$parameters$hessians[[k]])
        else
          matrix(0, M_minus_1 * D_distal, M_minus_1 * D_distal)

        for (m in seq_len(M_minus_1)) {
          cat(sprintf("  Outcome Category %d (vs Category 1) ON\n", m + 1L))
          for (v in seq_len(D_distal)) {
            est   <- distal_betas[k, m, v]
            idx   <- (m - 1L) * D_distal + v
            se    <- sqrt(max(0, Sigma[idx, idx]))
            z_val <- if (se > 0) est / se else NA_real_
            p_val <- if (!is.na(z_val)) 2 * (1 - pnorm(abs(z_val))) else NA_real_
            if (se > 0) {
              cat(sprintf("    %-13s %7.3f  [%6.3f, %6.3f]  %s\n",
                          var_names[v], exp(est),
                          exp(est - 1.96 * se), exp(est + 1.96 * se),
                          .fmt_pval(p_val)))
            } else {
              cat(sprintf("    %-13s %7.3f  [   N/A,    N/A]       N/A\n",
                          var_names[v], exp(est)))
            }
          }
        }
      }
    }
  }

  # D. Continuous distal (means)
  cont_sub <- NULL
  if (inherits(object$sm, "distal_continuous")) cont_sub <- object$sm
  if (inherits(object$sm, "nested") && "distal" %in% names(object$sm$models) &&
      inherits(object$sm$models$distal, "distal_continuous"))
    cont_sub <- object$sm$models$distal

  if (!is.null(cont_sub) && !is.null(cont_sub$parameters$means)) {
    cat("\nCONTINUOUS DISTAL OUTCOME (MEANS)\n")
    cat("---------------------------------------------------------\n")

    means    <- as.vector(cont_sub$parameters$means)
    ses      <- as.vector(cont_sub$parameters$ses)
    Sigma_mu <- cont_sub$parameters$Sigma_mu   # NULL for non-BCH models

    # Omnibus Wald test for equality of class means
    omni <- .wald_omnibus_means(means, ses, K, Sigma_mu = Sigma_mu)
    if (!is.na(omni$stat)) {
      cat(sprintf(
        "\nOmnibus test (class differences): Wald chi^2(%d) = %.2f, p%s\n",
        omni$df, omni$stat, .fmt_pval(omni$p)))
    }

    cat("\n                 Mean       [95% CI]        SE\n")
    for (k in seq_len(K)) {
      mu <- means[k]
      se <- ses[k]
      cat(sprintf("  Class %d      %7.3f  [%6.3f, %6.3f]   %7.3f\n",
                  k, mu, mu - 1.96 * se, mu + 1.96 * se, se))
    }
  }

  # E. Continuous distal regression
  cont_reg_sub <- NULL
  if (inherits(object$sm, "distal_continuous_regression")) cont_reg_sub <- object$sm
  if (inherits(object$sm, "nested") && "distal" %in% names(object$sm$models) &&
      inherits(object$sm$models$distal, "distal_continuous_regression"))
    cont_reg_sub <- object$sm$models$distal

  if (!is.null(cont_reg_sub) && !is.null(cont_reg_sub$parameters$betas)) {
    cat("\nCONTINUOUS DISTAL REGRESSION (Y ~ Z * Class)\n")
    cat("---------------------------------------------------------\n")

    betas     <- cont_reg_sub$parameters$betas
    ses       <- cont_reg_sub$parameters$ses
    D         <- ncol(betas)
    var_names <- if (!is.null(colnames(betas))) colnames(betas) else
      c("Intercept", paste0("Z", seq_len(D - 1L)))

    # Omnibus Wald test on class intercepts (class-specific means at Z = 0).
    # Uses the model-based SE for each intercept (sigma^2 * B_inv_k[1,1]),
    # assuming classes are independent - matching LatentGOLD's Wald(=) test.
    intercepts <- betas[, 1L]
    int_ses    <- ses[, 1L]   # already model-based if fitted with BCH v2
    prec_int   <- 1 / pmax(int_ses^2, 1e-15)
    mu_bar_int <- sum(prec_int * intercepts) / sum(prec_int)
    W_stat_int <- sum(prec_int * (intercepts - mu_bar_int)^2)
    omni <- list(stat = W_stat_int, df = K - 1L,
                 p = pchisq(W_stat_int, df = K - 1L, lower.tail = FALSE))
    if (!is.na(omni$stat)) {
      cov_note <- if (D > 1L) " (at covariate zero)" else ""
      cat(sprintf(
        "\nOmnibus test (class differences%s): Wald chi^2(%d) = %.2f, p%s\n",
        cov_note, omni$df, omni$stat, .fmt_pval(omni$p)))
    }

    cat("\n")
    for (k in seq_len(K)) {
      cat(sprintf("Class %d:\n", k))
      cat("                 Estimate   [95% CI]        P-Value\n")
      for (v in seq_len(D)) {
        est   <- betas[k, v]
        se    <- ses[k, v]
        z_val <- if (se > 0) est / se else NA_real_
        p_val <- if (!is.na(z_val)) 2 * (1 - pnorm(abs(z_val))) else NA_real_
        cat(sprintf("  %-13s %7.3f  [%6.3f, %6.3f]  %s\n",
                    var_names[v], est, est - 1.96 * se, est + 1.96 * se,
                    .fmt_pval(p_val)))
      }
      cat("\n")
    }

    #  Per-covariate Wald(=) tests: H0: slope_k equal across all classes
    # Uses the diagonal independence approximation (separate per-class
    # regressions), matching LatentGOLD's Wald(=) column.
    # Contrast matrix R = [-1 | I_{K-1}], df = K-1.
    if (K > 1L && D > 1L) {
      cat("---------------------------------------------------------\n")
      cat("Wald tests (equality of slopes across classes):\n")
      cat(sprintf("  %-13s   Wald(%s)%s  P-Value\n",
                  "", paste0("chi^2(", K - 1L, ")"), ""))
      R_eq <- cbind(-1, diag(K - 1L))
      for (v in 2L:D) {
        theta_v <- betas[, v]
        var_v   <- ses[, v]^2
        V_c     <- R_eq %*% diag(var_v) %*% t(R_eq)
        th_c    <- R_eq %*% theta_v
        W_v     <- tryCatch(
          as.numeric(t(th_c) %*% solve(V_c) %*% th_c),
          error = function(e) NA_real_)
        p_v     <- if (!is.na(W_v))
          pchisq(W_v, df = K - 1L, lower.tail = FALSE) else NA_real_
        cat(sprintf("  %-13s   %8.2f          %s\n",
                    var_names[v], W_v, .fmt_pval(p_v)))
      }
    }
  }

  # F. Continuous distal pooled
  cont_pool_sub <- NULL
  if (inherits(object$sm, "distal_continuous_pooled")) cont_pool_sub <- object$sm
  if (inherits(object$sm, "nested") && "distal" %in% names(object$sm$models) &&
      inherits(object$sm$models$distal, "distal_continuous_pooled"))
    cont_pool_sub <- object$sm$models$distal

  if (!is.null(cont_pool_sub) && !is.null(cont_pool_sub$parameters$beta_pooled)) {
    cat("\nCONTINUOUS DISTAL POOLED REGRESSION (Main Effects)\n")
    cat("---------------------------------------------------------\n")

    theta     <- as.vector(cont_pool_sub$parameters$beta_pooled)
    ses       <- as.vector(cont_pool_sub$parameters$ses)
    K_distal  <- K
    D_cov     <- length(theta) - K_distal
    # Use stored column names when available; fall back to Z1, Z2, ...
    stored_names <- colnames(cont_pool_sub$parameters$beta_pooled)
    var_names <- if (!is.null(stored_names) && length(stored_names) > K_distal)
      stored_names[(K_distal + 1L):length(stored_names)]
    else if (D_cov > 0) paste0("Z", seq_len(D_cov))
    else character(0)

    # Omnibus Wald test on class intercepts
    # When cov_theta is available (BCH step stored it), use the full
    # model-based contrast Wald: H0: int_k = int_1 for all k != 1.
    # This matches LatentGOLD's omnibus test and accounts for the
    # covariance between intercept estimates.
    cov_theta  <- cont_pool_sub$parameters$cov_theta
    intercepts <- theta[seq_len(K_distal)]
    int_ses    <- ses[seq_len(K_distal)]

    if (!is.null(cov_theta) && all(is.finite(cov_theta))) {
      cov_int   <- cov_theta[seq_len(K_distal), seq_len(K_distal)]
      R_int     <- cbind(-1, diag(K_distal - 1L))
      V_contr   <- R_int %*% cov_int %*% t(R_int)
      theta_c   <- R_int %*% intercepts
      W_stat    <- tryCatch(
        as.numeric(t(theta_c) %*% solve(V_contr) %*% theta_c),
        error = function(e) NA_real_
      )
      omni <- list(stat = W_stat, df = K_distal - 1L,
                   p = if (is.na(W_stat)) NA_real_
                   else pchisq(W_stat, df = K_distal - 1L, lower.tail = FALSE))
    } else {
      omni <- .wald_omnibus_means(intercepts, int_ses, K_distal)
    }

    if (!is.na(omni$stat)) {
      cov_note <- if (D_cov > 0) " (at covariate zero)" else ""
      cat(sprintf(
        "\nOmnibus test (class differences%s): Wald chi^2(%d) = %.2f, p%s\n",
        cov_note, omni$df, omni$stat, .fmt_pval(omni$p)))
    }

    cat("\n  Latent Class (Intercepts):\n")
    cat("                 Estimate   [95% CI]        P-Value\n")
    for (k in seq_len(K_distal)) {
      est   <- theta[k]
      se    <- ses[k]
      z_val <- if (se > 0) est / se else NA_real_
      p_val <- if (!is.na(z_val)) 2 * (1 - pnorm(abs(z_val))) else NA_real_
      cat(sprintf("    Class %d      %7.3f  [%6.3f, %6.3f]  %s\n",
                  k, est, est - 1.96 * se, est + 1.96 * se,
                  .fmt_pval(p_val)))
    }

    if (D_cov > 0) {
      cat("\n  Covariates (Pooled Slopes):\n")
      cat("                 Estimate   [95% CI]        P-Value\n")
      for (v in seq_len(D_cov)) {
        est   <- theta[K_distal + v]
        se    <- ses[K_distal + v]
        z_val <- if (se > 0) est / se else NA_real_
        p_val <- if (!is.na(z_val)) 2 * (1 - pnorm(abs(z_val))) else NA_real_
        cat(sprintf("    %-11s %7.3f  [%6.3f, %6.3f]  %s\n",
                    var_names[v], est, est - 1.96 * se, est + 1.96 * se,
                    .fmt_pval(p_val)))
      }
    }
  }

  cat("=========================================================\n")
}

#' Fit a Latent Mixture Model (LCA / LPA)
#'
#' @description
#' The core estimation function. Fits a latent class analysis (LCA) or latent
#' profile analysis (LPA) model using the EM algorithm. Optionally fits a
#' structural model (covariates or distal outcomes) using 1-, 2-, or 3-step
#' estimation with optional bias correction.
#'
#' @param X A numeric matrix or data frame of indicator variables for the
#'   measurement model. Rows are observations; columns are items or variables.
#' @param Y Optional numeric matrix or data frame of outcome or covariate
#'   variables for the structural model. Must be provided when
#'   \code{structural} is not \code{NULL}.
#' @param n_components Positive integer. Number of latent classes (or profiles)
#'   to estimate. Default is \code{2}.
#' @param measurement Character string or named list specifying the measurement
#'   model type. Accepted strings: \code{"binary"} / \code{"bernoulli"},
#'   \code{"categorical"} / \code{"multinoulli"},
#'   \code{"continuous"} / \code{"gaussian_diag"},
#'   \code{"gaussian"} / \code{"gaussian_unit"}.
#'   Missing values are handled automatically: any indicator column containing
#'   \code{NA} is estimated with a full-information (FIML) variant that masks the
#'   missing cells under a missing-at-random assumption, while complete columns
#'   use the faster complete-data estimator. A single specification (e.g.
#'   \code{"binary"}) therefore covers both complete and incomplete data, and the
#'   fitted object reports any missingness it found. The explicit \code{"*_nan"}
#'   forms (e.g. \code{"binary_nan"}, \code{"continuous_nan"}) remain accepted as
#'   aliases that force the missing-data variant.
#'   Pass a named list to specify a mixed (nested) measurement model with
#'   different variable types; each block's missing-data handling is resolved
#'   from the columns it governs. Default is \code{"binary"}.
#' @param structural Character string specifying the structural model type.
#'   One of \code{"covariate"}, \code{"distal_regression"},
#'   \code{"distal_pooled"}, \code{"distal_continuous"},
#'   \code{"distal_continuous_regression"}. Requires \code{Y}. Default is
#'   \code{NULL} (measurement model only).
#' @param n_steps Integer. Estimation approach: \code{1} for simultaneous
#'   1-step, \code{2} for 2-step, or \code{3} for bias-corrected 3-step.
#'   Default is \code{1}.
#' @param correction Character. Bias correction for 3-step estimation.
#'   One of \code{"none"}, \code{"BCH"}, or \code{"ML"}. Ignored when
#'   \code{n_steps} is not \code{3}. Default is \code{"none"}.
#' @param n_init Positive integer. Number of random restarts. The solution
#'   with the highest log-likelihood is retained. Default is \code{1}.
#' @param max_iter Positive integer. Maximum EM iterations per restart.
#'   Default is \code{1000}.
#' @param random_state Optional integer seed for reproducibility. Default is
#'   \code{NULL}.
#' @param order_by_size Logical. If \code{TRUE} (default), classes are sorted
#'   from largest to smallest after fitting.
#' @param weights Optional numeric vector of length \code{nrow(X)} for survey
#'   or case weights. Default is \code{NULL} (equal weights of 1).
#' @param strata Optional vector of stratum identifiers for complex survey designs.
#' @param cluster Optional vector of cluster identifiers for complex survey designs.
#' @param refine Logical. If \code{TRUE} (default), applies L-BFGS refinement
#'   after EM convergence to optimize the penalized maximum likelihood.
#' @param ... Additional arguments passed to the measurement or structural
#'   model constructors (e.g., \code{max_val} for multinoulli models).
#'
#' @return An object of class `mixture_model`, a list with:
#'   * `n_components` Number of latent classes.
#'   * `weights` Numeric vector of estimated class proportions.
#'   * `mm` Fitted measurement model state object.
#'   * `sm` Fitted structural model state object, or `NULL`.
#'   * `metrics` Named list: `ll` (log-likelihood), `aic`, `bic`, `sabic`,
#'     `n_params`, and `entropy` (relative entropy, 0-1 scale).
#'   * `log_resp` Matrix of log posterior class probabilities (n x K).
#'     Use `exp(fit$log_resp)` to obtain posterior probabilities.
#'   * `converged` Logical. Whether the EM algorithm converged.
#'   * `n_iter` Integer. Number of EM iterations run.
#'   * `step1_metrics` Named list of Step-1 fit indices (only when
#'     `n_steps = 3`).
#'
#' @examples
#' # Binary LCA with 3 classes and 5 random restarts
#' set.seed(1)
#' X <- matrix(rbinom(500, 1, 0.5), nrow = 100)
#' fit <- fit_mixture(X, n_components = 3, measurement = "binary", n_init = 5)
#' print(fit)
#' summary(fit)
#' measurement_summary(fit)
#'
#' # Continuous LPA (2 classes)
#' X_cont <- matrix(rnorm(300), nrow = 100)
#' fit_lpa <- fit_mixture(X_cont, n_components = 2, measurement = "continuous")
#'
#' # 3-step LCA with a covariate and ML correction
#' Z <- matrix(rnorm(100), nrow = 100)
#' fit_cov <- fit_mixture(X, Y = Z, n_components = 2, measurement = "binary",
#'                        structural = "covariate",
#'                        n_steps = 3, correction = "ML", n_init = 5)
#' summary(fit_cov)
#'
#' @export
#' @importFrom stats complete.cases cov dnorm optim pchisq plogis pnorm qlogis qnorm rbinom rnorm runif sd var
#' @importFrom utils setTxtProgressBar txtProgressBar
fit_mixture_internal <- function(X, Y = NULL, n_components = 2,
                                 measurement = "binary", structural = NULL,
                                 n_steps = 1, correction = "none", n_init = 1,
                                 max_iter = 1000, random_state = NULL,
                                 order_by_size = TRUE, weights = NULL,
                                 strata = NULL, cluster = NULL,
                                 refine = TRUE, ...) {

  if (is.data.frame(X)) X <- as.matrix(X)
  # Convert Y through prepare_covariates() so that:
  #   - numeric columns are passed through unchanged
  #   - factor / character columns are dummy-coded (first level = reference)
  #   - column names are always preserved for display in summary()
  if (!is.null(Y)) Y <- prepare_covariates(Y)

  n_samples <- nrow(X)

  # --- Input Validation ---

  # n_components must be a positive integer. Values of 0 or below would either
  # produce degenerate output silently or crash with uninformative messages.
  if (!is.numeric(n_components) || length(n_components) != 1 ||
      !is.finite(n_components) || n_components < 1L)
    stop(sprintf(
      "n_components must be a positive integer (>= 1). Got: %s",
      n_components
    ))

  # n_steps must be 1, 2, or 3
  if (!n_steps %in% c(1L, 2L, 3L))
    stop(sprintf("n_steps must be 1, 2, or 3. Got: %d", n_steps))

  # correction must be one of the three supported values.
  valid_corrections <- c("none", "ML", "BCH")
  if (!correction %in% valid_corrections)
    stop(sprintf(
      "correction '%s' not recognized. Choose from: %s",
      correction, paste(valid_corrections, collapse = ", ")
    ))

  # Missing values in the indicator matrix are handled automatically. The
  # measurement descriptor is resolved against the data so complete columns keep
  # the fast complete-data estimator while columns containing NA switch to the
  # FIML variant that masks missing cells (MAR assumption). The user-facing
  # specification (e.g. "binary") therefore covers both complete and incomplete
  # data; explicit "*_nan" strings remain accepted as aliases.
  measurement_engine <- .resolve_emission_descriptor(measurement, X)

  # Summarise missingness so it can be reported in the fitted object's print and
  # measurement summaries, and so the estimator used is explicit downstream.
  item_missing <- colSums(is.na(X))

  # Structural-side missingness (covariates / distal outcomes in Y). Covariates
  # are completed under the class-invariant Gaussian marginal
  # (endogenous-constrained-x; Sterba, 2014); a missing distal outcome is
  # handled by FIML inside the structural likelihood.
  y_any_missing <- !is.null(Y) && anyNA(Y)
  y_n_missing   <- if (!is.null(Y)) sum(is.na(Y)) else 0L
  struct_handled <- if (y_any_missing) {
    if (!is.null(structural) && structural %in%
        c("covariate", "predict_class"))
      "endogenous-constrained-x (Sterba, 2014)"
    else
      "endogenous-constrained-x covariates; FIML outcome"
  } else NA_character_

  missing_data <- list(
    any_missing      = anyNA(X) || y_any_missing,
    n_missing        = sum(is.na(X)),
    n_cells          = length(X),
    prop_missing     = if (length(X) > 0) mean(is.na(X)) else 0,
    per_item         = item_missing,
    n_items_affected = sum(item_missing > 0),
    handled_by       = if (anyNA(X)) "FIML (MAR assumption)" else NA_character_,
    y_any_missing    = y_any_missing,
    y_n_missing      = y_n_missing,
    structural_handled_by = struct_handled
  )

  # Validate binary data when a Bernoulli family is requested.
  if (is.character(measurement) &&
      measurement %in% c("binary", "bernoulli", "binary_nan", "bernoulli_nan")) {
    valid_vals <- X[!is.na(X)]
    if (length(valid_vals) > 0 && !all(valid_vals %in% c(0, 1)))
      stop(sprintf(
        paste0("measurement = '%s' requires X values in {0, 1}. ",
               "Found values outside this set: %s"),
        measurement,
        paste(sort(unique(valid_vals[!valid_vals %in% c(0, 1)]))[1:min(5, sum(!valid_vals %in% c(0,1)))],
              collapse = ", ")
      ))
  }

  if (is.null(weights)) {
    weights <- rep(1, n_samples)
  } else {
    if (length(weights) != n_samples)
      stop("Length of weights must match rows of X.")

    weights <- (weights / sum(weights, na.rm = TRUE)) * n_samples
  }

  # Survey design variables, when supplied, must align with the rows of X.
  # A design is considered present if either strata or cluster is given; the
  # other defaults so that every observation forms its own PSU or single
  # stratum, which leaves the linearization variance well defined.
  if (!is.null(strata) && length(strata) != n_samples)
    stop("Length of strata must match rows of X.")
  if (!is.null(cluster) && length(cluster) != n_samples)
    stop("Length of cluster must match rows of X.")
  has_survey_design <- !is.null(strata) || !is.null(cluster)

  # Structural model requires Y. Without this guard the SM is built but never
  # fitted (m_step_core gates on !is.null(Y)), so parameters$beta stays NULL
  # and every downstream function (coef, confint, wald tests) crashes with a
  # cryptic error rather than pointing here.
  if (!is.null(structural) && is.null(Y))
    stop(paste(
      "A structural model was specified but Y is NULL.",
      "Provide a Y matrix containing the outcome/covariate data,",
      "or set structural = NULL for a measurement-only model."
    ))

  model_state <- list(
    n_components          = n_components,
    weights               = rep(1 / n_components, n_components),
    mm                    = build_emission(measurement_engine, n_components = n_components, ...),
    sm                    = if (!is.null(structural))
      build_emission(structural, n_components = n_components, ...)
    else NULL,
    n_steps               = n_steps,
    correction            = correction,
    sample_weights        = weights,
    strata                = if (is.null(strata)) rep(1L, n_samples) else strata,
    cluster               = if (is.null(cluster)) seq_len(n_samples) else cluster,
    has_survey_design     = has_survey_design,
    # Retain the indicator matrix so plot() can scale continuous indicators
    # against their observed range (copy-on-write keeps this cheap).
    data                  = X,
    # Store the original descriptor so bootstrap.R can re-fit replicates
    # using the same measurement specification. Missing-data resolution is
    # re-applied per replicate, so the stored value is the user's spec, not the
    # resolved "*_nan" form.
    measurement_descriptor = measurement,
    # Record where and how missing data were handled (NA-free fits store a
    # summary with any_missing = FALSE).
    missing_data           = missing_data
  )
  class(model_state) <- "mixture_model"

  # Mirror the design onto the structural sub-model so that variance code
  # running inside m_step methods (which only receive the sub-model) can
  # reach the strata and cluster vectors. These are kept row-aligned with the
  # data the sub-model is fit on by any caller that subsets rows.
  if (!is.null(model_state$sm)) {
    model_state$sm$strata            <- model_state$strata
    model_state$sm$cluster           <- model_state$cluster
    model_state$sm$has_survey_design <- has_survey_design
  }

  if (n_steps == 1) {
    model_state <- fit_em(model_state, X, Y, n_init, max_iter, random_state,
                          refine = refine)

  } else if (n_steps == 2) {
    model_state <- fit_em(model_state, X, NULL, n_init, max_iter, random_state,
                          refine = refine)
    if (!is.null(Y) && !is.null(model_state$sm)) {
      resp <- exp(model_state$log_resp)
      model_state$sm <- init_params(model_state$sm, Y, resp)
      model_state$sm <- m_step(model_state$sm, Y, resp)
    }

  } else if (n_steps == 3) {
    model_state <- fit_em(model_state, X, NULL, n_init, max_iter, random_state,
                          refine = refine)

    # Step 1 metrics (measurement model only)
    n_params_s1  <- n_parameters(model_state$mm) + (model_state$n_components - 1)
    ll_s1        <- sum(model_state$sample_weights * model_state$lower_bound)
    resp_s1      <- exp(model_state$log_resp)
    abs_ent_s1   <- sum(model_state$sample_weights *
                          (-resp_s1 * log(resp_s1 + 1e-15)))
    # Use relative_entropy() to handle the K=1 edge case cleanly.
    rel_ent_s1   <- relative_entropy(abs_ent_s1,
                                     sum(model_state$sample_weights),
                                     model_state$n_components)

    model_state$step1_metrics <- list(
      ll       = ll_s1,
      n_params = n_params_s1,
      aic      = -2 * ll_s1 + 2 * n_params_s1,
      bic      = -2 * ll_s1 + log(n_samples) * n_params_s1,
      sabic    = -2 * ll_s1 + log((n_samples + 2) / 24) * n_params_s1,
      entropy  = rel_ent_s1
    )

    if (!is.null(Y) && !is.null(model_state$sm)) {
      if (correction == "ML") {
        model_state <- fit_ml(model_state, X, Y, max_iter = max_iter)
      } else if (correction == "BCH") {
        model_state <- fit_bch(model_state, X, Y)
      } else {
        # correction = "none": plain 2-step update on the structural model.
        # The measurement model is already frozen at this point; the SM is fit on the
        # posterior responsibilities from step 1 without any bias correction.
        resp <- exp(model_state$log_resp)
        model_state$sm <- init_params(model_state$sm, Y, resp)
        model_state$sm <- m_step(model_state$sm, Y, resp)
      }
    }
  }

  if (order_by_size) model_state <- sort_model_classes(model_state)

  # Attach column names.
  # Guard: only assign colnames to pis when dimensions match (Bernoulli: J cols;
  # Multinoulli: J*M cols - colnames(X) has length J so would misfire there).
  if (!is.null(colnames(X)) && !is.null(model_state$mm$parameters$pis) &&
      ncol(model_state$mm$parameters$pis) == ncol(X))
    colnames(model_state$mm$parameters$pis) <- colnames(X)

  # Attach column names to betas for distal_continuous_regression.
  if (!is.null(Y) && !is.null(model_state$sm) &&
      inherits(model_state$sm, "distal_continuous_regression") &&
      !is.null(model_state$sm$parameters$betas)) {
    K_br    <- model_state$n_components
    y_names <- if (!is.null(colnames(Y))) colnames(Y) else
      paste0("V", seq_len(ncol(Y)))
    cov_names_br <- if (ncol(Y) > 1L) y_names[-1L] else character(0L)
    br_names     <- c("Intercept", cov_names_br)
    if (length(br_names) == ncol(model_state$sm$parameters$betas))
      colnames(model_state$sm$parameters$betas) <- br_names
  }

  # Attach column names to beta_pooled for distal_continuous_pooled so that
  # summary() displays real variable names instead of Z1, Z2, etc.
  if (!is.null(Y) && !is.null(model_state$sm) &&
      inherits(model_state$sm, "distal_continuous_pooled") &&
      !is.null(model_state$sm$parameters$beta_pooled)) {
    K_bp    <- model_state$n_components
    y_names <- if (!is.null(colnames(Y))) colnames(Y) else
      paste0("V", seq_len(ncol(Y)))
    # First column of Y is the outcome; remaining are covariates
    cov_names_bp <- if (ncol(Y) > 1L) y_names[-1L] else character(0L)
    bp_names     <- c(paste0("Class_", seq_len(K_bp)), cov_names_bp)
    if (length(bp_names) == ncol(model_state$sm$parameters$beta_pooled))
      colnames(model_state$sm$parameters$beta_pooled) <- bp_names
  }

  # Attach covariate names to beta only when the SM has been initialised.
  # Guards against NULL beta on paths where the SM was not fit.
  if (!is.null(Y) && !is.null(model_state$sm) &&
      inherits(model_state$sm, "covariate") &&
      !is.null(model_state$sm$parameters$beta)) {
    intercept_flag <- isTRUE(model_state$sm$intercept)
    expected_D     <- ncol(model_state$sm$parameters$beta)

    y_col_names <- if (!is.null(colnames(Y))) colnames(Y)
    else paste0("V", seq_len(ncol(Y)))
    cov_names   <- if (intercept_flag) c("Intercept", y_col_names)
    else y_col_names
    if (!is.null(cov_names) && !is.null(expected_D) &&
        length(cov_names) == expected_D)
      colnames(model_state$sm$parameters$beta) <- cov_names
  }

  # Final metrics
  n_params <- n_parameters(model_state$mm) + (model_state$n_components - 1)
  if (!is.null(model_state$sm)) n_params <- n_params + n_parameters(model_state$sm)
  ll       <- sum(model_state$sample_weights * model_state$lower_bound)
  resp     <- exp(model_state$log_resp)
  abs_ent  <- sum(model_state$sample_weights * (-resp * log(resp + 1e-15)))
  # Use relative_entropy() to handle the K=1 edge case cleanly.
  ent      <- relative_entropy(abs_ent,
                               sum(model_state$sample_weights),
                               model_state$n_components)

  model_state$metrics <- list(
    ll       = ll,
    n_params = n_params,
    aic      = -2 * ll + 2 * n_params,
    bic      = -2 * ll + log(n_samples) * n_params,
    sabic    = -2 * ll + log((n_samples + 2) / 24) * n_params,
    entropy  = ent
  )

  return(model_state)
}

# ==============================================================================
# User-facing front-end for fit_mixture()
# ==============================================================================

# Measurement families whose complete-data descriptor has a missing-data (FIML)
# counterpart, mapping each base descriptor to the variant that masks NA during
# estimation. Descriptors absent from this table (e.g. structural families) have
# no missing-data variant and pass through resolution unchanged.
.nan_variant <- c(
  bernoulli     = "bernoulli_nan",
  binary        = "binary_nan",
  multinoulli   = "multinoulli_nan",
  categorical   = "categorical_nan",
  gaussian_diag = "gaussian_diag_nan",
  continuous    = "continuous_nan",
  gaussian_unit = "gaussian_unit_nan",
  gaussian      = "gaussian_nan"
)

# Resolve a measurement descriptor against its data so the estimator matches the
# data: complete columns keep the fast complete-data form, columns containing NA
# switch to the FIML variant that masks missing cells. Descriptors that already
# name a missing-data variant, or that have no variant, are returned unchanged.
#
# For a nested (mixed) measurement model each block is resolved against the
# columns it governs. Blocks are stored in order with consecutive column counts
# (.normalize_measurement() groups them this way), so a running offset maps each
# block to its columns.
.resolve_emission_descriptor <- function(descriptor, X) {
  nan_strings <- unname(.nan_variant)

  upgrade_one <- function(model, cols_have_na) {
    if (!cols_have_na) return(model)
    if (model %in% nan_strings) return(model)          # already a _nan variant
    variant <- .nan_variant[model]
    if (is.na(variant)) return(model)                  # no missing-data variant
    unname(variant)
  }

  if (is.character(descriptor) && length(descriptor) == 1L)
    return(upgrade_one(descriptor, anyNA(X)))

  if (is.list(descriptor)) {
    offset <- 0L
    for (name in names(descriptor)) {
      n_cols <- descriptor[[name]]$n_columns
      cols   <- seq.int(offset + 1L, offset + n_cols)
      descriptor[[name]]$model <-
        upgrade_one(descriptor[[name]]$model,
                    anyNA(X[, cols, drop = FALSE]))
      offset <- offset + n_cols
    }
    return(descriptor)
  }

  descriptor
}

# Translate a user measurement specification into the descriptor the fitting
# engine consumes, returning the (possibly re-grouped) indicator matrix.
#
# Accepts either a single type string (every indicator shares that type) or a
# named list/vector mapping a measurement type to the indicators it governs,
# for mixed-type models. Indicators may be referenced by column name or index:
#
#   measurement = "binary"
#   measurement = list(binary = c("q1", "q2"), continuous = c("score1"))
#   measurement = list(binary = 1:3, continuous = 4:5)
#
# Columns are grouped in the order given so the engine's block structure lines
# up; column names are preserved for display.
.normalize_measurement <- function(measurement, indicators) {
  indicators <- if (is.data.frame(indicators)) data.matrix(indicators)
  else as.matrix(indicators)

  if (is.character(measurement) && length(measurement) == 1L)
    return(list(descriptor = measurement, indicators = indicators))

  if (!is.list(measurement))
    stop("`measurement` must be a single type string or a named list mapping ",
         "measurement types to indicator columns.", call. = FALSE)

  parts <- as.list(measurement)
  if (is.null(names(parts)) || any(names(parts) == ""))
    stop("For a mixed measurement model, `measurement` must be a named list ",
         "whose names are measurement types (e.g. \"binary\", \"continuous\").",
         call. = FALSE)

  col_names <- colnames(indicators)
  resolve_cols <- function(sel) {
    if (is.character(sel)) {
      if (is.null(col_names))
        stop("Indicator columns were referenced by name, but `indicators` has ",
             "no column names.", call. = FALSE)
      idx <- match(sel, col_names)
      if (anyNA(idx))
        stop("Unknown indicator column name(s): ",
             paste(sel[is.na(idx)], collapse = ", "), call. = FALSE)
      idx
    } else {
      idx <- as.integer(sel)
      if (anyNA(idx) || any(idx < 1L) || any(idx > ncol(indicators)))
        stop("Indicator column indices in `measurement` are out of range.",
             call. = FALSE)
      idx
    }
  }

  keys        <- make.unique(names(parts), sep = "_")
  ordered_idx <- integer(0)
  descriptor  <- list()
  for (i in seq_along(parts)) {
    cols <- resolve_cols(parts[[i]])
    if (any(cols %in% ordered_idx))
      stop("An indicator column was assigned to more than one measurement type.",
           call. = FALSE)
    ordered_idx        <- c(ordered_idx, cols)
    descriptor[[keys[i]]] <- list(model = names(parts)[i], n_columns = length(cols))
  }

  unassigned <- setdiff(seq_len(ncol(indicators)), ordered_idx)
  if (length(unassigned) > 0)
    stop("Every indicator column must be assigned a measurement type. ",
         "Unassigned column(s): ", paste(unassigned, collapse = ", "),
         call. = FALSE)

  list(descriptor = descriptor,
       indicators  = indicators[, ordered_idx, drop = FALSE])
}

# Decide whether a distal outcome is continuous or categorical when the user
# leaves outcome_type = "auto". Factors, characters, and integer-valued numerics
# with few distinct values are treated as categorical.
.resolve_outcome_type <- function(outcome, outcome_type) {
  if (outcome_type != "auto") return(outcome_type)
  if (is.factor(outcome) || is.character(outcome)) return("categorical")
  v <- as.numeric(outcome)
  v <- v[!is.na(v)]
  if (length(v) == 0L) return("continuous")
  if (length(unique(v)) <= 10L && all(abs(v - round(v)) < 1e-8))
    return("categorical")
  "continuous"
}

# Pick a display label for the outcome column.
.outcome_label <- function(outcome) {
  if (!is.null(dim(outcome)) && !is.null(colnames(outcome)))
    return(colnames(outcome)[1])
  "outcome"
}

# Best-guess a variable name from the expression a user wrote, covering the
# common ways a single column is referenced: a bare symbol (age), extraction
# with `$` (data$age) or `[[` (data[["age"]]), and single-bracket indexing with
# a character column (data[, "age"]). Returns NULL when no name can be read off
# the expression (e.g. a positional index or a computed vector).
.derive_name <- function(expr) {
  if (is.symbol(expr)) return(as.character(expr))
  if (is.call(expr)) {
    op <- as.character(expr[[1]])
    if (op == "$")  return(as.character(expr[[3]]))
    if (op == "[[" && is.character(expr[[3]])) return(expr[[3]])
    if (op == "[") {
      args <- as.list(expr)
      for (i in seq_along(args)) {
        if (i <= 2L) next                 # skip `[` and the indexed object
        if (is.character(args[[i]])) return(args[[i]])
      }
    }
  }
  NULL
}

# Normalise a predictors / outcome_covariates argument so that prepare_covariates
# always receives something with proper column names. Matrices and data frames
# pass through unchanged (their names are kept); a single bare vector or factor
# is wrapped into a one-column data frame named from the user's expression, with
# `fallback` used only when no name can be inferred.
.as_named_covariates <- function(value, expr, fallback) {
  if (is.null(value) || !is.null(dim(value))) return(value)
  nm  <- .derive_name(expr)
  if (is.null(nm) || !nzchar(nm)) nm <- fallback
  out <- data.frame(value, check.names = FALSE, stringsAsFactors = FALSE)
  names(out) <- nm
  out
}

#' Fit a Latent Class or Latent Profile Mixture Model
#'
#' @description
#' Fits a finite mixture (latent class / latent profile) model. The latent
#' classes are defined by a set of measurement \code{indicators}; optionally, a
#' structural model relates the classes to external variables — either as
#' \code{predictors} of class membership, or as a distal \code{outcome} caused
#' by the classes.
#'
#' @param indicators Matrix or data frame of measurement items that define the
#'   latent classes (rows are observations, columns are items).
#' @param n_classes Number of latent classes/profiles to estimate.
#' @param measurement Either a single type string applied to every indicator
#'   (\code{"binary"}, \code{"categorical"}, \code{"continuous"},
#'   \code{"gaussian"}, and \code{"_nan"} missing-data variants), or, for a
#'   mixed-type model, a named list mapping each type to the indicator columns
#'   it governs by name or index, e.g.
#'   \code{list(binary = c("q1","q2"), continuous = "score")}.
#' @param predictors Optional covariates that predict latent class membership.
#'   Supplying this fits a class-membership regression (the "predict class"
#'   structural model). Mutually exclusive with \code{outcome}.
#' @param outcome Optional distal outcome caused by the latent classes.
#'   Mutually exclusive with \code{predictors}.
#' @param outcome_covariates Optional covariates that adjust the distal
#'   \code{outcome}.
#' @param outcome_type One of \code{"auto"}, \code{"continuous"}, or
#'   \code{"categorical"}. With \code{"auto"} (default) the type is inferred
#'   from \code{outcome}.
#' @param slopes When \code{outcome_covariates} are supplied, whether their
#'   effect is \code{"pooled"} (one slope shared across classes) or
#'   \code{"class_specific"} (a separate slope per class).
#' @param n_steps Estimation strategy: 1 (simultaneous), 2, or 3 (recommended
#'   when a structural model is present). Defaults to 3 when \code{predictors}
#'   or \code{outcome} is supplied and left unset, otherwise 1.
#' @param correction Bias correction for 3-step estimation: \code{"none"},
#'   \code{"ML"}, or \code{"BCH"}. When left unset for a 3-step structural
#'   model, a recommended default is chosen (ML for predictors and categorical
#'   outcomes, BCH for continuous outcomes).
#' @param weights,strata,cluster Optional survey design: sampling
#'   \code{weights}, and \code{strata}/\code{cluster} identifiers enabling
#'   design-based (linearization) standard errors.
#' @param n_init,max_iter,random_state,order_by_size,refine Estimation
#'   controls: number of random starts, maximum EM iterations, RNG seed,
#'   whether to order classes by size, and whether to run L-BFGS refinement.
#' @param X,Y,n_components,structural Deprecated legacy arguments retained for
#'   backward compatibility; prefer \code{indicators}, \code{n_classes},
#'   \code{predictors}, and \code{outcome}.
#' @param ... Passed through to the measurement-model constructors.
#'
#' @return A fitted \code{mixture_model} object.
#'
#' @examples
#' set.seed(1)
#' X <- matrix(rbinom(500, 1, 0.5), nrow = 100)
#' fit <- fit_mixture(X, n_classes = 2, measurement = "binary")
#'
#' \dontrun{
#' # Class membership predicted by a covariate (3-step, ML by default)
#' fit_mixture(X, n_classes = 2, predictors = age)
#'
#' # Distal outcome with a class-specific covariate slope
#' fit_mixture(X, n_classes = 2, outcome = bmi,
#'             outcome_covariates = age, slopes = "class_specific")
#'
#' # Mixed-type indicators
#' fit_mixture(items, n_classes = 3,
#'             measurement = list(binary = 1:5, continuous = 6:8))
#' }
#'
#' @export
fit_mixture <- function(indicators = NULL,
                        n_classes = 2,
                        measurement = "binary",
                        predictors = NULL,
                        outcome = NULL,
                        outcome_covariates = NULL,
                        outcome_type = c("auto", "continuous", "categorical"),
                        slopes = c("pooled", "class_specific"),
                        n_steps = 1,
                        correction = "none",
                        n_init = 1,
                        max_iter = 1000,
                        random_state = NULL,
                        order_by_size = TRUE,
                        weights = NULL,
                        strata = NULL,
                        cluster = NULL,
                        refine = TRUE,
                        X = NULL, Y = NULL, n_components = NULL, structural = NULL,
                        ...) {

  outcome_type <- match.arg(outcome_type)
  slopes       <- match.arg(slopes)
  steps_set    <- !missing(n_steps)
  corr_set     <- !missing(correction)

  # Capture the unevaluated expressions the user supplied so that a single
  # covariate passed as a bare vector (e.g. data$age or data[, "age"]) can be
  # given an informative name; see .as_named_covariates() below.
  predictors_expr <- substitute(predictors)
  outcome_cov_expr <- substitute(outcome_covariates)

  # --- Legacy interface bridge ------------------------------------------------
  legacy <- !is.null(X) || !is.null(Y) || !is.null(n_components) ||
    !is.null(structural)
  if (!is.null(X) && is.null(indicators)) indicators <- X
  if (!is.null(n_components))              n_classes  <- n_components

  if (legacy) {
    message("Note: `X`, `Y`, `n_components`, and `structural` are the legacy ",
            "interface. The current arguments are `indicators`, `n_classes`, ",
            "`predictors`, and `outcome` / `outcome_covariates`.")
    return(fit_mixture_internal(
      X = indicators, Y = Y, n_components = n_classes,
      measurement = measurement, structural = structural,
      n_steps = n_steps, correction = correction, n_init = n_init,
      max_iter = max_iter, random_state = random_state,
      order_by_size = order_by_size, weights = weights,
      strata = strata, cluster = cluster, refine = refine, ...))
  }

  if (is.null(indicators))
    stop("`indicators` is required: the matrix of measurement items that ",
         "define the latent classes.", call. = FALSE)

  if (!is.null(predictors) && !is.null(outcome))
    stop("Specify either `predictors` (to model class membership) or ",
         "`outcome` (a distal outcome), not both in one model.", call. = FALSE)
  if (!is.null(outcome_covariates) && is.null(outcome))
    stop("`outcome_covariates` requires an `outcome`.", call. = FALSE)

  # --- Measurement model (single-type or mixed) -------------------------------
  mm                 <- .normalize_measurement(measurement, indicators)
  X_use              <- mm$indicators
  measurement_engine <- mm$descriptor

  # --- Structural model -------------------------------------------------------
  structural_engine <- NULL
  Y_use             <- NULL

  if (!is.null(predictors)) {
    structural_engine <- "predict_class"
    Y_use             <- prepare_covariates(
      .as_named_covariates(predictors, predictors_expr, "predictor"))

  } else if (!is.null(outcome)) {
    otype <- .resolve_outcome_type(outcome, outcome_type)
    if (outcome_type == "auto")
      message(sprintf("Outcome treated as %s (set `outcome_type` to override).",
                      otype))

    has_cov <- !is.null(outcome_covariates)
    structural_engine <- if (otype == "continuous") {
      if (!has_cov)                  "continuous_outcome"
      else if (slopes == "pooled")   "continuous_outcome_adjusted"
      else                           "continuous_outcome_moderated"
    } else {
      if (!has_cov)                  "categorical_outcome"
      else if (slopes == "pooled")   "categorical_outcome_adjusted"
      else                           "categorical_outcome_moderated"
    }

    # Column 1 is always the outcome; covariates follow. The outcome is coerced
    # to a plain numeric column (categorical outcomes to 1-indexed integer
    # codes) so it is never dummy-coded as though it were a covariate.
    if (otype == "categorical") {
      out_col <- as.integer(as.factor(outcome))
    } else {
      out_col <- suppressWarnings(as.numeric(outcome))
      if (anyNA(out_col) && !anyNA(outcome))
        stop("A continuous `outcome` must be numeric.", call. = FALSE)
    }
    out_mat <- matrix(out_col, ncol = 1L,
                      dimnames = list(NULL, .outcome_label(outcome)))

    Y_use <- if (has_cov) cbind(out_mat, prepare_covariates(
      .as_named_covariates(outcome_covariates, outcome_cov_expr, "covariate")))
    else out_mat
  }

  # --- Friendly defaults when a structural model is present -------------------
  if (!is.null(structural_engine) && !steps_set) {
    n_steps <- 3L
    message("Using 3-step estimation (set `n_steps` to override).")
  }
  if (!is.null(structural_engine) && n_steps == 3L && !corr_set) {
    correction <- if (identical(structural_engine, "predict_class")) "ML"
    else if (startsWith(structural_engine, "categorical")) "ML"
    else "BCH"
    message(sprintf("Using '%s' bias correction (set `correction` to override).",
                    correction))
  }

  fit_mixture_internal(
    X = X_use, Y = Y_use, n_components = n_classes,
    measurement = measurement_engine, structural = structural_engine,
    n_steps = n_steps, correction = correction, n_init = n_init,
    max_iter = max_iter, random_state = random_state,
    order_by_size = order_by_size, weights = weights,
    strata = strata, cluster = cluster, refine = refine, ...)
}

#' Print a Brief Summary of a Fitted Mixture Model
#'
#' @description
#' Prints a compact overview of the fitted model including: number of classes,
#' estimation method, convergence status, log-likelihood, relative entropy,
#' and estimated class proportions. For full parameter tables, use
#' \code{\link{summary.mixture_model}} (structural parameters) or
#' \code{\link{measurement_summary}} (item parameters).
#'
#' @param x A fitted \code{mixture_model} object returned by
#'   \code{\link{fit_mixture}}.
#' @param ... Currently unused. Present for S3 method compatibility.
#'
#' @return Invisibly returns \code{x}. Called for its printed side-effect.
#'
#' @examples
#' set.seed(1)
#' X <- matrix(rbinom(300, 1, 0.5), nrow = 100)
#' fit <- fit_mixture(X, n_components = 2, measurement = "binary")
#' print(fit)
#' # or equivalently:
#' fit
#'
#' @export
print.mixture_model <- function(x, ...) {
  cat("=========================================================\n")
  cat("                  LATENT MIXTURE MODEL                   \n")
  cat("=========================================================\n")
  cat(sprintf("Classes Estimated  : %d\n", x$n_components))
  cat(sprintf("Estimation Method  : %d-step\n", x$n_steps))
  if (x$n_steps == 3) cat(sprintf("Correction Method  : %s\n", x$correction))
  cat(sprintf("Converged          : %s (in %d iterations)\n", x$converged, x$n_iter))
  if (!is.null(x$missing_data) && isTRUE(x$missing_data$any_missing)) {
    md <- x$missing_data
    cat(sprintf("Missing Data       : %d / %d cells (%.1f%%) in %d item%s \u2014 %s\n",
                md$n_missing, md$n_cells, 100 * md$prop_missing,
                md$n_items_affected, if (md$n_items_affected == 1L) "" else "s",
                md$handled_by))
  }
  cat("---------------------------------------------------------\n")
  if (!is.null(x$step1_metrics)) {
    cat(sprintf("  Log-Likelihood (Step 1) : %.2f\n", x$step1_metrics$ll))
    cat(sprintf("  Rel. Entropy   (Step 1) : %.4f\n", x$step1_metrics$entropy))
  } else {
    cat(sprintf("  Log-Likelihood : %.2f\n", x$metrics$ll))
    cat(sprintf("  Rel. Entropy   : %.4f\n", x$metrics$entropy))
  }
  cat("---------------------------------------------------------\n")
  cat("Class Weights (Sizes):\n")
  for (i in seq_along(x$weights))
    cat(sprintf("  Class %d: %.2f%%\n", i, x$weights[i] * 100))
  cat("=========================================================\n")
  cat("Type summary(model) for structural parameters or measurement_summary(model) for item parameters.\n")
}

#' Compare Mixture Models Across a Range of Class Numbers
#'
#' @description
#' Fits a sequence of measurement-only mixture models, one for each value of
#' \code{k} in \code{k_range}, and returns a table of fit indices to guide
#' class enumeration. The best model according to BIC is identified
#' automatically.
#'
#' @param X A numeric matrix or data frame of indicator variables.
#' @param k_range Integer vector of class numbers to fit. All values must be >= 1. Default is \code{1:5}.
#' @param measurement Character string specifying the measurement model type.
#'   See \code{\link{fit_mixture}} for accepted values. Default is
#'   \code{"binary"}.
#' @param n_init Positive integer. Number of random restarts per model.
#'   Default is \code{10}.
#' @param n_steps Integer. Estimation method: \code{1}, \code{2}, or \code{3}.
#'   Default is \code{1}.
#' @param ... Additional arguments passed to \code{\link{fit_mixture}}.
#'
#' @return A named list with three elements:
#'   * `fit_table` Data frame with one row per K and columns `Classes`, `LL`,
#'     `Params`, `AIC`, `BIC`, `SABIC`, `Entropy`.
#'   * `models` Named list of fitted `mixture_model` objects, one per K
#'     (names are `"K1"`, `"K2"`, etc.).
#'   * `best_k` Integer. The value of K with the lowest BIC.
#'
#' @examples
#' set.seed(1)
#' X <- matrix(rbinom(500, 1, 0.5), nrow = 100)
#' result <- compare_mixtures(X, k_range = 1:4, measurement = "binary",
#'                            n_init = 5)
#' result$fit_table
#' result$best_k
#'
#' @export
compare_mixtures <- function(X, k_range = 1:5, measurement = "binary",
                             n_init = 10, n_steps = 1, ...) {
  # k=0 would silently fit a degenerate model with LL=-Inf; negative
  # k values crash deep in initialisation with a cryptic error.
  if (any(k_range < 1L))
    stop(sprintf(
      "All values in k_range must be >= 1. Got invalid values: %s",
      paste(sort(unique(k_range[k_range < 1L])), collapse = ", ")
    ))
  cat(sprintf("Running Model Selection across K = %d to %d...\n\n",
              min(k_range), max(k_range)))

  # Resolve a single-type string or a mixed-type named list once, up front, so
  # every K is fit on the same (possibly re-grouped) indicators and descriptor.
  mm          <- .normalize_measurement(measurement, X)
  X           <- mm$indicators
  measurement <- mm$descriptor

  results <- list()
  models  <- list()
  for (k in k_range) {
    cat(sprintf("Fitting %d-class model...\n", k))
    fit <- fit_mixture_internal(X = X, Y = NULL, n_components = k,
                                measurement = measurement,
                                n_steps = n_steps, n_init = n_init, ...)
    models[[paste0("K", k)]] <- fit
    results[[k]] <- data.frame(
      Classes = k, LL = fit$metrics$ll, Params = fit$metrics$n_params,
      AIC = fit$metrics$aic, BIC = fit$metrics$bic,
      SABIC = fit$metrics$sabic, Entropy = fit$metrics$entropy
    )
  }
  fit_table   <- do.call(rbind, results)
  best_bic_k  <- fit_table$Classes[which.min(fit_table$BIC)]
  cat("\n=== Model Selection Summary ===\n")
  print(round(fit_table, 3))
  cat(sprintf("\n-> Best model according to BIC: %d classes\n", best_bic_k))
  return(list(fit_table = fit_table, models = models, best_k = best_bic_k))
}

#' Extract Covariate Odds Ratios from a Fitted Mixture Model
#'
#' @description
#' Extracts the logistic regression coefficients from a covariate structural
#' model and returns them as a matrix of odds ratios, centered on a reference
#' class. Only available when the model was fitted with
#' \code{structural = "covariate"}.
#'
#' @param object A fitted \code{mixture_model} object with a covariate
#'   structural model.
#' @param ref_class Integer. The reference class for centering. All other
#'   class odds ratios are expressed relative to this class. Default is
#'   \code{1}.
#' @param covariate_names Optional character vector of predictor names to
#'   override the column names stored in the model. Default is \code{NULL}.
#' @param ... Currently unused. Present for S3 method compatibility.
#'
#' @return A K x D numeric matrix of odds ratios, where rows are latent
#'   classes and columns are predictors (including the intercept). The
#'   reference class row will always show \code{1.000}.
#'
#' @examples
#' set.seed(1)
#' X <- matrix(rbinom(500, 1, 0.5), nrow = 100)
#' Z <- matrix(rnorm(100), nrow = 100)
#' colnames(Z) <- "age"
#' fit <- fit_mixture(X, Y = Z, n_components = 2, measurement = "binary",
#'                    structural = "covariate",
#'                    n_steps = 3, correction = "ML", n_init = 5)
#' coef(fit)
#'
#' @export
coef.mixture_model <- function(object, ref_class = 1, covariate_names = NULL, ...) {
  if (is.null(object$sm) || !inherits(object$sm, "covariate"))
    stop("No covariate model.")
  K     <- object$n_components
  betas <- object$sm$parameters$beta
  if (!is.null(covariate_names))
    colnames(betas) <- c("Intercept", covariate_names)
  betas_ref   <- sweep(betas, 2, betas[ref_class, ], "-")
  odds_ratios <- exp(betas_ref)
  rownames(odds_ratios) <- paste("Class", 1:K)
  rownames(odds_ratios)[ref_class] <- paste("Class", ref_class, "(Ref)")
  return(odds_ratios)
}
