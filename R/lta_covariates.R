# ==============================================================================
# Latent transition analysis - multiple groups and covariates
# ==============================================================================
#
# Collins & Lanza (2010, ch. 8) extend the latent transition model in two ways:
# a covariate may predict the latent status at Time 1 (sec. 8.10.1) or the
# transitions between statuses (sec. 8.10.2), and a grouping variable may shift
# the status prevalences and the transition probabilities (sec. 8.2-8.3). Their
# sec. 8.14 shows the two are equivalent when the measurement model is invariant
# across groups, so multiple-group LTA is implemented here as a set of dummy
# predictors, saturated over the transition rows to give each group its own free
# transition matrix. Group differences (sec. 8.6-8.8) are then tested by fitting
# with and without the group in the relevant regression.
#
# Every one of these regressions is a multinomial logit over the destination
# status, fitted by the .fit_mnl() that already serves the cross-sectional
# covariate model. The responsibilities weighting them come from the E-step:
# gamma at Time 1 for the initial-status regression, and the pairwise posteriors
# xi for the transition regressions.
#
# A saturated covariate model on transitions carries K(K-1) coefficients per
# covariate per occasion. Real transition tables have many near-empty cells, so
# it often fails to converge to anything interpretable. The default,
# transition_effects = "common", instead fits one slope per covariate shared
# across origin statuses (Wang & Wang, 2020, eq. 6.28).

# ------------------------------------------------------------------------------
# Design matrices
# ------------------------------------------------------------------------------

# Convert a user-supplied predictor argument into a numeric design matrix with
# an intercept, reusing prepare_covariates() so that factors are dummy-coded and
# column names survive for display.
.lta_design <- function(predictors, n, label) {
  if (is.null(predictors)) return(NULL)
  Z <- prepare_covariates(.as_named_covariates(predictors, NULL, label))
  Z <- complete_covariates(as.matrix(Z))
  if (nrow(Z) != n)
    stop(sprintf("`%s` must have one row per case.", label), call. = FALSE)
  cbind(Intercept = 1, Z)
}

# Dummy code a grouping variable, dropping the first level as the reference.
.lta_group_design <- function(group, n) {
  g <- factor(group)
  if (length(g) != n)
    stop("`group` must have one entry per case.", call. = FALSE)
  if (nlevels(g) < 2L)
    stop("`group` must have at least two levels.", call. = FALSE)
  Z <- stats::model.matrix(~ g)[, -1L, drop = FALSE]
  colnames(Z) <- paste0("group:", levels(g)[-1L])
  list(design = Z, levels = levels(g), factor = g)
}

# ------------------------------------------------------------------------------
# Parameters -> probabilities
# ------------------------------------------------------------------------------

# Log initial-status probabilities: a length-K vector, or an n x K matrix when a
# covariate predicts the initial status.
.lta_log_delta <- function(state) {
  if (is.null(state$delta_beta)) return(log(pmax(state$delta, 1e-300)))
  log(pmax(softmax_rows(state$Z_delta %*% t(state$delta_beta)), 1e-300))
}

# Coefficient matrix governing transitions out of origin status k for the m-th
# estimated transition matrix. Under "by_origin" each origin has its own fit;
# under "common" a single fit covers them all and the origin enters the design.
.lta_tau_beta <- function(state, m, k) {
  B <- state$tau_beta[[m]]
  if (identical(state$transition_effects, "by_origin")) B[[k]] else B
}

# Log transition probabilities for every pair of adjacent occasions, in the
# shape .lta_forward_backward() expects: a K x K matrix per occasion, or a list
# of K matrices of size n x K when a covariate predicts transitions.
.lta_log_tau <- function(state) {
  if (is.null(state$tau_beta))
    return(lapply(state$tau, function(m) log(pmax(m, 1e-300))))

  Tn <- state$n_times
  lapply(seq_len(Tn - 1L), function(t) {
    m <- if (isTRUE(state$tau_homogeneous)) 1L else t
    lapply(seq_len(state$n_statuses), function(k)
      log(pmax(softmax_rows(.lta_tau_design(state, k) %*%
                              t(.lta_tau_beta(state, m, k))), 1e-300)))
  })
}

# Design matrix used by the transition regression for origin status k.
#   "common"    : [intercept, origin dummies, covariates] - one slope per
#                 covariate shared across origin statuses (Wang & Wang eq. 6.28)
#   "by_origin" : [intercept, covariates] fitted separately per origin status
.lta_tau_design <- function(state, k) {
  Z <- state$Z_tau
  K <- state$n_statuses
  if (identical(state$transition_effects, "by_origin")) return(Z)
  dummies <- matrix(0, nrow(Z), K - 1L)
  if (k < K) dummies[, k] <- 1
  colnames(dummies) <- paste0("from:", seq_len(K - 1L))
  cbind(Z[, 1L, drop = FALSE], dummies, Z[, -1L, drop = FALSE])
}

# ------------------------------------------------------------------------------
# M-step for the regressions
# ------------------------------------------------------------------------------

# Initial status: an ordinary weighted multinomial logit on the Time 1
# posteriors, which already sum to one across statuses.
.lta_mstep_delta_cov <- function(state, gamma1) {
  fit <- .fit_mnl(state$Z_delta, gamma1, weights = state$weights_vec,
                  start = state$delta_beta)
  state$delta_beta    <- fit$beta
  state$delta_hessian <- fit$hessian
  state$delta         <- colSums(fit$prob * state$weights_vec) /
    sum(state$weights_vec)
  state
}

# Transitions. For every (case, occasion, origin) triple the pairwise posterior
# xi gives the expected destination distribution; its total mass is gamma, the
# probability of being in that origin status at all. .fit_mnl() expects rows
# summing to one, so the mass is normalised out of the responsibilities and
# folded into the case weight instead - algebraically identical, and it keeps
# the shared fitter free of special cases.
.lta_mstep_tau_cov <- function(state, fb) {
  K  <- state$n_statuses
  Tn <- state$n_times
  w  <- state$weights_vec
  n  <- length(w)

  mats <- if (isTRUE(state$tau_homogeneous)) 1L else Tn - 1L
  betas <- vector("list", mats)
  hess  <- vector("list", mats)

  for (m in seq_len(mats)) {
    ts <- if (isTRUE(state$tau_homogeneous)) seq_len(Tn - 1L) else m

    if (identical(state$transition_effects, "by_origin")) {
      bk <- vector("list", K); hk <- vector("list", K)
      for (k in seq_len(K)) {
        resp <- matrix(0, n * length(ts), K)
        wt   <- numeric(n * length(ts))
        Zs   <- matrix(0, n * length(ts), ncol(state$Z_tau))
        for (i_t in seq_along(ts)) {
          rows <- (i_t - 1L) * n + seq_len(n)
          pk   <- fb$pairwise[[ts[i_t]]][[k]]
          mass <- pmax(rowSums(pk), 1e-12)
          resp[rows, ] <- pk / mass
          wt[rows]     <- w * mass
          Zs[rows, ]   <- state$Z_tau
        }
        f <- .fit_mnl(Zs, resp, weights = wt,
                      start = state$tau_beta[[m]][[k]])
        bk[[k]] <- f$beta; hk[[k]] <- f$hessian
      }
      betas[[m]] <- bk; hess[[m]] <- hk

    } else {
      # The stacked design is the same at every EM iteration - only the
      # responsibilities and weights change - so it is built once and cached.
      key <- as.character(length(ts))
      if (is.null(state$.tau_design_cache[[key]])) {
        state$.tau_design_cache[[key]] <- do.call(
          rbind, rep(lapply(seq_len(K), function(k) .lta_tau_design(state, k)),
                     length(ts)))
      }
      Zs <- state$.tau_design_cache[[key]]

      n_rows <- n * K * length(ts)
      resp <- matrix(0, n_rows, K)
      wt   <- numeric(n_rows)
      r0 <- 0L
      for (i_t in seq_along(ts)) for (k in seq_len(K)) {
        rows <- r0 + seq_len(n); r0 <- r0 + n
        pk   <- fb$pairwise[[ts[i_t]]][[k]]
        mass <- pmax(rowSums(pk), 1e-12)
        resp[rows, ] <- pk / mass
        wt[rows]     <- w * mass
      }
      f <- .fit_mnl(Zs, resp, weights = wt, start = state$tau_beta[[m]])
      betas[[m]] <- f$beta; hess[[m]] <- f$hessian
    }
  }

  state$tau_beta    <- betas
  state$tau_hessian <- hess
  state$tau_n_params <- .lta_tau_cov_n_params(state)

  # Report the average transition matrix implied by the fitted regression, so
  # that transition_matrix() keeps returning something interpretable.
  lt <- .lta_log_tau(state)
  state$tau <- lapply(lt, function(per_k)
    t(vapply(per_k, function(M) colSums(exp(M) * w) / sum(w), numeric(K))))
  state
}

.lta_tau_cov_n_params <- function(state) {
  K <- state$n_statuses
  mats <- if (isTRUE(state$tau_homogeneous)) 1L else state$n_times - 1L
  D <- ncol(.lta_tau_design(state, 1L))
  if (identical(state$transition_effects, "by_origin"))
    mats * K * (K - 1L) * D
  else
    mats * (K - 1L) * D
}

# ------------------------------------------------------------------------------
# Reporting
# ------------------------------------------------------------------------------

#' Covariate Effects in a Latent Transition Model
#'
#' @description
#' Prints the multinomial-logit coefficients for the covariates predicting the
#' initial latent status and, where fitted, the transitions between statuses.
#' Coefficients are contrasts against the last latent status, which is the
#' reference category; exponentiating gives an odds ratio.
#'
#' For transitions fitted with `transition_effects = "common"` the design
#' contains the origin-status dummies (labelled `from:k`) that supply the
#' row-specific intercepts, followed by the covariate slopes, which are shared
#' across origin statuses. With `"by_origin"` a separate table is printed per
#' origin status.
#'
#' @param object A model fitted by [`fit_lta()`] with covariates or a group.
#' @param digits Number of digits to print.
#' @return `object`, invisibly.
#' @export
lta_covariate_summary <- function(object, digits = 3) {
  if (!inherits(object, "lta_model"))
    stop("`object` must be a fitted latent transition model.", call. = FALSE)
  K <- object$n_statuses

  if (is.null(object$delta_beta) && is.null(object$tau_beta)) {
    message("This model has no covariates or grouping variable.")
    return(invisible(object))
  }

  cat("\n=========================================================\n")
  cat("   LATENT TRANSITION MODEL - COVARIATE EFFECTS (logits)\n")
  cat("=========================================================\n")
  cat(sprintf("Reference status: %d. Coefficients are contrasts against it;\n",
              K))
  cat("exp(coefficient) is an odds ratio.\n")

  # Coefficient table with Wald tests from the multinomial-logit information.
  # The Hessian covers the free classes only, packed row-major by class.
  show_block <- function(beta, hessian, cov_names, row_lab) {
    free <- K - 1L
    D    <- ncol(beta)
    V    <- if (is.null(hessian)) NULL else tryCatch(pinv(-hessian),
                                                    error = function(e) NULL)
    se   <- if (is.null(V)) rep(NA_real_, free * D) else
      sqrt(pmax(diag(V), 0))
    est  <- as.vector(t(beta[seq_len(free), , drop = FALSE]))
    z    <- est / se
    out  <- data.frame(
      Status      = rep(row_lab, each = D),
      Term        = rep(cov_names, times = free),
      Estimate    = round(est, digits),
      SE          = round(se, digits),
      z           = round(z, 2),
      p           = format.pval(2 * stats::pnorm(-abs(z)), digits = 3, eps = 1e-16),
      OR          = round(exp(est), digits)
    )
    print(out, row.names = FALSE)
  }

  if (!is.null(object$delta_beta)) {
    cat("\nPREDICTING LATENT STATUS AT THE FIRST OCCASION\n")
    show_block(object$delta_beta, object$delta_hessian,
               colnames(object$Z_delta),
               paste0("Status ", seq_len(K - 1L)))
  }

  if (!is.null(object$tau_beta)) {
    cat("\nPREDICTING TRANSITIONS\n")
    by_origin <- identical(object$transition_effects, "by_origin")
    for (m in seq_along(object$tau_beta)) {
      tag <- if (length(object$tau_beta) == 1L) "all occasions" else
        sprintf("occasion %d -> %d", m, m + 1L)
      cat(sprintf("\n  [%s]\n", tag))
      if (by_origin) {
        for (k in seq_len(K)) {
          cat(sprintf("\n    from Status %d:\n", k))
          show_block(object$tau_beta[[m]][[k]], object$tau_hessian[[m]][[k]],
                     colnames(object$Z_tau),
                     paste0("to Status ", seq_len(K - 1L)))
        }
      } else {
        show_block(object$tau_beta[[m]], object$tau_hessian[[m]],
                   colnames(.lta_tau_design(object, 1L)),
                   paste0("to Status ", seq_len(K - 1L)))
      }
    }
  }
  cat("\n=========================================================\n")
  invisible(object)
}
