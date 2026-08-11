#' Confidence Intervals for Odds Ratios in a Mixture Model
#'
#' @description
#' Computes confidence intervals for the odds ratios of covariate effects on
#' latent class membership, using whichever analytical variance the model was
#' fitted with — see the \code{se} argument of \code{\link{fit_mixture}} and
#' \code{\link{covariate_se}}. The estimator used is named in the printed
#' output and carried in the \code{method} attribute of the result.
#' Bootstrapped standard errors from \code{\link{bootstrap_covariates}} can be
#' supplied instead, which is worth doing in small samples, when classes
#' overlap heavily, or as a check on the analytical intervals.
#'
#' @param object A fitted \code{mixture_model} object with a covariate
#'   structural model (fitted with \code{structural = "covariate"}).
#' @param parm A specification of which parameters are to be given confidence
#'   intervals. Currently ignored (intervals are returned for all covariate parameters).
#' @param level Numeric. The confidence level. Default is \code{0.95}.
#' @param boot_results Optional. A list returned by
#'   \code{\link{bootstrap_covariates}}. When provided, bootstrapped standard
#'   errors are used. When \code{NULL} (default), the analytical variance
#'   stored on the fitted model is used.
#' @param ref_class Integer. The reference latent class. Odds ratios for all
#'   other classes are expressed relative to this class. Default is \code{1}.
#' @param ... Currently unused. Present for S3 method compatibility.
#'
#' @return A named list with one data frame per predictor variable (including
#'   the intercept). Each data frame has columns `OR` (odds ratio), `Lower`,
#'   and `Upper` (confidence bounds), with one row per latent class. Values
#'   are rounded to three decimal places.
#'
#' @examples
#' set.seed(1)
#' X <- matrix(rbinom(500, 1, 0.5), nrow = 100)
#' Z <- matrix(rnorm(100), nrow = 100)
#' colnames(Z) <- "age"
#' fit <- fit_mixture(X, Y = Z, n_components = 2, measurement = "binary",
#'                    structural = "covariate",
#'                    n_steps = 3, correction = "ML", n_init = 5)
#'
#' # Analytical CIs (from Hessian)
#' confint(fit)
#'
#' \dontrun{
#' # Bootstrapped CIs
#' boot <- bootstrap_covariates(fit, X, Z, n_reps = 200)
#' confint(fit, boot_results = boot)
#' }
#'
#' @export
confint.mixture_model <- function(object, parm = NULL, level = 0.95,
                                  boot_results = NULL, ref_class = 1, ...) {
  if (is.null(object$sm) || !inherits(object$sm, "covariate"))
    stop("No covariate model.")

  K <- object$n_components

  if (!is.numeric(ref_class) || length(ref_class) != 1 ||
      ref_class < 1 || ref_class > K)
    stop(sprintf(
      "ref_class must be an integer between 1 and %d. Got: %s",
      K, ref_class
    ))

  D      <- ncol(object$sm$parameters$beta)
  z_crit <- qnorm(1 - (1 - level) / 2)

  # 1. Extract betas centered on ref_class
  betas     <- object$sm$parameters$beta
  betas_ref <- sweep(betas, 2, betas[ref_class, ], "-")

  # 2. Determine standard errors
  if (!is.null(boot_results)) {
    method <- "Bootstrap"

    # boot_betas is n_reps x K x D, centered on whatever ref_class was used
    # when bootstrap_covariates() was called.  Re-center on the requested
    # ref_class so SEs match the contrasts in betas_ref.
    boot_betas <- boot_results$boot_betas   # n_reps x K x D
    ref_mat    <- boot_betas[, ref_class, ] # n_reps x D

    boot_recentered <- boot_betas
    for (k in seq_len(K))
      boot_recentered[, k, ] <- boot_betas[, k, ] - ref_mat

    se <- apply(boot_recentered, c(2, 3), sd)

  } else {
    V_robust <- object$sm$parameters$V_robust
    if (!is.null(V_robust)) {
      method <- object$sm$parameters$V_method %||% "Survey-robust (linearization)"
      Sigma  <- V_robust
    } else {
      method <- object$sm$parameters$V_method %||% "Q-function Hessian"
      H <- object$sm$parameters$hessian
      if (is.null(H)) stop("Hessian missing. Refit model.")
      Sigma <- pinv(-H)
    }
    se      <- matrix(0, nrow = K, ncol = D)
    clamped <- 0L
    for (c in seq_len(K)) {
      for (v in seq_len(D)) {
        idx_c   <- (c - 1L) * D + v
        idx_ref <- (ref_class - 1L) * D + v
        var_diff <- Sigma[idx_c, idx_c] + Sigma[idx_ref, idx_ref] -
          2 * Sigma[idx_c, idx_ref]
        # The clamp keeps an NaN out of the printed table, but a negative
        # contrast variance means the covariance matrix is not positive
        # semi-definite, and the zero-width interval it produces looks like a
        # precise estimate rather than a failure. Say so.
        scale <- max(1, abs(Sigma[idx_c, idx_c]), abs(Sigma[idx_ref, idx_ref]))
        if (is.finite(var_diff) && var_diff < -1e-8 * scale) clamped <- clamped + 1L
        se[c, v] <- sqrt(max(0, var_diff))
      }
    }
    if (clamped > 0L)
      warning(sprintf(paste("%d of %d standard errors came from a negative",
                            "variance and were set to zero, giving intervals of",
                            "zero width. The covariance matrix is not positive",
                            "semi-definite; refit with se = \"robust\" for",
                            "usable intervals."),
                      clamped, K * D))
  }

  # 3. Compute bounds on the log-OR scale, then exponentiate
  lower_beta <- betas_ref - z_crit * se
  upper_beta <- betas_ref + z_crit * se

  OR <- exp(betas_ref)
  LB <- exp(lower_beta)
  UB <- exp(upper_beta)

  # 4. Build result list (one data frame per predictor)
  cov_names <- colnames(betas)
  if (is.null(cov_names)) cov_names <- paste0("V", seq_len(D))

  row_labels              <- paste("Class", seq_len(K))
  row_labels[ref_class]   <- paste("Class", ref_class, "(Ref)")

  res_list <- vector("list", D)
  names(res_list) <- cov_names

  for (v in seq_len(D)) {
    df <- data.frame(OR    = round(OR[, v], 3),
                     Lower = round(LB[, v], 3),
                     Upper = round(UB[, v], 3),
                     row.names = row_labels)
    res_list[[v]] <- df
  }

  # Attach metadata for the print method
  attr(res_list, "ref_class") <- ref_class
  attr(res_list, "level")     <- level
  attr(res_list, "method")    <- method
  class(res_list) <- c("mixture_confint", "list")

  return(res_list)
}

#' @export
print.mixture_confint <- function(x, ...) {
  ref_class <- attr(x, "ref_class")
  level     <- attr(x, "level")
  method    <- attr(x, "method")
  pct       <- paste0(round(level * 100), "%")

  cat("=========================================================\n")
  cat("        CONFIDENCE INTERVALS FOR ODDS RATIOS             \n")
  cat("=========================================================\n")
  cat(sprintf("Reference Class : %d\n", ref_class))
  cat(sprintf("Level           : %s   Method: %s\n", pct, method))
  cat("---------------------------------------------------------\n")
  cat(sprintf("  %-16s  %7s  %7s  %7s\n", "", "OR", "Lower", "Upper"))

  for (nm in names(x)) {
    cat(sprintf("\n%s\n", nm))
    df <- x[[nm]]
    for (i in seq_len(nrow(df))) {
      is_ref <- grepl("\\(Ref\\)", rownames(df)[i])
      if (is_ref) {
        cat(sprintf("  %-16s  %7.3f  %7s  %7s\n",
                    rownames(df)[i], df$OR[i], "-", "-"))
      } else {
        cat(sprintf("  %-16s  %7.3f  %7.3f  %7.3f\n",
                    rownames(df)[i], df$OR[i], df$Lower[i], df$Upper[i]))
      }
    }
  }
  cat("=========================================================\n")
  invisible(x)
}
