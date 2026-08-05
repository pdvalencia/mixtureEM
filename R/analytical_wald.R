# ==============================================================================
# Analytical Wald Omnibus Test
# ==============================================================================

# Columns of beta that belong to one covariate.
#
# The recorded term grouping (see prepare_covariates()) is authoritative when
# present: it is the only thing that knows the two dummies of a three-level
# factor are one variable. Substring matching remains as a fallback for models
# fitted before the grouping was recorded, or built by hand, but it is exactly
# what the grouping exists to replace — `grep("Age", ...)` also captures
# `Age_Decades` and would silently test two covariates as one — so it warns
# when it lands on columns from more than one apparent variable.
.wald_term_columns <- function(params, term_name) {
  col_names <- colnames(params$beta)
  terms     <- params$terms

  if (!is.null(terms) && length(terms) == ncol(params$beta)) {
    idx <- which(terms == term_name)
    if (length(idx)) return(idx)
    # A user may name a single dummy column ("Marital.Single") rather than the
    # variable; honour that, but keep it to the columns actually named.
    idx <- which(col_names == term_name)
    if (length(idx)) return(idx)
    stop(sprintf("Covariate '%s' not found. Available: %s",
                 term_name, paste(unique(setdiff(terms, "Intercept")),
                                  collapse = ", ")), call. = FALSE)
  }

  idx <- grep(term_name, col_names, fixed = TRUE)
  if (length(idx) == 0L) stop("Variable name not found.")
  if (length(idx) > 1L)
    warning(sprintf(paste0(
      "Matched %d columns (%s) by name for covariate '%s'; this model carries ",
      "no term grouping, so the match may span more than one variable. Refit ",
      "to record it."),
      length(idx), paste(col_names[idx], collapse = ", "), term_name),
      call. = FALSE)
  idx
}

#' Analytical Wald Test for a Single Covariate
#'
#' @description
#' Performs an omnibus Wald chi-squared test for the effect of a single
#' covariate on latent class membership. The null hypothesis is that the
#' covariate has no effect on class membership across all non-reference classes
#' simultaneously; the degrees of freedom are
#' (classes - 1) x (covariate columns), the conventional df for this test.
#'
#' The test uses whichever analytical variance the model was fitted with (see
#' \code{\link{covariate_se}}) and names it in the returned \code{Method}
#' column. For small samples or poorly conditioned Hessians, consider
#' \code{\link{wald_omnibus_test}} instead, which uses bootstrapped standard
#' errors from \code{\link{bootstrap_covariates}}.
#'
#' A caveat that applies to any Wald test of a strong effect: under the
#' Hauck-Donner phenomenon a coefficient large enough to nearly separate a
#' class drives the Wald statistic back towards zero. A non-significant omnibus
#' test standing beside large coefficients is therefore a warning, not a green
#' light; a likelihood-ratio or bootstrap test settles it.
#'
#' @param model A fitted \code{mixture_model} object with a covariate
#'   structural model (fitted with \code{structural = "covariate"}).
#' @param term_name Character string. The name of the covariate to test, as
#'   supplied to the model. For a factor this is the variable, not one of its
#'   dummy columns: naming a three-level \code{Marital} tests both of its
#'   contrasts jointly, which is the point of the omnibus test. Naming a single
#'   column (\code{"Marital.Single"}) tests that column alone. Models fitted
#'   before the term grouping was recorded fall back to substring matching
#'   against the column names, with a warning when that is ambiguous.
#' @param ref_class Integer. The reference latent class. Default is \code{1}.
#'
#' @return A single-row data frame with columns `Covariate`, `Wald_Chi2`,
#'   `df`, `p_value`, and `Method`, the name of the variance estimator behind
#'   the statistic.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' X <- matrix(rbinom(500, 1, 0.5), nrow = 100)
#' Z <- matrix(rnorm(100), nrow = 100)
#' colnames(Z) <- "age"
#' fit <- fit_mixture(X, Y = Z, n_components = 3, measurement = "binary",
#'                    structural = "covariate",
#'                    n_steps = 3, correction = "ML", n_init = 5)
#' analytical_wald_test(fit, term_name = "age")
#' }
#'
#' @export
analytical_wald_test <- function(model, term_name, ref_class = 1) {
  if (is.null(model$sm) || !inherits(model$sm, "covariate")) stop("No covariate model.")

  K <- model$n_components

  if (!is.numeric(ref_class) || length(ref_class) != 1 ||
      ref_class < 1 || ref_class > K)
    stop(sprintf(
      "ref_class must be an integer between 1 and %d. Got: %s",
      K, ref_class
    ))

  # K=1 means no contrasts can be formed (there is only one class).
  if (K == 1)
    stop("analytical_wald_test requires at least 2 classes. Got n_components = 1.")

  target_col_idx <- .wald_term_columns(model$sm$parameters, term_name)
  stat <- .wald_omnibus_core(model$sm$parameters, K, ref_class, target_col_idx)

  result <- data.frame(
    Covariate = term_name,
    Wald_Chi2 = round(stat$W, 3),
    df        = stat$df,
    p_value   = stat$p_value,
    Method    = stat$method
  )
  class(result) <- c("mixture_wald", "data.frame")
  return(result)
}

# Wald statistic for the joint null that a set of beta columns has no effect on
# class membership, i.e. that every contrast against the reference class is zero
# for every one of those columns. df is (K-1) x (columns),
# which for a k-level factor is (K-1) x (k-1).
#
# The variance is whichever one the model was fitted with, so the test inherits
# the step-3 correction automatically (see covariate_se()); it is named in the
# return value so a printed statistic can say which one produced it.
.wald_omnibus_core <- function(params, K, ref_class, cols) {
  V_robust <- params$V_robust
  if (!is.null(V_robust)) {
    Sigma_full  <- V_robust
    test_method <- params$V_method %||% "Survey-robust"
  } else {
    H <- params$hessian
    if (is.null(H) || all(H == 0)) stop("Hessian matrix is missing. Refit the model.")
    Sigma_full  <- pinv(-H)
    test_method <- params$V_method %||% "Q-function Hessian"
  }

  D            <- ncol(params$beta)
  total_params <- K * D
  test_classes <- setdiff(1:K, ref_class)
  num_tests    <- length(test_classes) * length(cols)

  C   <- matrix(0, nrow = num_tests, ncol = total_params)
  row <- 1
  for (c in test_classes) {
    for (v in cols) {
      C[row, (c - 1) * D + v]         <- 1
      C[row, (ref_class - 1) * D + v] <- -1
      row <- row + 1
    }
  }

  diff_beta <- C %*% as.vector(t(params$beta))
  V_cov     <- C %*% Sigma_full %*% t(C)

  # Safe matrix inversion
  W  <- as.numeric(t(diff_beta) %*% pinv(V_cov) %*% diff_beta)
  list(W = W, df = num_tests,
       p_value = pchisq(W, df = num_tests, lower.tail = FALSE),
       method = test_method)
}

#' @export
print.mixture_wald <- function(x, ...) {
  cat("=========================================================\n")
  cat("                 WALD TEST (COVARIATE)                   \n")
  cat("=========================================================\n")
  cat(sprintf("  Covariate : %s\n", x$Covariate))
  cat(sprintf("  Method    : %s\n", x$Method))
  cat("---------------------------------------------------------\n")
  cat(sprintf("  Wald \u03c7\u00b2(%d) = %.3f,  p%s\n",
              x$df, x$Wald_Chi2,
              if (x$p_value < 0.001) " < .001"
              else sprintf(" = %.3f", x$p_value)))
  cat("=========================================================\n")
  invisible(x)
}
