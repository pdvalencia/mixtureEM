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
  # Poisson LCA: class-specific rates are indexed by class like pis and means
  if (!is.null(model_state$mm$parameters[["rates"]])) {
    model_state$mm$parameters$rates <-
      model_state$mm$parameters$rates[new_order, , drop = FALSE]
  }
  # LCGA: growth coefficients are one row per class, and the gaussian family's
  # residual variance is one entry per class.
  if (!is.null(model_state$mm$parameters[["coefs"]])) {
    model_state$mm$parameters$coefs <-
      model_state$mm$parameters$coefs[new_order, , drop = FALSE]
    model_state$mm$parameters$dispersion <-
      model_state$mm$parameters$dispersion[new_order]
  }
  # GMM: growth-factor means are one row per class, residual variances one row
  # per class, and the growth-factor covariance one matrix per class. All are
  # stored per class even where a constraint makes the classes share a value, so
  # all are permuted unconditionally. The growth-factor regressions on
  # covariates are stored the same way and join them when there are any.
  if (!is.null(model_state$mm$parameters[["alpha"]])) {
    model_state$mm$parameters$alpha <-
      model_state$mm$parameters$alpha[new_order, , drop = FALSE]
    model_state$mm$parameters$theta <-
      model_state$mm$parameters$theta[new_order, , drop = FALSE]
    model_state$mm$parameters$psi <-
      model_state$mm$parameters$psi[new_order]
    if (!is.null(model_state$mm$parameters[["gamma"]]))
      model_state$mm$parameters$gamma <-
        model_state$mm$parameters$gamma[new_order]
  }

  # --- Sort nested measurement model sub-model parameters ---
  # The flat-parameter block above only touches model_state$mm$parameters, which
  # is empty for nested models.  Sub-model parameters live one level deeper at
  # model_state$mm$models[[name]]$parameters and must be sorted independently.
  if (inherits(model_state$mm, c("nested", "blocks"))) {
    for (name in names(model_state$mm$models)) {
      sub <- model_state$mm$models[[name]]
      if (!is.null(sub$parameters[["pis"]]))
        sub$parameters$pis <- sub$parameters$pis[new_order, , drop = FALSE]
      if (!is.null(sub$parameters[["means"]]))
        sub$parameters$means <- sub$parameters$means[new_order, , drop = FALSE]
      if (!is.null(sub$parameters[["covariances"]]))
        sub$parameters$covariances <-
          sub$parameters$covariances[new_order, , drop = FALSE]
      if (!is.null(sub$parameters[["rates"]]))
        sub$parameters$rates <- sub$parameters$rates[new_order, , drop = FALSE]
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
#' @param ... Passed to methods.
#'
#' @return Invisibly, a data frame in long format with one row per item,
#'   response category (polytomous items only, \code{NA} otherwise), and
#'   class: columns \code{block} (sub-model name for mixed measurement
#'   models, \code{NA} otherwise), \code{parameter} (\code{"probability"},
#'   \code{"mean"}, or \code{"rate"}), \code{item}, \code{category},
#'   \code{class}, and \code{estimate}. The same numbers are printed as
#'   formatted tables.
#'
#' @examples
#' set.seed(1)
#' X <- matrix(rbinom(500, 1, 0.5), nrow = 100)
#' fit <- fit_mixture(X, n_components = 2, measurement = "binary")
#' measurement_summary(fit)
#' params <- measurement_summary(fit)   # reuse the table programmatically
#'
#' @export
measurement_summary <- function(object, ...) UseMethod("measurement_summary")

#' @rdname measurement_summary
#' @export
measurement_summary.default <- function(object, ...) {
  K <- object$n_components
  cat("=========================================================\n")
  cat("             MEASUREMENT MODEL PARAMETERS                \n")
  cat("=========================================================\n")

  collected <- list()

  print_item_matrix <- function(mat, title, sub_model = NULL,
                                parameter = "estimate",
                                block = NA_character_) {
    cat(sprintf("\n%s\n", title))
    item_names <- colnames(mat)
    base_items <- NULL
    categories <- NULL

    if (is.null(item_names)) {
      if (!is.null(sub_model) && !is.null(sub_model$max_val)) {
        M       <- sub_model$max_val
        n_items <- ncol(mat) / M
        base    <- sub_model$item_names
        if (is.null(base) || length(base) != n_items)
          base <- paste0("Poly_Item_", seq_len(n_items))
        item_names <- paste0(rep(base, each = M),
                             " (Cat ", rep(seq_len(M), times = n_items), ")")
        base_items <- rep(base, each = M)
        categories <- rep(seq_len(M), times = n_items)
      } else {
        item_names <- paste0("Item_", 1:ncol(mat))
      }
    }
    if (is.null(base_items)) base_items <- item_names
    if (is.null(categories)) categories <- rep(NA_integer_, length(item_names))

    # Long indicator names widen the table only up to the shortening cap;
    # anything longer is abbreviated with a key printed under the table. The
    # returned data frame keeps the full names.
    disp    <- .shorten_labels(item_names, width = 30L)
    label_w <- max(20L, max(nchar(disp)))
    cat(sprintf("%-*s", label_w, "Indicator"))
    for (k in 1:K) cat(sprintf(" | Class %d", k))
    cat("\n")
    cat(paste0(rep("-", label_w + K * 10), collapse = ""), "\n")

    for (j in 1:ncol(mat)) {
      cat(sprintf("%-*s", label_w, disp[j]))
      for (k in 1:K) cat(sprintf(" | %7.3f", mat[k, j]))
      cat("\n")
    }
    .cat_label_legend(disp, indent = "")

    # Long-format rows for the returned data frame. as.vector(mat) walks the
    # K x J matrix column by column, so class varies fastest within each item.
    J <- ncol(mat)
    collected[[length(collected) + 1L]] <<- data.frame(
      block     = rep(block, K * J),
      parameter = rep(parameter, K * J),
      item      = rep(base_items, each = K),
      category  = rep(categories, each = K),
      class     = rep(seq_len(K), times = J),
      estimate  = as.vector(mat),
      stringsAsFactors = FALSE)
  }

  mm <- object$mm
  if (inherits(mm, c("nested", "blocks"))) {
    for (name in names(mm$models)) {
      sub_mm <- mm$models[[name]]
      if (!is.null(sub_mm$parameters$pis))
        print_item_matrix(sub_mm$parameters$pis,
                          paste("Categorical Probabilities:", toupper(name)),
                          sub_mm, "probability", name)
      if (!is.null(sub_mm$parameters$means))
        print_item_matrix(sub_mm$parameters$means,
                          paste("Continuous Means:", toupper(name)),
                          sub_mm, "mean", name)
      if (!is.null(sub_mm$parameters$rates))
        print_item_matrix(sub_mm$parameters$rates,
                          paste("Count Rates:", toupper(name)),
                          sub_mm, "rate", name)
    }
  } else {
    if (!is.null(mm$parameters$pis))
      print_item_matrix(mm$parameters$pis, "CATEGORICAL PROBABILITIES", mm,
                        "probability")
    if (!is.null(mm$parameters$means))
      print_item_matrix(mm$parameters$means, "CONTINUOUS MEANS", mm, "mean")
    if (!is.null(mm$parameters$rates))
      print_item_matrix(mm$parameters$rates, "COUNT RATES (lambda)", mm,
                        "rate")
  }
  if (!is.null(object$missing_data) && isTRUE(object$missing_data$any_missing)) {
    md <- object$missing_data
    cat(sprintf("\nMissing data: %d of %d cells (%.1f%%) across %d item%s, handled via %s.\n",
                md$n_missing, md$n_cells, 100 * md$prop_missing,
                md$n_items_affected, if (md$n_items_affected == 1L) "" else "s",
                md$handled_by))
  }
  cat("=========================================================\n")
  invisible(do.call(rbind, collected))
}

#' Print Classification Diagnostics
#'
#' @description
#' Prints the two tables that describe how cleanly the model assigns cases.
#'
#' The **Average Posterior Probability (AvePP)** matrix has one row per set of
#' observations modally assigned to a class and one column per class, holding
#' the mean posterior probability. High values on the diagonal, low values off
#' it, indicate well-separated classes.
#'
#' The **classification table** cross-classifies the probabilistic memberships
#' against the modal assignment, and yields the classification error: the
#' proportion of cases the modal rule is expected to place in the wrong class.
#' See [`classification_table()`] for the details, and [`absolute_fit()`] and
#' [`bivariate_residuals()`] for the fit of the model itself rather than the
#' quality of its assignments.
#'
#' Both tables use the case weights when the model was fitted with any.
#'
#' @param object A fitted \code{mixture_model} object returned by
#'   \code{\link{fit_mixture}}.
#' @param ... Passed to methods.
#'
#' @return Invisibly, a list with `ave_pp` (the K x K matrix), `table` (the
#'   classification table) and `error` (the classification error). All are also
#'   printed to the console.
#'
#' @examples
#' set.seed(1)
#' X <- matrix(rbinom(500, 1, 0.5), nrow = 100)
#' fit <- fit_mixture(X, n_components = 2, measurement = "binary")
#' classification_diagnostics(fit)
#'
#' @export
classification_diagnostics <- function(object, ...)
  UseMethod("classification_diagnostics")

#' @rdname classification_diagnostics
#' @export
classification_diagnostics.default <- function(object, ...) {
  resp <- exp(object$log_resp)
  K    <- object$n_components
  w    <- object$sample_weights %||% rep(1, nrow(resp))

  ave_pp <- .ave_pp(resp, w, K)
  rownames(ave_pp) <- paste("Assigned Class", 1:K)
  colnames(ave_pp) <- paste("Prob C", 1:K)

  cat("=========================================================\n")
  cat("          AVERAGE POSTERIOR PROBABILITIES (AvePP)        \n")
  cat("=========================================================\n")
  cat("Rows: Modal Assignment | Columns: Mean Probability\n\n")
  print(round(ave_pp, 3))
  cat("=========================================================\n")

  cat("\n")
  tab <- .classification_table(resp, w, K)
  print(tab)

  invisible(list(ave_pp = ave_pp, table = tab, error = attr(tab, "error")))
}

#' Class Sizes of a Fitted Mixture Model
#'
#' @description
#' Returns the estimated size of each latent class in the three forms applied
#' papers report: the model's class proportion, the expected number of cases,
#' and the number of cases modally assigned to the class. Case weights are
#' used when the model was fitted with any.
#'
#' @param object A fitted \code{mixture_model} object returned by
#'   \code{\link{fit_mixture}}.
#' @param ... Passed to methods.
#'
#' @return A data frame with one row per class: \code{class},
#'   \code{proportion} (model-estimated class weight), \code{n_expected}
#'   (proportion times the analysed sample size), and \code{n_modal} (cases
#'   assigned by highest posterior probability).
#'
#' @examples
#' set.seed(1)
#' X <- matrix(rbinom(500, 1, 0.5), nrow = 100)
#' fit <- fit_mixture(X, n_classes = 2)
#' class_sizes(fit)
#'
#' @export
class_sizes <- function(object, ...) UseMethod("class_sizes")

#' @rdname class_sizes
#' @export
class_sizes.mixture_model <- function(object, ...) {
  K     <- object$n_components
  resp  <- exp(object$log_resp)
  w     <- object$sample_weights %||% rep(1, nrow(resp))
  modal <- max.col(resp, ties.method = "first")
  n_tot <- sum(w)

  data.frame(
    class      = seq_len(K),
    proportion = as.vector(object$weights),
    n_expected = as.vector(object$weights) * n_tot,
    n_modal    = vapply(seq_len(K), function(k) sum(w[modal == k]), numeric(1))
  )
}

# Weighted average posterior probability by modal class. The weights matter
# whenever a case stands for more than one observation, which is the norm for
# the frequency-weighted pattern files these models are often fitted to.
.ave_pp <- function(resp, w, K) {
  modal  <- max.col(resp, ties.method = "first")
  ave_pp <- matrix(NA_real_, nrow = K, ncol = K)
  for (k in seq_len(K)) {
    idx <- modal == k
    if (any(idx))
      ave_pp[k, ] <- colSums(resp[idx, , drop = FALSE] * w[idx]) / sum(w[idx])
  }
  ave_pp
}

# ==============================================================================
# Internal helpers for summary.mixture_model
# ==============================================================================

# Bind collected table rows into one clean data frame (NULL when empty).
.rbind_tidy <- function(rows) {
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (!length(rows)) return(NULL)
  df <- do.call(rbind, rows)
  rownames(df) <- NULL
  df
}

# Format p-values to publication conventions.
# Values below .001 are shown as "< .001"; NaN/NA appear as a dash.
.fmt_pval <- function(p) {
  if (is.na(p) || is.nan(p)) return("       -")
  if (p < 0.001)              return("  < .001")
  sprintf("   %5.3f", p)
}

# Per-covariate omnibus Wald tests, printed under the class-predictor table.
#
# The coefficient table above answers "does this covariate distinguish class c
# from the reference class?" once per contrast. Two things are missing from it,
# and this table supplies both:
#
#   * Multiplicity. A 3-class model with five covariates prints twelve
#     coefficient tests, and some of them will be significant.
#   * For a covariate with more than two categories the table contains no test
#     of the covariate at all — a three-level marital status becomes two dummies
#     times two class contrasts, and none of those four rows answers "does
#     marital status predict class membership?", which is the question that was
#     asked. The omnibus test is the one that does.
#
# Nothing is printed when every covariate is a single column and there are only
# two classes, since each omnibus test is then the square of a z already shown.
.print_covariate_omnibus <- function(object, sm_sub, ref_class) {
  params <- sm_sub$parameters
  if (is.null(params$beta) || ncol(params$beta) == 0L) return(invisible(NULL))
  K <- object$n_components
  if (K < 2L) return(invisible(NULL))

  terms <- params$terms
  if (is.null(terms) || length(terms) != ncol(params$beta))
    terms <- colnames(params$beta) %||% paste0("V", seq_len(ncol(params$beta)))

  cov_terms <- unique(setdiff(terms, "Intercept"))
  if (!length(cov_terms)) return(invisible(NULL))
  if (K == 2L && length(cov_terms) == length(terms) - sum(terms == "Intercept"))
    return(invisible(NULL))

  rows <- lapply(cov_terms, function(tm) {
    cols <- which(terms == tm)
    st   <- tryCatch(.wald_omnibus_core(params, K, ref_class, cols),
                     error = function(e) NULL)
    if (is.null(st)) return(NULL)
    data.frame(term = tm, chi2 = st$W, df = st$df, p = st$p_value,
               stringsAsFactors = FALSE)
  })
  rows <- do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
  if (is.null(rows) || !nrow(rows)) return(invisible(NULL))

  cat("\nOMNIBUS TEST PER COVARIATE (effect across all classes)\n")
  cat("---------------------------------------------------------\n")
  # Same shortening policy as the coefficient table above, so a long variable
  # name reads identically in both places.
  disp    <- .shorten_labels(rows$term)
  label_w <- .label_width(disp, min = 20L)
  cat(sprintf("  %-*s %11s %4s  %s\n", label_w, "", "Wald Chi2", "df", "P-Value"))
  for (i in seq_len(nrow(rows)))
    cat(sprintf("  %-*s %11.3f %4d  %s\n",
                label_w, disp[i], rows$chi2[i], rows$df[i], .fmt_pval(rows$p[i])))
  .cat_label_legend(disp)
  # The Wald statistic is non-monotone in the effect size: a covariate strong
  # enough to nearly separate a class drives it back towards zero
  # (Hauck & Donner, 1977). The test gates in one direction only.
  cat("  Note: a non-significant test beside large coefficients can be the\n")
  cat("        Hauck-Donner effect; confirm with wald_omnibus_test().\n")
  invisible(rows)
}

# Omnibus Wald test for equality of K means (continuous distal outcomes).
#
# When Sigma_mu (the full K x K sandwich variance-covariance matrix of the
# means) is available, uses the full-covariance formulation:
#   W = c^T V^{-1} c,  where c = R * mu,  V = R * Sigma_mu * R^T
#   R = contrast matrix [class k vs class 1, k = 2..K]  (df = K-1)
#
# This is the robust Wald formulation of Bakk, Oberski and Vermunt
# (2014), accounting for the cross-class covariance induced by BCH weights.
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
#' @return Invisibly, a list holding the printed numbers in tidy form, ready
#'   for further use: \code{$coefficients} (one row per class contrast and
#'   covariate: estimate, SE, p, odds ratio and its confidence limits),
#'   \code{$omnibus} (the per-covariate omnibus Wald tests), and, when a
#'   distal outcome is present, \code{$outcome} (predicted probabilities or
#'   class means/estimates with their tests). Returns \code{NULL} when the
#'   model has no structural part.
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

  # Repeated ahead of everything else: a structural result read off a fit whose
  # measurement model has a collapsed class is not worth interpreting, so the
  # note has to come before the coefficients rather than after them.
  .print_degenerate_note(object)

  if (is.null(object$sm)) {
    cat("Notice: No structural model found. Use measurement_summary() for item parameters.\n")
    return(invisible())
  }

  # Everything printed below is also collected here and returned invisibly,
  # so vignettes and downstream code can use the numbers without re-deriving
  # them from the model internals.
  out <- list(ref_class  = ref_class,
              n_classes  = K,
              n_steps    = object$n_steps,
              correction = object$correction)

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
    # A covariate standard error can come from any of several estimators that
    # differ by up to a factor of two, so the printed table names the one used
    # rather than leaving the reader to guess (Bakk, Oberski & Vermunt, 2014).
    cat(sprintf("Standard errors: %s\n",
                sm_sub$parameters$V_method %||%
                  if (!is.null(sm_sub$parameters$V_robust))
                    "Survey-robust (linearization)" else "Q-function Hessian"))
    cat("---------------------------------------------------------\n")

    betas     <- sm_sub$parameters$beta
    D         <- ncol(betas)
    var_names <- if (!is.null(colnames(betas))) colnames(betas) else paste0("V", 1:D)
    # The label column sizes itself to the longest (shortened) predictor name:
    # a dummy-coded predictor name runs past fifteen characters routinely, and
    # names past the shortening cap are abbreviated with a key printed under
    # the table. The data frame returned by summary() keeps the full names.
    disp      <- .shorten_labels(var_names)
    label_w   <- .label_width(disp, min = 20L)
    Sigma     <- if (!is.null(sm_sub$parameters$V_robust))
      sm_sub$parameters$V_robust else pinv(-sm_sub$parameters$hessian)

    cat(sprintf("  %-*s %9s  %s  %s\n", label_w,
                "", "OR", "       [95% CI]       ", "P-Value"))

    cov_rows <- vector("list", (K - 1L) * D)
    row_i    <- 0L
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
        cat(sprintf("  %-*s %9.3f  [%9.3f, %9.3f]  %s\n",
                    label_w, disp[v], exp(est),
                    exp(est - 1.96 * se), exp(est + 1.96 * se),
                    .fmt_pval(p_val)))
        row_i <- row_i + 1L
        cov_rows[[row_i]] <- data.frame(
          class = c, term = var_names[v], estimate = est, se = se,
          z = z_val, p = p_val, OR = exp(est),
          OR_lower = exp(est - 1.96 * se), OR_upper = exp(est + 1.96 * se),
          stringsAsFactors = FALSE)
      }
    }
    out$coefficients <- .rbind_tidy(cov_rows[seq_len(row_i)])
    .cat_label_legend(disp)

    out$omnibus <- .print_covariate_omnibus(object, sm_sub, ref_class)
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
      pp <- matrix(NA_real_, K_distal, M,
                   dimnames = list(paste("Class", seq_len(K_distal)),
                                   paste("Cat", seq_len(M))))
      cat("\nPredicted Probabilities:\n")
      cat(sprintf("  %-12s", ""))
      for (m in seq_len(M)) cat(sprintf(" Cat %-4d", m))
      cat("\n")
      for (k in seq_len(K_distal)) {
        probs <- .pred_probs_pooled(beta_mat, K_distal, D_cov = 0L, k)
        pp[k, ] <- probs
        cat(sprintf("  Class %-6d", k))
        for (m in seq_len(M)) cat(sprintf("  %6.3f ", probs[m]))
        cat("\n")
      }

      # --- Pairwise OR table ---
      or_rows <- list()
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
          or_rows[[length(or_rows) + 1L]] <- data.frame(
            category = m + 1L, class = c, estimate = est, se = se, p = p_val,
            OR = exp(est), OR_lower = exp(est - 1.96 * se),
            OR_upper = exp(est + 1.96 * se), stringsAsFactors = FALSE)
        }
      }
      out$outcome <- list(
        type    = "categorical",
        omnibus = data.frame(chi2 = omni$stat, df = omni$df, p = omni$p),
        predicted_probabilities = pp,
        odds_ratios = .rbind_tidy(or_rows))
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
      pp <- matrix(NA_real_, K_distal, M,
                   dimnames = list(paste("Class", seq_len(K_distal)),
                                   paste("Cat", seq_len(M))))
      cov_note <- if (D_cov > 0) " (covariates held at zero)" else ""
      cat(sprintf("\nPredicted Probabilities%s:\n", cov_note))
      cat(sprintf("  %-12s", ""))
      for (m in seq_len(M)) cat(sprintf(" Cat %-4d", m))
      cat("\n")
      for (k in seq_len(K_distal)) {
        probs <- .pred_probs_pooled(beta_mat, K_distal, D_cov, k)
        pp[k, ] <- probs
        cat(sprintf("  Class %-6d", k))
        for (m in seq_len(M)) cat(sprintf("  %6.3f ", probs[m]))
        cat("\n")
      }

      # --- Pairwise OR table (secondary) ---
      or_rows  <- list()
      cov_rows <- list()
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
          or_rows[[length(or_rows) + 1L]] <- data.frame(
            category = m + 1L, class = c, estimate = est, se = se, p = p_val,
            OR = exp(est), OR_lower = exp(est - 1.96 * se),
            OR_upper = exp(est + 1.96 * se), stringsAsFactors = FALSE)
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
            cov_rows[[length(cov_rows) + 1L]] <- data.frame(
              category = m + 1L, term = var_names[v], estimate = est, se = se,
              p = p_val, OR = exp(est), OR_lower = exp(est - 1.96 * se),
              OR_upper = exp(est + 1.96 * se), stringsAsFactors = FALSE)
          }
        }
      }
      out$outcome <- list(
        type    = "categorical",
        omnibus = data.frame(chi2 = omni$stat, df = omni$df, p = omni$p),
        predicted_probabilities = pp,
        odds_ratios = .rbind_tidy(or_rows),
        covariate_effects = .rbind_tidy(cov_rows))
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
      pp <- matrix(NA_real_, K_distal, M,
                   dimnames = list(paste("Class", seq_len(K_distal)),
                                   paste("Cat", seq_len(M))))
      cat(sprintf("\nPredicted Probabilities%s:\n", cov_note))
      cat(sprintf("  %-12s", ""))
      for (m in seq_len(M)) cat(sprintf(" Cat %-4d", m))
      cat("\n")
      for (k in seq_len(K_distal)) {
        betas_k <- matrix(distal_betas[k, , ], nrow = M_minus_1, ncol = D_distal)
        probs   <- .pred_probs_reg(betas_k, D_distal)
        pp[k, ] <- probs
        cat(sprintf("  Class %-6d", k))
        for (m in seq_len(M)) cat(sprintf("  %6.3f ", probs[m]))
        cat("\n")
      }

      # --- Class-specific OR tables (secondary) ---
      est_rows <- list()
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
            est_rows[[length(est_rows) + 1L]] <- data.frame(
              class = k, category = m + 1L, term = var_names[v],
              estimate = est, se = if (se > 0) se else NA_real_, p = p_val,
              OR = exp(est), stringsAsFactors = FALSE)
          }
        }
      }
      out$outcome <- list(
        type = "categorical",
        predicted_probabilities = pp,
        estimates = .rbind_tidy(est_rows))
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
    out$outcome <- list(
      type    = "continuous",
      omnibus = data.frame(chi2 = omni$stat, df = omni$df, p = omni$p),
      means   = data.frame(class = seq_len(K), mean = means, se = ses,
                           lower = means - 1.96 * ses,
                           upper = means + 1.96 * ses))
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
    disp      <- .shorten_labels(var_names)
    label_w   <- .label_width(disp, min = 13L)
    hdr_pad   <- strrep(" ", 17L + label_w - 13L)

    # Omnibus Wald test on class intercepts (class-specific means at Z = 0).
    # Uses the model-based SE for each intercept (sigma^2 * B_inv_k[1,1]),
    # assuming classes are independent (a Wald test of equality).
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

    est_rows <- list()
    cat("\n")
    for (k in seq_len(K)) {
      cat(sprintf("Class %d:\n", k))
      cat(hdr_pad, "Estimate   [95% CI]        P-Value\n", sep = "")
      for (v in seq_len(D)) {
        est   <- betas[k, v]
        se    <- ses[k, v]
        z_val <- if (se > 0) est / se else NA_real_
        p_val <- if (!is.na(z_val)) 2 * (1 - pnorm(abs(z_val))) else NA_real_
        cat(sprintf("  %-*s %7.3f  [%6.3f, %6.3f]  %s\n",
                    label_w, disp[v], est, est - 1.96 * se, est + 1.96 * se,
                    .fmt_pval(p_val)))
        est_rows[[length(est_rows) + 1L]] <- data.frame(
          class = k, term = var_names[v], estimate = est, se = se, p = p_val,
          lower = est - 1.96 * se, upper = est + 1.96 * se,
          stringsAsFactors = FALSE)
      }
      cat("\n")
    }

    #  Per-covariate Wald(=) tests: H0: slope_k equal across all classes
    # Uses the diagonal independence approximation (separate per-class
    # regressions).
    # Contrast matrix R = [-1 | I_{K-1}], df = K-1.
    eq_rows <- list()
    if (K > 1L && D > 1L) {
      cat("---------------------------------------------------------\n")
      cat("Wald tests (equality of slopes across classes):\n")
      cat(sprintf("  %-*s   Wald(%s)%s  P-Value\n",
                  label_w, "", paste0("chi^2(", K - 1L, ")"), ""))
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
        cat(sprintf("  %-*s   %8.2f          %s\n",
                    label_w, disp[v], W_v, .fmt_pval(p_v)))
        eq_rows[[length(eq_rows) + 1L]] <- data.frame(
          term = var_names[v], chi2 = W_v, df = K - 1L, p = p_v,
          stringsAsFactors = FALSE)
      }
    }
    .cat_label_legend(disp)
    out$outcome <- list(
      type      = "continuous",
      omnibus   = data.frame(chi2 = omni$stat, df = omni$df, p = omni$p),
      estimates = .rbind_tidy(est_rows),
      slope_equality = .rbind_tidy(eq_rows))
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
    # This omnibus formulation accounts for the
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

    int_rows <- list()
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
      int_rows[[length(int_rows) + 1L]] <- data.frame(
        class = k, estimate = est, se = se, p = p_val,
        lower = est - 1.96 * se, upper = est + 1.96 * se,
        stringsAsFactors = FALSE)
    }

    slope_rows <- list()
    if (D_cov > 0) {
      disp    <- .shorten_labels(var_names)
      label_w <- .label_width(disp, min = 11L)
      cat("\n  Covariates (Pooled Slopes):\n")
      cat(strrep(" ", 17L + label_w - 11L),
          "Estimate   [95% CI]        P-Value\n", sep = "")
      for (v in seq_len(D_cov)) {
        est   <- theta[K_distal + v]
        se    <- ses[K_distal + v]
        z_val <- if (se > 0) est / se else NA_real_
        p_val <- if (!is.na(z_val)) 2 * (1 - pnorm(abs(z_val))) else NA_real_
        cat(sprintf("    %-*s %7.3f  [%6.3f, %6.3f]  %s\n",
                    label_w, disp[v], est, est - 1.96 * se, est + 1.96 * se,
                    .fmt_pval(p_val)))
        slope_rows[[length(slope_rows) + 1L]] <- data.frame(
          term = var_names[v], estimate = est, se = se, p = p_val,
          lower = est - 1.96 * se, upper = est + 1.96 * se,
          stringsAsFactors = FALSE)
      }
      .cat_label_legend(disp)
    }
    out$outcome <- list(
      type       = "continuous",
      omnibus    = data.frame(chi2 = omni$stat, df = omni$df, p = omni$p),
      intercepts = .rbind_tidy(int_rows),
      covariate_effects = .rbind_tidy(slope_rows))
  }

  cat("=========================================================\n")
  invisible(out)
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
#'   \code{"gaussian"} / \code{"gaussian_unit"},
#'   \code{"count"} / \code{"poisson"}.
#'   Missing values are handled automatically: any indicator column containing
#'   \code{NA} is estimated with a full-information (FIML) variant that masks the
#'   missing cells under a missing-at-random assumption, while complete columns
#'   use the faster complete-data estimator. A single specification (e.g.
#'   \code{"binary"}) therefore covers both complete and incomplete data, and the
#'   fitted object reports any missingness it found. Cases missing on
#'   \emph{every} indicator are the exception: they carry no information about
#'   class membership, so they are deleted before estimation and
#'   reported by \code{print()} and in \code{$missing_data$n_empty_rows}. The explicit \code{"*_nan"}
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
#' @param weight_type Either \code{"sampling"} (default; weights are rescaled to
#'   sum to the number of cases) or \code{"frequency"} (weights are counts of
#'   identical cases and set the effective sample size).
#' @param strata Optional vector of stratum identifiers for complex survey designs.
#' @param cluster Optional vector of cluster identifiers for complex survey designs.
#' @param bayes_constants Optional named list of prior strengths
#'   (\code{latent}, \code{categorical}, \code{poisson}, \code{variances}), each
#'   defaulting to \code{1}. See \code{\link{fit_mixture}}.
#' @param refine Logical. If \code{TRUE} (default), applies L-BFGS refinement
#'   after EM convergence to optimize the penalized maximum likelihood.
#' @param se Character. How standard errors for a covariate (class-prediction)
#'   structural model are computed when \code{n_steps} is \code{2} or \code{3}.
#'   \code{"corrected"} (default) is the first-order corrected estimator of
#'   Bakk et al. (2014): the step-3 sandwich plus the variance
#'   propagated from step 1. \code{"robust"} keeps only the sandwich.
#'   \code{"hessian"} inverts the
#'   step-3 observed information alone. See \code{\link{covariate_se}} for the
#'   differences and when they matter. Ignored for other structural models and
#'   for \code{n_steps = 1}.
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
#' @importFrom stats complete.cases cov dnorm dpois optim pchisq plogis pnorm qlogis qnorm rbinom rnorm rpois runif sd var
#' @importFrom utils setTxtProgressBar txtProgressBar
fit_mixture_internal <- function(X, Y = NULL, n_components = 2,
                                 measurement = "binary", structural = NULL,
                                 n_steps = 1, correction = "none", n_init = 1,
                                 max_iter = 1000, random_state = NULL,
                                 order_by_size = TRUE, weights = NULL,
                                 weight_type = c("sampling", "frequency"),
                                 strata = NULL, cluster = NULL,
                                 refine = TRUE,
                                 bayes_constants = NULL,
                                 warm_start = NULL,
                                 se = c("corrected", "robust", "hessian"), ...) {

  weight_type <- match.arg(weight_type)
  se          <- match.arg(se)
  bayes_constants <- .resolve_bayes_constants(bayes_constants)

  if (is.data.frame(X)) X <- as.matrix(X)
  # Convert Y through prepare_covariates() so that:
  #   - numeric columns are passed through unchanged
  #   - factor / character columns are dummy-coded (first level = reference)
  #   - column names are always preserved for display in summary()
  if (!is.null(Y)) Y <- prepare_covariates(Y)

  # --- Cases with no observed indicator data ----------------------------------
  # Dropped before anything else is computed, so that the sample size, the
  # missingness summary, the entropy, and the modal class counts all describe
  # the cases actually analysed. See .empty_rows() for why these cases distort
  # those quantities while contributing nothing to the likelihood.
  empty_rows <- .empty_rows(X)
  n_input_rows <- nrow(X)

  if (length(empty_rows) > 0L) {
    if (length(empty_rows) == n_input_rows)
      stop("Every case is missing on all indicators, so there are no data to ",
           "fit. This usually means the indicator columns were not read ",
           "correctly.", call. = FALSE)

    # Row-aligned arguments are length-checked against the input before any
    # subsetting: a mis-specified vector must still produce its own clear error
    # rather than being silently truncated to the retained rows.
    if (!is.null(weights) && length(weights) != n_input_rows)
      stop(sprintf("`weights` must have one entry per case (expected %d, got %d).",
                   n_input_rows, length(weights)), call. = FALSE)
    if (!is.null(strata) && length(strata) != n_input_rows)
      stop("Length of strata must match rows of X.", call. = FALSE)
    if (!is.null(cluster) && length(cluster) != n_input_rows)
      stop("Length of cluster must match rows of X.", call. = FALSE)

    keep <- setdiff(seq_len(n_input_rows), empty_rows)
    X    <- X[keep, , drop = FALSE]
    if (!is.null(Y))       Y       <- Y[keep, , drop = FALSE]
    if (!is.null(weights)) weights <- weights[keep]
    if (!is.null(strata))  strata  <- strata[keep]
    if (!is.null(cluster)) cluster <- cluster[keep]

    warning(sprintf(
      paste0("%d case%s had no observed value on any indicator and %s removed ",
             "before estimation (n = %d analysed). Rows: %s."),
      length(empty_rows), if (length(empty_rows) == 1L) "" else "s",
      if (length(empty_rows) == 1L) "was" else "were", length(keep),
      .abbreviate_indices(empty_rows)), call. = FALSE)
  }

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
    structural_handled_by = struct_handled,
    # Cases removed for having no observed indicator at all. Kept separate from
    # the FIML summary above, which describes only the analysed cases.
    n_empty_rows     = length(empty_rows),
    empty_rows       = empty_rows,
    n_input_rows     = n_input_rows
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

  # Validate count data when a Poisson family is requested. Negative or
  # fractional values give dpois() a probability of zero for every class, so the
  # fit would fail with a log-likelihood of -Inf rather than a usable message.
  if (is.character(measurement) &&
      measurement %in% c("count", "poisson", "count_nan", "poisson_nan")) {
    valid_vals <- X[!is.na(X)]
    bad <- valid_vals[valid_vals < 0 | abs(valid_vals - round(valid_vals)) > 1e-8]
    if (length(bad) > 0)
      stop(sprintf(
        paste0("measurement = '%s' requires non-negative integer counts. ",
               "Found values outside this set: %s"),
        measurement,
        paste(utils::head(sort(unique(bad)), 5), collapse = ", ")
      ), call. = FALSE)
  }

  # Weights that are whole numbers adding up to far more than the number of rows
  # are almost always frequency counts, and treating them as sampling weights
  # would report a sample size of 31 patterns where the study had 631 people.
  # Say so rather than silently producing a BIC on the wrong scale.
  if (weight_type == "sampling" && .looks_like_frequencies(weights, n_samples))
    message("These weights look like frequency counts (whole numbers summing ",
            "to ", format(sum(weights)), " across ", n_samples, " rows). ",
            "They are being treated as sampling weights; use ",
            "weight_type = \"frequency\" if each row stands for that many cases.")

  wt <- .resolve_weights(weights, n_samples, weight_type)
  weights <- wt$weights
  n_eff   <- wt$n_eff

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
    weight_type           = weight_type,
    n_eff                 = n_eff,
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
    missing_data           = missing_data,
    # Prior strengths (see R/bayes_constants.R). Stored on the fit so that
    # print()/summary() can report a non-default choice and bootstrap replicates
    # inherit it.
    bayes_constants        = bayes_constants
  )
  class(model_state) <- "mixture_model"

  # Push the constants down onto the emissions, recursively, so that every
  # M-step and refine_lbfgs() read one object rather than each holding its own
  # copy of a default.
  model_state$mm <- .attach_bayes_constants(model_state$mm, bayes_constants)
  model_state$sm <- .attach_bayes_constants(model_state$sm, bayes_constants)

  # Mirror the survey design onto the structural sub-model (see
  # .mirror_design_onto_sm in R/stepwise.R).
  model_state <- .mirror_design_onto_sm(model_state)

  if (n_steps == 1) {
    model_state <- fit_em(model_state, X, Y, n_init, max_iter, random_state,
                          refine = refine, warm_start = warm_start)

  } else {
    model_state <- fit_em(model_state, X, NULL, n_init, max_iter, random_state,
                          refine = refine, warm_start = warm_start)

    # Step 1 metrics (measurement model only)
    if (n_steps == 3)
      model_state$step1_metrics <- .step1_metrics(model_state)

    model_state <- .apply_structural_steps(model_state, X, Y, n_steps,
                                           correction, max_iter, se)
  }

  # Class sorting, display names, and combined-model metrics (see
  # .finalize_model_state in R/stepwise.R).
  model_state <- .finalize_model_state(model_state, X, Y, order_by_size)

  # Collapsed-variance check, after sorting so the class numbers it reports are
  # the ones the user will see. See R/gaussian_boundary.R.
  model_state <- .check_gaussian_degeneracy(model_state, X)

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
  gaussian      = "gaussian_nan",
  poisson       = "poisson_nan",
  count         = "count_nan"
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
#'   \code{"gaussian"}, \code{"count"}, and \code{"_nan"} missing-data
#'   variants), or, for a mixed-type model, a named list mapping each type to
#'   the indicator columns it governs by name or index, e.g.
#'   \code{list(binary = c("q1","q2"), continuous = "score")}.
#'
#'   \code{"count"} fits a Poisson measurement model: each item is
#'   Poisson-distributed within class with its own rate, reported by
#'   \code{\link{measurement_summary}} as a rate (lambda) per item and class.
#'   It requires non-negative integer indicators.
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
#' @param group Optional observed grouping variable for a multiple-group
#'   model (Collins & Lanza, 2010, sec. 5.7-5.12), e.g. grade or gender.
#'   Unlike \code{predictors}, which only ever shifts class membership, a
#'   group can also be allowed to shift the item-response probabilities
#'   themselves; see \code{group_effects}.
#' @param group_effects Which parameters \code{group} is allowed to shift.
#'   \code{"both"} (default) frees both the item-response probabilities and
#'   the class prevalences across groups, i.e. fits each group's own
#'   model. \code{"measurement"} frees only the item-response probabilities
#'   (prevalences stay pooled). \code{"prevalence"} frees only the class
#'   prevalences, by entering \code{group} as a class-membership covariate
#'   exactly like \code{predictors} (item-response probabilities stay
#'   invariant across groups) — this is the same equivalence
#'   [`fit_rmlca()`] documents for its own \code{predictors}, sec. 6.10.2.
#'   \code{"none"} ignores \code{group} for estimation. Fit the same data
#'   under two of these and compare them with [`longitudinal_lrt()`] to get
#'   Collins & Lanza's measurement-invariance test (sec. 5.8, comparing
#'   \code{"both"} against \code{"prevalence"}) and prevalence-equivalence
#'   test (sec. 5.11, comparing \code{"prevalence"} against \code{"none"}).
#'   Because a multiple-group model like Collins & Lanza's is fit as one
#'   simultaneous model rather than through the auxiliary-variable 3-step
#'   approximation, pass \code{n_steps = 1} explicitly to reproduce their
#'   numbers exactly.
#'
#'   With \code{"both"} or \code{"measurement"}, each group's item-response
#'   probabilities are estimated from that group's cases alone, tied to the
#'   other groups only by shared initialization — unlike occasions in
#'   [`fit_rmlca()`], which every case informs simultaneously. A likelihood-
#'   ratio comparison via [`longitudinal_lrt()`] is unaffected (it only
#'   compares total log-likelihoods), but the *class labels* a configural fit
#'   assigns are not guaranteed to line up across groups: "Class 1" in one
#'   group's profile need not be the same kind of class as "Class 1" in
#'   another's. Use a larger \code{n_init} and a fixed \code{random_state}
#'   for configural fits, and when reading per-group profiles, match classes
#'   by their item-response pattern rather than by position.
#' @param group_invariant_items Item indices or names held equal across
#'   groups even when \code{group_effects} frees the measurement model
#'   (Collins & Lanza's partial-invariance models, sec. 5.9). \code{NULL}
#'   (the default) leaves every item free, i.e. a fully configural model.
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
#' @param weight_type What the numbers in \code{weights} mean.
#'   \code{"sampling"} (the default) treats them as survey or probability
#'   weights, saying how much of the population each case represents; only their
#'   relative sizes matter, so they are rescaled to sum to the number of cases.
#'   \code{"frequency"} treats them as counts of identical cases, as in a
#'   response-pattern table where one row stands for many respondents; the
#'   sample size behind AIC and BIC is then the sum of the counts. Getting this
#'   wrong changes BIC but not the parameter estimates.
#' @param n_init,max_iter,random_state,order_by_size,refine Estimation
#'   controls: number of random starts, maximum EM iterations, RNG seed,
#'   whether to order classes by size, and whether to run L-BFGS refinement.
#' @param bayes_constants Optional list adjusting the strength of the weak
#'   priors the estimator places on each block of parameters. Named
#'   \code{latent} (class weights and, in the transition models, the initial
#'   and transition probabilities), \code{categorical} (item-response
#'   probabilities), \code{poisson} (count rates), and \code{variances}
#'   (class-specific variances of continuous indicators); all default to
#'   \code{1}. Each is a number of pseudo-observations spread over the classes,
#'   so its influence shrinks as the sample grows.
#'
#'   Two uses. \strong{Reproducing an unregularized fit:} setting a constant to
#'   \code{0} removes that prior and gives plain maximum likelihood for that
#'   block. This is an escape hatch for matching a reference analysis, not a
#'   recommended setting — the unpenalised mixture likelihood is unbounded, and
#'   maximum likelihood for a mixture of normals is known to be inconsistent
#'   without some such restriction (Kiefer & Wolfowitz, 1956).
#'   \strong{Rescuing a collapsed fit:} if a fit warns that a class variance has
#'   collapsed, \code{bayes_constants = list(variances = 5)} is the recommended
#'   remedy. That value is empirical — calibrated on continuous indicators
#'   scored on a five-point scale — rather than derived, so check that it has
#'   not moved the parameters you care about.
#'
#'   This is not a tuning menu. The defaults are the intended settings, and
#'   \code{n_init} with \code{random_state} remains the way to search harder for
#'   a solution.
#' @param se How standard errors for \code{predictors} are computed in a 2- or
#'   3-step model. \code{"corrected"} (the default) adds the variance carried
#'   over from step 1 to the step-3 sandwich, following Bakk et al.
#'   (2014); \code{"robust"} reports the sandwich alone;
#'   \code{"hessian"} inverts the step-3 observed
#'   information only. See \code{\link{covariate_se}}.
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
                        group = NULL,
                        group_effects = c("both", "measurement", "prevalence", "none"),
                        group_invariant_items = NULL,
                        n_steps = 1,
                        correction = "none",
                        n_init = 1,
                        max_iter = 1000,
                        random_state = NULL,
                        order_by_size = TRUE,
                        weights = NULL,
                        weight_type = c("sampling", "frequency"),
                        strata = NULL,
                        cluster = NULL,
                        refine = TRUE,
                        bayes_constants = NULL,
                        se = c("corrected", "robust", "hessian"),
                        X = NULL, Y = NULL, n_components = NULL, structural = NULL,
                        ...) {

  se            <- match.arg(se)
  outcome_type  <- match.arg(outcome_type)
  slopes        <- match.arg(slopes)
  group_effects <- match.arg(group_effects)
  weight_type   <- match.arg(weight_type)
  steps_set     <- !missing(n_steps)
  corr_set      <- !missing(correction)

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
      weight_type = weight_type,
      strata = strata, cluster = cluster, refine = refine,
      bayes_constants = bayes_constants, se = se, ...))
  }

  if (is.null(indicators))
    stop("`indicators` is required: the matrix of measurement items that ",
         "define the latent classes.", call. = FALSE)

  if (!is.null(predictors) && !is.null(outcome))
    stop("Specify either `predictors` (to model class membership) or ",
         "`outcome` (a distal outcome), not both in one model. To run both ",
         "analyses from one solution, fit the unconditional model and use ",
         "add_covariates() and add_outcome() on it.", call. = FALSE)
  if (!is.null(outcome_covariates) && is.null(outcome))
    stop("`outcome_covariates` requires an `outcome`.", call. = FALSE)
  if (!is.null(group) && group_effects %in% c("both", "prevalence") &&
      !is.null(outcome))
    stop("`group` with a prevalence effect (`group_effects = \"both\"` or ",
         '"prevalence") uses the same structural-model slot as `outcome`; ',
         "combine `group` with `predictors` instead, or set `group_effects` ",
         'to "measurement" or "none".', call. = FALSE)

  # --- Measurement model (single-type or mixed) -------------------------------
  mm                 <- .normalize_measurement(measurement, indicators)
  X_use              <- mm$indicators
  measurement_engine <- mm$descriptor

  # --- Grouping variable: measurement effect (multiple-group model) ----------
  # Reuses the time-blocks trick across a group axis instead of a time axis
  # (R/group_blocks.R): each group gets its own copy of the item block, and a
  # case's other-group blocks are structurally missing, which FIML already
  # handles correctly. The prevalence effect is handled below, alongside
  # `predictors`, since both use the same class-membership regression.
  group_info       <- NULL
  group_extra_args <- list()
  group_warm_start <- NULL
  if (!is.null(group)) {
    group_info <- .lta_group_design(group, nrow(X_use))

    if (group_effects %in% c("both", "measurement")) {
      # Warm start (see R/group_blocks.R). The pooled measurement model is fitted
      # first, on the unpadded data, and replicated into every group block as one
      # extra restart of the group-varying search. It is the same model the
      # invariance test compares against, so this is what makes the comparison a
      # test rather than a lower bound. Fitted before X_use is padded and
      # `measurement_engine` is overwritten below, since both are needed as they
      # stand here.
      .sub_fit <- function(rows) {
        args <- c(list(X = if (is.null(rows)) X_use else X_use[rows, , drop = FALSE],
                       Y = NULL, n_components = n_classes,
                       measurement = measurement_engine, structural = NULL,
                       n_steps = 1, correction = "none", n_init = n_init,
                       max_iter = max_iter, random_state = random_state,
                       order_by_size = order_by_size,
                       weights = if (is.null(rows)) weights else weights[rows],
                       weight_type = weight_type,
                       strata  = if (is.null(rows)) strata  else strata[rows],
                       cluster = if (is.null(rows)) cluster else cluster[rows],
                       refine = refine, bayes_constants = bayes_constants,
                       se = se),
                  list(...))
        out <- try(suppressWarnings(suppressMessages(
          do.call(fit_mixture_internal, args))), silent = TRUE)
        if (inherits(out, "try-error")) NULL else out
      }

      pooled_fit <- .sub_fit(NULL)
      # A failure here costs the warm start and nothing else: the group model is
      # still fitted from random starts exactly as before.
      if (!is.null(pooled_fit)) {
        g_idx     <- as.integer(group_info$factor)
        group_mms <- lapply(seq_len(nlevels(group_info$factor)), function(g) {
          rows <- which(g_idx == g)
          # A group too small to identify K classes on its own is left to the
          # pooled parameters rather than fitted badly.
          if (length(rows) < 2L * n_classes) return(NULL)
          one <- .sub_fit(rows)
          if (is.null(one)) NULL else .align_to_pooled(one$mm, pooled_fit$mm)
        })
        group_warm_start <- .group_blocks_warm_start(pooled_fit, group_mms)
      }

      item_names <- colnames(X_use) %||% paste0("Item", seq_len(ncol(X_use)))
      grp_spec <- .resolve_invariance(
        if (is.null(group_invariant_items)) "none" else "partial",
        group_invariant_items, item_names, measurement)
      n_groups <- nlevels(group_info$factor)
      X_grp  <- .pad_group_blocks(X_use, group_info$factor)
      engine <- .longitudinal_measurement_spec(measurement, X_grp,
                                               n_items = ncol(X_use),
                                               n_times = n_groups)
      X_use              <- X_grp
      measurement_engine <- "group_blocks"
      group_extra_args <- list(
        n_items         = length(item_names),
        n_groups        = n_groups,
        sub_model       = engine$sub_model,
        invariant_items = grp_spec$invariant_items,
        max_val         = engine$max_val
      )
    }
  }

  # --- Structural model -------------------------------------------------------
  structural_engine <- NULL
  Y_use             <- NULL

  if (!is.null(predictors)) {
    structural_engine <- "predict_class"
    Y_use             <- prepare_covariates(
      .as_named_covariates(predictors, predictors_expr, "predictor"))

  } else if (!is.null(outcome)) {
    outcome_spec      <- .build_outcome_spec(outcome, outcome_covariates,
                                             outcome_type, slopes,
                                             outcome_cov_expr)
    structural_engine <- outcome_spec$engine
    Y_use             <- outcome_spec$Y
  }

  # --- Grouping variable: prevalence effect -----------------------------------
  # Enters `group` as a class-membership covariate, exactly like `predictors`
  # (and combined with it, if both are given): each group gets its own free
  # prevalences while the measurement model above is whatever the
  # `group_effects` "measurement" branch left it (pooled, unless that branch
  # also ran).
  if (!is.null(group) && group_effects %in% c("both", "prevalence")) {
    structural_engine <- "predict_class"
    Y_use <- if (is.null(Y_use)) group_info$design
             else .cbind_covariates(Y_use, group_info$design)
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

  dots <- if (length(group_extra_args)) utils::modifyList(list(...), group_extra_args)
          else list(...)

  fit <- do.call(fit_mixture_internal, c(list(
    X = X_use, Y = Y_use, n_components = n_classes,
    measurement = measurement_engine, structural = structural_engine,
    n_steps = n_steps, correction = correction, n_init = n_init,
    max_iter = max_iter, random_state = random_state,
    order_by_size = order_by_size, weights = weights,
    weight_type = weight_type,
    strata = strata, cluster = cluster, refine = refine,
    bayes_constants = bayes_constants, warm_start = group_warm_start,
    se = se), dots))

  if (!is.null(group)) {
    fit$group_info    <- group_info
    fit$group_effects <- group_effects
  }

  # Non-convergence used to be near-impossible to hit and was reported only by
  # print(), which a user working from summary() or the coefficients never sees.
  # It became reachable when the emissions L-BFGS does not polish were given a
  # stopping rule tight enough to be trusted: those models genuinely need
  # hundreds of iterations, and the default cap is 1000. Silently returning the
  # iterate EM happened to be on at the cap is the one outcome worth
  # interrupting for, so say it here rather than leaving it to be noticed.
  if (isFALSE(fit$converged))
    warning(sprintf(
      paste0("EM did not converge within max_iter = %d iterations. The ",
             "estimates are wherever the algorithm had reached, which need ",
             "not be a maximum. Refit with a larger `max_iter`; if that does ",
             "not help, the model is probably weakly identified at this ",
             "number of classes."),
      max_iter), call. = FALSE)

  fit
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
  # Deleted cases are reported on their own line, ahead of the FIML summary, so
  # the printed sample size can always be reconciled with the input data.
  if (isTRUE(x$missing_data$n_empty_rows > 0L))
    cat(sprintf("Cases Removed      : %d of %d with no observed indicator (n = %d analysed)\n",
                x$missing_data$n_empty_rows, x$missing_data$n_input_rows,
                x$missing_data$n_input_rows - x$missing_data$n_empty_rows))
  if (!is.null(x$missing_data) && isTRUE(x$missing_data$any_missing)) {
    md <- x$missing_data
    cat(sprintf("Missing Data       : %d / %d cells (%.1f%%) in %d item%s \u2014 %s\n",
                md$n_missing, md$n_cells, 100 * md$prop_missing,
                md$n_items_affected, if (md$n_items_affected == 1L) "" else "s",
                md$handled_by))
  }
  if (identical(x$weight_type, "frequency"))
    cat(sprintf("Case Weights       : frequency counts (%s cases in %d rows)\n",
                format(x$n_eff), length(x$sample_weights)))
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
  # A warning is transient; someone opening a saved fit months later should
  # still see that its variances collapsed. See R/gaussian_boundary.R.
  .print_degenerate_note(x)
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
