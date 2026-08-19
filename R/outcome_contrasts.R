# ==============================================================================
# Class-vs-class contrasts on a distal outcome
# ==============================================================================
# The omnibus test printed by summary() answers "do the classes differ on this
# outcome at all", which is the question a reviewer asks first and almost never
# the question the paper is about. What gets written up is which classes differ:
# the high-risk class scores half a standard deviation above the normative one,
# the two intermediate classes are indistinguishable. Those are contrasts, and
# every one of them is a linear combination of parameters the fit already
# carries, together with a covariance matrix it already estimated.
#
# The covariance is the whole reason this belongs in the package rather than in
# a user's script. Two class means from the same fit are not independent - the
# same posteriors, and under BCH the same inverted classification-error matrix,
# enter both - so the standard error of their difference is not the root of the
# sum of their squared standard errors. Off the diagonal of Sigma_mu the two
# can differ by a good deal, and always in the direction that matters: positive
# covariance between the two means makes the naive interval too wide and the
# test too conservative, so the difference a user computes by hand is the one
# they are most likely to miss.

# Find a distal sub-model of a given class, whether the structural model is that
# model or a nested one holding it alongside a covariate model. Three copies of
# this unwrapping already existed in summary(); this is the one the new code
# uses, and `exclude` covers the one case where the classes overlap
# (distal_categorical inherits from distal_pooled).
.distal_submodel <- function(sm, what, exclude = NULL) {
  ok <- function(m) {
    if (is.null(m) || !inherits(m, what)) return(FALSE)
    if (!is.null(exclude) && inherits(m, exclude)) return(FALSE)
    TRUE
  }
  if (ok(sm)) return(sm)
  if (inherits(sm, "nested") && "distal" %in% names(sm$models) &&
      ok(sm$models$distal))
    return(sm$models$distal)
  NULL
}

# The (class, reference) pairs to report, in printing order.
#
# With no reference given, every unordered pair, each reported once as the
# higher-numbered class against the lower. With one given, that class is held
# fixed on the right-hand side, which is the layout the categorical outcome
# tables already use.
.contrast_pairs <- function(K, ref = NULL) {
  if (is.null(ref)) {
    pairs <- utils::combn(K, 2L)
    return(data.frame(class = pairs[2L, ], reference = pairs[1L, ]))
  }
  data.frame(class = setdiff(seq_len(K), ref), reference = ref)
}

# One contrast: the difference between two entries of theta, and its standard
# error from the corresponding 2x2 block of V.
#
# Var(theta_a - theta_b) = V[a,a] + V[b,b] - 2 V[a,b]. The cross term is the
# point; dropping it is the mistake this function exists to prevent.
.contrast_row <- function(theta, V, a, b) {
  est <- theta[a] - theta[b]
  var <- V[a, a] + V[b, b] - 2 * V[a, b]
  se  <- sqrt(max(0, var))
  z   <- if (is.finite(se) && se > 0) est / se else NA_real_
  p   <- if (is.na(z)) NA_real_ else 2 * stats::pnorm(-abs(z))
  list(estimate = est, se = se, z = z, p = p)
}

#' Class-vs-Class Contrasts on a Distal Outcome
#'
#' @description
#' Tests which latent classes differ on a distal outcome attached by
#' [add_outcome()], rather than whether any of them do. The omnibus Wald test
#' printed by [summary()] is a single joint statistic; this is the set of
#' contrasts behind it, each with a standard error, a confidence interval and a
#' p-value.
#'
#' Each contrast is reported as `class` minus `reference`, so a positive
#' estimate means the class scores above the class it is compared with. With
#' `ref = NULL` every pair of classes is reported once; naming a class in `ref`
#' holds it fixed as the comparison.
#'
#' The standard errors account for the covariance between the class parameters,
#' which is not optional here. Two class means from one fit are estimated from
#' the same posteriors, and under the BCH correction from the same inverted
#' classification-error matrix, so they are correlated; the standard error of
#' their difference is not the root of the sum of their squared standard errors,
#' and the naive version is conservative exactly when the correlation is
#' positive. The full sandwich covariance the fit already carries is used
#' wherever it is available.
#'
#' @param fit A fitted model with a distal outcome attached, as returned by
#'   [add_outcome()].
#' @param ref Optional reference class, as an integer between 1 and the number
#'   of classes. `NULL` (the default) reports all pairs.
#' @param adjust Multiplicity adjustment across the reported contrasts, passed
#'   to [stats::p.adjust()]: `"none"` (the default), `"holm"` or
#'   `"bonferroni"`. With `"none"` the p-values are the per-contrast ones, which
#'   is what to report when the contrasts were named in advance; an adjustment
#'   belongs on the all-pairs table read as a family.
#' @param level Confidence level for the intervals. Default `0.95`.
#' @param ... Currently unused.
#'
#' @return A data frame of class `outcome_contrasts`, one row per contrast, with
#'   columns `category` (the outcome category the contrast is on, `NA` for a
#'   continuous outcome), `class`, `reference`, `estimate`, `se`, `lower`,
#'   `upper`, `z` and `p`, plus `p_adj` when `adjust` is not `"none"` and
#'   `OR`, `OR_lower`, `OR_upper` for a categorical outcome, whose estimates are
#'   log odds. The `method` attribute names the covariance the standard errors
#'   came from.
#'
#' @seealso [add_outcome()] for fitting the outcome model and [summary()] for
#'   the omnibus test; [wald_omnibus_test()] and [bootstrap_covariates()] for
#'   the bootstrap route, which is what to use for an outcome whose parameters
#'   are estimated one class at a time.
#'
#' @examples
#' set.seed(1)
#' items <- matrix(rbinom(300, 1, 0.5), nrow = 100)
#' bmi   <- rnorm(100, mean = 25)
#' fit   <- fit_mixture(items, n_classes = 2, measurement = "binary")
#' fit_out <- add_outcome(fit, bmi)
#' outcome_contrasts(fit_out)
#'
#' @export
outcome_contrasts <- function(fit, ref = NULL,
                              adjust = c("none", "holm", "bonferroni"),
                              level = 0.95, ...) {
  adjust <- match.arg(adjust)
  if (!inherits(fit, "mixture_model"))
    stop("`fit` must be a fitted mixture model.", call. = FALSE)

  K <- fit$n_components
  if (is.null(K) || K < 2L)
    stop("Class contrasts need at least two classes.", call. = FALSE)
  if (!is.null(ref)) {
    if (!is.numeric(ref) || length(ref) != 1L || ref < 1 || ref > K)
      stop(sprintf("`ref` must be a class index between 1 and %d.", K),
           call. = FALSE)
    ref <- as.integer(ref)
  }
  if (!is.numeric(level) || length(level) != 1L || level <= 0 || level >= 1)
    stop("`level` must be a confidence level strictly between 0 and 1.",
         call. = FALSE)

  built <- .outcome_contrast_parts(fit, K)
  pairs <- .contrast_pairs(K, ref)

  rows <- list()
  for (part in built$parts) {
    for (i in seq_len(nrow(pairs))) {
      a <- pairs$class[i]; b <- pairs$reference[i]
      r <- .contrast_row(part$theta, part$V, part$index[a], part$index[b])
      rows[[length(rows) + 1L]] <- data.frame(
        category  = part$category,
        class     = a,
        reference = b,
        estimate  = r$estimate, se = r$se, z = r$z, p = r$p,
        stringsAsFactors = FALSE)
    }
  }

  out <- .rbind_tidy(rows)
  crit <- stats::qnorm(1 - (1 - level) / 2)
  out$lower <- out$estimate - crit * out$se
  out$upper <- out$estimate + crit * out$se

  if (adjust != "none") out$p_adj <- stats::p.adjust(out$p, method = adjust)

  cols <- c("category", "class", "reference", "estimate", "se",
            "lower", "upper", "z", "p")
  if (adjust != "none") cols <- c(cols, "p_adj")
  out <- out[, cols, drop = FALSE]

  if (identical(built$type, "categorical")) {
    out$OR       <- exp(out$estimate)
    out$OR_lower <- exp(out$lower)
    out$OR_upper <- exp(out$upper)
  }

  rownames(out) <- NULL
  attr(out, "method")       <- built$method
  attr(out, "outcome_type") <- built$type
  attr(out, "level")        <- level
  attr(out, "adjust")       <- adjust
  class(out) <- c("outcome_contrasts", "data.frame")
  out
}

# Reduce whichever distal model is attached to the pieces a contrast needs: a
# parameter vector, its covariance, and the positions in that covariance the
# classes occupy. One "part" per outcome category, so a polytomous outcome is a
# stack of the same contrasts on each of its non-reference categories.
.outcome_contrast_parts <- function(fit, K) {

  # --- Continuous outcome: class means -----------------------------------
  cont <- .distal_submodel(fit$sm, "distal_continuous")
  if (!is.null(cont) && !is.null(cont$parameters$means)) {
    means <- as.vector(cont$parameters$means)
    Sigma <- cont$parameters$Sigma_mu

    # Sigma_mu is stored only by the BCH correction. Without it the class means
    # are treated as independent, which is what the omnibus test falls back to
    # as well; the attribute says so, because a user comparing the two needs to
    # know which one they got.
    if (is.null(Sigma) || !all(is.finite(Sigma))) {
      ses    <- as.vector(cont$parameters$ses)
      Sigma  <- diag(ses^2, nrow = K)
      method <- "independent class means (no cross-class covariance stored)"
    } else {
      method <- "sandwich covariance of the class means"
    }

    return(list(
      type   = "continuous",
      method = method,
      parts  = list(list(category = NA_integer_, theta = means, V = Sigma,
                         index = seq_len(K)))))
  }

  # --- Categorical outcome: class intercepts, one set per category --------
  #
  # beta_pooled is (M-1) x L with L = K + D_cov, and is vectorised row-major
  # into the Hessian, so category m's intercept for class k sits at
  # (m-1)*L + k. That indexing is fixed by .wald_omnibus_pooled() and by the
  # odds-ratio tables in summary(); it is repeated rather than re-derived.
  pooled <- .distal_submodel(fit$sm, "distal_pooled")
  if (!is.null(pooled) && !is.null(pooled$parameters$beta_pooled) &&
      nrow(pooled$parameters$beta_pooled) > 0L) {
    beta <- pooled$parameters$beta_pooled
    if (is.null(pooled$parameters$hessian))
      stop("This outcome model carries no covariance for its parameters, so ",
           "contrasts cannot be formed. Refit it with `se = \"corrected\"` or ",
           "`se = \"robust\"`.", call. = FALSE)

    V <- pinv(-pooled$parameters$hessian)
    L <- ncol(beta)

    parts <- lapply(seq_len(nrow(beta)), function(m) {
      list(category = m + 1L,
           theta    = as.vector(t(beta))[((m - 1L) * L + 1L):(m * L)],
           V        = V[((m - 1L) * L + 1L):(m * L),
                        ((m - 1L) * L + 1L):(m * L), drop = FALSE],
           index    = seq_len(K))
    })

    return(list(type = "categorical",
                method = "robust sandwich covariance of the outcome model",
                parts = parts))
  }

  # --- Everything else ----------------------------------------------------
  if (!is.null(.distal_submodel(fit$sm, "distal_regression")) ||
      !is.null(.distal_submodel(fit$sm, "distal_continuous_regression")))
    stop("This outcome model estimates a separate parameter block per class ",
         "and stores no covariance between the blocks, so class-vs-class ",
         "contrasts cannot be formed from it. Fit the outcome with ",
         "`slopes = \"pooled\"`, or use bootstrap_covariates() with ",
         "wald_omnibus_test().", call. = FALSE)

  stop("No distal outcome found on this fit. Attach one with add_outcome() ",
       "first.", call. = FALSE)
}

#' @export
print.outcome_contrasts <- function(x, ...) {
  cat("Class contrasts on the distal outcome\n")
  cat(sprintf("Standard errors: %s\n", attr(x, "method")))
  if (!identical(attr(x, "adjust"), "none"))
    cat(sprintf("P-values adjusted across %d contrasts (%s)\n",
                nrow(x), attr(x, "adjust")))
  cat("\n")

  poly <- any(!is.na(x$category))
  lvl  <- sprintf("[%g%% CI]", 100 * attr(x, "level"))
  # The widths here track the value format below exactly; the interval field is
  # 18 characters wide ("[" + 7 + ", " + 7 + "]"), so its heading is padded to
  # the same and the columns line up.
  cat(sprintf("  %-13s %9s %9s  %-18s  %s\n",
              if (poly) "Cat  Contrast" else "Contrast",
              "Estimate", "SE", lvl, "P-Value"))

  p_col <- if ("p_adj" %in% names(x)) x$p_adj else x$p
  for (i in seq_len(nrow(x))) {
    lab <- sprintf("%d vs %d", x$class[i], x$reference[i])
    if (poly) lab <- sprintf("%-4d %s", x$category[i], lab)
    cat(sprintf("  %-13s %9.3f %9.3f  [%7.3f, %7.3f]  %s\n",
                lab, x$estimate[i], x$se[i], x$lower[i], x$upper[i],
                .fmt_pval(p_col[i])))
  }
  invisible(x)
}
