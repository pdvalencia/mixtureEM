# ------------------------------------------------------------------------------
# Standard errors for the measurement parameters of a fit_mixture() model
# ------------------------------------------------------------------------------

#' Measurement Parameter Estimates with Standard Errors
#'
#' @description
#' Returns the measurement parameters of a [`fit_mixture()`] fit -- the item
#' probabilities, the class means and variances, the count rates -- and the
#' class-size logits, each with a standard error, a z statistic and a two-sided
#' p-value. It is the standard-error companion to [`measurement_summary()`],
#' which prints the same estimates without them.
#'
#' @details
#' The standard errors come from the log-likelihood of the fitted model,
#' differentiated with respect to the parameters on the scale the fit
#' estimates them on (logits for probabilities, logs for rates and standard
#' deviations, means as they are) and carried to the reported scale by the
#' delta method. The three estimators are those of [`growth_estimates()`]:
#'
#' * `"robust"` (the default) is the sandwich \eqn{H^{-1} B H^{-1}} (White,
#'   1982), which stays valid when the within-class model is only
#'   approximately right -- for instance, normal indicators whose real
#'   within-class distribution is skewed, or variances held equal across
#'   classes when they are not.
#' * `"hessian"` is the inverse observed information, the textbook
#'   maximum-likelihood variance when the model is correct.
#' * `"opg"` is the inverse outer product of the case scores.
#'
#' A gap between `"robust"` and `"hessian"` is a sign of misspecification.
#' Sample `weights` and a `cluster` variable are honoured as in
#' [`growth_estimates()`].
#'
#' A parameter held equal across classes -- a variance under
#' `variances_equal = TRUE` -- is repeated once per class, with the same
#' standard error, as in [`measurement_summary()`]. The prior that keeps
#' estimates off the edge of the parameter space (`bayes_constants` in
#' [`fit_mixture()`]) is not part of the variance: these are the standard
#' errors of the likelihood. With covariates on class membership fitted in one
#' step, the class probabilities belong to the covariate model and no
#' class-size logits are reported (see [`coef()`]); under `n_steps = 3` the
#' standard errors are those of the step-1 model, which is where the
#' measurement parameters are estimated.
#'
#' @section Probabilities at the boundary:
#' An item probability that has run to 0 or 1 (within 0.001) is flagged in
#' the `boundary` column, and a message says how many there are. Its standard
#' error is not interpretable: the likelihood is flat in that direction, so
#' any number would measure how far the optimiser stopped from the edge
#' rather than how precisely the data pin the probability down. A probability
#' whose logit has run past \eqn{\pm 15} is therefore held where it stopped,
#' as a fixed value, and its standard error is `NA`; the other parameters'
#' standard errors are those of the model with it fixed.
#'
#' @param fit A model returned by [`fit_mixture()`].
#' @param se Which standard-error estimator: `"robust"` (default), `"hessian"`
#'   or `"opg"`. See Details.
#'
#' @return A data frame with one row per parameter and class: `block` (the
#'   sub-model for a mixed or blocked measurement model, `NA` otherwise),
#'   `parameter` (`"probability"`, `"mean"`, `"variance"`, `"rate"` or
#'   `"class_logit"`), `item`, `category` (polytomous items only), `class`,
#'   `estimate`, `se`, `z`, `p_value` and `boundary`. The covariance matrix
#'   of the reported parameters is attached as the attribute `"vcov"`.
#'
#' @references
#' White, H. (1982). Maximum likelihood estimation of misspecified models.
#' *Econometrica*, *50*(1), 1-25. \doi{10.2307/1912526}
#'
#' Louis, T. A. (1982). Finding the observed information matrix when using the
#' EM algorithm. *Journal of the Royal Statistical Society: Series B*, *44*(2),
#' 226-233. \doi{10.1111/j.2517-6161.1982.tb01203.x}
#'
#' @seealso [`measurement_summary()`], [`growth_estimates()`] for growth
#'   models, and [`plot.mixture_model()`], whose `error_bars` draw these
#'   standard errors on the profile bar chart.
#'
#' @examples
#' set.seed(1)
#' X <- rbind(matrix(rnorm(300, 0), 100), matrix(rnorm(300, 2), 100))
#' fit <- fit_mixture(X, n_classes = 2, measurement = "continuous",
#'                    n_init = 5)
#' measurement_estimates(fit)
#' @export
measurement_estimates <- function(fit, se = c("robust", "hessian", "opg")) {
  se <- match.arg(se)
  if (inherits(fit, c("gmm", "lcga")))
    stop("This is a growth model; its parameters and standard errors are ",
         "reported by growth_estimates().", call. = FALSE)
  if (!inherits(fit, "mixture_model") || is.null(fit$mm))
    stop("`fit` must be a model returned by fit_mixture().", call. = FALSE)
  if (is.null(fit$data))
    stop("Raw indicators are not stored on this fit, so the likelihood ",
         "cannot be differentiated. Refit with the current version.",
         call. = FALSE)

  v   <- .reported_vcov(fit, .measurement_reported, "measurement model",
                        drop_walls = TRUE)
  V   <- v[[se]]
  out <- v$table
  est <- out$estimate
  s   <- if (is.null(V)) rep(NA_real_, length(est)) else sqrt(pmax(diag(V), 0))
  s[!is.na(s) & s < 1e-12 * pmax(1, abs(est))] <- NA_real_
  out$se       <- s
  out$z        <- est / s
  out$p_value  <- 2 * stats::pnorm(-abs(out$z))
  # Same edge as the note measurement_summary() prints (.print_boundary_note()).
  out$boundary <- out$parameter == "probability" &
    (est <= 1e-3 | est >= 1 - 1e-3)
  rownames(out) <- NULL
  if (any(out$boundary))
    message(sprintf(paste0(
      "%d probabilit%s sit%s at the boundary (column `boundary`); %s ",
      "standard error%s not interpretable. See ?measurement_estimates."),
      sum(out$boundary), if (sum(out$boundary) == 1L) "y" else "ies",
      if (sum(out$boundary) == 1L) "s" else "",
      if (sum(out$boundary) == 1L) "its" else "their",
      if (sum(out$boundary) == 1L) " is" else "s are"))
  attr(out, "vcov") <- V
  out
}

# The measurement parameters as reported, from a model state: every block's
# item parameters on their natural scale, then the class-size logits.
.measurement_reported <- function(fit, st, logits) {
  mm <- st$mm
  K  <- fit$n_components
  tab <- if (inherits(mm, c("nested", "blocks")))
    do.call(rbind, lapply(names(mm$models), function(nm)
      .measurement_rows(mm$models[[nm]], K, nm)))
  else .measurement_rows(mm, K, NA_character_)
  if (logits && K > 1L) {
    w <- pmax(st$weights, 1e-300)
    tab <- rbind(tab, data.frame(
      block = NA_character_, parameter = "class_logit",
      item = sprintf("log(p%d/p%d)", seq_len(K - 1L), K),
      category = NA_integer_, class = seq_len(K - 1L),
      estimate = log(w[-K] / w[K]), stringsAsFactors = FALSE))
  }
  tab
}

# One flat emission's rows, class varying fastest within each item, in the
# column layout measurement_summary() returns.
.measurement_rows <- function(emis, K, block) {
  p <- emis$parameters
  rows <- function(mat, parameter, item, category = NA_integer_) {
    J <- ncol(mat)
    data.frame(block = block, parameter = parameter,
               item = rep(item, each = K),
               category = rep(rep_len(category, J), each = K),
               class = rep(seq_len(K), times = J),
               estimate = as.vector(mat), stringsAsFactors = FALSE)
  }
  names_or <- function(mat, prefix)
    colnames(mat) %||% paste0(prefix, seq_len(ncol(mat)))
  out <- list()
  if (!is.null(p$pis)) {
    if (!is.null(emis$max_val)) {
      cats <- .multinoulli_cats(emis)
      base <- emis$item_names
      if (is.null(base) || length(base) != length(cats))
        base <- paste0("Poly_Item_", seq_along(cats))
      out$pis <- rows(p$pis[, .multinoulli_real_cols(emis), drop = FALSE],
                      "probability", rep(base, cats),
                      unlist(lapply(cats, seq_len)))
    } else {
      out$pis <- rows(p$pis, "probability", names_or(p$pis, "Item_"))
    }
  }
  if (!is.null(p$means)) {
    nm <- names_or(p$means, "Cont_")
    out$means <- rows(p$means, "mean", nm)
    if (!is.null(p$covariances))
      out$var <- rows(p$covariances, "variance", nm)
  }
  if (!is.null(p$rates))
    out$rates <- rows(p$rates, "rate", names_or(p$rates, "Count_"))
  do.call(rbind, unname(out))
}
