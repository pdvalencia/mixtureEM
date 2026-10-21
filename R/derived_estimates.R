# ==============================================================================
# Derived parameters with delta-method standard errors
# ==============================================================================

#' Derived Parameters with Standard Errors and Intervals
#'
#' @description
#' Estimates functions of a fitted model's parameters -- a sum of two
#' coefficients, a difference, a probability built from several of them --
#' each with a standard error, a z test against zero and a confidence
#' interval. The parameters are named by the labels [`wald_test()`] uses.
#'
#' @details
#' Each derived parameter is given as one named argument in `...`, either
#' \itemize{
#'   \item a named numeric vector of weights, `c(a = 1, b = 1, c = -1)`,
#'     meaning \eqn{a + b - c}; its variance is \eqn{c^\top V c}, exact; or
#'   \item a function of the named vector of estimates, returning one number,
#'     such as `function(p) plogis(p[["a"]] + p[["b"]])`. Its standard error
#'     is the delta method's, \eqn{\sqrt{g^\top V g}} with \eqn{g} the
#'     gradient of the function at the estimates, taken by central
#'     differences.
#' }
#' \eqn{V} is the covariance matrix the fit reports (for a growth model, the
#' one `se` selects). The interval is the estimate plus or minus the normal
#' quantile times the standard error, on the scale of the derived parameter.
#' For a probability near 0 or 1 an interval built on the logit scale and
#' transformed back stays inside the unit interval; derive the logit, and
#' transform its bounds with `plogis()`.
#'
#' For [`fit_joint_classes()`] the labels are those of
#' [`joint_coefficients()`]: `"CS#j Intercept"`, `"CS#j ON CA#k"`, `"CS#j ON
#' x"` and `"CS#j ON x | CA#k"`. A covariate `x` acting on `cs` by `ca` class
#' has slope `"CS#j ON x"` in the last `ca` class and that plus `"CS#j ON x |
#' CA#k"` in class k, so the association among cases with `x = 1` is
#' `c("CS#j ON CA#k" = 1, "CS#j ON x | CA#k" = 1)`.
#' [`joint_probabilities()`] builds the conditional class probabilities this
#' way.
#'
#' @param fit A model [`wald_test()`] accepts.
#' @param ... Named derived parameters: weight vectors or functions (see
#'   Details).
#' @param level Confidence level of the intervals.
#' @param se For growth models, which covariance estimator to use, as in
#'   [`wald_test()`]. Ignored otherwise.
#'
#' @return A data frame with one row per derived parameter: `parameter`,
#'   `estimate`, `se`, `z`, `p`, `lower`, `upper`.
#'
#' @references
#' Oehlert, G. W. (1992). A note on the delta method. *The American
#' Statistician*, *46*(1), 27-29.
#'
#' @seealso [`wald_test()`] to test several restrictions at once,
#'   [`joint_probabilities()`].
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' n  <- 800
#' w  <- rbinom(n, 1, 0.5)
#' ca <- sample(1:2, n, TRUE)
#' cs <- ifelse(runif(n) < plogis(-1 + 2 * (ca == 1) + w * (ca == 1)), 1L, 2L)
#' A  <- sapply(1:4, function(j) rbinom(n, 1, ifelse(ca == 1, 0.85, 0.15)))
#' S  <- sapply(1:4, function(j) rbinom(n, 1, ifelse(cs == 1, 0.85, 0.15)))
#' fit <- fit_joint_classes(A, S, predictors_cs = data.frame(w = w),
#'                          cs_effects = "by_ca", n_init = 5, random_state = 1)
#' # The association among cases with w = 1, and the probability of cs class 1
#' # for ca class 1 at w = 1.
#' derived_estimates(fit,
#'   assoc_w1 = c("CS#1 ON CA#1" = 1, "CS#1 ON w | CA#1" = 1),
#'   p_w1 = function(p) plogis(p[["CS#1 Intercept"]] + p[["CS#1 ON CA#1"]] +
#'                             p[["CS#1 ON w"]] + p[["CS#1 ON w | CA#1"]]))
#' }
#' @export
derived_estimates <- function(fit, ..., level = 0.95,
                              se = c("robust", "hessian", "opg")) {
  se   <- match.arg(se)
  defs <- list(...)
  if (!length(defs) || is.null(names(defs)) || any(!nzchar(names(defs))))
    stop("Give each derived parameter as a named argument, e.g. ",
         "derived_estimates(fit, diff = c(a = 1, b = -1)).", call. = FALSE)
  if (!is.numeric(level) || length(level) != 1L || level <= 0 || level >= 1)
    stop("`level` must be a single number between 0 and 1.", call. = FALSE)
  pe  <- .wald_parameters(fit, se)
  th  <- pe$estimate
  V   <- pe$vcov
  lab <- names(th)

  rows <- lapply(names(defs), function(nm) {
    d <- defs[[nm]]
    if (is.function(d)) {
      f0 <- d(th)
      if (!is.numeric(f0) || length(f0) != 1L || !is.finite(f0))
        stop(sprintf("`%s` must return one finite number.", nm), call. = FALSE)
      g <- vapply(seq_along(th), function(i) {
        h <- 1e-5 * max(1, abs(th[[i]]))
        up <- th; up[i] <- up[i] + h
        dn <- th; dn[i] <- dn[i] - h
        (d(up) - d(dn)) / (2 * h)
      }, numeric(1))
      est <- f0
    } else {
      if (!is.numeric(d) || is.null(names(d)) || any(!nzchar(names(d))))
        stop(sprintf(paste0("`%s` must be a named numeric vector of weights or ",
                            "a function of the named estimates."), nm),
             call. = FALSE)
      g <- .wald_constraint_matrix(list(d), lab)[1L, ]
      est <- sum(g * th)
    }
    used <- g != 0
    if (anyNA(V[used, used]))
      stop(sprintf("`%s` uses a parameter with no standard error.", nm),
           call. = FALSE)
    s <- sqrt(max(drop(crossprod(g[used], V[used, used, drop = FALSE] %*% g[used])), 0))
    data.frame(parameter = nm, estimate = est, se = s, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  q <- stats::qnorm(1 - (1 - level) / 2)
  out$z <- out$estimate / out$se
  out$p <- 2 * stats::pnorm(-abs(out$z))
  out$lower <- out$estimate - q * out$se
  out$upper <- out$estimate + q * out$se
  attr(out, "method") <- pe$method
  out
}

#' Conditional Class Probabilities of a Two-Latent-Variable Model
#'
#' @description
#' P(cs = j | ca = k) at given covariate values, from a
#' [`fit_joint_classes()`] fit, with delta-method standard errors and
#' intervals ([`derived_estimates()`]). The fit's own `conditional` table
#' averages over the sample's covariates; this one fixes them, which is how a
#' moderated association is described: the probabilities among men and among
#' women, say.
#'
#' @param fit A [`fit_joint_classes()`] fit.
#' @param at A named list or vector of covariate values for the covariates on
#'   `cs`. Covariates not named are held at 0. `NULL` (the default) holds all
#'   at 0.
#' @param level Confidence level of the intervals.
#'
#' @return A data frame with columns `ca`, `cs`, `estimate`, `se`, `lower`,
#'   `upper`, one row per pair of classes.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' n  <- 800
#' w  <- rbinom(n, 1, 0.5)
#' ca <- sample(1:2, n, TRUE)
#' cs <- ifelse(runif(n) < plogis(-1 + 2 * (ca == 1) + w * (ca == 1)), 1L, 2L)
#' A  <- sapply(1:4, function(j) rbinom(n, 1, ifelse(ca == 1, 0.85, 0.15)))
#' S  <- sapply(1:4, function(j) rbinom(n, 1, ifelse(cs == 1, 0.85, 0.15)))
#' fit <- fit_joint_classes(A, S, predictors_cs = data.frame(w = w),
#'                          cs_effects = "by_ca", n_init = 5, random_state = 1)
#' joint_probabilities(fit, at = list(w = 1))
#' }
#' @export
joint_probabilities <- function(fit, at = NULL, level = 0.95) {
  stopifnot(inherits(fit, "joint_classes"))
  K1 <- fit$n_ca; K2 <- fit$n_cs
  a  <- fit$association
  covs <- setdiff(unique(a$term), c("Intercept", sprintf("CA#%d", seq_len(K1))))
  at <- unlist(at)
  if (length(at) && (is.null(names(at)) || any(!names(at) %in% covs)))
    stop("`at` must name covariates of `cs`: ",
         if (length(covs)) paste(covs, collapse = ", ") else "this fit has none",
         ".", call. = FALSE)
  lab <- .joint_labels(a)
  get <- function(p, nm) if (nm %in% lab) p[[nm]] else 0
  # The logit of cs class j (< K2) against the last, for ca class k.
  logit <- function(p, k, j) {
    v <- get(p, sprintf("CS#%d Intercept", j))
    if (k < K1) v <- v + get(p, sprintf("CS#%d ON CA#%d", j, k))
    for (x in names(at))
      v <- v + at[[x]] * (get(p, sprintf("CS#%d ON %s", j, x)) +
                          if (k < K1) get(p, sprintf("CS#%d ON %s | CA#%d", j, x, k))
                          else 0)
    v
  }
  prob <- function(k, j) function(p) {
    e <- exp(c(vapply(seq_len(K2 - 1L), function(i) logit(p, k, i), numeric(1)), 0))
    e[j] / sum(e)
  }
  grid <- expand.grid(cs = seq_len(K2), ca = seq_len(K1))
  defs <- stats::setNames(
    lapply(seq_len(nrow(grid)), function(r) prob(grid$ca[r], grid$cs[r])),
    sprintf("P(CS#%d | CA#%d)", grid$cs, grid$ca))
  d <- do.call(derived_estimates, c(list(fit), defs, list(level = level)))
  data.frame(ca = sprintf("CA#%d", grid$ca), cs = sprintf("CS#%d", grid$cs),
             estimate = d$estimate, se = d$se, lower = d$lower,
             upper = d$upper, stringsAsFactors = FALSE)
}
