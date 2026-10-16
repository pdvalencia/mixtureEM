# ------------------------------------------------------------------------------
# Standard errors for a growth mixture model's parameters
# ------------------------------------------------------------------------------

#' Growth Mixture Parameter Estimates with Standard Errors
#'
#' @description
#' Returns every parameter of a [`fit_gmm()`] fit -- the growth-factor means
#' (intercepts, when there are `growth_predictors`), the growth-factor
#' regressions, the growth-factor variances and covariances, the residual
#' variances, the fitted trajectory and the class-size logits -- with a
#' standard error, a z statistic and a two-sided p-value. For a [`fit_lcga()`]
#' fit it returns the growth coefficients of each class (on the link scale:
#' logits for `"binomial"`, log rates for `"poisson"`), the residual variance
#' (`"gaussian"` only), the fitted trajectory on the response scale and the
#' class-size logits.
#'
#' @details
#' The standard errors are computed from the log-likelihood of the fitted model
#' itself, by differentiating it with respect to the parameters the fit
#' estimated, and are then carried to the reported scale by the delta method.
#' Three estimators are offered:
#'
#' * `"robust"` (the default) is the sandwich \eqn{H^{-1} B H^{-1}}, where
#'   \eqn{H} is the observed information (minus the matrix of second
#'   derivatives of the log-likelihood) and \eqn{B} the sum over cases of the
#'   outer product of each case's score (White, 1982). It stays a valid
#'   variance when the within-class normal model is only approximately right,
#'   which for a growth mixture is the usual situation: the classes are partly a
#'   device for describing a non-normal distribution of trajectories (Bauer &
#'   Curran, 2003).
#' * `"hessian"` is \eqn{H^{-1}}, the inverse observed information, the
#'   textbook maximum-likelihood variance when the model is correct.
#' * `"opg"` is \eqn{B^{-1}}, the inverse outer product of the case scores,
#'   which estimates the same quantity as `"hessian"` under a correct model
#'   from first derivatives alone.
#'
#' When the model is right the three agree in large samples; a marked gap
#' between `"robust"` and the other two is a sign of misspecification.
#' With sample `weights` every case's contribution is weighted, and with a
#' `cluster` variable the scores are summed within clusters before \eqn{B} is
#' formed.
#'
#' Constrained parameters are repeated across the classes they are shared by,
#' with the same standard error, as in [`measurement_summary()`]. The prior
#' that keeps the variances away from zero during estimation is not part of
#' the variance: these are the standard errors of the likelihood.
#'
#' A growth-factor covariance matrix outside the positive semi-definite cone
#' still gets standard errors (see `psi` in [`fit_gmm()`]); the point is a
#' maximum of the likelihood like any other. Where the information matrix is
#' singular -- a parameter the data do not identify -- the standard errors are
#' `NA` and a warning says so.
#'
#' The class-size logits are \eqn{\log(\pi_k / \pi_K)}. When the class sizes
#' are modelled by `predictors` they are not reported here; see [`coef()`].
#' Under `n_steps = 3` the standard errors here are those of the step-1 model,
#' which is where the growth parameters are estimated.
#'
#' @param fit A model returned by [`fit_gmm()`] or [`fit_lcga()`].
#' @param se Which standard-error estimator: `"robust"` (default), `"hessian"`
#'   or `"opg"`. See Details.
#'
#' @return A data frame with one row per parameter and class: `parameter`
#'   (`"growth_mean"`, `"growth_regression"`, `"growth_variance"`,
#'   `"growth_covariance"`, `"residual_variance"`, `"fitted"`,
#'   `"class_logit"`, or, for a continuous `outcome` fitted in one step,
#'   `"outcome_mean"` and `"outcome_variance"`), `item`, `class`, `estimate`, `se`, `z` and `p_value`.
#'   The covariance matrix of the reported parameters is attached as the
#'   attribute `"vcov"`.
#'
#' @references
#' White, H. (1982). Maximum likelihood estimation of misspecified models.
#' *Econometrica*, *50*(1), 1-25. \doi{10.2307/1912526}
#'
#' Louis, T. A. (1982). Finding the observed information matrix when using the
#' EM algorithm. *Journal of the Royal Statistical Society: Series B*, *44*(2),
#' 226-233. \doi{10.1111/j.2517-6161.1982.tb01203.x}
#'
#' Bauer, D. J., & Curran, P. J. (2003). Distributional assumptions of growth
#' mixture models: Implications for overextraction of latent trajectory
#' classes. *Psychological Methods*, *8*(3), 338-363.
#' \doi{10.1037/1082-989X.8.3.338}
#'
#' @seealso [`fit_gmm()`], [`fit_lcga()`], [`measurement_summary()`].
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' n <- 400
#' cls <- rbinom(n, 1, 0.5)
#' icept <- ifelse(cls == 1, 3, 1) + rnorm(n, 0, 1)
#' slope <- ifelse(cls == 1, 1.0, 0.2) + rnorm(n, 0, 0.4)
#' y <- outer(icept, rep(1, 4)) + outer(slope, 0:3) + matrix(rnorm(n * 4, 0, 0.7), n, 4)
#' fit <- fit_gmm(y, times = 4, n_classes = 2, n_init = 10, random_state = 1)
#' growth_estimates(fit)
#' }
#' @export
growth_estimates <- function(fit, se = c("robust", "hessian", "opg")) {
  se <- match.arg(se)
  if (!inherits(fit, c("gmm", "lcga")) ||
      !isTRUE(fit$growth$model %in% c("gmm", "lcga")))
    stop("`fit` must be a growth model from fit_gmm() or fit_lcga().",
         call. = FALSE)
  v <- .gmm_vcov(fit)
  V <- v[[se]]
  est <- v$table$estimate
  s   <- if (is.null(V)) rep(NA_real_, length(est)) else sqrt(pmax(diag(V), 0))
  # A parameter with no sampling variance (a fixed value) has no z.
  s[!is.na(s) & s < 1e-12 * pmax(1, abs(est))] <- NA_real_
  out <- v$table
  out$se      <- s
  out$z       <- est / s
  out$p_value <- 2 * stats::pnorm(-abs(out$z))
  rownames(out) <- NULL
  attr(out, "vcov") <- V
  out
}

# The growth parameters as reported, from a model state: the rows of
# measurement_summary.gmm()'s table, then the class-size logits, then a
# one-step continuous outcome's means and variances.
.gmm_reported <- function(fit, st, logits) {
  rbind(.gmm_reported_growth(fit, st, logits), .gmm_outcome_rows(fit, st))
}

# A shared variance is repeated across the classes, as any shared parameter is.
.gmm_outcome_rows <- function(fit, st) {
  sm <- st$sm
  if (is.null(sm) || !identical(as.integer(fit$n_steps), 1L) ||
      !inherits(sm, "distal_continuous")) return(NULL)
  K  <- fit$n_components
  nm <- colnames(fit$Y)[1L] %||% "outcome"
  data.frame(parameter = rep(c("outcome_mean", "outcome_variance"), each = K),
             item = nm, class = rep(seq_len(K), 2L),
             estimate = c(sm$parameters$means[, 1L], sm$parameters$covariances[, 1L]),
             stringsAsFactors = FALSE)
}

.gmm_reported_growth <- function(fit, st, logits) {
  g  <- fit$growth
  mm <- st$mm
  K  <- fit$n_components
  if (identical(g$model, "lcga")) {
    g$coefficients <- mm$parameters$coefs
    g$fitted <- .lcga_fitted(mm)
    if (identical(g$family, "gaussian"))
      g$residual_variance <- mm$parameters$dispersion
    return(.growth_reported_table(g, K, st, logits))
  }
  g$means <- mm$parameters$alpha
  if (!is.null(g$coefficients))
    g$coefficients <- .gmm_name_coefficients(mm, g$covariate_names)
  g$psi <- mm$parameters$psi[seq_len(if (isTRUE(mm$psi_equal)) 1L else K)]
  g$residual_variance <- mm$parameters$theta
  g$fitted <- .sn_fitted(mm)
  .growth_reported_table(g, K, st, logits)
}

.growth_reported_table <- function(g, K, st, logits) {
  blocks <- .growth_blocks(list(growth = g, n_components = K))
  tab <- do.call(rbind, lapply(blocks, function(b)
    .growth_rows(b$mat, b$parameter, K)))
  tab <- tab[, c("parameter", "item", "class", "estimate")]
  if (logits && K > 1L) {
    w <- pmax(st$weights, 1e-300)
    tab <- rbind(tab, data.frame(parameter = "class_logit",
                                 item = sprintf("log(p%d/p%d)", seq_len(K - 1L), K),
                                 class = seq_len(K - 1L),
                                 estimate = log(w[-K] / w[K]),
                                 stringsAsFactors = FALSE))
  }
  tab
}

# The three covariance matrices of the reported parameters.
#
# The parameters are the ones the finish works on (.qn_pack_sn(), R/qn_pack.R):
# Psi as its elements, so a matrix outside the cone is an ordinary point, and
# Theta as log standard deviations. The scores are central differences of each
# case's log-likelihood (.qn_fd_scores(), one-sided where the other side has no
# likelihood), the observed information a central difference of their sum, and
# the delta method carries all three to the reported scale.
.gmm_vcov <- function(fit) {
  one_step <- !is.null(fit$sm) && identical(as.integer(fit$n_steps), 1L)
  st <- fit
  st$frozen <- NULL
  if (!one_step) st$sm <- NULL
  Y  <- if (one_step) fit$Y else NULL
  pr <- .qn_mixture_problem(st, fit$data, Y)
  if (is.null(pr))
    stop("This growth model has no parameter packing to differentiate.",
         call. = FALSE)
  par <- pr$par
  p   <- length(par)
  n   <- nrow(fit$data)
  w   <- if (is.null(pr$w)) rep(1, n) else pr$w
  idx <- seq_len(p)

  S <- .qn_fd_scores(pr$case_ll, par, idx)
  grad <- function(x) colSums(w * .qn_fd_scores(pr$case_ll, x, idx))
  h <- 1e-4 * pmax(1, abs(par))
  H <- vapply(idx, function(j) {
    a <- par; a[j] <- a[j] + h[j]
    b <- par; b[j] <- b[j] - h[j]
    (grad(a) - grad(b)) / (2 * h[j])
  }, numeric(p))
  info <- -(H + t(H)) / 2

  cl <- fit$cluster
  if (is.null(cl) || length(cl) != n) cl <- seq_len(n)
  Sc   <- rowsum(w * S, cl, reorder = FALSE)
  meat <- crossprod(Sc)
  opg  <- crossprod(sqrt(w) * S)

  inv <- function(M) {
    e <- eigen(M, symmetric = TRUE)
    if (min(e$values) <= 1e-8 * max(abs(e$values))) return(NULL)
    e$vectors %*% (t(e$vectors) / e$values)
  }
  Hi <- inv(info)
  Bi <- inv(opg)
  if (is.null(Hi) || is.null(Bi))
    warning(paste0(
      "The information matrix of this growth model is singular, so some ",
      "parameter is not identified by the data and its standard error is ",
      "reported as NA. A class with almost no cases, or a variance at the ",
      "edge of what the model allows, is the usual reason."), call. = FALSE)

  logits <- !is.null(pr$lay$idx$weights)
  rep_at <- function(x) .gmm_reported(fit, pr$unpack(x), logits)$estimate
  tab <- .gmm_reported(fit, pr$unpack(par), logits)
  J <- vapply(idx, function(j) {
    a <- par; a[j] <- a[j] + h[j]
    b <- par; b[j] <- b[j] - h[j]
    (rep_at(a) - rep_at(b)) / (2 * h[j])
  }, numeric(nrow(tab)))
  J <- matrix(J, nrow(tab))
  map <- function(V) if (is.null(V)) NULL else J %*% V %*% t(J)

  list(table   = tab,
       hessian = map(Hi),
       opg     = map(Bi),
       robust  = if (is.null(Hi)) NULL else map(Hi %*% meat %*% Hi))
}
