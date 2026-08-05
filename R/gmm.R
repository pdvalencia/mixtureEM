# ==============================================================================
# Growth Mixture Modeling (GMM)
# ==============================================================================
#
# The user-facing layer over the structured multivariate-normal emission in
# R/structured_normal.R. A growth mixture model is a latent class growth model
# plus within-class random effects: classes still differ in the shape of their
# average trajectory, but cases are no longer assumed to sit exactly on their
# class's curve. That is the whole difference, and it is usually a large one -
# in benchmark runs the two random effects are worth ~474 log-likelihood
# units for two parameters.
#
# Everything else the package offers applies unchanged, because the mixture
# E-step never learns that the emission has an internal structure: model
# selection, the bootstrap likelihood-ratio test, predictors of class
# membership, distal outcomes, survey designs, and FIML for missing occasions.

# Whether the covariance structure asked for can be told apart from the data.
#
# A class supplies T(T+1)/2 distinct covariance elements. The structured model
# spends q(q+1)/2 of them on the growth-factor covariance and one per free
# residual variance; asking for more than the data contain gives an unidentified
# model that will still "converge", to an arbitrary point on a flat ridge. The
# check is cheap and the failure it prevents is silent, which is the argument
# for making it an error rather than a warning.
.gmm_check_identified <- function(times, q, theta_shared_occasions) {
  available <- times * (times + 1L) / 2L
  spent <- q * (q + 1L) / 2L + (if (theta_shared_occasions) 1L else times)
  if (spent > available)
    stop(sprintf(
      paste0("The covariance structure asked for is not identified: %d ",
             "parameters (%d for the growth factors, %d residual variance%s) ",
             "from the %d distinct covariance elements %d occasions supply. ",
             "Use fewer random effects, or residual = \"constant\"."),
      spent, q * (q + 1L) / 2L,
      if (theta_shared_occasions) 1L else times,
      if (theta_shared_occasions) "" else "s",
      available, times), call. = FALSE)
  invisible(NULL)
}

# Covariates on the growth factors, checked and put on the emission's terms.
#
# Two things are settled here rather than deeper down. First the covariates are
# required to be complete: a case missing x has no model-implied mean at all, so
# unlike a missing occasion -- which the multivariate normal integrates out
# exactly -- there is nothing to condition on. Listwise deletion is the usual
# treatment elsewhere; this package refuses instead, because deleting them here would leave
# `weights`, `strata` and `cluster` (which arrive through `...`) pointing at the
# wrong rows, and a silent misalignment of a survey design is a far worse
# outcome than an error the user resolves with complete.cases().
#
# Second, the rows are put in the frame the emission will see. The emission
# holds the covariate matrix in its own state, so it must line up with the X
# that fit_mixture_internal() actually estimates on -- and that function deletes
# cases with no observed occasion before building the emission. Deleting the
# same rows here keeps the two aligned; the deletion is reported there, as it
# already is for every other model.
.gmm_growth_covariates <- function(value, expr, X) {
  if (is.null(value)) return(NULL)

  gx <- prepare_covariates(.as_named_covariates(value, expr, "growth predictor"))
  if (!is.matrix(gx)) gx <- as.matrix(gx)

  if (nrow(gx) != nrow(X))
    stop(sprintf(
      "`growth_predictors` must have one row per case (%d given, %d needed).",
      nrow(gx), nrow(X)), call. = FALSE)

  if (anyNA(gx)) {
    bad <- which(!stats::complete.cases(gx))
    stop(sprintf(
      paste0("`growth_predictors` must be complete, but %d case%s missing on ",
             "at least one: %s. A covariate on the growth factors is what the ",
             "class mean is built from, so a case missing it has no fitted ",
             "trajectory; missing occasions are handled by FIML, missing ",
             "covariates are not. Drop or impute those cases first."),
      length(bad), if (length(bad) == 1L) " is" else "s are",
      .abbreviate_indices(bad)), call. = FALSE)
  }

  drop <- .empty_rows(X)
  if (length(drop)) gx <- gx[-drop, , drop = FALSE]
  gx
}

# Growth factor by covariate, one matrix per class, for the fitted object.
.gmm_name_coefficients <- function(mm, cov_names) {
  out <- lapply(mm$parameters$gamma, function(G) {
    dimnames(G) <- list(colnames(mm$design), cov_names)
    G
  })
  if (isTRUE(mm$gamma_equal)) out[1L] else out
}

#' Growth Mixture Modeling
#'
#' @description
#' Fits class-specific growth curves to a continuous outcome measured
#' repeatedly, allowing cases to vary about their own class's trajectory through
#' random growth factors. This is the growth mixture model of Muthen & Shedden
#' (1999) and, with `random_effects = "none"`, latent
#' class growth analysis for a continuous outcome.
#'
#' The model is
#' \deqn{y_i \mid c_i = k \sim N(\Lambda (\alpha_k + \Gamma_k x_i),\; \Lambda_r \Psi_k \Lambda_r' + \Theta_k),}
#' where \eqn{\Lambda} is the polynomial design in time, \eqn{\alpha_k} the
#' growth-factor means of class \eqn{k}, \eqn{\Psi_k} the covariance of the
#' growth factors that are allowed to vary within a class, and \eqn{\Theta_k}
#' the diagonal matrix of residual variances. Because the outcome is continuous
#' the random effects integrate out in closed form, so no numerical integration
#' is involved and the estimator is the ordinary EM algorithm the rest of the
#' package uses.
#'
#' \eqn{\Gamma_k} is present only when `growth_predictors` are supplied, and is
#' the within-class growth-factor regression: covariates shift a case's own growth factors *within* its
#' class, answering "who, in this trajectory group, starts higher or grows
#' faster?". It is a different question from `predictors`, which is `c ON x` and
#' asks who is *in* the group; the two can be asked together, and are then
#' estimated in one pass with `n_steps = 1`. With covariates
#' present \eqn{\alpha_k} is the growth-factor *intercept* rather than its mean,
#' and the reported trajectory is evaluated at the sample mean of the
#' covariates.
#'
#' Contrast [`fit_lcga()`], which fixes the growth-factor variances at zero, so
#' every case in a class sits on that class's curve up to occasion-level noise.
#' LCGA is the more parsimonious model and, in criminology especially, often the
#' preferred one; GMM is the more realistic and is what developmental psychology
#' and epidemiology usually publish. The two are nested, so the choice can be
#' made on BIC.
#'
#' @param indicator The repeated outcome. Either a wide matrix or data frame
#'   with one column per occasion, a three-dimensional array with dimensions n
#'   by 1 by times, or a long data frame together with `id` and `time`. Exactly
#'   one variable may be modelled.
#' @param n_classes Integer. Number of trajectory classes.
#' @param times Integer. Number of occasions. Required for wide input; inferred
#'   otherwise.
#' @param degree Degree of the polynomial in time: `1` for a linear trajectory
#'   (intercept and slope), `2` for a quadratic, and so on.
#' @param random_effects Which growth factors vary within a class:
#'   * `"intercept_slope"` (default) — both the intercept and the linear slope,
#'     the standard growth mixture model. Cases in a class differ both in where
#'     they start and in how fast they change.
#'   * `"intercept"` — the intercept only. Cases in a class start at different
#'     levels but change in parallel.
#'   * `"none"` — no random effects, which is latent class growth analysis; the
#'     same model as `fit_lcga(family = "gaussian")`, but with the residual
#'     variances free over occasions by default.
#'   * `"all"` — every growth factor, including a quadratic or higher term.
#' @param psi Whether the growth-factor covariance \eqn{\Psi} is held `"equal"`
#'   across classes (the default) or estimated `"free"`ly in each.
#'   Equality is the better default in practice as well as by convention: it is
#'   what BIC typically selects, and a class-specific \eqn{\Psi} is
#'   where growth mixture models most often produce a negative variance.
#' @param residual Whether the residual variances are `"occasion"`-specific (the
#'   default) or `"constant"` across occasions.
#' @param residual_equal Logical. Hold the residual variances equal across
#'   classes (the conventional default).
#' @param time_scores Numeric values of time used in the polynomial, one per
#'   occasion. Defaults to `0, 1, ..., times - 1`, which makes the intercept
#'   growth factor the level at the first occasion. Supply the actual
#'   measurement times when the occasions are unequally spaced.
#' @param layout For wide input with a three-dimensional array or several
#'   columns per occasion, whether columns run `"time_major"` or `"item_major"`.
#' @param id,time For long input, the case and occasion identifiers, given
#'   either as column names or as vectors.
#' @param item For long input, the column holding the outcome.
#' @param time_labels Optional display labels for the occasions.
#' @param predictors Optional predictors of class membership (`c ON x`), passed
#'   to [`fit_mixture()`]'s three-step machinery.
#' @param growth_predictors Optional covariates regressing the growth factors
#'   themselves (`i s ON x`): a vector, matrix or data frame with one row per
#'   case, factors dummy-coded. These enter the measurement model rather than
#'   the structural one, so they are estimated jointly with the trajectories in
#'   a single pass regardless of `n_steps`. They must be complete: a case
#'   missing a covariate has no model-implied mean, and unlike a missing
#'   occasion there is nothing here to integrate it out, so remove or impute
#'   those cases first.
#' @param growth_predictors_equal Logical. Hold the growth-factor regressions
#'   equal across classes (the conventional default). Set `FALSE` to let each class have its
#'   own, which asks whether a covariate matters differently in different
#'   trajectory groups.
#' @param ... Further arguments passed to [`fit_mixture()`], such as `outcome`
#'   and its companions, `n_init`, `random_state`, `weights`, `strata` or
#'   `cluster`.
#'
#' @return An object of class `c("gmm", "mixture_model")`. In addition to the
#'   usual fields it carries `$growth`, holding the design matrix, time scores,
#'   the per-class growth-factor means, the fitted class trajectories, the
#'   growth-factor covariance of each class, the residual variances and, when
#'   `growth_predictors` were given, `$growth$coefficients`: one
#'   growth-factor-by-covariate matrix per class.
#'
#' @references
#' Muthen, B., & Shedden, K. (1999). Finite mixture modeling with mixture
#' outcomes using the EM algorithm. *Biometrics*, *55*(2), 463-469.
#'
#' @seealso [`fit_lcga()`] for the no-random-effects version and
#'   [`fit_rmlca()`] for trajectory classes with no growth curve at all.
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
#' fit
#' }
#' @export
fit_gmm <- function(indicator,
                    n_classes = 2,
                    times = NULL,
                    degree = 1,
                    random_effects = c("intercept_slope", "intercept", "none",
                                       "all"),
                    psi = c("equal", "free"),
                    residual = c("occasion", "constant"),
                    residual_equal = TRUE,
                    time_scores = NULL,
                    layout = c("time_major", "item_major"),
                    id = NULL, time = NULL, item = NULL,
                    time_labels = NULL,
                    predictors = NULL,
                    growth_predictors = NULL,
                    growth_predictors_equal = TRUE,
                    ...) {

  growth_predictors_expr <- substitute(growth_predictors)
  random_effects <- match.arg(random_effects)
  psi            <- match.arg(psi)
  residual       <- match.arg(residual)
  layout         <- match.arg(layout)

  prep <- .prepare_longitudinal(indicator, times = times, items = item,
                                layout = layout, id = id, time = time,
                                time_labels = time_labels)

  if (prep$n_items != 1L)
    stop(sprintf(
      paste0("A growth mixture model has one repeated outcome, but %d were ",
             "found at each of the %d occasions. Select a single variable."),
      prep$n_items, prep$n_times), call. = FALSE)

  Tn <- prep$n_times
  if (is.null(time_scores)) time_scores <- seq.int(0L, Tn - 1L)
  if (length(time_scores) != Tn)
    stop(sprintf("`time_scores` must have one value per occasion (%d given, %d needed).",
                 length(time_scores), Tn), call. = FALSE)

  degree <- as.integer(degree)
  if (is.na(degree) || degree < 0L)
    stop("`degree` must be a non-negative integer.", call. = FALSE)
  if (degree + 1L >= Tn)
    stop(sprintf(
      paste0("A degree-%d trajectory needs at least %d occasions but %d were ",
             "given. With degree + 1 = %d coefficients per class the curve ",
             "would reproduce the occasion means exactly."),
      degree, degree + 2L, Tn, degree + 1L), call. = FALSE)

  .lcga_check_outcome(prep$X, "gaussian")

  design <- .lcga_design(time_scores, degree)
  q <- length(structured_normal_model(1L, design,
                                      random_effects = random_effects)$r_cols)
  .gmm_check_identified(Tn, q, identical(residual, "constant"))

  if (!is.null(predictors))
    predictors <- .as_named_covariates(predictors, substitute(predictors),
                                       "predictor")

  gx <- .gmm_growth_covariates(growth_predictors, growth_predictors_expr,
                               prep$X)

  fit <- fit_mixture(
    indicators              = prep$X,
    n_classes               = n_classes,
    measurement             = "structured_normal",
    predictors              = predictors,
    design                  = design,
    random_effects          = random_effects,
    psi                     = psi,
    residual                = residual,
    residual_equal          = residual_equal,
    growth_covariates       = gx,
    growth_covariates_equal = growth_predictors_equal,
    ...
  )

  mm <- fit$mm
  fit$growth <- list(
    model             = "gmm",
    family            = "gaussian",
    degree            = degree,
    random_effects    = random_effects,
    time_scores       = time_scores,
    time_labels       = prep$time_labels,
    design            = design,
    factor_names      = colnames(design),
    means             = mm$parameters$alpha,
    coefficients      = if (is.null(gx)) NULL else
      .gmm_name_coefficients(mm, colnames(gx)),
    covariate_names   = colnames(gx),
    covariate_equal   = isTRUE(mm$gamma_equal),
    covariate_means   = if (is.null(gx)) NULL else colMeans(mm$xmat),
    fitted            = .sn_fitted(mm),
    psi               = mm$parameters$psi[seq_len(if (mm$psi_equal) 1L else
                                                    n_classes)],
    psi_equal         = mm$psi_equal,
    residual_variance = mm$parameters$theta,
    residual_equal    = mm$theta_equal,
    residual_constant = mm$theta_shared_occasions,
    random_names      = colnames(design)[mm$r_cols],
    wave_missing      = prep$wave_missing
  )
  fit$growth$boundary <- .gmm_boundary(mm)
  if (length(fit$growth$boundary))
    warning(sprintf(
      paste0("The solution is on the boundary of the parameter space: %s. ",
             "Growth mixture models reach this most often when the ",
             "growth-factor covariance or the residual variances are free ",
             "across classes and the data do not support that many; the ",
             "estimates at a boundary are not interpretable as variances and ",
             "their standard errors are not valid. Consider psi = \"equal\", ",
             "residual_equal = TRUE, or fewer classes."),
      paste(fit$growth$boundary, collapse = "; ")), call. = FALSE)

  class(fit) <- c("gmm", class(fit))
  fit
}

# Variance parameters sitting at the floor the M-step imposes.
#
# A growth mixture model whose covariance structure the data do not support does
# not fail; it produces a variance that wants to be negative. Some programs
# let it go negative and print it. This
# package floors instead, which keeps the likelihood a likelihood, but a floored
# estimate is the same diagnosis and would otherwise be reported as though it
# were an ordinary small variance. Hence the warning: an inadmissible solution
# should be visible without the user checking for it.
.gmm_boundary <- function(mm) {
  msg <- character(0)

  at_floor <- which(mm$parameters$theta <= .sn_theta_floor * 10,
                    arr.ind = TRUE)
  if (nrow(at_floor))
    msg <- c(msg, sprintf(
      "residual variance at zero for %s",
      paste(sprintf("class %d, occasion %d", at_floor[, 1], at_floor[, 2]),
            collapse = ", ")))

  if (length(mm$r_cols)) {
    for (k in seq_along(mm$parameters$psi)) {
      P <- mm$parameters$psi[[k]]
      if (!length(P)) next
      if (min(eigen(P, symmetric = TRUE, only.values = TRUE)$values) <= 1e-7) {
        msg <- c(msg, sprintf(
          "growth-factor covariance singular in class %d (no within-class variation left in the growth factors, which is the latent class growth model fit_lcga() estimates directly)",
          k))
        if (mm$psi_equal) break
      }
    }
  }
  msg
}

# ------------------------------------------------------------------------------
# Observed data on a trajectory plot
# ------------------------------------------------------------------------------

# The signature LCGA/GMM figure is the estimated curves *with the data behind
# them*: without it the reader has no way to tell a class that describes a real
# subgroup from one that has split a single cloud in half. Shared by
# plot.lcga() and plot.gmm() because the two figures differ only in what the
# curve is fitted to.
#
#   "means" draws each class's observed mean at each occasion, over the cases
#           modally assigned to it. The comparison it invites is the right one
#           for a trajectory model: does the fitted curve pass through the
#           observed means of the people the model says belong to that class?
#   "cases" draws individual observed trajectories, thinly and translucently.
#           Honest about within-class spread, which is exactly what a GMM
#           models and an LCGA assumes away, and unreadable above a few hundred
#           cases, hence not the default.
.growth_observed <- function(x, kind) {
  if (identical(kind, "none")) return(NULL)
  Y <- x$data
  if (is.null(Y) || !is.matrix(Y) || ncol(Y) != length(x$growth$time_scores))
    return(NULL)
  list(Y = Y, modal = max.col(exp(x$log_resp)), kind = kind)
}

.growth_draw_observed <- function(obs, time_scores, cols) {
  if (is.null(obs)) return(invisible(NULL))
  if (identical(obs$kind, "cases")) {
    for (k in sort(unique(obs$modal))) {
      Yk <- obs$Y[obs$modal == k, , drop = FALSE]
      graphics::matlines(time_scores, t(Yk), lty = 1, lwd = 0.5,
                         col = grDevices::adjustcolor(cols[k], alpha.f = 0.12))
    }
  } else {
    for (k in sort(unique(obs$modal))) {
      Yk <- obs$Y[obs$modal == k, , drop = FALSE]
      graphics::lines(time_scores, colMeans(Yk, na.rm = TRUE), lty = 3, lwd = 2,
                      col = cols[k])
      graphics::points(time_scores, colMeans(Yk, na.rm = TRUE), pch = 1, cex = 1.2,
                       col = cols[k])
    }
  }
  invisible(NULL)
}

# Range that has to fit in the panel: the fitted curves, plus whatever of the
# observed data is being drawn on top of them.
.growth_observed_range <- function(obs) {
  if (is.null(obs)) return(NULL)
  if (identical(obs$kind, "cases")) return(range(obs$Y, na.rm = TRUE))
  means <- vapply(sort(unique(obs$modal)),
                  function(k) range(colMeans(obs$Y[obs$modal == k, , drop = FALSE],
                                             na.rm = TRUE)),
                  numeric(2))
  range(means, na.rm = TRUE)
}

#' Trajectory Plot for a Growth Mixture Model
#'
#' @description
#' Draws the estimated mean trajectory of each class, with the observed data of
#' the cases assigned to it behind the curves.
#'
#' @param x An object returned by [`fit_gmm()`].
#' @param observed What of the observed data to draw: `"means"` (the observed
#'   mean at each occasion among cases modally assigned to the class, dotted),
#'   `"cases"` (individual trajectories, translucent) or `"none"`.
#' @param main Plot title.
#' @param class_labels Optional class labels for the legend.
#' @param colors Optional colour vector, recycled across classes.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @importFrom graphics matplot matlines lines points axis legend
#' @importFrom grDevices adjustcolor
#' @export
plot.gmm <- function(x, observed = c("means", "cases", "none"), main = NULL,
                     class_labels = NULL, colors = NULL, ...) {
  observed <- match.arg(observed)
  g  <- x$growth
  K  <- x$n_components
  ft <- g$fitted

  cols   <- if (is.null(colors)) rep(.okabe_ito, length.out = K) else
    rep(colors, length.out = K)
  shapes <- rep(15:20, length.out = K)
  base   <- if (is.null(class_labels)) paste("Class", seq_len(K)) else class_labels
  labels <- .class_plot_labels(base, x$weights)

  obs <- .growth_observed(x, observed)
  rng <- range(c(ft, .growth_observed_range(obs)), na.rm = TRUE)
  pad <- if (diff(rng) > 0) 0.05 * diff(rng) else max(abs(rng[1]), 1) * 0.05
  ylim <- rng + c(-pad, pad)

  matplot(g$time_scores, t(ft), type = "n", ylim = ylim, xaxt = "n", las = 1,
          bty = "l", xlab = "Occasion", ylab = "Mean",
          main = main %||% "Estimated class trajectories")
  .growth_draw_observed(obs, g$time_scores, cols)
  matlines(g$time_scores, t(ft), type = "b", pch = shapes, lty = 1, lwd = 2.5,
           col = cols)
  axis(1, at = g$time_scores, labels = g$time_labels)
  legend("topleft", legend = labels, col = cols, pch = shapes, lty = 1,
         lwd = 2, bty = "n", cex = 0.85)

  invisible(x)
}

#' Print a Fitted Growth Mixture Model
#'
#' @param x An object returned by [`fit_gmm()`].
#' @param ... Passed to the next method.
#' @return `x`, invisibly.
#' @export
print.gmm <- function(x, ...) {
  g <- x$growth
  K <- x$n_components
  shape <- c("intercept only", "linear", "quadratic", "cubic")[
    min(g$degree + 1L, 4L)]

  cat("\n")
  cat("=========================================================\n")
  cat("             GROWTH MIXTURE MODEL\n")
  cat("=========================================================\n")
  cat(sprintf("Occasions          : %d (time scores %s)\n",
              length(g$time_scores),
              paste(format(g$time_scores, trim = TRUE), collapse = ", ")))
  cat(sprintf("Trajectory         : %s\n", shape))
  cat(sprintf("Random effects     : %s\n",
              if (!length(g$random_names)) "none (latent class growth analysis)"
              else paste(g$random_names, collapse = ", ")))
  if (any(g$wave_missing))
    cat(sprintf("Wave attrition     : %d case-occasions unobserved\n",
                sum(g$wave_missing)))

  cat(sprintf("\nGROWTH FACTOR %s\n",
              if (is.null(g$coefficients)) "MEANS" else "INTERCEPTS"))
  mu <- g$means
  dimnames(mu) <- list(paste("Class", seq_len(K)), g$factor_names)
  print(round(mu, 3))

  if (!is.null(g$coefficients)) {
    cat(sprintf("\nGROWTH FACTORS ON COVARIATES%s\n",
                if (g$covariate_equal) " (held equal across classes)" else ""))
    for (i in seq_along(g$coefficients)) {
      if (!g$covariate_equal) cat(sprintf("  Class %d\n", i))
      print(round(g$coefficients[[i]], 3))
    }
  }

  if (length(g$random_names)) {
    cat(sprintf("\nGROWTH FACTOR (CO)VARIANCE%s\n",
                if (g$psi_equal) " (held equal across classes)" else ""))
    for (i in seq_along(g$psi)) {
      P <- g$psi[[i]]
      dimnames(P) <- list(g$random_names, g$random_names)
      if (!g$psi_equal) cat(sprintf("  Class %d\n", i))
      print(round(P, 3))
    }
  }

  cat(sprintf("\nRESIDUAL VARIANCE%s\n",
              if (g$residual_equal) " (held equal across classes)" else ""))
  rv <- g$residual_variance[seq_len(if (g$residual_equal) 1L else K), ,
                            drop = FALSE]
  dimnames(rv) <- list(
    if (g$residual_equal) "Variance" else paste("Class", seq_len(nrow(rv))),
    g$time_labels)
  print(round(rv, 3))

  cat(sprintf("\nFITTED TRAJECTORY (%s)\n",
              if (is.null(g$coefficients)) "mean" else
                "at the mean of the covariates"))
  ft <- g$fitted
  dimnames(ft) <- list(paste("Class", seq_len(K)), g$time_labels)
  print(round(ft, 3))

  NextMethod()
}
