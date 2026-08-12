# ==============================================================================
# Latent Class Growth Analysis (LCGA)
# ==============================================================================
#
# One outcome observed at T occasions. Each class carries its own fixed-effect
# growth curve - a polynomial in time on the link scale - and there are no
# within-class random effects, so conditional on class the occasions are
# independent and the likelihood factorises exactly as in an ordinary latent
# class model. This is Nagin's group-based trajectory model: a growth curve
# with no growth-factor variances.
#
# Parameterisation. A class is K x p coefficients on the design matrix
#
#     Lambda = [1, t, t^2, ...]      (T x p, p = degree + 1)
#
# so the class trajectory on the link scale is Lambda %*% beta_k. Threshold-parameterised
# software writes the same model with a threshold held equal across classes and the intercept
# growth factor fixed at zero in the last class; that is a reparameterisation of
# free per-class intercepts, not a restriction, and both forms have p*K free
# growth parameters. Estimates therefore translate as
#
#     intercept_k = mean(I)_k - threshold,     slope_k = mean(S)_k
#
# which is how the ex8.9 validation compares the two.

# Polynomial design matrix in the time scores.
#
# Raw rather than orthogonal polynomials: an applied reader expects the
# intercept to be the fitted value at time zero and the linear coefficient to be
# change per unit of time, which is the conventional reporting scale.
.lcga_design <- function(time_scores, degree) {
  design <- outer(time_scores, seq.int(0L, degree), "^")
  colnames(design) <- c("intercept", "linear", "quadratic", "cubic",
                        paste0("degree", seq_len(max(degree, 4L))))[
                          seq_len(degree + 1L)]
  design
}

#' Constructor for latent class growth models
#'
#' @description
#' Sets up the emission state for a latent class growth model, in which each
#' class follows its own fixed-effect polynomial trajectory over the occasions
#' and there is no within-class random effect.
#'
#' @param n_components Integer. The number of latent classes to estimate.
#' @param design Numeric matrix with one row per occasion and one column per
#'   growth coefficient, as built by the polynomial in the time scores.
#' @param family Character. `"binomial"`, `"gaussian"` or `"poisson"`.
#' @param ... Additional arguments, ignored.
#'
#' @return A list object of class \code{c("lcga", "emission")}.
#' @export
lcga_model <- function(n_components, design, family = "binomial", ...) {
  state <- list(
    n_components = n_components,
    design       = design,
    family       = family,
    fam          = .wglm_family(family),
    parameters   = list()
    # No `em_tol` here: refine_lbfgs() does not cover this emission, and
    # fit_single_init() gives every unpolished emission the tighter stopping
    # rule automatically. On a published benchmark the package default stopped 0.09
    # short of the maximum with the linear coefficient 6% out; the shared rule
    # lands within 0.011. See `.em_tol_unpolished` in R/em_core.R.
  )
  class(state) <- c("lcga", "emission")
  state
}

# The occasion-stacked ("long") view of the n x T outcome matrix, which is what
# the weighted GLM consumes: element (i, t) of X sits at row (t - 1) * n + i, so
# the design row for that element is row t of Lambda and a case's responsibility
# repeats across the T occasions. Built once per M-step rather than per class.
.lcga_long <- function(model_state, X) {
  n  <- nrow(X)
  Tn <- ncol(X)
  list(
    n = n, Tn = Tn,
    y = as.vector(X),
    D = model_state$design[rep(seq_len(Tn), each = n), , drop = FALSE],
    observed = !is.na(as.vector(X))
  )
}

# Observed marginal response at each occasion, which centres the prior that
# keeps a class from separating. Weighted when a survey design is in play, for
# the same reason refine_lbfgs() weights its Bernoulli marginal: the same data
# supplied as weighted patterns and as expanded rows must give one answer.
.lcga_marginal <- function(X, weights = NULL) {
  if (is.null(weights)) return(colMeans(X, na.rm = TRUE))
  obs <- !is.na(X)
  X0  <- replace(X, !obs, 0)
  colSums(X0 * weights) / pmax(colSums(obs * weights), 1e-12)
}

# Starting values: the pooled trajectory, perturbed once per class.
#
# A pooled GLM puts every class on the data's own scale, which matters because
# nothing else here is scale-free - a count outcome averaging 0.3 events and one
# averaging 300 need different starting coefficients. The perturbation is added
# to the fitted linear predictor at each occasion and then projected back onto
# the design, so a class starts on a smooth polynomial rather than a jagged one,
# and both level and shape are perturbed.
#
# The perturbation scale is one unit on the link scale for binomial and poisson,
# where logits and log rates are absolute, and the outcome's own standard
# deviation for gaussian, where the identity link means the linear predictor is
# measured in the units of y.
#' @exportS3Method
init_params.lcga <- function(model_state, X, resp, random_state = NULL, ...) {
  if (!is.null(random_state)) set.seed(random_state)

  K   <- model_state$n_components
  fam <- model_state$fam
  lng <- .lcga_long(model_state, X)
  p   <- ncol(model_state$design)

  pooled <- .wglm_fit(lng$D, lng$y, as.numeric(lng$observed), fam)
  eta0   <- as.vector(model_state$design %*% pooled$coefficients)

  s <- if (fam$name == "gaussian") {
    sd_y <- stats::sd(lng$y[lng$observed])
    if (is.finite(sd_y) && sd_y > 0) sd_y else 1
  } else 1

  coefs <- matrix(0, nrow = K, ncol = p)
  for (k in seq_len(K)) {
    eta_k <- eta0 + stats::rnorm(length(eta0), 0, s)
    coefs[k, ] <- qr.solve(model_state$design, eta_k)
  }

  model_state$parameters$coefs      <- coefs
  model_state$parameters$dispersion <- rep(pooled$dispersion, K)
  model_state
}

# M-step: one posterior-weighted GLM per class.
#' @exportS3Method
# The separation-guarding prior here is deliberately NOT wired to
# `bayes_constants`: it is a prior on GLM coefficients on a link scale, shared by
# the binomial, Poisson and Gaussian families, and none of the four named
# constants describes it. Giving it one of their names would let a user who
# turned off, say, the categorical prior silently change a Gaussian trajectory
# model too. It stays at 1 until it gets a name of its own.
m_step.lcga <- function(model_state, X, resp, weights = NULL, alpha = 1.0, ...) {
  if (!is.null(weights)) resp <- sweep(resp, 1, weights, "*")

  K   <- model_state$n_components
  fam <- model_state$fam
  lng <- .lcga_long(model_state, X)

  prior <- .wglm_prior_rows(model_state$design,
                            .lcga_marginal(X, weights),
                            alpha / K)

  y_obs <- replace(lng$y, !lng$observed, 0)

  coefs      <- model_state$parameters$coefs
  dispersion <- model_state$parameters$dispersion

  for (k in seq_len(K)) {
    w_k <- rep(resp[, k], lng$Tn)
    w_k[!lng$observed] <- 0     # FIML: a missing occasion informs nothing

    fit <- .wglm_fit(D     = rbind(lng$D, prior$D),
                     y     = c(y_obs,    prior$y),
                     w     = c(w_k,      prior$w),
                     fam   = fam,
                     start = coefs[k, ])
    coefs[k, ]    <- fit$coefficients
    dispersion[k] <- fit$dispersion
  }

  model_state$parameters$coefs      <- coefs
  model_state$parameters$dispersion <- dispersion
  model_state
}

#' @exportS3Method
log_likelihood.lcga <- function(model_state, X, ...) {
  n   <- nrow(X)
  Tn  <- ncol(X)
  K   <- model_state$n_components
  fam <- model_state$fam

  eta <- model_state$design %*% t(model_state$parameters$coefs)   # T x K
  log_eps <- matrix(0, nrow = n, ncol = K)

  for (k in seq_len(K)) {
    mu <- fam$linkinv(eta[, k])
    if (fam$name == "binomial") mu <- pmin(pmax(mu, 1e-15), 1 - 1e-15)
    if (fam$name == "poisson")  mu <- pmax(mu, 1e-15)
    ll <- fam$log_dens(X, matrix(mu, nrow = n, ncol = Tn, byrow = TRUE),
                       model_state$parameters$dispersion[k])
    # FIML: an unobserved occasion drops out of the case's sum. Also guards the
    # -Inf a saturated class would contribute at a probability pinned to 0 or 1.
    ll[!is.finite(ll)] <- 0
    log_eps[, k] <- rowSums(ll)
  }
  log_eps
}

#' @exportS3Method
n_parameters.lcga <- function(model_state, ...) {
  np <- length(model_state$parameters$coefs)
  if (isTRUE(model_state$fam$has_dispersion))
    np <- np + length(model_state$parameters$dispersion)
  np
}

# What each family will accept as an outcome.
#
# Checked here rather than left to the M-step because the failure a wrong family
# produces downstream is uninformative: glm.fit() reports "y values must be 0 <=
# y <= 1" for a binomial fit, and a Poisson fit of non-integer data simply
# returns a quasi-likelihood answer with no complaint at all. Only the values
# actually observed are checked; a missing occasion is missing, not invalid.
.lcga_check_outcome <- function(X, family) {
  vals <- X[!is.na(X)]
  if (!length(vals)) return(invisible(NULL))

  if (family == "binomial" && !all(vals %in% c(0, 1)))
    stop('family = "binomial" requires outcome values in {0, 1}.', call. = FALSE)

  if (family == "poisson") {
    if (any(vals < 0))
      stop('family = "poisson" requires non-negative counts, but negative ',
           'values are present.', call. = FALSE)
    # Counts read from a text file arrive as doubles, so what matters is that a
    # value is a whole number, not that it is stored as an integer.
    if (any(abs(vals - round(vals)) > 1e-8))
      stop('family = "poisson" requires whole-number counts, but fractional ',
           'values are present. For a continuous outcome use ',
           'family = "gaussian".', call. = FALSE)
  }

  if (family == "gaussian" && !all(is.finite(vals)))
    stop('family = "gaussian" requires finite outcome values; Inf or NaN are ',
         'present. Code unobserved occasions as NA, which FIML handles.',
         call. = FALSE)

  invisible(NULL)
}

# Fitted trajectory on the response scale: K x T, the quantity a trajectory plot
# draws and print() reports.
.lcga_fitted <- function(mm) {
  eta <- mm$design %*% t(mm$parameters$coefs)     # T x K
  t(mm$fam$linkinv(eta))
}

# ------------------------------------------------------------------------------
# User wrapper
# ------------------------------------------------------------------------------

#' Latent Class Growth Analysis
#'
#' @description
#' Fits class-specific growth curves to a single outcome measured repeatedly:
#' each latent class follows its own polynomial trajectory over time, with no
#' within-class random effects. This is Nagin's group-based trajectory modelling
#' (Nagin, 2005). Because there are no random
#' effects, the occasions are conditionally independent given class, so
#' everything the package already offers - model selection, the bootstrap
#' likelihood-ratio test, predictors of class membership, distal outcomes,
#' survey designs and FIML for missing occasions - applies unchanged.
#'
#' The trajectory is modelled on the link scale,
#' \deqn{g(E[y_{it} \mid c_i = k]) = \beta_{k0} + \beta_{k1} t + \beta_{k2} t^2 + \dots,}
#' with the logit link for a binary outcome. Contrast [`fit_rmlca()`], which
#' also produces trajectory classes but leaves each occasion's response
#' parameter free rather than constraining it to a smooth curve; LCGA is the
#' more parsimonious model and the one that extrapolates.
#'
#' **What the zero-variance assumption costs when you enumerate classes.** In an
#' LCGA the information criteria can keep improving as classes are added,
#' because the within-class variation the model refuses to estimate has to be
#' absorbed somewhere, and more classes is where it goes. Berlin et al. (2014,
#' p. 196) hit this in their own analysis — "increasing the number of latent
#' classes resulted in increasingly better (i.e., smaller) AICs, BICs, and
#' SSA-BICs, without any detriment to entropy" — and read it as a symptom
#' rather than a result: "this is probably a poorly specified model … The
#' assumption of zero variance within classes (as modeled here) is not likely
#' tenable and might account for each successive model … seeming to improve the
#' fit." So before reading a monotone BIC as evidence for many classes, fit a
#' [`fit_gmm()`] at the same K and compare.
#'
#' That is not an argument against the model. Where within-class homogeneity
#' genuinely holds, LCGAs "often allow for more straightforward
#' interpretations" (Ram & Grimm, 2009, p. 574), and the constraint "may be
#' particularly helpful when working with smaller sample sizes or when more
#' complex models fail due to nonconvergence, out of range estimates, or other
#' statistical problems, or as an initial modeling step prior to specifying a
#' GMM model" (Berlin et al., 2014, p. 191) — which is also how Jung and
#' Wickrama (2008, p. 304) recommend starting.
#'
#' @param indicator The repeated outcome. Either a wide matrix or data frame
#'   with one column per occasion, a three-dimensional array with dimensions n
#'   by 1 by times, or a long data frame together with `id` and `time`. Exactly
#'   one variable may be modelled; for several at once see [`fit_rmlca()`].
#'   Named for consistency with [`fit_rmlca()`]'s `indicators`, and so that
#'   `outcome=` stays free for a distal outcome caused by the trajectory class.
#' @param n_classes Integer. Number of trajectory classes.
#' @param times Integer. Number of occasions. Required for wide input; inferred
#'   otherwise.
#' @param degree Degree of the polynomial in time: `1` for a linear trajectory
#'   (intercept and slope), `2` for a quadratic, and so on. Must leave at least
#'   one occasion beyond the coefficients being estimated, so `degree` cannot
#'   exceed `times - 2`. In practice that means three time points support a
#'   linear pattern, four a quadratic as well, and five a cubic (Berlin et al.,
#'   2014, p. 191).
#' @param family Distribution of the outcome given class, which also sets the
#'   link the trajectory is linear on:
#'   * `"binomial"` — a binary outcome in `{0, 1}`, logit link. The trajectory
#'     is a curve in the probability of the event.
#'   * `"gaussian"` — a continuous outcome, identity link. The trajectory is a
#'     curve in the mean, and each class also carries a residual variance,
#'     constant across occasions. This is LCGA in the sense of Nagin's
#'     censored-normal model without the censoring, and the no-random-effects
#'     special case of a growth mixture model.
#'   * `"poisson"` — a count outcome, log link. The trajectory is a curve in the
#'     event rate. Counts must be whole and non-negative; over-dispersed or
#'     zero-inflated counts are not yet modelled as such.
#' @param time_scores Numeric values of time used in the polynomial, one per
#'   occasion. Defaults to `0, 1, ..., times - 1`, which makes the intercept the
#'   fitted value at the first occasion. Supply the actual measurement times
#'   when the occasions are unequally spaced.
#' @param layout For wide input with a three-dimensional array or several
#'   columns per occasion, whether columns run `"time_major"` or `"item_major"`.
#' @param id,time For long input, the case and occasion identifiers, given
#'   either as column names or as vectors.
#' @param item For long input, the column holding the outcome.
#' @param time_labels Optional display labels for the occasions.
#' @param predictors Optional predictors of class membership, passed to
#'   [`fit_mixture()`]'s three-step machinery.
#' @param ... Further arguments passed to [`fit_mixture()`], such as `outcome`
#'   and its companions, `n_init`, `random_state`, `weights`, `strata` or
#'   `cluster`.
#'
#' @return An object of class `c("lcga", "mixture_model")`. In addition to the
#'   usual fields it carries `$growth`, holding the design matrix, time scores,
#'   family, per-class coefficients, the fitted trajectories, and — for
#'   `family = "gaussian"` — the per-class residual variance.
#'
#' @seealso [`fit_rmlca()`] for unconstrained trajectory classes and
#'   [`fit_lta()`] for a model in which class membership itself changes.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' n <- 400
#' cls <- rbinom(n, 1, 0.5)
#' eta <- outer(ifelse(cls == 1, 0.0, -0.8), rep(1, 4)) +
#'   outer(ifelse(cls == 1, 1.1, -0.1), 0:3)
#' y <- matrix(rbinom(n * 4, 1, plogis(eta)), n, 4)
#' fit <- fit_lcga(y, times = 4, n_classes = 2, n_init = 5, random_state = 1)
#' fit
#' }
#' @references
#' Nagin, D. S. (1999). Analyzing developmental trajectories: A semiparametric,
#' group-based approach. \emph{Psychological Methods}, \emph{4}(2), 139-157.
#' \doi{10.1037/1082-989X.4.2.139}
#'
#' Nagin, D. S., & Odgers, C. L. (2010). Group-based trajectory modeling in
#' clinical research. \emph{Annual Review of Clinical Psychology}, \emph{6},
#' 109-138. \doi{10.1146/annurev.clinpsy.121208.131413}
#'
#' Berlin, K. S., Parra, G. R., & Williams, N. A. (2014). An introduction to
#' latent variable mixture modeling (part 2): Longitudinal latent class growth
#' analysis and growth mixture models. \emph{Journal of Pediatric Psychology},
#' \emph{39}(2), 188-203. \doi{10.1093/jpepsy/jst085}
#'
#' Jung, T., & Wickrama, K. A. S. (2008). An introduction to latent class growth
#' analysis and growth mixture modeling. \emph{Social and Personality Psychology
#' Compass}, \emph{2}(1), 302-317. \doi{10.1111/j.1751-9004.2007.00054.x}
#'
#' Lee, T. K., Wickrama, K. A. S., & O'Neal, C. W. (2023). An introduction to
#' growth mixture models (GMM). In \emph{International Encyclopedia of
#' Education} (4th ed., Vol. 14, pp. 646-655). Elsevier.
#' \doi{10.1016/B978-0-12-818630-5.10076-4}
#'
#' Ram, N., & Grimm, K. J. (2009). Growth mixture modeling: A method for
#' identifying differences in longitudinal change among unobserved groups.
#' \emph{International Journal of Behavioral Development}, \emph{33}(6),
#' 565-576. \doi{10.1177/0165025409343765}
#' @export
fit_lcga <- function(indicator,
                     n_classes = 2,
                     times = NULL,
                     degree = 1,
                     family = c("binomial", "gaussian", "poisson"),
                     time_scores = NULL,
                     layout = c("time_major", "item_major"),
                     id = NULL, time = NULL, item = NULL,
                     time_labels = NULL,
                     predictors = NULL,
                     ...) {

  family <- match.arg(family)
  layout <- match.arg(layout)

  prep <- .prepare_longitudinal(indicator, times = times, items = item,
                                layout = layout, id = id, time = time,
                                time_labels = time_labels)

  if (prep$n_items != 1L)
    stop(sprintf(
      paste0("LCGA models one repeated outcome, but %d were found at each of ",
             "the %d occasions. Select a single variable, or use fit_rmlca() ",
             "for trajectory classes over several indicators."),
      prep$n_items, prep$n_times), call. = FALSE)

  Tn <- prep$n_times
  if (is.null(time_scores)) time_scores <- seq.int(0L, Tn - 1L)
  if (length(time_scores) != Tn)
    stop(sprintf("`time_scores` must have one value per occasion (%d given, %d needed).",
                 length(time_scores), Tn), call. = FALSE)

  degree <- as.integer(degree)
  if (is.na(degree) || degree < 0L)
    stop("`degree` must be a non-negative integer.", call. = FALSE)
  # p = degree + 1 coefficients from T occasions. Allowing p == T would fit the
  # occasion means exactly, making the growth constraint vacuous and the model
  # an unconstrained RMLCA in disguise; requiring one spare occasion keeps the
  # trajectory a restriction that can be tested.
  if (degree + 1L >= Tn)
    stop(sprintf(
      paste0("A degree-%d trajectory needs at least %d occasions but %d were ",
             "given. With degree + 1 = %d coefficients per class the curve ",
             "would reproduce the occasion means exactly."),
      degree, degree + 2L, Tn, degree + 1L), call. = FALSE)

  .lcga_check_outcome(prep$X, family)

  design <- .lcga_design(time_scores, degree)

  if (!is.null(predictors))
    predictors <- .as_named_covariates(predictors, substitute(predictors),
                                       "predictor")

  fit <- fit_mixture(
    indicators  = prep$X,
    n_classes   = n_classes,
    measurement = "lcga",
    predictors  = predictors,
    design      = design,
    family      = family,
    ...
  )

  fit$growth <- list(
    model       = "lcga",
    family      = family,
    degree      = degree,
    time_scores = time_scores,
    time_labels = prep$time_labels,
    design      = design,
    coefficients = fit$mm$parameters$coefs,
    fitted      = .lcga_fitted(fit$mm),
    wave_missing = prep$wave_missing
  )
  # Carried only where it means something: the other two families have no free
  # dispersion, and the emission holds a placeholder 1 for them that would read
  # as an estimate if it were reported.
  if (family == "gaussian")
    fit$growth$residual_variance <- fit$mm$parameters$dispersion
  class(fit) <- c("lcga", class(fit))
  fit
}

#' Trajectory Plot for a Latent Class Growth Model
#'
#' @description
#' Draws the estimated trajectory of each class on the response scale, with the
#' occasions on the x-axis — the figure an LCGA is read from, since the classes
#' *are* the trajectories — and the observed data of the cases assigned to each
#' class behind them.
#'
#' @param x An object returned by [`fit_lcga()`].
#' @param observed What of the observed data to draw: `"means"` (the observed
#'   mean at each occasion among cases modally assigned to the class, dotted),
#'   `"cases"` (individual trajectories, translucent — informative for a
#'   continuous or count outcome, much less so for a binary one) or `"none"`.
#' @param main Plot title.
#' @param class_labels Optional class labels for the legend.
#' @param colors Optional colour vector, recycled across classes.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @importFrom graphics matplot matlines axis legend
#' @export
plot.lcga <- function(x, observed = c("means", "cases", "none"), main = NULL,
                      class_labels = NULL, colors = NULL, ...) {
  observed <- match.arg(observed)
  g  <- x$growth
  K  <- x$n_components
  ft <- g$fitted

  cols   <- if (is.null(colors)) rep(.okabe_ito, length.out = K)
            else rep(colors, length.out = K)
  shapes <- rep(15:20, length.out = K)
  base   <- if (is.null(class_labels)) paste("Class", seq_len(K)) else class_labels
  labels <- .class_plot_labels(base, x$weights)

  obs  <- .growth_observed(x, observed)
  orng <- .growth_observed_range(obs)

  # A probability panel is always the full 0-1 range, so classes stay comparable
  # across figures. A rate is anchored at zero, because how close a class sits to
  # "no events" is the thing being read. A mean is anchored at nothing: forcing
  # zero in would squash every curve into a corner of an outcome measured around,
  # say, 100, so the panel is the range of the curves with a margin.
  ylim <- switch(g$family,
                 binomial = c(0, 1),
                 poisson  = range(c(0, ft, orng), na.rm = TRUE),
                 gaussian = {
                   r <- range(c(ft, orng), na.rm = TRUE)
                   # A single flat class gives a zero-width range, which is not
                   # a drawable panel.
                   pad <- if (diff(r) > 0) 0.05 * diff(r) else max(abs(r[1]), 1) * 0.05
                   r + c(-pad, pad)
                 })

  matplot(g$time_scores, t(ft), type = "n", ylim = ylim, xaxt = "n", las = 1,
          bty = "l", xlab = "Occasion", ylab = x$mm$fam$scale_label,
          main = main %||% "Estimated class trajectories")
  .growth_draw_observed(obs, g$time_scores, cols)
  matlines(g$time_scores, t(ft), type = "b", pch = shapes, lty = 1, lwd = 2.5,
           col = cols)
  axis(1, at = g$time_scores, labels = g$time_labels)
  legend("topleft", legend = labels, col = cols, pch = shapes, lty = 1,
         lwd = 2, bty = "n", cex = 0.85)

  invisible(x)
}

#' Print a Fitted Latent Class Growth Model
#'
#' @param x An object returned by [`fit_lcga()`].
#' @param ... Passed to the next method.
#' @return `x`, invisibly.
#' @export
print.lcga <- function(x, ...) {
  g <- x$growth
  shape <- c("intercept only", "linear", "quadratic", "cubic")[
    min(g$degree + 1L, 4L)]
  cat("\n")
  cat("=========================================================\n")
  cat("           LATENT CLASS GROWTH ANALYSIS\n")
  cat("=========================================================\n")
  cat(sprintf("Occasions          : %d (time scores %s)\n",
              length(g$time_scores),
              paste(format(g$time_scores, trim = TRUE), collapse = ", ")))
  cat(sprintf("Trajectory         : %s, %s link\n", shape,
              switch(g$family, binomial = "logit", poisson = "log",
                     gaussian = "identity")))
  if (any(g$wave_missing))
    cat(sprintf("Wave attrition     : %d case-occasions unobserved\n",
                sum(g$wave_missing)))

  cat("\nGROWTH COEFFICIENTS (link scale)\n")
  co <- g$coefficients
  dimnames(co) <- list(paste("Class", seq_len(nrow(co))), colnames(g$design))
  print(round(co, 3))

  cat(sprintf("\nFITTED TRAJECTORY (%s)\n", tolower(x$mm$fam$scale_label)))
  ft <- g$fitted
  dimnames(ft) <- list(paste("Class", seq_len(nrow(ft))), g$time_labels)
  print(round(ft, 3))

  # The residual variance is a reported quantity in a continuous LCGA: it is
  # what the trajectory does not explain, and the ratio between classes says
  # whether one is tighter around its curve than another.
  if (!is.null(g$residual_variance)) {
    cat("\nRESIDUAL VARIANCE (within class, constant over occasions)\n")
    rv <- matrix(g$residual_variance, nrow = 1,
                 dimnames = list("Variance",
                                 paste("Class", seq_along(g$residual_variance))))
    print(round(rv, 3))
  }

  NextMethod()
}
