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
#' Applied papers usually name a growth mixture model by what the classes are
#' allowed to differ in, and Ram and Grimm (2009, pp. 569-570) give the three
#' levels. *Means* — classes differ in mean change only — is `psi = "equal"`,
#' the default here. *Means + Covs* — classes also differ in how much
#' interindividual variation there is within them — is `psi = "free"`; the
#' residual variances are a separate switch, `residual_equal = FALSE`, and
#' should not be merged with it. *Means + Covs + Pattern*, where classes differ
#' in the shape of change as well, would need class-specific time scores;
#' `fit_gmm()` builds one design for every class, so that level is not
#' available.
#'
#' Two cautions worth carrying into any such analysis. Classes will be found
#' whether or not there are groups to find: a skewed or otherwise non-normal
#' outcome can be fitted by extra classes that are not subgroups of anyone (Jung
#' & Wickrama, 2008, p. 305; Lee et al., 2023, p. 652, both citing Bauer &
#' Curran, 2003). And as Ram and Grimm put it (p. 574), "groups and differences
#' among groups will be found; but whether they represent true processes that
#' generated the data is unknown" — the remedies they name are replication on
#' new data and checking that class membership relates to other measured
#' variables in the ways theory predicts.
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
#'   (intercept and slope), `2` for a quadratic, and so on. Each degree needs
#'   occasions to identify it — with three time points a linear pattern can be
#'   modelled, with four a quadratic as well, and with five a cubic (Berlin et
#'   al., 2014, p. 191). Asking for more than the occasions support is an error
#'   rather than a warning.
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
#'
#'   It is not free, though, and the cost runs the other way. Simulation work
#'   reported by Lee et al. (2023, p. 651) finds that equality restrictions on
#'   the growth-factor variances and covariances "could result in the
#'   over-extraction of latent classes and biased parameter estimates" — the
#'   within-class heterogeneity the constraint refuses to let differ between
#'   classes has to go somewhere, and an extra class is where it goes. So
#'   `psi = "equal"` buys stability and bounds the likelihood, and it can buy an
#'   extra class that is an artefact of the constraint. Where the number of
#'   classes is itself the finding, fit both and report which was used.
#'
#'   This pulls against the collapsed-variance warning, which fires more often
#'   under `psi = "free"`. That is not a contradiction: they are the two sides
#'   of one trade-off between a model flexible enough to be realistic and one
#'   constrained enough to be estimable.
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
#' Berlin, K. S., Parra, G. R., & Williams, N. A. (2014). An introduction to
#' latent variable mixture modeling (part 2): Longitudinal latent class growth
#' analysis and growth mixture models. *Journal of Pediatric Psychology*,
#' *39*(2), 188-203. \doi{10.1093/jpepsy/jst085}
#'
#' Jung, T., & Wickrama, K. A. S. (2008). An introduction to latent class growth
#' analysis and growth mixture modeling. *Social and Personality Psychology
#' Compass*, *2*(1), 302-317. \doi{10.1111/j.1751-9004.2007.00054.x}
#'
#' Lee, T. K., Wickrama, K. A. S., & O'Neal, C. W. (2023). An introduction to
#' growth mixture models (GMM). In *International Encyclopedia of Education*
#' (4th ed., Vol. 14, pp. 646-655). Elsevier.
#' \doi{10.1016/B978-0-12-818630-5.10076-4}
#'
#' Ram, N., & Grimm, K. J. (2009). Growth mixture modeling: A method for
#' identifying differences in longitudinal change among unobserved groups.
#' *International Journal of Behavioral Development*, *33*(6), 565-576.
#' \doi{10.1177/0165025409343765}
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
  fit$growth$boundary <- .gmm_boundary(mm, fit$data, fit$sample_weights)
  if (length(fit$growth$boundary))
    warning(.gmm_boundary_message(fit$growth$boundary, mm), call. = FALSE)

  # Issued here rather than in fit_mixture(), which returned before the boundary
  # flag above existed: the two warnings give opposite advice about n_init, so
  # the one that reads the flag has to run after it is set.
  .check_replication(fit)

  class(fit) <- c("gmm", class(fit))
  fit
}

# Variance parameters that have collapsed relative to the data, not merely to
# the numerical floor.
#
# A growth mixture model whose covariance structure the data do not support does
# not fail; it produces a variance that wants to be negative. Some programs let
# it go negative and print it. This package floors instead, which keeps the
# likelihood a likelihood, but a floored estimate is the same diagnosis and
# would otherwise be reported as though it were an ordinary small variance.
#
# The floor tests this function used to run on their own only fire once the
# M-step has pinned a parameter at 1e-6 or an eigenvalue at the 1e-8 clip. That
# is the last stage of a collapse, not the diagnostic one: on a reference
# class-varying fit the shared residual variance at the final occasion reached
# 1.19e-4 -- 0.08% of that occasion's observed variance, three orders of
# magnitude below every other occasion -- and neither test fired, while the fit
# carried the lowest BIC of the whole model set. The ratio to the observed
# marginal is the statistic `.gaussian_boundary()` already uses for the same
# failure in a latent profile model, at the same 1% threshold and for the same
# reason (R/gaussian_boundary.R:72): it means the same thing whatever the
# outcome is measured in. The floor tests are kept as an unconditional second
# trigger, since a parameter genuinely pinned at the floor must still be
# reported even if the marginal is unavailable.
#
# For Psi the ratio is taken on the random effects' *contribution to the implied
# variance of the outcome*, max_t diag(Lambda_r Psi_k Lambda_r')_t / s2_t, not on
# Psi's own entries. Psi's entries are on the growth factors' scale, which for a
# slope is the outcome's scale divided by time squared, so a bare comparison
# would flag every slope variance ever estimated. The eigenvalue *ratio* is not
# usable either: it is 0.004 on a collapsed Psi and 0.004 on a healthy one.
#
# The share is taken over the whole of Psi_k, per class, rather than per growth
# factor. A slope variance going to zero under a healthy intercept variance is
# not a degeneracy: it is the random_effects = "intercept" model, which this
# package offers and which is a reportable finding. What is a degeneracy is the
# whole of Psi_k going to zero -- that class is an LCGA class -- and a residual
# variance going to zero, which is the unbounded-likelihood spike.
.gmm_boundary <- function(mm, X = NULL, weights = NULL, threshold = 0.01) {
  msg <- character(0)

  s2 <- .gmm_marginal(mm, X, weights)
  Tn <- nrow(mm$design)

  theta      <- mm$parameters$theta
  theta_rows <- if (isTRUE(mm$theta_equal)) 1L else nrow(theta)
  # With residual = "constant" the columns of theta are one parameter, not T of
  # them, so it is judged once. The marginal it is compared with is the smallest
  # of the occasions', which is the conservative choice: a single variance
  # shared over occasions of different spread is only a collapse if it is small
  # relative to even the least variable of them.
  shared_t <- isTRUE(mm$theta_shared_occasions)
  occasions <- if (shared_t) 1L else seq_len(Tn)
  for (k in seq_len(theta_rows)) {
    for (t in occasions) {
      v  <- theta[k, t]
      s2t <- if (is.null(s2)) NA_real_ else if (shared_t) min(s2) else s2[t]
      at_floor <- v <= .sn_theta_floor * 10
      ratio    <- v / s2t
      if (!at_floor && !(is.finite(ratio) && ratio < threshold)) next
      where <- paste0(
        if (theta_rows == 1L) "" else sprintf("in class %d ", k),
        if (shared_t) "(constant over occasions)" else sprintf("at occasion %d", t))
      msg <- c(msg, if (is.finite(ratio))
        sprintf("residual variance %s (%.3g vs %.3g for %s, %s)",
                where, v, s2t,
                if (shared_t) "the least variable occasion"
                else "that occasion overall",
                .gmm_pct(ratio))
      else sprintf("residual variance %s at the estimation floor (%.3g)",
                   where, v))
    }
  }

  if (length(mm$r_cols)) {
    Lr <- mm$design[, mm$r_cols, drop = FALSE]
    psi_rows <- if (isTRUE(mm$psi_equal)) 1L else length(mm$parameters$psi)
    for (k in seq_len(psi_rows)) {
      P <- mm$parameters$psi[[k]]
      if (is.null(P) || !length(P)) next
      implied  <- diag(Lr %*% P %*% t(Lr))
      singular <- min(eigen(P, symmetric = TRUE, only.values = TRUE)$values) <=
        1e-7
      share <- if (is.null(s2)) NA_real_ else max(implied / s2)
      if (!singular && !(is.finite(share) && share < threshold)) next
      where <- if (psi_rows == 1L) "the growth factors"
               else sprintf("the growth factors in class %d", k)
      detail <- if (is.finite(share))
        sprintf(" (contributing at most %s of the outcome's variance)",
                .gmm_pct(share))
      else ""
      msg <- c(msg, sprintf(
        paste0("no within-class variation left in %s%s, which is the latent ",
               "class growth model fit_lcga() estimates directly"),
        where, detail))
    }
  }
  msg
}

# Observed variance of each occasion, or NULL when the data are not available.
#
# Wrapped rather than inlined because a caller that hands over a matrix of the
# wrong width -- a mis-shaped fit, or an internal call made before the data are
# resolved -- should fall back to the floor tests rather than compare a variance
# with an unrelated column's.
.gmm_marginal <- function(mm, X, weights = NULL) {
  if (is.null(X) || !is.matrix(X) || ncol(X) != nrow(mm$design)) return(NULL)
  if (!is.null(weights) && length(weights) != nrow(X)) weights <- NULL
  s2 <- .marginal_var(X, weights)
  if (any(!is.finite(s2)) || any(s2 <= 0)) return(NULL)
  s2
}

# A ratio as a percentage, at a precision that survives three orders of
# magnitude: "0.08%" says what "0%" would hide.
.gmm_pct <- function(ratio) {
  pct <- 100 * ratio
  sprintf("%s%%", format(signif(pct, 2), scientific = FALSE, trim = TRUE))
}

# The warning text. Two things it must do that naming the parameter does not.
#
# It must say that BIC cannot be compared across a degenerate fit and a clean
# one, because a collapsed variance *wins* on BIC -- on the reference data the
# degenerate solution beat every admissible model in the set by 40 points -- so
# a user following ordinary model-selection practice is led straight to it.
#
# And it must say that raising n_init is the wrong response. This is not a
# convergence failure: on that same fit the spike was found by 1 start in 50 and
# the clean solution by 21 in 30, so more starts means more chances to find the
# spike. `.check_gaussian_degeneracy()` already says this for the latent profile
# case, for the same reason, and the two warnings are meant to read as one
# family.
#
# The remedies are named in order of how much structure they impose, and only
# when they are not already in force -- advice a user has already taken reads as
# noise and buries the advice they have not.
.gmm_boundary_message <- function(lines, mm = NULL) {
  fixes <- character(0)
  if (is.null(mm) || !isTRUE(mm$psi_equal))
    fixes <- c(fixes, "hold the growth-factor covariance equal across classes with psi = \"equal\"")
  if (is.null(mm) || !isTRUE(mm$theta_equal))
    fixes <- c(fixes, "hold the residual variances equal across classes with residual_equal = TRUE")
  if (is.null(mm) || length(mm$r_cols) > 1L)
    fixes <- c(fixes, "let only the intercept vary within a class with random_effects = \"intercept\"")
  fixes <- c(fixes, "fit fewer classes, since a class the data cannot support usually means there are too many")

  sprintf(
    paste0("A variance has collapsed towards zero: %s. The likelihood of a ",
           "mixture of normals is unbounded in this direction, so this ",
           "solution can score better than any meaningful one while ",
           "describing a handful of near-identical cases rather than a ",
           "subgroup. Do not interpret these estimates as variances -- their ",
           "standard errors are not valid either -- and do not compare this ",
           "fit's BIC with a clean one's, since it is inflated by the spike. ",
           "This is not a convergence failure, so raising n_init will not fix ",
           "it and can make it worse. Ways out, to choose between on ",
           "substantive grounds: %s. A class with no growth-factor variation ",
           "left is a latent class growth class, which fit_lcga() estimates ",
           "directly and without this failure mode."),
    paste(lines, collapse = "; "),
    paste(sprintf("(%d) %s", seq_along(fixes), fixes), collapse = "; "))
}

# Repeat the flag in print(). Someone opening a saved fit months later should
# still see that its variances collapsed; the warning at fit time is gone by
# then. `.print_degenerate_note()` (R/gaussian_boundary.R) is the model.
.print_gmm_boundary_note <- function(x) {
  lines <- x$growth$boundary
  if (!length(lines)) return(invisible(NULL))
  cat("\nWARNING - collapsed variance:\n")
  for (line in lines) cat("  ", line, "\n", sep = "")
  cat("  These estimates are not interpretable, and this fit's BIC cannot be\n")
  cat("  compared with a clean one's. See ?fit_gmm for what to do.\n")
  invisible(NULL)
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

  .print_gmm_boundary_note(x)

  NextMethod()
}

# ------------------------------------------------------------------------------
# The growth parameters as a table
# ------------------------------------------------------------------------------

# One printed block: a K x J matrix with the classes across the top, in the
# layout measurement_summary.default() uses for item parameters, so a growth
# model's table reads like every other model's.
.growth_print_block <- function(title, mat, K) {
  cat(sprintf("\n%s\n", title))
  disp    <- .shorten_labels(colnames(mat), width = 30L)
  label_w <- max(20L, max(nchar(disp)))
  cat(sprintf("%-*s", label_w, "Parameter"))
  for (k in seq_len(K)) cat(sprintf(" | Class %d", k))
  cat("\n")
  cat(paste0(rep("-", label_w + K * 10), collapse = ""), "\n")
  for (j in seq_len(ncol(mat))) {
    cat(sprintf("%-*s", label_w, disp[j]))
    for (k in seq_len(K)) cat(sprintf(" | %7.3f", mat[k, j]))
    cat("\n")
  }
  .cat_label_legend(disp, indent = "")
  invisible(NULL)
}

# Long-format rows for one block, in the schema measurement_summary() returns.
.growth_rows <- function(mat, parameter, K) {
  J <- ncol(mat)
  data.frame(block     = rep(NA_character_, K * J),
             parameter = rep(parameter, K * J),
             item      = rep(colnames(mat), each = K),
             category  = rep(NA_integer_, K * J),
             class     = rep(seq_len(K), times = J),
             estimate  = as.vector(mat),
             stringsAsFactors = FALSE)
}

# A constrained parameter is repeated across the classes rather than reported
# once. The alternative -- one row with class = NA -- would make the table
# unjoinable to class_sizes() or to the posterior, which is most of what a tidy
# table is for, and would quietly drop the parameter out of any per-class
# summary a user writes. The constraint is stated in the printed heading
# instead, exactly as print.gmm() states it.
.growth_expand <- function(mat, K) {
  if (nrow(mat) == K) return(mat)
  out <- matrix(rep(mat[1L, ], each = K), nrow = K,
                dimnames = list(NULL, colnames(mat)))
  out
}

# Every measurement parameter of a growth model, as one K-by-parameter matrix
# per family, in print order.
.growth_blocks <- function(x) {
  g <- x$growth
  K <- x$n_components
  out <- list()

  add <- function(title, parameter, mat) {
    if (is.null(mat) || !length(mat)) return(invisible(NULL))
    out[[length(out) + 1L]] <<- list(title = title, parameter = parameter,
                                     mat = .growth_expand(mat, K))
  }

  if (identical(g$model, "lcga")) {
    means <- g$coefficients
    colnames(means) <- colnames(g$design)
    add(sprintf("GROWTH COEFFICIENTS (%s scale)",
                if (identical(g$family, "gaussian")) "response" else "link"),
        "growth_mean", means)
  } else {
    means <- g$means
    colnames(means) <- g$factor_names
    add(if (is.null(g$coefficients)) "GROWTH FACTOR MEANS"
        else "GROWTH FACTOR INTERCEPTS", "growth_mean", means)

    # The covariate regressions are measurement parameters here -- they enter
    # the emission, not the structural model -- so a table that left them out
    # would be missing estimated parameters on exactly the models where the
    # growth means are no longer means.
    if (!is.null(g$coefficients)) {
      cf <- lapply(g$coefficients, function(B) {
        v <- as.vector(B)
        names(v) <- as.vector(outer(rownames(B), colnames(B),
                                    function(a, b) paste(a, "ON", b)))
        v
      })
      M <- do.call(rbind, cf)
      colnames(M) <- names(cf[[1L]])
      add(sprintf("GROWTH FACTORS ON COVARIATES%s",
                  if (g$covariate_equal) " (held equal across classes)" else ""),
          "growth_regression", M)
    }

    if (length(g$random_names)) {
      nm  <- g$random_names
      q   <- length(nm)
      idx <- which(upper.tri(diag(q), diag = TRUE), arr.ind = TRUE)
      lab <- ifelse(idx[, 1] == idx[, 2], nm[idx[, 1]],
                    paste(nm[idx[, 1]], nm[idx[, 2]], sep = " with "))
      V <- do.call(rbind, lapply(g$psi, function(P) P[idx]))
      colnames(V) <- lab
      diag_cols <- idx[, 1] == idx[, 2]
      add(sprintf("GROWTH FACTOR VARIANCES%s",
                  if (g$psi_equal) " (held equal across classes)" else ""),
          "growth_variance", V[, diag_cols, drop = FALSE])
      if (any(!diag_cols))
        add(sprintf("GROWTH FACTOR COVARIANCES%s",
                    if (g$psi_equal) " (held equal across classes)" else ""),
            "growth_covariance", V[, !diag_cols, drop = FALSE])
    }
  }

  rv <- g$residual_variance
  if (!is.null(rv)) {
    if (!is.matrix(rv)) rv <- matrix(rv, ncol = 1L)
    colnames(rv) <- if (ncol(rv) == length(g$time_labels)) g$time_labels
                    else "Variance"
    add(sprintf("RESIDUAL VARIANCE%s",
                if (isTRUE(g$residual_equal)) " (held equal across classes)"
                else ""),
        "residual_variance", rv)
  }

  ft <- g$fitted
  colnames(ft) <- g$time_labels
  add("FITTED TRAJECTORY", "fitted", ft)

  out
}

#' @rdname measurement_summary
#' @export
measurement_summary.gmm <- function(object, ...) {
  K <- object$n_components
  cat("=========================================================\n")
  cat("             MEASUREMENT MODEL PARAMETERS                \n")
  cat("=========================================================\n")

  blocks <- .growth_blocks(object)
  for (b in blocks) .growth_print_block(b$title, b$mat, K)

  if (!is.null(object$missing_data) && isTRUE(object$missing_data$any_missing)) {
    md <- object$missing_data
    cat(sprintf(paste0("\nMissing data: %d of %d cells (%.1f%%) across %d ",
                       "occasion%s, handled via %s.\n"),
                md$n_missing, md$n_cells, 100 * md$prop_missing,
                md$n_items_affected, if (md$n_items_affected == 1L) "" else "s",
                md$handled_by))
  }
  .print_gmm_boundary_note(object)
  cat("=========================================================\n")

  invisible(do.call(rbind, lapply(blocks, function(b)
    .growth_rows(b$mat, b$parameter, K))))
}

#' @rdname measurement_summary
#' @export
measurement_summary.lcga <- measurement_summary.gmm
