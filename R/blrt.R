# ==============================================================================
# Bootstrap Likelihood Ratio Test (BLRT) - DYNAMIC MEASUREMENT
# ==============================================================================

# ------------------------------------------------------------------------------
# HELPER: Dynamically generate synthetic data from a fitted measurement model
# (Uses "Duck Typing" to perfectly handle any S3 class name)
# ------------------------------------------------------------------------------
generate_synthetic_data <- function(mm, classes, N) {
  # A time-block model is a column-partitioned model just like `nested`: each
  # occasion's sub-model generates its own block and the blocks are laid out
  # time-major, which reproduces the original column order exactly.
  #
  # Deliberately NOT widened to also cover `group_blocks`: unlike time blocks,
  # a group-blocks row is structurally missing everywhere but its own group's
  # block, and this generator has no group label to reproduce that pattern, so
  # naively generating fully-observed data for every group would not resemble
  # the real design the null distribution is supposed to approximate. A
  # `group_blocks` model falls through to the `stop()` below until that's
  # built; use `blrt()`/`compare_mixtures()` on the pooled (non-grouped) data
  # to choose the number of classes, as Collins & Lanza do (sec. 5.7.1),
  # before fitting the multiple-group model.
  if (inherits(mm, c("nested", "time_blocks"))) {
    res <- list()
    for (name in names(mm$models)) {
      res[[name]] <- generate_synthetic_data(mm$models[[name]], classes, N)
    }
    return(do.call(cbind, res))
  }

  is_cont   <- !is.null(mm$parameters$means)
  is_poly   <- !is.null(mm$max_val)
  is_bin    <- !is.null(mm$parameters$pis) && is.null(mm$max_val)
  is_count  <- !is.null(mm$parameters$rates)
  is_growth <- !is.null(mm$parameters$coefs)
  is_gmm    <- !is.null(mm$parameters$alpha)

  # A growth mixture model draws a whole T-vector at once rather than one value
  # per occasion independently: the within-class random effects are exactly what
  # makes a case's occasions correlated, so generating them occasion by occasion
  # would produce a null distribution for a model nobody fitted. The class's
  # implied covariance is the structured Sigma, and a multivariate normal draw
  # from it is a Cholesky factor times a standard normal vector.
  if (is_gmm) {
    Tn    <- nrow(mm$design)
    X_gen <- matrix(0, nrow = N, ncol = Tn)
    chols <- lapply(seq_len(mm$n_components),
                    function(k) .sn_chol(.sn_sigma(mm, k)))

    # With covariates on the growth factors the class mean is a case-level
    # quantity, so the replicate is generated *conditional on the observed x* --
    # the covariates are exogenous and are not resampled, which is what makes
    # the null distribution the one for the model that was fitted. That requires
    # the replicate to have as many cases as the original.
    has_x <- !is.null(mm$xmat)
    if (has_x && N != nrow(mm$xmat))
      stop(sprintf(
        paste0("A growth mixture model with covariates on the growth factors ",
               "generates replicates conditional on the observed covariates, ",
               "so a replicate must have one case per original case (%d ",
               "requested, %d available)."), N, nrow(mm$xmat)), call. = FALSE)
    mu <- if (has_x) lapply(seq_len(mm$n_components), function(k) .sn_mu(mm, k))
          else NULL
    fitted <- if (has_x) NULL else .sn_fitted(mm)      # K x T

    for (i in seq_len(N)) {
      k   <- classes[i]
      m_i <- if (has_x) mu[[k]][i, ] else fitted[k, ]
      X_gen[i, ] <- m_i + as.vector(crossprod(chols[[k]], rnorm(Tn)))
    }
    return(X_gen)
  }

  # A growth model's class parameters are coefficients rather than one value per
  # column, so the trajectory has to be evaluated before anything can be drawn;
  # after that it generates exactly like the flat model of the same family, one
  # draw per occasion. Class enumeration matters more here than almost anywhere
  # else in the package — how many trajectory groups there are is usually the
  # research question — so the BLRT has to reach this model.
  if (is_growth) {
    mu <- .lcga_fitted(mm)                    # K x T on the response scale
    Tn <- ncol(mu)
    X_gen <- matrix(0, nrow = N, ncol = Tn)
    for (i in seq_len(N)) {
      m_i <- mu[classes[i], ]
      X_gen[i, ] <- switch(
        mm$fam$name,
        binomial = rbinom(Tn, 1, prob = m_i),
        poisson  = rpois(Tn, lambda = m_i),
        gaussian = rnorm(Tn, mean = m_i,
                         sd = sqrt(mm$parameters$dispersion[classes[i]])))
    }
    return(X_gen)

  } else if (is_count) {
    rates <- mm$parameters$rates
    D     <- ncol(rates)
    X_gen <- matrix(0L, nrow = N, ncol = D)
    for (i in 1:N)
      X_gen[i, ] <- rpois(D, lambda = rates[classes[i], ])
    return(X_gen)

  } else if (is_cont) {
    D <- ncol(mm$parameters$means)
    X_gen <- matrix(0, nrow = N, ncol = D)
    for (i in 1:N) {
      vars <- if (!is.null(mm$parameters$covariances)) mm$parameters$covariances[classes[i], ] else rep(1, D)
      X_gen[i, ] <- rnorm(D, mean = mm$parameters$means[classes[i], ], sd = sqrt(vars))
    }
    return(X_gen)

  } else if (is_poly) {
    D <- ncol(mm$parameters$pis)
    M <- mm$max_val
    n_items <- D / M
    X_gen <- matrix(0L, nrow = N, ncol = n_items)
    for (i in 1:N) {
      probs <- matrix(mm$parameters$pis[classes[i], ], nrow = n_items, ncol = M, byrow = TRUE)
      for (j in 1:n_items) {
        X_gen[i, j] <- sample(1:M, 1, prob = probs[j, ])
      }
    }
    return(X_gen)

  } else if (is_bin) {
    D <- ncol(mm$parameters$pis)
    X_gen <- matrix(0, nrow = N, ncol = D)
    for (i in 1:N) {
      X_gen[i, ] <- rbinom(D, 1, prob = mm$parameters$pis[classes[i], ])
    }
    return(X_gen)

  } else {
    stop(sprintf("Unsupported measurement model parameters. S3 Class: %s", paste(class(mm), collapse=", ")))
  }
}

# The descriptor and engine arguments that reproduce a fitted growth model at a
# different number of classes.
#
# Everything here is already on the fitted object; the only thing missing was a
# door. Building it by hand meant calling .lcga_design(), which is internal, and
# knowing which of the emission's fields are engine arguments -- so in practice
# the BLRT was unreachable for the models whose class enumeration is most often
# the research question.
#
# The covariates on the growth factors come across too. They are exogenous and
# are not resampled (see generate_synthetic_data()), so a replicate has one case
# per original case, which is what the BLRT does anyway.
.blrt_growth_spec <- function(fit) {
  if (!inherits(fit, c("gmm", "lcga")))
    stop(sprintf(
      paste0("`from_fit` takes a model fitted by fit_gmm() or fit_lcga(); ",
             "this is a %s. For every other model the specification is the ",
             "`measurement` argument, which blrt() already takes."),
      paste(class(fit), collapse = "/")), call. = FALSE)

  g <- fit$growth
  if (is.null(fit$data))
    stop("`from_fit` needs the fitted model's data, which this object does ",
         "not carry.", call. = FALSE)

  args <- if (inherits(fit, "lcga"))
    list(design = g$design, family = g$family)
  else
    list(design                  = g$design,
         random_effects          = g$random_effects,
         psi                     = if (isTRUE(g$psi_equal)) "equal" else "free",
         residual                = if (isTRUE(g$residual_constant)) "constant"
                                   else "occasion",
         residual_equal          = isTRUE(g$residual_equal),
         growth_covariates       = fit$mm$xmat,
         growth_covariates_equal = isTRUE(g$covariate_equal))

  list(indicators  = fit$data,
       measurement = if (inherits(fit, "lcga")) "lcga" else "structured_normal",
       args        = args)
}

# ------------------------------------------------------------------------------
# MAIN BLRT FUNCTION
# ------------------------------------------------------------------------------

#' Bootstrap Likelihood Ratio Test (BLRT) for Class Enumeration
#'
#' @description
#' Tests whether a model with more classes fits significantly better than one
#' with fewer, using a parametric bootstrap to approximate the null
#' distribution of the likelihood-ratio statistic. This avoids the known
#' violation of standard chi-squared regularity conditions in mixture models,
#' where the null places a parameter on the boundary of the parameter space.
#'
#' Both \code{blrt()} and \code{calc_blrt()} fit the smaller- and larger-class
#' models on the observed data, compute the observed likelihood-ratio
#' statistic, then generate \code{n_reps} synthetic datasets under the smaller
#' model to build the reference distribution. \code{blrt()} is the preferred
#' name; \code{calc_blrt()} is retained for backward compatibility.
#'
#' A growth model is specified by a design matrix in time and a covariance
#' structure rather than by a measurement string, which is more than a
#' \code{measurement =} argument can carry. Pass the fitted model itself with
#' \code{from_fit =} instead: the data, the time design, the random effects and
#' the covariance constraints are all read off it, so the null and alternative
#' models differ from the fit in the number of classes and in nothing else.
#' That is the condition the test requires: likelihood-ratio tests "compare
#' models that differ only in the number of classes ... but are not appropriate
#' for comparing models that allow for different types of between-class
#' differences" (Ram & Grimm, 2009, p. 571). Passing the fit guarantees it by
#' construction.
#'
#' @param indicators Matrix or data frame of measurement items. (\code{X} is
#'   accepted as a deprecated alias.) Not needed when \code{from_fit} is given.
#' @param from_fit A model fitted by \code{\link{fit_gmm}} or
#'   \code{\link{fit_lcga}}, whose data and specification are used for both
#'   models and every replicate. \code{k_small} and \code{k_large} still say
#'   which class counts to compare; everything else comes from the fit.
#' @param k_small Number of classes in the smaller (null) model.
#' @param k_large Number of classes in the larger (alternative) model. Must be
#'   strictly greater than \code{k_small}.
#' @param measurement Measurement specification, as in \code{\link{fit_mixture}}
#'   (a single type string or a named list for mixed-type indicators).
#' @param n_reps Number of bootstrap replications. Default \code{100},
#'   following Dziak et al. (2014) and the general advice in Davison and
#'   Hinkley (1997, p. 143) that the number of replicates be at least 99. It is
#'   a floor rather than a target: the attainable p-values are
#'   \eqn{1/(B+1), 2/(B+1), \ldots}, so a decision that turns on the third
#'   decimal needs \code{n_reps = 999}, which recovers about 0.95 of the power
#'   of the full test where 99 draws recover about 0.83 (and only about 0.60 at
#'   \eqn{\alpha = .01}). See \code{vignette("estimation")}.
#' @param n_init_base Random restarts when fitting the observed-data models.
#'   Default \code{20}, as elsewhere in the package.
#' @param n_init_boot Random restarts per bootstrap replicate. Default
#'   \code{10}. This is a compute compromise rather than a recommended value:
#'   the two models are refitted \code{2 * n_reps} times, so the replicate
#'   search is where the cost of the test lives. Dziak et al. (2014) used 50
#'   and note that too few restarts under the alternative can make the
#'   likelihood ratio come out negative. \code{blrt()} counts those draws and
#'   warns when there are any; if it does, raise this to \code{50}.
#' @param verbose Logical; print progress while bootstrapping. Default
#'   \code{TRUE}.
#' @param ... Additional arguments passed to the fitting engine.
#' @param X Deprecated alias for \code{indicators}.
#'
#' @return An object of class \code{blrt_test}: a list with \code{p_value},
#'   \code{obs_diff} (the observed \eqn{2\,\Delta\ell} statistic),
#'   \code{null_dist} (the bootstrap null distribution), \code{n_negative} (how
#'   many draws produced a negative statistic), and the compared class counts.
#'   It has \code{print} and \code{plot} methods.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' X <- matrix(rbinom(500, 1, 0.5), nrow = 100)
#' res <- blrt(X, k_small = 2, k_large = 3, measurement = "binary", n_reps = 50)
#' res                 # clean printed summary
#' plot(res)           # null distribution with the observed statistic marked
#' res$p_value
#' }
#'
#' @references
#' Davison, A. C., & Hinkley, D. V. (1997). \emph{Bootstrap Methods and Their
#' Application} (ch. 4). Cambridge University Press.
#'
#' Dziak, J. J., Lanza, S. T., & Tan, X. (2014). Effect size, statistical power
#' and sample size requirements for the bootstrap likelihood ratio test in
#' latent class analysis. \emph{Structural Equation Modeling}, \emph{21}(4),
#' 534-552. \doi{10.1080/10705511.2014.919819}
#'
#' McLachlan, G. J. (1987). On bootstrapping the likelihood ratio test statistic
#' for the number of components in a normal mixture. \emph{Applied Statistics},
#' \emph{36}(3), 318-324. \doi{10.2307/2347790}
#'
#' Ram, N., & Grimm, K. J. (2009). Growth mixture modeling: A method for
#' identifying differences in longitudinal change among unobserved groups.
#' \emph{International Journal of Behavioral Development}, \emph{33}(6),
#' 565-576. \doi{10.1177/0165025409343765}
#'
#' Nylund, K. L., Asparouhov, T., & Muthen, B. O. (2007). Deciding on the number
#' of classes in latent class analysis and growth mixture modeling: A Monte
#' Carlo simulation study. \emph{Structural Equation Modeling}, \emph{14}(4),
#' 535-569. \doi{10.1080/10705510701575396}
#'
#' @export
blrt <- function(indicators, k_small, k_large, measurement = "binary",
                 n_reps = 100, n_init_base = 20, n_init_boot = 10,
                 verbose = TRUE, ..., from_fit = NULL, X = NULL) {

  supplied <- !missing(indicators)
  if (!is.null(X) && !supplied) {
    indicators <- X
    supplied   <- TRUE
  }

  if (k_small >= k_large)
    stop(sprintf(
      "k_large (%d) must be strictly greater than k_small (%d) for the BLRT.",
      k_large, k_small), call. = FALSE)

  # A growth model's specification arrives as an object rather than a string.
  # The engine arguments are prepended to `...` rather than appended so that
  # anything the caller passes explicitly still wins.
  extra <- list(...)
  if (!is.null(from_fit)) {
    spec        <- .blrt_growth_spec(from_fit)
    if (!supplied) indicators <- spec$indicators
    measurement <- spec$measurement
    extra       <- utils::modifyList(spec$args, extra)
  }

  # Resolve single-type or mixed-type measurement once, so the observed fits and
  # every bootstrap replicate share the same indicators and descriptor.
  mm          <- .normalize_measurement(measurement, indicators)
  Xd          <- mm$indicators
  measurement <- mm$descriptor

  # One place that calls the engine, so the observed fits and the replicates
  # cannot drift apart in what they were given.
  fit_engine <- function(data, k, ...)
    do.call(fit_mixture_internal,
            c(list(X = data, n_components = k, measurement = measurement),
              list(...), extra))

  if (verbose)
    message(sprintf("BLRT: comparing %d vs %d classes with %d bootstrap draws...",
                    k_small, k_large, n_reps))

  null_model <- fit_engine(Xd, k_small, n_init = n_init_base)
  alt_model  <- fit_engine(Xd, k_large, n_init = n_init_base)

  obs_diff  <- 2 * (alt_model$metrics$ll - null_model$metrics$ll)
  null_dist <- numeric(n_reps)
  N         <- nrow(Xd)

  for (i in seq_len(n_reps)) {
    classes <- sample(seq_len(k_small), size = N, replace = TRUE,
                      prob = null_model$weights)
    X_gen   <- generate_synthetic_data(null_model$mm, classes, N)

    # refine = FALSE: replicates only need the likelihood ratio, not polished
    # estimates, which makes each draw far cheaper without affecting the p-value.
    m_null_gen <- fit_engine(X_gen, k_small, n_init = n_init_boot,
                             refine = FALSE)
    m_alt_gen  <- fit_engine(X_gen, k_large, n_init = n_init_boot,
                             refine = FALSE)

    null_dist[i] <- 2 * (m_alt_gen$metrics$ll - m_null_gen$metrics$ll)

    if (verbose && (i %% max(1L, n_reps %/% 10L) == 0L))
      message(sprintf("  %d / %d draws complete", i, n_reps))
  }

  p_val <- (sum(null_dist >= obs_diff) + 1) / (n_reps + 1)

  .blrt_check_granularity(p_val, n_reps)
  n_negative <- sum(null_dist < 0)
  .blrt_check_negative(n_negative, n_reps, k_small, k_large)

  result <- list(
    p_value    = p_val,
    obs_diff   = obs_diff,
    null_dist  = null_dist,
    k_small    = k_small,
    k_large    = k_large,
    n_reps     = n_reps,
    n_negative = n_negative,
    ll_small   = null_model$metrics$ll,
    ll_large   = alt_model$metrics$ll
  )
  class(result) <- c("blrt_test", "list")
  result
}

# A Monte Carlo p-value can only take the values 1/(B+1), 2/(B+1), ..., 1, so at
# 100 draws it moves in steps of about 0.01 and cannot separate .04 from .06.
# The granularity only matters where the decision turns on it, which is within
# one step of the conventional level.
#
# Silent from 999 draws upward, where the remedy would be the number already in
# use: a warning whose advice is "do what you did" is noise.
.blrt_check_granularity <- function(p_val, n_reps) {
  if (n_reps >= 999L) return(invisible(NULL))
  if (abs(p_val - 0.05) > 1 / (n_reps + 1)) return(invisible(NULL))
  warning(sprintf(paste0(
    "The bootstrap p-value (%.4f) is within one Monte Carlo step of .05, and ",
    "with %d draws the attainable p-values are 1/%d, 2/%d, ... - the test ",
    "cannot tell .04 from .06 here. Refit with `n_reps = 999`, which recovers ",
    "about 0.95 of the power of the full test where 99 draws recover about ",
    "0.83."),
    p_val, n_reps, n_reps + 1L, n_reps + 1L), call. = FALSE)
  invisible(NULL)
}

# A replicate whose k-class fit lands below the (k-1)-class model nested inside
# it is not a feature of the null distribution; it is a search that stopped
# short on the larger model. The statistic is biased downward on that draw, and
# the p-value with it.
.blrt_check_negative <- function(n_negative, n_reps, k_small, k_large) {
  if (n_negative <= 0L) return(invisible(NULL))
  warning(sprintf(paste0(
    "%d of %d bootstrap draws produced a negative likelihood-ratio statistic, ",
    "which means the %d-class replicate fitted worse than the %d-class model ",
    "nested inside it. That is a local maximum in the replicate search rather ",
    "than a null distribution: refit with `n_init_boot = 50`. Until then the ",
    "p-value is biased and should not be reported."),
    n_negative, n_reps, k_large, k_small), call. = FALSE)
  invisible(NULL)
}

#' @rdname blrt
#' @export
calc_blrt <- function(X, k_small, k_large, measurement = "binary",
                      n_reps = 100, n_init_base = 20, n_init_boot = 10,
                      verbose = TRUE, ...) {
  blrt(indicators = X, k_small = k_small, k_large = k_large,
       measurement = measurement, n_reps = n_reps,
       n_init_base = n_init_base, n_init_boot = n_init_boot,
       verbose = verbose, ...)
}

#' @export
print.blrt_test <- function(x, ...) {
  decision <- if (x$p_value < 0.05)
    sprintf("the %d-class model fits significantly better", x$k_large)
  else
    sprintf("no significant improvement over the %d-class model", x$k_small)

  p_fmt <- if (x$p_value < 0.001) "< .001" else sprintf("= %.3f", x$p_value)

  cat("=========================================================\n")
  cat("       BOOTSTRAP LIKELIHOOD RATIO TEST (BLRT)            \n")
  cat("=========================================================\n")
  cat(sprintf("  Comparison       : %d vs %d classes\n", x$k_small, x$k_large))
  cat(sprintf("  Log-Likelihoods  : %.2f (%d-class) -> %.2f (%d-class)\n",
              x$ll_small, x$k_small, x$ll_large, x$k_large))
  cat(sprintf("  LR statistic     : %.2f\n", x$obs_diff))
  cat(sprintf("  Bootstrap draws  : %d\n", x$n_reps))
  # Only when there are any: on a clean run this line would be a zero the reader
  # has to interpret.
  if (isTRUE(x$n_negative > 0L))
    cat(sprintf("  Negative draws   : %d (replicate search too shallow)\n",
                x$n_negative))
  cat("---------------------------------------------------------\n")
  cat(sprintf("  Bootstrap p      : p %s\n", p_fmt))
  cat(sprintf("  Conclusion       : %s.\n", decision))
  cat("=========================================================\n")
  invisible(x)
}

#' @importFrom graphics hist abline legend
#' @export
plot.blrt_test <- function(x, ...) {
  hist(x$null_dist,
       breaks = "FD",
       main   = sprintf("BLRT null distribution (%d vs %d classes)",
                        x$k_small, x$k_large),
       xlab   = "2 * (LL_large - LL_small) under the null",
       col    = "grey85", border = "white",
       xlim   = range(c(x$null_dist, x$obs_diff)))
  abline(v = x$obs_diff, col = "red", lwd = 2)
  legend("topright",
         legend = c(sprintf("Observed = %.2f", x$obs_diff),
                    sprintf("p %s", if (x$p_value < 0.001) "< .001"
                            else sprintf("= %.3f", x$p_value))),
         bty = "n")
  invisible(x)
}
