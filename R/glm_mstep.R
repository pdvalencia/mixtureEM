# ==============================================================================
# Weighted-GLM M-step (IRLS)
# ==============================================================================
#
# Shared infrastructure for every model whose M-step is "fit one regression per
# class, weighting each case by its posterior responsibility". Latent class
# growth analysis is the first consumer; mixture regression and count indicators
# with covariates use the same component.
#
# The fitting itself is delegated to stats::glm.fit(), which is the same
# Fortran/C IRLS that glm() runs: step-halving on divergence, a pivoted QR that
# survives rank deficiency, and convergence bookkeeping. Re-implementing IRLS in
# R would be slower and less robust, and would have to re-derive behaviour that
# glm.fit() already gets right. What this file adds is the part glm.fit() does
# not know about: responsibility weights, missing cells, empty classes, and the
# prior that keeps a class from separating.

# The three families, each carrying the two things the EM loop needs outside the
# M-step: a log density (for the E-step) and a link inverse (for reporting and
# for the trajectory a class implies).
#
# `dispersion` is the residual variance for the gaussian family and is ignored
# by the other two, which have none.
.wglm_family <- function(family = c("gaussian", "binomial", "poisson")) {
  family <- match.arg(family)
  switch(
    family,
    gaussian = list(
      name       = "gaussian",
      glm_family = stats::gaussian(),
      linkinv    = function(eta) eta,
      log_dens   = function(y, mu, dispersion)
        stats::dnorm(y, mean = mu, sd = sqrt(dispersion), log = TRUE),
      has_dispersion = TRUE,
      # What a class's fitted value is called when printed.
      scale_label = "Mean"
    ),
    binomial = list(
      name       = "binomial",
      glm_family = stats::binomial(),
      linkinv    = function(eta) stats::plogis(eta),
      log_dens   = function(y, mu, dispersion)
        stats::dbinom(y, size = 1, prob = mu, log = TRUE),
      has_dispersion = FALSE,
      scale_label = "Probability"
    ),
    poisson = list(
      name       = "poisson",
      glm_family = stats::poisson(),
      linkinv    = function(eta) exp(eta),
      log_dens   = function(y, mu, dispersion)
        stats::dpois(y, lambda = mu, log = TRUE),
      has_dispersion = FALSE,
      scale_label = "Event rate"
    )
  )
}

# Fit a single posterior-weighted GLM.
#
#   D      n x p design matrix (an intercept column, if wanted, is already in D)
#   y      response, length n; NA entries are dropped by zeroing their weight
#   w      non-negative case weights (responsibility x sampling weight)
#   fam    a .wglm_family() list
#   start  optional starting coefficients, normally the previous EM iterate
#
# Returns the coefficient vector, the dispersion (gaussian only), and whether
# IRLS converged. A class that has emptied out - no case still carries weight -
# keeps the coefficients it came in with rather than producing NaNs that would
# then poison the E-step for every remaining class.
.wglm_fit <- function(D, y, w, fam, start = NULL) {
  ok <- is.finite(y) & is.finite(w) & w > 0
  if (!any(ok))
    return(list(coefficients = start %||% rep(0, ncol(D)),
                dispersion = 1, converged = FALSE))

  Dk <- D[ok, , drop = FALSE]
  yk <- y[ok]
  wk <- w[ok]

  # glm.fit() warns about two things that are expected here rather than
  # symptomatic: fractional "successes", which is what a responsibility weight
  # on a 0/1 response always produces, and fitted probabilities at 0 or 1, which
  # a well-separated latent class legitimately has. Both are suppressed;
  # non-convergence is not suppressed but read back off the fitted object below,
  # since it is the one signal worth acting on.
  run <- function(start)
    suppressWarnings(tryCatch(
      stats::glm.fit(x = Dk, y = yk, weights = wk, start = start,
                     family = fam$glm_family,
                     control = list(epsilon = 1e-10, maxit = 50, trace = FALSE),
                     intercept = TRUE),
      error = function(e) NULL))

  fit <- run(start)
  # The previous EM iterate is normally an excellent start, but glm.fit()
  # refuses a start whose fitted values fall outside the family's domain, which
  # a class that moved a long way in the last iteration can produce. Its own
  # automatic start does not have that problem, so fall back to it rather than
  # abandoning the class.
  if (is.null(fit) && !is.null(start)) fit <- run(NULL)

  if (is.null(fit))
    return(list(coefficients = start %||% rep(0, ncol(D)),
                dispersion = 1, converged = FALSE))

  beta <- fit$coefficients
  # A rank-deficient design (a time score repeated, or a class with observations
  # at a single occasion) aliases columns to NA. Zero is the neutral value: the
  # aliased term drops out of the linear predictor instead of propagating NA.
  beta[!is.finite(beta)] <- 0
  names(beta) <- NULL

  dispersion <- 1
  if (isTRUE(fam$has_dispersion)) {
    mu  <- fam$linkinv(as.vector(Dk %*% beta))
    ssw <- sum(wk * (yk - mu)^2)
    # Floored for the same reason the gaussian emission floors its variances: a
    # class that fits its members exactly would otherwise claim infinite density.
    dispersion <- max(ssw / sum(wk), 1e-6)
  }

  list(coefficients = beta,
       dispersion   = dispersion,
       converged    = isTRUE(fit$converged))
}

# Pseudo-observations that keep a class from separating.
#
# The package's Bernoulli and Poisson emissions each carry a conjugate prior of
# alpha/K pseudo-observations per item, centred on that item's observed
# marginal. A GLM has no conjugate prior in that sense, but the same idea
# expresses itself exactly: append one row per design point, holding the pooled
# marginal response there, carrying alpha/K weight. A class with plenty of
# members barely notices it; a class drifting towards a fitted probability of 0
# or 1 - where the coefficient it implies is infinite and no later iteration can
# come back - is pulled towards the pooled trajectory instead.
#
# Returned in the same (design, y, w) form as the real data so the caller just
# rbind()s them.
.wglm_prior_rows <- function(design, marginal, prior_obs) {
  keep <- is.finite(marginal)
  list(D = design[keep, , drop = FALSE],
       y = marginal[keep],
       w = rep(prior_obs, sum(keep)))
}
