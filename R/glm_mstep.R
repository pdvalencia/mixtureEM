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

# The five families, each carrying the two things the EM loop needs outside the
# M-step: a log density (for the E-step) and a link inverse (for reporting and
# for the trajectory a class implies).
#
# `dispersion` is the residual variance for the gaussian family and is ignored
# by the others, which have none.
.wglm_family <- function(family = c("gaussian", "binomial", "poisson", "zip",
                                    "negbin")) {
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
    ),
    # Zero-inflated Poisson: a zero is either structural, with probability
    # `pi`, or a Poisson zero. The rate is fitted by the Poisson GLM and the
    # inflation by a logistic one (m_step.lcga()), so `glm_family` is the
    # count part's. `mu` is the Poisson rate, not the mean (1 - pi) mu.
    zip = list(
      name       = "zip",
      glm_family = stats::poisson(),
      linkinv    = function(eta) exp(eta),
      log_dens   = function(y, mu, dispersion, pi) {
        out <- log1p(-pi) + stats::dpois(y, lambda = mu, log = TRUE)
        z <- !is.na(y) & y == 0
        out[z] <- log(pi[z] + (1 - pi[z]) * exp(-mu[z]))
        out
      },
      has_dispersion = FALSE,
      has_inflation  = TRUE,
      scale_label = "Event rate"
    ),
    # Negative binomial (NB2): var = mu + alpha mu^2. The rate is fitted by
    # Fisher scoring at fixed alpha (.nb_rate_fit()) and alpha by a 1-D search
    # (.nb_alpha_fit()), both in m_step.lcga(), so `glm_family` is unused.
    negbin = list(
      name       = "negbin",
      glm_family = stats::poisson(),
      linkinv    = function(eta) exp(eta),
      log_dens   = function(y, mu, dispersion, alpha) .nb_log_dens(y, mu, alpha),
      has_dispersion     = FALSE,
      has_overdispersion = TRUE,
      scale_label = "Event rate"
    )
  )
}

# The NB2 log density, cell by cell, keeping the shape of `y`. An alpha at the
# Poisson boundary is the Poisson density itself, not a size of 1e300.
#
# The density is evaluated once per distinct (count, mean, alpha) cell and
# looked up for the rest: an E-step on thousands of records has only a few
# distinct means per class and a few dozen distinct counts, and evaluating
# dnbinom() cell by cell was half of a fit's time.
.nb_log_dens <- function(y, mu, alpha) {
  alpha <- rep_len(alpha, length(y))
  mu    <- rep_len(mu, length(y))
  out <- y
  ok  <- !is.na(y)
  if (!any(ok)) return(out)
  yo  <- y[ok]; mo <- mu[ok]; ao <- alpha[ok]
  um  <- unique(mo); ua <- unique(ao)
  key <- ((match(mo, um) - 1) * length(ua) + match(ao, ua) - 1) *
    (max(yo) + 1) + yo
  u   <- !duplicated(key)
  yu  <- yo[u]; mu_u <- mo[u]; au <- ao[u]
  du  <- yu
  po  <- au < 1e-12
  du[po]  <- stats::dpois(yu[po], lambda = mu_u[po], log = TRUE)
  du[!po] <- stats::dnbinom(yu[!po], size = 1 / au[!po], mu = mu_u[!po],
                            log = TRUE)
  out[ok] <- du[match(key, key[u])]
  out
}

# The NB2 rate regression at fixed alpha: Fisher scoring with step halving on
# the weighted kernel y log mu - (y + 1/alpha) log(1 + alpha mu). A row with
# alpha = 0 is a Poisson row, which is how the Poisson-form prior rows of
# .lcga_prior_rows() enter. Score D'w(y - mu)/(1 + alpha mu), information
# D' diag(w mu/(1 + alpha mu)) D.
.nb_rate_fit <- function(D, y, w, alpha, start) {
  ok <- is.finite(y) & is.finite(w) & w > 0
  if (!any(ok)) return(start)
  D <- D[ok, , drop = FALSE]; y <- y[ok]; w <- w[ok]; a <- alpha[ok]
  po <- a <= 0
  kern <- function(b) {
    eta <- as.vector(D %*% b)
    mu  <- exp(eta)
    sum(w * ifelse(po, y * eta - mu, y * eta - (y + 1 / a) * log1p(a * mu)))
  }
  beta <- start
  cur  <- kern(beta)
  for (it in seq_len(100)) {
    mu   <- exp(as.vector(D %*% beta))
    v    <- 1 + a * mu
    sc   <- crossprod(D, w * (y - mu) / v)
    info <- crossprod(D, D * (w * mu / v))
    step <- tryCatch(as.vector(solve(info, sc)), error = function(e) NULL)
    if (is.null(step) || !all(is.finite(step))) break
    for (h in 0:30) {
      cand <- beta + step
      new  <- kern(cand)
      if (is.finite(new) && new >= cur - 1e-12) break
      step <- step / 2
    }
    if (!is.finite(new) || new < cur - 1e-12) break
    beta <- cand
    done <- max(abs(step)) < 1e-10 || new - cur < 1e-12
    cur  <- new
    if (done) break
  }
  beta
}

# One alpha shared by a set of cells: the weighted NB2 log-likelihood
# maximised on the log scale. The lower end is the Poisson boundary.
.nb_alpha_fit <- function(y, mu, w) {
  ok <- is.finite(y) & w > 0
  if (!any(ok)) return(1e-8)
  y <- y[ok]; mu <- mu[ok]; w <- w[ok]
  # The objective is a weighted sum over (count, mean) pairs, and without
  # covariates a growth model has only K x T distinct means and few distinct
  # counts, so the pairs are pooled first: the search below evaluates the
  # density some thirty times per M-step, and on thousands of records that
  # was most of the fit's time.
  um  <- unique(mu)
  key <- (match(mu, um) - 1) * (max(y) + 1) + y
  w   <- rowsum(w, key, reorder = FALSE)[, 1]
  u   <- !duplicated(key)
  y   <- y[u]; mu <- mu[u]
  f <- function(la) sum(w * stats::dnbinom(y, size = exp(-la), mu = mu, log = TRUE))
  exp(stats::optimize(f, c(log(1e-8), log(1e4)), maximum = TRUE,
                      tol = 1e-10)$maximum)
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

# One safeguarded Newton step of the weighted binomial GLM from `start`,
# instead of the full IRLS solve above: a generalised-EM M-step in the sense
# of Dempster, Laird and Rubin (1977), which only has to *improve* the
# complete-data log-likelihood, not maximise it. Used by the random-intercept
# search when `state$ri$gem` is set. The difference matters from a wide random
# start: solving each item's GLM to convergence there drives the loadings to
# 1e12 within three iterations (RECORDS.md, "R12"), while one step follows a
# damped path that the E-step can correct on the next iteration. Step-halving
# on the aggregated log-likelihood keeps the EM guarantee that the objective
# never falls.
.wglm_newton_step <- function(D, y, w, start) {
  ok <- is.finite(y) & is.finite(w) & w > 0
  if (!any(ok)) return(list(coefficients = start, dispersion = 1, converged = FALSE))
  Dk <- D[ok, , drop = FALSE]
  yk <- y[ok]
  wk <- w[ok]
  wll <- function(b) {
    eta <- as.vector(Dk %*% b)
    sum(wk * (yk * eta - log1p(exp(eta))))
  }
  beta <- start
  cur  <- wll(beta)
  mu   <- stats::plogis(as.vector(Dk %*% beta))
  sc   <- crossprod(Dk, wk * (yk - mu))
  info <- crossprod(Dk, Dk * (wk * mu * (1 - mu)))
  step <- tryCatch(solve(info + diag(1e-8, ncol(Dk)), sc), error = function(e) NULL)
  if (!is.null(step) && all(is.finite(step))) {
    step <- as.vector(step)
    for (h in 0:12) {
      cand <- beta + step
      new  <- wll(cand)
      if (is.finite(new) && new >= cur - 1e-10) { beta <- cand; break }
      step <- step / 2
    }
  }
  names(beta) <- NULL
  list(coefficients = beta, dispersion = 1, converged = FALSE)
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
# rbind()s them. `prior_obs` is one weight for every row, or one per row of
# `design` when the rows should not carry the prior equally (a continuous
# random intercept's quadrature nodes, .lta_ri_prior_obs()).
.wglm_prior_rows <- function(design, marginal, prior_obs) {
  keep <- is.finite(marginal)
  list(D = design[keep, , drop = FALSE],
       y = marginal[keep],
       w = rep_len(prior_obs, length(marginal))[keep])
}
