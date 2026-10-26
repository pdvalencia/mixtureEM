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
# so the class trajectory on the link scale is Lambda %*% beta_k. The
# threshold parameterisation writes the same model with a threshold held equal
# across classes and the intercept growth factor fixed at zero in the last
# class; that is a reparameterisation of free per-class intercepts, not a
# restriction, and both forms have p*K free growth parameters. Estimates
# translate as
#
#     intercept_k = mean(I)_k - threshold,     slope_k = mean(S)_k

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
#' @param family Character. `"binomial"`, `"gaussian"`, `"poisson"`, `"zip"`
#'   or `"negbin"`.
#' @param dispersion_structure Character. For `family = "negbin"`, whether the
#'   overdispersion is one value per occasion shared by the classes
#'   (`"occasion"`) or one per class constant over occasions (`"class"`).
#'   Ignored by the other families.
#' @param inflation_design For `family = "zip"`, the design of the
#'   zero-inflation logit, one row per occasion. `NULL`, the default, is a
#'   single constant column. Ignored by the other families.
#' @param ... Additional arguments, ignored.
#'
#' @return A list object of class \code{c("lcga", "emission")}.
#' @export
lcga_model <- function(n_components, design, family = "binomial",
                       dispersion_structure = "occasion",
                       inflation_design = NULL, ...) {
  state <- list(
    n_components = n_components,
    design       = design,
    family       = family,
    fam          = .wglm_family(family),
    # The zero-inflation logit's own design, one row per occasion: a constant
    # unless fit_lcga() passes a polynomial, so by default each class has one
    # probability of a structural zero at every occasion.
    inflation_design = if (identical(family, "zip"))
      inflation_design %||%
        matrix(1, nrow(design), 1L, dimnames = list(NULL, "intercept")),
    # Which cells of the K x T overdispersion matrix share one free value.
    dispersion_structure = if (identical(family, "negbin"))
      match.arg(dispersion_structure, c("occasion", "class")),
    # The Gaussian family's EM is accelerated (.em_squarem(), R/em_core.R) when
    # nothing but the measurement model is estimated: a start that wanders to
    # a poor maximum took hundreds of iterations to get there. On simulated
    # stress sets it reaches plain EM's maximum and parameters (to 1e-8); a
    # class regression was not among them, so that model keeps plain EM.
    # Binomial and Poisson keep it too: on their stress sets the accelerated
    # search reached the same maxima but replicated fewer starts where a
    # class sits near probability 0 or 1, along the flat ridge such a class
    # has.
    em_squarem   = if (identical(family, "gaussian")) "unconditional",
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

# What the priors are centred on, pooled over every observed record: the mean
# response `m` and its variance `s2` (divisor N, not N - 1), plus `occ`, the
# occasions observed at least once, which are the design points the binomial
# and Poisson pseudo-rows sit on. Weighted when a survey design is in play, for
# the same reason refine_lbfgs() weights its Bernoulli marginal: the same data
# supplied as weighted patterns and as expanded rows must give one answer.
.lcga_marginal <- function(X, weights = NULL) {
  obs <- !is.na(X)
  w   <- if (is.null(weights)) obs * 1 else obs * weights
  X0  <- replace(X, !obs, 0)
  sw  <- max(sum(w), 1e-12)
  m   <- sum(X0 * w) / sw
  list(m = m, s2 = sum(w * (X0 - m)^2 * obs) / sw, occ = colSums(obs) > 0)
}

# The binomial and Poisson priors as GLM pseudo-rows: one row per observed
# occasion, holding the pooled mean m. A binomial row of weight alpha / (K U)
# adds (alpha / (K U)) [m log mu + (1 - m) log(1 - mu)] to a class's objective,
# the Dirichlet prior of Galindo-Garre & Vermunt (2006) spread evenly over the U
# design points, so a class gains alpha / K pseudo-cases in all, whatever the
# number of occasions. A Poisson row of weight alpha / (K U m) adds
# (alpha / (K U)) [log lambda - lambda / m], the gamma prior that adds alpha / K
# events to each class at the pooled rate. The Gaussian family puts no prior on
# its coefficients (its prior is on the residual variance, m_step.lcga()), and a
# constant of zero switches either prior off.
.lcga_prior_rows <- function(mm, cen) {
  D <- mm$design[cen$occ, , drop = FALSE]
  U <- nrow(D)
  # A zero-inflated or overdispersed count carries the Poisson prior on its
  # rate.
  count <- mm$fam$name %in% c("poisson", "zip", "negbin")
  a <- switch(mm$fam$name,
              binomial = .bayes_alpha(mm, "categorical"),
              poisson  = , zip = , negbin = .bayes_alpha(mm, "poisson"),
              0)
  if (a <= 0 || U == 0L || (count && cen$m <= 0))
    return(list(D = D[0L, , drop = FALSE], y = numeric(0), w = numeric(0)))
  w <- a / (mm$n_components * U)
  if (count) w <- w / cen$m
  list(D = D, y = rep(cen$m, U), w = rep(w, U))
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
  # On the pattern table each row stands for its count of cases
  # (`row_counts`, set where fit_mixture() collapsed the data), so every
  # pooled quantity below counts it that often.
  rc  <- model_state$row_counts
  rc  <- if (length(rc) == nrow(X)) rep(rc, ncol(X)) else rep(1, length(lng$y))

  pooled <- .wglm_fit(lng$D, lng$y, as.numeric(lng$observed) * rc, fam)
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

  # The inflation starts at the zeros the pooled Poisson cannot explain: the
  # observed share of zeros less the Poisson's own, as a share of the cases
  # that are not Poisson zeros, kept off 0 and 1 and perturbed per class.
  yo <- lng$y[lng$observed]
  wo <- rc[lng$observed]
  if (isTRUE(fam$has_inflation)) {
    p0  <- exp(-stats::weighted.mean(yo, wo))
    exc <- (stats::weighted.mean(yo == 0, wo) - p0) / (1 - p0)
    g0  <- stats::qlogis(min(max(exc, 0.05), 0.95))
    q   <- ncol(model_state$inflation_design)
    model_state$parameters$inflation <-
      cbind(g0 + stats::rnorm(K), matrix(0, K, q - 1L))
  }

  # The overdispersion starts at the pooled moment estimate (var - mean) /
  # mean^2, which the mixture of class means inflates; kept off the Poisson
  # boundary, from which the one-dimensional search could still leave.
  if (isTRUE(fam$has_overdispersion)) {
    mo <- stats::weighted.mean(yo, wo)
    m  <- max(mo, 1e-2)
    v  <- sum(wo * (yo - mo)^2) / (sum(wo) - 1)
    a0 <- min(max((v - m) / m^2, 0.05), 10)
    model_state$parameters$overdispersion <- matrix(a0, K, lng$Tn)
  }
  model_state
}

# M-step: one posterior-weighted GLM per class, with each family's prior taken
# from its own `bayes_constants` entry: `categorical` for binomial, `poisson`
# for Poisson (both as pseudo-rows, .lcga_prior_rows()), and `variances` for the
# Gaussian residual variance. The variance prior is the inverse-gamma term
# -(1/2)(alpha / K)[log sigma2 + s2 / sigma2], s2 the pooled variance of the
# outcome: alpha / K pseudo-observations one pooled standard deviation from the
# class mean, so the update is (RSS + (alpha / K) s2) / (N + alpha / K) and a
# class cannot shrink its variance to zero.
#' @exportS3Method
m_step.lcga <- function(model_state, X, resp, weights = NULL, ...) {
  if (!is.null(weights)) resp <- sweep(resp, 1, weights, "*")

  K   <- model_state$n_components
  fam <- model_state$fam
  if (fam$name == "gaussian")
    return(.lcga_gaussian_m_step(model_state, X, resp, weights))
  if (fam$name %in% c("binomial", "poisson", "zip"))
    return(.lcga_count_m_step(model_state, X, resp, weights))
  # What is left is the negative binomial: ECM. Each class's rate at the
  # current alpha, then alpha at the new rates once every class has moved. The
  # rate's kernel y eta - (y + 1/alpha) log(1 + alpha mu) is linear in the
  # response, and every record at an occasion shares that occasion's alpha
  # and design row, so the rate is fitted to each occasion's total weight and
  # weighted mean, as the other families' are. The prior rows are Poisson
  # rows (alpha 0), so the count prior is the Poisson family's.
  lng   <- .lcga_long(model_state, X)
  prior <- .lcga_prior_rows(model_state, .lcga_marginal(X, weights))
  y_obs <- replace(lng$y, !lng$observed, 0)
  coefs <- model_state$parameters$coefs
  od    <- model_state$parameters$overdispersion
  obs   <- !is.na(X)
  W     <- crossprod(obs * 1, resp)                 # T x K weight by occasion
  S     <- crossprod(replace(X, !obs, 0), resp)     # T x K weighted sums
  D     <- model_state$design

  for (k in seq_len(K))
    coefs[k, ] <- .nb_rate_fit(D     = rbind(D, prior$D),
                               y     = c(S[, k] / W[, k], prior$y),
                               w     = c(W[, k], prior$w),
                               alpha = c(od[k, ], rep(0, length(prior$y))),
                               start = coefs[k, ])

  # One alpha per occasion pools that occasion's records over the classes; one
  # per class pools the class's records over the occasions.
  occ <- rep(seq_len(lng$Tn), each = lng$n)
  MU  <- exp(lng$D %*% t(coefs))                          # (n T) x K
  W   <- resp[rep(seq_len(lng$n), lng$Tn), , drop = FALSE] * lng$observed
  if (identical(model_state$dispersion_structure, "class")) {
    for (k in seq_len(K))
      od[k, ] <- .nb_alpha_fit(y_obs, MU[, k], W[, k])
  } else {
    for (t in seq_len(lng$Tn)) {
      r <- occ == t
      od[, t] <- .nb_alpha_fit(rep(y_obs[r], K), as.vector(MU[r, ]),
                               as.vector(W[r, ]))
    }
  }

  model_state$parameters$coefs      <- coefs
  model_state$parameters$dispersion <- rep(1, K)
  model_state$parameters$overdispersion <- od
  model_state
}

# The Gaussian M-step of m_step.lcga(), on the n x T matrix itself rather than
# its occasion-stacked records. Every case shares the design row of its
# occasion, so the weighted least-squares normal equations D'W(y - D b) = 0 of a
# class depend on the data only through its total weight and weighted response
# sum at each occasion: two T x K cross-products for all the classes at once,
# then a solve on T rows. The trajectory is the one the records give (to
# rounding), and the cost no longer grows with n beyond the cross-products. The
# residual variance needs the spread within each occasion, which those sums do
# not carry, so it is read off the cases. `resp` arrives already multiplied by
# the case weights (m_step.lcga()); `weights` is read only for the prior's
# pooled variance.
.lcga_gaussian_m_step <- function(model_state, X, resp, weights = NULL) {
  K     <- model_state$n_components
  n     <- nrow(X)
  D     <- model_state$design
  a_var <- .bayes_alpha(model_state, "variances") / K
  cen   <- if (a_var > 0) .lcga_marginal(X, weights)

  obs <- !is.na(X)
  X0  <- replace(X, !obs, 0)
  W   <- crossprod(obs * 1, resp)                   # T x K weight by occasion
  S   <- crossprod(X0, resp)                        # T x K weighted sums

  coefs      <- model_state$parameters$coefs
  dispersion <- model_state$parameters$dispersion
  for (k in seq_len(K)) {
    ok <- is.finite(W[, k]) & W[, k] > 0
    # .wglm_fit()'s answer for a class with no weight anywhere: unchanged.
    if (!any(ok)) { dispersion[k] <- 1; next }
    # The rank tolerance glm.fit() uses at .wglm_fit()'s epsilon of 1e-10.
    beta <- stats::lm.wfit(D[ok, , drop = FALSE], S[ok, k] / W[ok, k],
                           W[ok, k], tol = 1e-13)$coefficients
    beta[!is.finite(beta)] <- 0
    names(beta) <- NULL
    coefs[k, ] <- beta

    E <- X0 - rep(as.vector(D %*% beta), each = n)
    E[!obs] <- 0
    rss <- sum(resp[, k] * rowSums(E * E))
    nk  <- sum(W[, k])
    dispersion[k] <- if (a_var > 0)
      max((rss + a_var * cen$s2) / (nk + a_var), 1e-6)
    else max(rss / nk, 1e-6)
  }
  model_state$parameters$coefs      <- coefs
  model_state$parameters$dispersion <- dispersion
  model_state
}

# The binomial, Poisson and zero-inflated Poisson M-step of m_step.lcga(), on
# the occasion sums for the reason the Gaussian one is: a class's weighted
# log-likelihood depends on the data only through its total weight W and
# weighted response sum S at each occasion, two T x K cross-products for all
# the classes at once. The prior's pseudo-rows (.lcga_prior_rows()) sit on the
# observed occasions with response m, so they add their weight to W and
# weight x m to S there. Each class is then solved by Newton on T rows
# (.lcga_count_newton()), with .wglm_fit() on the same rows when Newton cannot
# proceed.
#
# Zero inflation keeps the inner E-step of the per-record form: the share z
# of a zero that is structural depends only on its class and occasion, so with
# Z0 the weight on observed zeros at each occasion the rate is fitted on
# weight W - z Z0 (a count's events are all Poisson) and the inflation logit is
# a binomial with z Z0 successes in W trials, on the inflation design's rows.
.lcga_count_m_step <- function(model_state, X, resp, weights = NULL) {
  K   <- model_state$n_components
  fam <- model_state$fam
  D   <- model_state$design
  cen   <- .lcga_marginal(X, weights)
  prior <- .lcga_prior_rows(model_state, cen)

  obs <- !is.na(X)
  X0  <- replace(X, !obs, 0)
  W   <- crossprod(obs * 1, resp)                   # T x K weight by occasion
  S   <- crossprod(X0, resp)                        # T x K weighted sums
  Wp  <- matrix(0, nrow(W), K)                      # the prior's share of W
  if (length(prior$w)) {
    occ <- which(cen$occ)
    Wp[occ, ] <- prior$w
    S[occ, ] <- S[occ, ] + prior$w * prior$y
  }

  coefs <- model_state$parameters$coefs
  infl  <- model_state$parameters$inflation
  zip   <- isTRUE(fam$has_inflation)
  if (zip) {
    Z0 <- crossprod((obs & X == 0) * 1, resp)       # T x K weight on zeros
    Z0[is.na(Z0)] <- 0
    Dz <- model_state$inflation_design
  }
  for (k in seq_len(K)) {
    Wk <- W[, k]
    if (!any(is.finite(Wk) & Wk > 0)) next
    if (zip) {
      pi_k <- stats::plogis(as.vector(Dz %*% infl[k, ]))
      mu_k <- exp(as.vector(D %*% coefs[k, ]))
      z    <- pi_k / (pi_k + (1 - pi_k) * exp(-mu_k))
      zk   <- z * Z0[, k]
      g    <- .lcga_count_newton(Dz, zk, Wk, "binomial", infl[k, ])
      if (is.null(g)) {
        ok <- Wk > 0
        g  <- .wglm_fit(Dz[ok, , drop = FALSE], zk[ok] / Wk[ok], Wk[ok],
                        .wglm_family("binomial"), start = infl[k, ])$coefficients
      }
      infl[k, ] <- g
      Wk <- Wk - zk
    }
    Wk <- Wk + Wp[, k]
    beta <- .lcga_count_newton(D, S[, k], Wk, if (zip) "poisson" else fam$name,
                               coefs[k, ])
    if (is.null(beta)) {
      ok   <- Wk > 0
      beta <- .wglm_fit(D[ok, , drop = FALSE], S[ok, k] / Wk[ok], Wk[ok],
                        fam, start = coefs[k, ])$coefficients
    }
    coefs[k, ] <- beta
  }
  model_state$parameters$coefs      <- coefs
  model_state$parameters$dispersion <- rep(1, K)
  model_state$parameters$inflation  <- infl
  model_state
}

# The probability of a structural zero, K x T.
.lcga_inflation <- function(mm)
  t(stats::plogis(mm$inflation_design %*% t(mm$parameters$inflation)))

#' @exportS3Method
log_likelihood.lcga <- function(model_state, X, ...) {
  n   <- nrow(X)
  Tn  <- ncol(X)
  K   <- model_state$n_components
  fam <- model_state$fam

  eta <- model_state$design %*% t(model_state$parameters$coefs)   # T x K
  log_eps <- matrix(0, nrow = n, ncol = K)

  # Gaussian: the normal log density summed over a case's observed occasions,
  # written out rather than evaluated cell by cell through dnorm(); the same
  # value to rounding, at a fraction of the cost on a large n.
  if (fam$name == "gaussian") {
    obs   <- !is.na(X)
    X0    <- replace(X, !obs, 0)
    n_obs <- rowSums(obs)
    for (k in seq_len(K)) {
      s2 <- model_state$parameters$dispersion[k]
      E  <- X0 - rep(eta[, k], each = n)
      E[!obs] <- 0
      log_eps[, k] <- -0.5 * (n_obs * log(2 * base::pi * s2) + rowSums(E * E) / s2)
    }
    return(log_eps)
  }
  todo <- seq_len(K)
  # Binomial, Poisson and zero-inflated Poisson the same way: a case's sum
  # over its observed occasions is a few matrix products for all the classes,
  # y'log mu + (o - y)'log(1 - mu), or y'log mu - o'mu less the case's own
  # log y! terms (the zero-inflated split below).
  # The binomial keeps its two terms apart rather than writing y'(log mu -
  # log(1 - mu)) + o'log(1 - mu): near mu = 1 that difference is a large
  # number minus itself, and the rounding it leaves is what the numerical
  # Hessian of the standard errors then differentiates. A class with a term
  # that is not finite takes the cell-by-cell path below, whose guard
  # handles it.
  if (fam$name %in% c("binomial", "poisson", "zip")) {
    obs <- !is.na(X)
    X0  <- replace(X, !obs, 0)
    mu  <- fam$linkinv(eta)
    if (fam$name == "binomial") {
      mu <- pmin(pmax(mu, 1e-15), 1 - 1e-15)
      log_eps <- X0 %*% log(mu) + (obs - X0) %*% log1p(-mu)
    } else if (fam$name == "poisson") {
      mu <- pmax(mu, 1e-15)
      log_eps <- X0 %*% log(mu) - obs %*% mu - rowSums(lgamma(X0 + 1))
    } else {
      # Zero-inflated: a positive count is log(1 - pi) plus its Poisson term,
      # a zero log(pi + (1 - pi) e^-mu), with pi held off 0 and 1 as below.
      mu  <- pmax(mu, 1e-15)
      pz  <- t(pmin(pmax(.lcga_inflation(model_state), 1e-15), 1 - 1e-15))
      pos <- obs & X0 > 0
      log_eps <- X0 %*% log(mu) + pos %*% (log1p(-pz) - mu) -
        rowSums(lgamma(X0 + 1)) + (obs & X0 == 0) %*% log(pz + (1 - pz) * exp(-mu))
    }
    todo <- which(colSums(!is.finite(log_eps)) > 0)
    if (!length(todo)) return(unname(log_eps))
    log_eps[, todo] <- 0
  }
  # Held off 0 and 1 for the same reason as a binomial mean: a class whose
  # inflation reached 1 would otherwise give a positive count the -Inf that the
  # FIML guard below zeroes, and so claim a perfect fit to data it excludes.
  if (isTRUE(fam$has_inflation))
    pi <- pmin(pmax(.lcga_inflation(model_state), 1e-15), 1 - 1e-15)

  for (k in todo) {
    mu <- fam$linkinv(eta[, k])
    if (fam$name == "binomial") mu <- pmin(pmax(mu, 1e-15), 1 - 1e-15)
    if (fam$name %in% c("poisson", "zip", "negbin")) mu <- pmax(mu, 1e-15)
    as_cells <- function(v) matrix(v, nrow = n, ncol = Tn, byrow = TRUE)
    ll <- if (isTRUE(fam$has_inflation))
      fam$log_dens(X, as_cells(mu), 1, as_cells(pi[k, ]))
    else if (isTRUE(fam$has_overdispersion))
      fam$log_dens(X, as_cells(mu), 1,
                   as_cells(model_state$parameters$overdispersion[k, ]))
    else fam$log_dens(X, as_cells(mu), model_state$parameters$dispersion[k])
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
  np + length(model_state$parameters$inflation) + .lcga_n_overdispersion(model_state)
}

# The free overdispersion values: one per occasion or one per class.
.lcga_n_overdispersion <- function(mm) {
  od <- mm$parameters$overdispersion
  if (is.null(od)) return(0L)
  if (identical(mm$dispersion_structure, "class")) nrow(od) else ncol(od)
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

  if (family %in% c("poisson", "zip", "negbin")) {
    if (any(vals < 0))
      stop(sprintf('family = "%s" requires non-negative counts, but negative ',
                   family), 'values are present.', call. = FALSE)
    # Counts read from a text file arrive as doubles, so what matters is that a
    # value is a whole number, not that it is stored as an integer.
    if (any(abs(vals - round(vals)) > 1e-8))
      stop(sprintf('family = "%s" requires whole-number counts, but fractional ',
                   family),
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
#'     special case of a growth mixture model. A class with more than 90% of
#'     its scores, over all occasions, on one value (a floor of zeros, say) has
#'     no residual spread left; a start that ends there is set aside, as in
#'     every model the package fits, and a warning is given only if every start
#'     did.
#'   * `"poisson"` — a count outcome, log link. The trajectory is a curve in the
#'     event rate. Counts must be whole and non-negative. A Poisson count's
#'     variance equals its mean; for counts that spread more, use `"negbin"`.
#'   * `"zip"` — a zero-inflated count, log link. A zero has two sources: a
#'     structural zero, from a case that is not at risk of the event at that
#'     occasion, and a Poisson zero, from a case that is at risk but had no
#'     event (Lambert, 1992; Nagin & Land, 1993, call the first intermittency).
#'     Each class has its own trajectory in the rate of the at-risk part and
#'     its own probability of a structural zero, estimated on the logit scale:
#'     constant over occasions, or linear in time with `inflation_degree = 1`.
#'     Excess zeros fitted with `"poisson"`
#'     have nowhere to go but extra classes, so where zeros pile up beyond
#'     what the rate implies, compare the two families before reading the
#'     class count.
#'   * `"negbin"` — an overdispersed count, log link: the negative binomial
#'     in its NB2 form, variance mu + alpha mu^2 (Cameron & Trivedi, 2013,
#'     ch. 3), so alpha = 0 is the Poisson. Each class has its own trajectory
#'     in the mean count; how alpha is shared is set by `dispersion`. A
#'     Poisson fit can only absorb extra spread by adding classes, so on
#'     overdispersed counts compare the two families before reading the class
#'     count.
#' @param dispersion Character, `family = "negbin"` only. `"occasion"` (the
#'   default) estimates one alpha per occasion, shared by the classes;
#'   `"class"` estimates one alpha per class, constant over occasions.
#' @param inflation_degree `family = "zip"` only. Degree of the polynomial in
#'   time for the logit of the structural-zero probability, in each class: `0`
#'   (the default) holds the probability constant over occasions, `1` lets it
#'   rise or fall linearly on the logit scale, which is how a share of cases
#'   that stop or start being at risk with age is written (Nagin & Land,
#'   1993). Uses the same `time_scores` as the trajectory. Higher degrees are
#'   not offered: a structural zero and a Poisson zero look alike at every
#'   occasion, so the inflation curve is the weakest-identified part of the
#'   model.
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
#'   `family = "gaussian"` — the per-class residual variance, or — for
#'   `family = "zip"` — the per-class inflation logits (`inflation`, one
#'   column per inflation coefficient) and the probability of a structural
#'   zero they imply (`inflation_probability`: one per class, or a K x T
#'   matrix when `inflation_degree = 1`).
#'   The fitted trajectory of a `"zip"` fit is the at-risk rate; the expected
#'   count is that rate times one minus the inflation probability. A
#'   `"negbin"` fit carries `overdispersion`, the K x T matrix of alpha
#'   (repeated across the classes or the occasions it is shared by), and
#'   `dispersion`, the structure.
#'
#' @seealso [`growth_estimates()`] for the coefficients with standard errors,
#'   [`fit_rmlca()`] for unconstrained trajectory classes and
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
#' Cameron, A. C., & Trivedi, P. K. (2013). \emph{Regression analysis of
#' count data} (2nd ed.). Cambridge University Press.
#' \doi{10.1017/CBO9781139013567}
#'
#' Lambert, D. (1992). Zero-inflated Poisson regression, with an application
#' to defects in manufacturing. \emph{Technometrics}, \emph{34}(1), 1-14.
#' \doi{10.2307/1269547}
#'
#' Nagin, D. S., & Land, K. C. (1993). Age, criminal careers, and population
#' heterogeneity: Specification and estimation of a nonparametric, mixed
#' Poisson model. \emph{Criminology}, \emph{31}(3), 327-362.
#' \doi{10.1111/j.1745-9125.1993.tb01133.x}
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
                     family = c("binomial", "gaussian", "poisson", "zip",
                                "negbin"),
                     dispersion = c("occasion", "class"),
                     inflation_degree = 0,
                     time_scores = NULL,
                     layout = c("time_major", "item_major"),
                     id = NULL, time = NULL, item = NULL,
                     time_labels = NULL,
                     predictors = NULL,
                     ...) {

  family <- match.arg(family)
  dispersion <- match.arg(dispersion)
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

  inflation_degree <- as.integer(inflation_degree)
  if (length(inflation_degree) != 1L || is.na(inflation_degree) ||
      !inflation_degree %in% 0:1)
    stop("`inflation_degree` must be 0 (constant) or 1 (linear in time).",
         call. = FALSE)
  if (inflation_degree > 0L && family != "zip")
    stop('`inflation_degree` applies to family = "zip" only.', call. = FALSE)
  # The same spare occasion as the trajectory: one inflation value per occasion
  # would leave nothing for the rate and the inflation to be told apart by.
  if (inflation_degree + 1L >= Tn)
    stop(sprintf(
      "A linear inflation logit needs at least 3 occasions but %d were given.",
      Tn), call. = FALSE)

  .lcga_check_outcome(prep$X, family)

  design <- .lcga_design(time_scores, degree)
  inflation_design <- if (family == "zip")
    .lcga_design(time_scores, inflation_degree)

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
    dispersion_structure = dispersion,
    inflation_design = inflation_design,
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
  if (family == "zip") {
    fit$growth$inflation_degree <- inflation_degree
    fit$growth$inflation_design <- inflation_design
    fit$growth$inflation <- fit$mm$parameters$inflation
    # One probability per class when it is constant; K x T when it moves.
    fit$growth$inflation_probability <- if (inflation_degree == 0L)
      .lcga_inflation(fit$mm)[, 1L] else .lcga_inflation(fit$mm)
  }
  if (family == "negbin") {
    fit$growth$overdispersion <- fit$mm$parameters$overdispersion
    fit$growth$dispersion     <- dispersion
  }
  class(fit) <- c("lcga", class(fit))
  fit$call <- match.call()
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
                 poisson  = , zip = , negbin = range(c(0, ft, orng), na.rm = TRUE),
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
              switch(g$family, binomial = "logit", poisson = , zip = , negbin = "log",
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

  if (!is.null(g$inflation_probability)) {
    ip <- g$inflation_probability
    if (is.matrix(ip)) {
      cat("\nZERO INFLATION (probability of a structural zero, by occasion)\n")
      dimnames(ip) <- list(paste("Class", seq_len(nrow(ip))), g$time_labels)
    } else {
      cat("\nZERO INFLATION (probability of a structural zero, constant over occasions)\n")
      ip <- matrix(ip, nrow = 1,
                   dimnames = list("Probability", paste("Class", seq_along(ip))))
    }
    print(round(ip, 3))
  }

  if (!is.null(g$overdispersion)) {
    od <- g$overdispersion
    if (identical(g$dispersion, "class")) {
      cat("
OVERDISPERSION (alpha, variance = mean + alpha mean^2; per class)
")
      od <- matrix(od[, 1L], nrow = 1,
                   dimnames = list("Alpha", paste("Class", seq_len(nrow(od)))))
    } else {
      cat("
OVERDISPERSION (alpha, variance = mean + alpha mean^2; per occasion, shared by the classes)
")
      od <- matrix(od[1L, ], nrow = 1, dimnames = list("Alpha", g$time_labels))
    }
    print(round(od, 3))
  }

  NextMethod()
}
