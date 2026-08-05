# ==============================================================================
# Structured multivariate-normal emission
# ==============================================================================
#
# A class's implied distribution over the T repeated measures is multivariate
# normal with a *structured* covariance:
#
#     y_i | c_i = k  ~  N( Lambda (alpha_k + Gamma_k x_i) ,
#                          Lambda_r Psi_k Lambda_r' + Theta_k )
#
# where Lambda is the fixed T x p growth design (1, t, t^2, ...), Lambda_r is the
# subset of its columns carrying a random effect, Psi_k is the q x q covariance
# of those random effects and Theta_k is a diagonal matrix of residual
# variances. Gamma_k regresses the growth factors on case-level covariates
# x_i, and is absent (m = 0) in the unconditional model,
# where alpha_k is simply the vector of growth-factor means.
#
# The crucial fact, and the reason this is worth building: for a *continuous*
# outcome the random effects integrate out in closed form. There is no numerical
# integration anywhere below. The mixture E-step is untouched — it still asks
# each class for a per-case log density — and the M-step works on expected
# sufficient statistics obtained from the conditional distribution of the random
# effects given the data, which is itself normal.
#
# One emission therefore covers a family of models usually presented as quite
# different specifications:
#
#   q = 0                        latent class growth analysis, continuous
#                                (Theta free per occasion and/or per class)
#   q = 1, intercept only        random-intercept growth mixture model
#   q = 2, intercept and slope   the standard growth mixture model
#
# and the same component is what a factor mixture model will need later, with
# Lambda estimated rather than fixed.
#
# Constraints carried in the state rather than in the parameters:
#
#   r_cols                  which columns of Lambda are random (possibly none)
#   psi_equal               Psi held equal across classes (the usual default)
#   theta_equal             Theta held equal across classes (the usual default)
#   theta_shared_occasions  one residual variance for all occasions, rather
#                           than one per occasion
#   gamma_equal             the growth-factor regressions held equal across
#                           classes (the usual default)
#
# Parameters are always stored per class, including where a constraint makes the
# classes share a value; the constraint is imposed in the M-step and counted in
# n_parameters(). Keeping the storage rectangular means everything downstream —
# sorting, printing, plotting, the BLRT generator — sees one layout.

# ------------------------------------------------------------------------------
# Constructor
# ------------------------------------------------------------------------------

#' Constructor for structured multivariate-normal emissions
#'
#' @description
#' Sets up the emission state for a growth mixture model: each class follows its
#' own polynomial trajectory in the growth-factor means, and cases vary about
#' their class's trajectory through random effects with a class-specific or
#' class-common covariance.
#'
#' @param n_components Integer. The number of latent classes to estimate.
#' @param design Numeric matrix with one row per occasion and one column per
#'   growth coefficient, as built by the polynomial in the time scores.
#' @param random_effects Which growth factors vary within a class: `"none"`,
#'   `"intercept"`, `"intercept_slope"`, or `"all"`.
#' @param psi Whether the growth-factor covariance is held `"equal"` across
#'   classes or estimated `"free"`ly in each.
#' @param residual Whether residual variances are `"occasion"`-specific or
#'   `"constant"` across occasions.
#' @param residual_equal Logical. Hold the residual variances equal across
#'   classes.
#' @param growth_covariates Optional numeric matrix of case-level covariates,
#'   one row per case and one column per covariate, regressing the growth
#'   factors on external variables. Must be row-aligned with the outcome matrix
#'   the emission is fitted to, and complete.
#' @param growth_covariates_equal Logical. Hold those regressions equal across
#'   classes.
#' @param ... Additional arguments, ignored.
#'
#' @return A list object of class \code{c("structured_normal", "emission")}.
#' @export
structured_normal_model <- function(n_components, design,
                                    random_effects = "intercept_slope",
                                    psi = "equal",
                                    residual = "occasion",
                                    residual_equal = TRUE,
                                    growth_covariates = NULL,
                                    growth_covariates_equal = TRUE, ...) {
  p <- ncol(design)
  r_cols <- switch(
    random_effects,
    none            = integer(0),
    intercept       = 1L,
    intercept_slope = seq_len(min(2L, p)),
    all             = seq_len(p),
    stop("Unknown `random_effects` specification: ", random_effects,
         call. = FALSE))

  if (!is.null(growth_covariates)) {
    growth_covariates <- as.matrix(growth_covariates)
    if (!is.numeric(growth_covariates))
      stop("`growth_covariates` must be numeric.", call. = FALSE)
    if (anyNA(growth_covariates))
      stop("`growth_covariates` must be complete; see fit_gmm() for why.",
           call. = FALSE)
    if (ncol(growth_covariates) == 0L) growth_covariates <- NULL
  }

  state <- list(
    n_components           = n_components,
    design                 = design,
    r_cols                 = r_cols,
    random_effects         = random_effects,
    psi_equal              = identical(psi, "equal"),
    theta_shared_occasions = identical(residual, "constant"),
    theta_equal            = isTRUE(residual_equal),

    # The covariate matrix lives in the emission state rather than being passed
    # to each method, because the E-step contract hands an emission nothing but
    # the outcome matrix. It is therefore the caller's job to keep it aligned
    # with the rows of X; fit_gmm() does that, including across the empty-row
    # deletion fit_mixture_internal() performs.
    xmat                   = growth_covariates,
    gamma_equal            = isTRUE(growth_covariates_equal),
    parameters             = list(),

    # No `em_tol`: refine_lbfgs() does not cover this emission, so
    # fit_single_init() applies the tighter unpolished stopping rule for us.
    #
    # EM with random effects converges slowly, and measurably so: on a
    # benchmark dataset (n = 500, four occasions) the median start needs ~660
    # iterations to reach the shared abs = 1e-4 rule, and with Psi and Theta
    # free across classes most starts were still climbing at the package
    # default ceiling of 1000. Reported estimates from a run that stopped
    # mid-climb are not the maximum-likelihood ones, so the ceiling is raised
    # here rather than left to the user to discover from a warning.
    em_max_iter            = 5000L,

    # Raising the ceiling alone would make a 30-restart fit cost 150,000
    # iterations, nearly all of them spent climbing hills that lose. A two-stage
    # search costs a twentieth of that and converges only the starts that
    # matter; see fit_em().
    #
    # It is enabled only where it was measured to be free. With Psi and Theta
    # held equal across classes, a 20-iteration first stage ranks the start that
    # wins *first* out of 30 on the ex8.1 data, so the staged search returns the
    # same optimum as running everything to convergence, six times faster. With
    # them free across classes it ranks that start *last* — and lengthening the
    # first stage to 50, 100 or 200 iterations makes it 27th, then 30th, then
    # 30th. The ranking is anti-correlated there rather than merely noisy,
    # because the winning basin is the slow one creeping along the ridge towards
    # the boundary where a variance would go negative, and slow early progress
    # is exactly what a short first stage penalises. No first-stage length fixes
    # that, so the class-free model runs every start to convergence instead:
    # -3221.290 against the -3223.379 a staged search settles for.
    em_stage1              = if (identical(psi, "equal") &&
                                 isTRUE(residual_equal))
      list(iter = 20L, frac = 0.2, min_keep = 4L) else NULL
  )
  class(state) <- c("structured_normal", "emission")
  state
}

# ------------------------------------------------------------------------------
# Internals
# ------------------------------------------------------------------------------

# Cases grouped by which occasions they observe.
#
# Every quantity this emission needs — the marginal covariance to invert, the
# conditional mean and variance of the random effects, the design used in the
# growth-mean update — depends on the case only through its set of observed
# occasions. Grouping first means one Cholesky per (class, pattern) instead of
# one per (class, case), which is what makes FIML here no more expensive than
# complete data when attrition is monotone (T + 1 patterns, not n).
.sn_patterns <- function(X) {
  obs <- !is.na(X)
  key <- do.call(paste0, as.data.frame(obs + 0L))
  ukey <- unique(key)
  lapply(ukey, function(k) {
    rows <- which(key == k)
    list(rows = rows, cols = which(obs[rows[1L], ]))
  })
}

# Sigma_k = Lambda_r Psi_k Lambda_r' + Theta_k, the marginal covariance of a
# class's repeated measures.
.sn_sigma <- function(ms, k) {
  Tn <- nrow(ms$design)
  S  <- diag(ms$parameters$theta[k, ], nrow = Tn)
  if (length(ms$r_cols)) {
    Lr <- ms$design[, ms$r_cols, drop = FALSE]
    S  <- S + Lr %*% ms$parameters$psi[[k]] %*% t(Lr)
  }
  (S + t(S)) / 2
}

# Cholesky factor, with a ridge as a last resort.
#
# Theta is floored strictly above zero and Psi is a sum of outer products, so
# Sigma is positive definite by construction and this normally returns on the
# first try. The ridge exists so that a numerically borderline class degrades
# into a slightly over-dispersed one rather than aborting the whole fit.
.sn_chol <- function(S) {
  ch <- tryCatch(chol(S), error = function(e) NULL)
  if (!is.null(ch)) return(ch)
  d <- mean(diag(S))
  for (eps in c(1e-8, 1e-6, 1e-4, 1e-2)) {
    ch <- tryCatch(chol(S + diag(eps * d, nrow(S))), error = function(e) NULL)
    if (!is.null(ch)) return(ch)
  }
  NULL
}

# Where each class's p x (1 + m) coefficient block sits in the joint system,
# as a p x (1 + m) matrix of global indices in vec order.
#
# The intercepts alpha_1..alpha_K come first, one p-block per class, and the
# shared Gamma follows as a single p*m block that every class points into. That
# pointing is the constraint: a class contributes to the same equations the
# others do, which is what makes the solution the constrained MLE rather than an
# average of unconstrained ones.
.sn_coef_index <- function(K, p, m) {
  lapply(seq_len(K), function(k) {
    I <- matrix(0L, p, m + 1L)
    I[, 1L] <- (k - 1L) * p + seq_len(p)
    I[, -1L] <- K * p + seq_len(p * m)
    I
  })
}

# Smallest residual variance allowed. Same role as the 1e-6 floor the diagonal
# Gaussian emission applies: a class that fitted its members exactly would
# otherwise claim infinite density and swallow the sample.
.sn_theta_floor <- 1e-6

# Symmetrise Psi and clip its eigenvalues at zero.
#
# The M-step update is a weighted sum of outer products plus a posterior
# variance, so it is positive semi-definite in exact arithmetic; this removes
# the asymmetry and the small negative eigenvalues that accumulate in floating
# point. A tiny positive floor rather than zero, because a growth factor whose
# variance reaches exactly zero can never leave — it contributes nothing to the
# conditional mean of the random effects, so no later iteration can revive it.
.sn_psi_pd <- function(P) {
  if (!length(P)) return(P)
  P <- (P + t(P)) / 2
  ev <- eigen(P, symmetric = TRUE)
  if (min(ev$values) < 1e-8) {
    ev$values <- pmax(ev$values, 1e-8)
    P <- ev$vectors %*% diag(ev$values, nrow(P)) %*% t(ev$vectors)
    P <- (P + t(P)) / 2
  }
  P
}

# Number of growth-factor covariates (zero for the unconditional model).
.sn_m <- function(ms) if (is.null(ms$xmat)) 0L else ncol(ms$xmat)

# Covariates with the leading column of ones: the design of the growth-factor
# regression, n x (1 + m). A one-column matrix of ones when there are none, so
# that alpha is simply the intercept of a regression on nothing.
.sn_z <- function(ms, n) if (is.null(ms$xmat)) matrix(1, n, 1L) else
  cbind(1, ms$xmat)

# A class's mean trajectory. Without covariates this is the T-vector Lambda
# alpha_k and every case shares it; with them it is one row per case,
# Lambda (alpha_k + Gamma_k x_i). The two shapes are kept distinct rather than
# always returning the n x T matrix because the unconditional model is the
# common one and would otherwise pay n times the memory and arithmetic for a
# quantity that does not vary.
.sn_mu <- function(ms, k) {
  a <- ms$parameters$alpha[k, ]
  if (is.null(ms$xmat)) return(as.vector(ms$design %*% a))
  B <- cbind(a, ms$parameters$gamma[[k]])            # p x (1 + m)
  tcrossprod(.sn_z(ms, nrow(ms$xmat)) %*% t(B), ms$design)
}

# Residuals of the rows and columns a missingness pattern selects, against
# whichever of the two mean shapes the model has.
.sn_resid <- function(X, mu, rows, cols) {
  Y <- X[rows, cols, drop = FALSE]
  if (is.matrix(mu)) Y - mu[rows, cols, drop = FALSE] else
    sweep(Y, 2, mu[cols], "-")
}

# Fitted trajectory on the outcome scale: K x T, what the plot draws and print()
# reports. The random effects have mean zero within a class, so the class mean
# curve is Lambda alpha_k regardless of Psi.
#
# With covariates there is no single class curve -- every case has its own -- so
# the reported one is evaluated at the sample mean of the covariates. That is
# the adjusted trajectory: the curve of a case at the average of x, not the
# average curve of the class, which would also carry the class's own covariate
# composition and so confound the trajectory with who is in the class.
.sn_fitted <- function(ms) {
  A <- ms$parameters$alpha
  if (!is.null(ms$xmat)) {
    p    <- ncol(ms$design)
    xbar <- colMeans(ms$xmat)
    shift <- vapply(ms$parameters$gamma, function(G) as.vector(G %*% xbar),
                    numeric(p))
    A <- A + t(matrix(shift, nrow = p))
  }
  t(ms$design %*% t(A))
}

# ------------------------------------------------------------------------------
# Initialisation
# ------------------------------------------------------------------------------

# Per-case ordinary least squares growth factors, used to place the starting
# values on the data's own scale.
#
# Returns the factor estimates (NA for a case observing fewer occasions than
# there are coefficients) along with the pooled within-case residual variance,
# which is a much better starting Theta than any constant would be.
.sn_case_factors <- function(design, X) {
  n <- nrow(X); p <- ncol(design)
  fac <- matrix(NA_real_, n, p)
  ss <- 0; df <- 0

  for (pat in .sn_patterns(X)) {
    cols <- pat$cols; rows <- pat$rows
    if (length(cols) < p) next
    Lp <- design[cols, , drop = FALSE]
    Yp <- X[rows, cols, drop = FALSE]
    B  <- tryCatch(t(qr.solve(Lp, t(Yp))), error = function(e) NULL)
    if (is.null(B)) next
    fac[rows, ] <- B
    if (length(cols) > p) {
      resid <- Yp - B %*% t(Lp)
      ss <- ss + sum(resid^2)
      df <- df + length(rows) * (length(cols) - p)
    }
  }

  theta0 <- if (df > 0) ss / df else stats::var(as.vector(X), na.rm = TRUE) / 2
  list(fac = fac, theta0 = max(theta0, 1e-3))
}

# Starting values: a random case's own growth curve per class.
#
# The package seeds a Gaussian mixture by picking K rows of X as the starting
# means; the analogue here is to pick K cases and take their fitted growth
# factors, which puts every class on a trajectory the data actually contain and
# is scale-free in the way a fixed perturbation is not. Psi starts at the
# between-case covariance of those factor estimates with the sampling
# contribution removed, which is the classical method-of-moments estimate.
#' @exportS3Method
init_params.structured_normal <- function(model_state, X, resp,
                                          random_state = NULL, ...) {
  if (!is.null(random_state)) set.seed(random_state)

  K  <- model_state$n_components
  L  <- model_state$design
  Tn <- nrow(L); p <- ncol(L)
  q  <- length(model_state$r_cols)

  m      <- .sn_m(model_state)
  cf     <- .sn_case_factors(L, X)
  fac    <- cf$fac
  theta0 <- cf$theta0
  ok     <- which(stats::complete.cases(fac))

  # With covariates on the growth factors, start Gamma at the pooled regression
  # of the case-level OLS factors on x and work from then on with the residuals.
  # Two reasons, both about starting where the model will end up: alpha is now
  # an intercept rather than a mean, so seeding it from a raw case factor would
  # place it off by Gamma x_i; and the between-case spread that seeds Psi is the
  # spread of the *residualised* factors, since whatever x explains is no longer
  # random variation within a class.
  gamma0 <- matrix(0, p, m)
  if (m > 0L && length(ok) > m + 1L) {
    Zok <- cbind(1, model_state$xmat[ok, , drop = FALSE])
    B   <- tryCatch(qr.solve(Zok, fac[ok, , drop = FALSE]),
                    error = function(e) NULL)          # (1 + m) x p
    if (!is.null(B) && all(is.finite(B))) {
      gamma0 <- t(B[-1L, , drop = FALSE])              # p x m
      fac    <- fac - model_state$xmat %*% t(gamma0)
    }
  }

  pooled <- if (length(ok)) colMeans(fac[ok, , drop = FALSE]) else
    qr.solve(L, colMeans(X, na.rm = TRUE))
  spread <- if (length(ok) > 1L)
    apply(fac[ok, , drop = FALSE], 2, stats::sd) else rep(1, p)
  spread[!is.finite(spread) | spread <= 0] <- 1

  alpha <- matrix(0, K, p)
  seeds <- if (length(ok) >= K) sample(ok, K) else rep(NA_integer_, K)
  for (k in seq_len(K)) {
    alpha[k, ] <- if (!is.na(seeds[k])) fac[seeds[k], ] else
      pooled + stats::rnorm(p, 0, spread)
  }

  # Between-case covariance of the OLS factors overstates Psi by the sampling
  # variance of the estimates themselves, theta * (L'L)^{-1}; subtracting it is
  # the standard moment correction. Halved as a starting value because starting
  # too small converges up more gracefully than starting too large converges
  # down (an over-large Psi flattens the responsibilities in the first E-step).
  psi0 <- diag(0, q)
  if (q > 0) {
    S <- if (length(ok) > p) stats::cov(fac[ok, , drop = FALSE]) else
      diag(spread^2, p)
    S <- S - theta0 * tryCatch(solve(crossprod(L)),
                               error = function(e) diag(0, p))
    psi0 <- .sn_psi_pd(S[model_state$r_cols, model_state$r_cols,
                         drop = FALSE] / 2)
  }

  model_state$parameters$alpha <- alpha
  model_state$parameters$gamma <- rep(list(gamma0), K)
  model_state$parameters$psi   <- rep(list(psi0), K)
  model_state$parameters$theta <- matrix(theta0, K, Tn)
  model_state
}

# ------------------------------------------------------------------------------
# Likelihood
# ------------------------------------------------------------------------------

#' @exportS3Method
log_likelihood.structured_normal <- function(model_state, X, ...) {
  n <- nrow(X)
  K <- model_state$n_components
  L <- model_state$design

  pats    <- .sn_patterns(X)
  log_eps <- matrix(0, n, K)

  for (k in seq_len(K)) {
    Sig <- .sn_sigma(model_state, k)
    mu  <- .sn_mu(model_state, k)

    for (pat in pats) {
      cols <- pat$cols; rows <- pat$rows
      # A case observing nothing contributes a flat likelihood, which is what
      # keeps it in the sample with its class probabilities equal to the
      # marginal ones rather than dropping it.
      if (!length(cols)) next

      ch <- .sn_chol(Sig[cols, cols, drop = FALSE])
      if (is.null(ch)) { log_eps[rows, k] <- -1e10; next }

      R <- .sn_resid(X, mu, rows, cols)
      z <- backsolve(ch, t(R), transpose = TRUE)          # |cols| x n_rows
      log_eps[rows, k] <- -0.5 * (length(cols) * log(2 * pi) +
                                    2 * sum(log(diag(ch))) + colSums(z^2))
    }
  }
  log_eps
}

# ------------------------------------------------------------------------------
# M-step
# ------------------------------------------------------------------------------
#
# Two passes over the (class, pattern) grid.
#
# Pass 1 computes, for every case, the conditional mean b_ik and variance V of
# its random effects given its own data and a hypothetical class membership,
#
#     b_ik = Psi_k Lambda_r' Sigma_k^{-1} (y_i - Lambda alpha_k)
#     V_kP = Psi_k - Psi_k Lambda_r' Sigma_k^{-1} Lambda_r Psi_k
#
# and accumulates from them the growth-mean normal equations and the Psi update.
# Note that V depends on the missingness pattern but not on the case, which is
# what makes the pattern grouping exact rather than an approximation.
#
# Pass 2 recomputes the residuals at the *new* growth means and updates Theta.
# Splitting the two is a conditional maximisation step (ECM): alpha is maximised
# given the old Theta and Theta given the new alpha, which is monotone in the
# Q-function just as a joint maximisation would be, and avoids the fixed-point
# iteration that solving them simultaneously would require.
#
# Covariates on the growth factors change only the first of the two normal
# equations, and only in its bookkeeping. Writing B_k = [alpha_k, Gamma_k] and
# z_i = (1, x_i), the mean of case i in class k is Lambda B_k z_i =
# (z_i' kron Lambda) vec(B_k), so the weighted least-squares system in vec(B_k)
# accumulates
#
#     A   += (sum_i w_ik z_i z_i')  kron  (Lambda_P' W_P Lambda_P)
#     rhs += vec( Lambda_P' W_P (sum_i w_ik r_i z_i') )
#
# which is the unconditional system exactly when z_i = 1. The one thing that is
# genuinely new is the *constraint*: with Gamma held equal across classes but
# alpha free -- the conventional default -- the K systems no
# longer separate, so they are assembled into one system of K p + p m equations
# by scattering each class's block into a shared index map. With Gamma free
# across classes, or absent, the system is block diagonal again and each class
# is solved on its own exactly as it was before covariates existed.
#' @exportS3Method
m_step.structured_normal <- function(model_state, X, resp, weights = NULL, ...) {
  if (!is.null(weights)) resp <- sweep(resp, 1, weights, "*")

  K  <- model_state$n_components
  L  <- model_state$design
  Tn <- nrow(L); p <- ncol(L)
  q  <- length(model_state$r_cols)
  Lr <- if (q) L[, model_state$r_cols, drop = FALSE] else NULL

  m  <- .sn_m(model_state)
  d  <- m + 1L
  Z  <- .sn_z(model_state, nrow(X))

  pats  <- .sn_patterns(X)
  alpha <- model_state$parameters$alpha
  gamma <- model_state$parameters$gamma
  psi   <- model_state$parameters$psi
  theta <- model_state$parameters$theta

  # Posterior means of the random effects, kept for the second pass: one
  # n x q matrix per class.
  bhat  <- lapply(seq_len(K), function(k) matrix(0, nrow(X), q))
  # Posterior variances, one q x q matrix per (class, pattern).
  vhat  <- vector("list", K)

  psi_num <- lapply(seq_len(K), function(k) matrix(0, q, q))
  psi_den <- numeric(K)
  alpha_new <- alpha
  gamma_new <- gamma

  # Per-class blocks of the growth-coefficient system, in vec(B_k) order.
  A_k   <- lapply(seq_len(K), function(k) matrix(0, p * d, p * d))
  rhs_k <- lapply(seq_len(K), function(k) numeric(p * d))

  # ---- Pass 1: random-effect posteriors, growth coefficients, Psi ------------
  for (k in seq_len(K)) {
    vhat[[k]] <- vector("list", length(pats))
    mu  <- .sn_mu(model_state, k)
    Sig <- .sn_sigma(model_state, k)
    th  <- theta[k, ]

    for (ip in seq_along(pats)) {
      cols <- pats[[ip]]$cols; rows <- pats[[ip]]$rows
      if (!length(cols)) next

      wk <- resp[rows, k]
      sw <- sum(wk)

      ch <- .sn_chol(Sig[cols, cols, drop = FALSE])
      if (is.null(ch)) next
      Sinv <- chol2inv(ch)

      R  <- .sn_resid(X, mu, rows, cols)

      if (q) {
        LrP <- Lr[cols, , drop = FALSE]
        G   <- psi[[k]] %*% t(LrP) %*% Sinv            # q x |cols|
        B   <- R %*% t(G)                              # n_rows x q
        V   <- psi[[k]] - G %*% LrP %*% psi[[k]]
        vhat[[k]][[ip]] <- V
        bhat[[k]][rows, ] <- B
        psi_num[[k]] <- psi_num[[k]] + crossprod(B, wk * B) + sw * V
        adj <- B %*% t(LrP)                            # random part, removed
      } else {
        adj <- 0
      }
      psi_den[k] <- psi_den[k] + sw

      # Growth-coefficient normal equations, weighted by the inverse residual
      # variance of the occasions this pattern observes.
      LP   <- L[cols, , drop = FALSE]
      Winv <- 1 / th[cols]
      LWL  <- crossprod(LP * Winv, LP)                 # p x p
      Zp   <- Z[rows, , drop = FALSE]                  # n_rows x d
      Res  <- X[rows, cols, drop = FALSE] - adj

      A_k[[k]]   <- A_k[[k]] + kronecker(crossprod(Zp, wk * Zp), LWL)
      rhs_k[[k]] <- rhs_k[[k]] +
        as.vector(crossprod(LP * Winv, crossprod(Res, wk * Zp)))
    }
  }

  # ---- Solve for alpha (and Gamma), with the class-equality constraint -------
  #
  # A class that has emptied out, or whose remaining cases cannot identify the
  # trajectory, keeps the coefficients it came in with rather than producing
  # NaNs that would poison the E-step for every other class. Without covariates
  # that is a per-class decision because the systems are separate; with a shared
  # Gamma they are not, so the empty class's block is dropped from the joint
  # system instead and the rest is solved without it.
  if (m == 0L || !model_state$gamma_equal) {
    # Block diagonal: solve each class on its own, as before.
    for (k in seq_len(K)) {
      sol <- tryCatch(solve(A_k[[k]], rhs_k[[k]]), error = function(e) NULL)
      if (!is.null(sol) && all(is.finite(sol))) {
        Bk <- matrix(sol, p, d)
        alpha_new[k, ]  <- Bk[, 1L]
        if (m > 0L) gamma_new[[k]] <- Bk[, -1L, drop = FALSE]
      }
    }
  } else {
    n_free <- K * p + p * m
    idx    <- .sn_coef_index(K, p, m)
    A   <- matrix(0, n_free, n_free)
    rhs <- numeric(n_free)
    for (k in seq_len(K)) {
      ii <- idx[[k]]
      A[ii, ii] <- A[ii, ii] + A_k[[k]]
      rhs[ii]   <- rhs[ii] + rhs_k[[k]]
    }
    active <- setdiff(seq_len(n_free),
                      unlist(lapply(which(psi_den <= 1e-8),
                                    function(k) (k - 1L) * p + seq_len(p))))
    sub <- tryCatch(solve(A[active, active, drop = FALSE], rhs[active]),
                    error = function(e) NULL)
    if (!is.null(sub) && all(is.finite(sub))) {
      sol <- rep(NA_real_, n_free)
      sol[active] <- sub
      for (k in seq_len(K)) {
        Bk <- matrix(sol[idx[[k]]], p, d)
        if (!all(is.finite(Bk))) next
        alpha_new[k, ]  <- Bk[, 1L]
        gamma_new[[k]] <- Bk[, -1L, drop = FALSE]
      }
    }
  }

  # ---- Psi, with the class-equality constraint imposed by pooling ------------
  if (q) {
    if (model_state$psi_equal) {
      num <- Reduce(`+`, psi_num)
      den <- sum(psi_den)
      shared <- .sn_psi_pd(if (den > 1e-8) num / den else psi[[1L]])
      psi <- rep(list(shared), K)
    } else {
      for (k in seq_len(K))
        if (psi_den[k] > 1e-8)
          psi[[k]] <- .sn_psi_pd(psi_num[[k]] / psi_den[k])
    }
  }

  # ---- Pass 2: residual variances at the new growth means --------------------
  th_num <- matrix(0, K, Tn)
  th_den <- matrix(0, K, Tn)

  updated <- model_state
  updated$parameters$alpha <- alpha_new
  updated$parameters$gamma <- gamma_new

  for (k in seq_len(K)) {
    mu <- .sn_mu(updated, k)
    for (ip in seq_along(pats)) {
      cols <- pats[[ip]]$cols; rows <- pats[[ip]]$rows
      if (!length(cols)) next

      wk <- resp[rows, k]
      sw <- sum(wk)
      E  <- .sn_resid(X, mu, rows, cols)

      corr <- 0
      if (q) {
        LrP <- Lr[cols, , drop = FALSE]
        E   <- E - bhat[[k]][rows, , drop = FALSE] %*% t(LrP)
        V   <- vhat[[k]][[ip]]
        # The residual carries the uncertainty in the random effect as well as
        # the squared point residual: E[(y - Lambda_r u)^2] = e^2 + Lambda_r V
        # Lambda_r' on the diagonal. Omitting it is the classic mistake that
        # makes Theta collapse towards zero.
        if (!is.null(V)) corr <- rowSums((LrP %*% V) * LrP)
      }

      th_num[k, cols] <- th_num[k, cols] + colSums(wk * E^2) + sw * corr
      th_den[k, cols] <- th_den[k, cols] + sw
    }
  }

  if (model_state$theta_shared_occasions) {
    th_num[] <- rowSums(th_num)
    th_den[] <- rowSums(th_den)
  }
  if (model_state$theta_equal) {
    th_num[] <- rep(colSums(th_num), each = K)
    th_den[] <- rep(colSums(th_den), each = K)
  }

  keep  <- th_den > 1e-8
  theta[keep] <- pmax(th_num[keep] / th_den[keep], .sn_theta_floor)

  model_state$parameters$alpha <- alpha_new
  model_state$parameters$gamma <- gamma_new
  model_state$parameters$psi   <- psi
  model_state$parameters$theta <- theta
  model_state
}

# ------------------------------------------------------------------------------
# Free parameters
# ------------------------------------------------------------------------------

#' @exportS3Method
n_parameters.structured_normal <- function(model_state, ...) {
  K  <- model_state$n_components
  p  <- ncol(model_state$design)
  Tn <- nrow(model_state$design)
  q  <- length(model_state$r_cols)
  m  <- .sn_m(model_state)

  n_psi   <- q * (q + 1L) / 2L * (if (model_state$psi_equal) 1L else K)
  n_theta <- (if (model_state$theta_shared_occasions) 1L else Tn) *
    (if (model_state$theta_equal) 1L else K)
  n_gamma <- p * m * (if (model_state$gamma_equal) 1L else K)

  K * p + n_gamma + n_psi + n_theta
}
