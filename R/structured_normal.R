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
#   psi_type                "equal", "variances" or "free": Psi held equal
#                           across classes (the usual default), class-specific
#                           variances with the covariances held equal, or free
#   psi_equal               psi_type == "equal", kept for the readers that only
#                           ask whether Psi is one matrix
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
#'   classes, has class-specific `"variances"` with the covariances held equal,
#'   or is estimated `"free"`ly in each.
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
  if (!psi %in% c("equal", "variances", "free"))
    stop("Unknown `psi` specification: ", psi, call. = FALSE)
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
    psi_type               = psi,
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
    # wins *first* out of 30 on a growth-mixture benchmark, so the staged
    # search returns the same optimum as running everything to convergence,
    # six times faster. With
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

# Whether a class's Sigma_k has gone singular: the spike at which the
# likelihood of a mixture of normals is unbounded (Kiefer and Wolfowitz 1956;
# Day 1969). On outcomes with many tied values -- a scale floor, say -- a class
# can sit on the tied cases with no spread left in some direction, and its
# density there grows without limit. Such a point is not a maximum of the
# likelihood in any useful sense: it describes the ties, not a subgroup, and it
# outscores every meaningful solution. The search therefore treats a start that
# ends there as one that failed to converge (fit_em()), which is the standard
# handling of degenerate solutions in normal mixtures (McLachlan and Peel 2000,
# ch. 3). The floor on Theta only keeps such a start numerically alive.
#
# The test is the smallest eigenvalue of any Sigma_k relative to the smallest
# observed variance of an occasion, so it reads the same in any units. Measured
# in simulated growth data with a scale floor, spike end points sit at 3.4e-6
# and below, admissible solutions -- with a growth-factor
# covariance outside the cone included -- at 8.5e-3 and above.
# The threshold is the log-midpoint of that gap. Without the observed variances
# the test falls back to an eigenvalue at the Theta floor.
.sn_spike_ratio <- 1.7e-4
.sn_degenerate <- function(ms, X = NULL, weights = NULL,
                           s2 = .gmm_marginal(ms, X, weights)) {
  ev <- vapply(seq_len(ms$n_components), function(k)
    min(eigen(.sn_sigma(ms, k), symmetric = TRUE, only.values = TRUE)$values),
    numeric(1))
  if (any(!is.finite(ev))) return(TRUE)
  if (is.null(s2)) return(min(ev) <= 10 * .sn_theta_floor)
  min(ev) / min(s2) < .sn_spike_ratio
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

# Which Psi form a state carries. A state built before `psi_type` existed knew
# only the two forms `psi_equal` distinguishes.
.sn_psi_type <- function(ms)
  ms$psi_type %||% (if (isTRUE(ms$psi_equal)) "equal" else "free")

# The Psi update under psi = "variances": class-specific variances, one set of
# covariances shared by every class.
#
# The Psi part of the expected complete-data log-likelihood is
#
#     Q(Psi) = -1/2 sum_k n_k [ log|Psi_k| + tr(Psi_k^{-1} S_k) ],
#
# with S_k = psi_num_k / n_k the class's expected cross-product of the random
# effects. Free across classes it is maximised by Psi_k = S_k and held equal by
# the pooled S; with only the covariances shared there is no closed form, so the
# maximum is found numerically. Its gradient is n_k/2 (Psi^-1 S Psi^-1 - Psi^-1),
# taken element by element on each class's diagonal and summed over classes for
# the shared off-diagonal. The variances are searched on the log scale, and a
# point where some Psi_k is not positive definite scores as infinitely bad, so
# the update stays inside the cone: log|Psi_k| is a barrier there. The result is
# taken only if Q rises from where the iteration came in, which is what keeps EM
# monotone (a generalised EM step).
.sn_psi_variances <- function(psi_num, psi_den, psi_old) {
  K <- length(psi_old); q <- nrow(psi_old[[1L]])
  live <- psi_den > 1e-8
  if (!any(live)) return(psi_old)
  S  <- lapply(seq_len(K), function(k)
    if (live[k]) psi_num[[k]] / psi_den[k] else psi_old[[k]])
  lt <- which(lower.tri(diag(q)), arr.ind = TRUE)
  nc <- nrow(lt)

  build <- function(v) {
    D <- matrix(exp(v[seq_len(K * q)]), K, q)
    C <- v[K * q + seq_len(nc)]
    lapply(seq_len(K), function(k) {
      P <- diag(D[k, ], q)
      P[lt] <- C; P[lt[, 2:1, drop = FALSE]] <- C
      P
    })
  }
  negQ <- function(v) {
    out <- 0
    P <- build(v)
    for (k in which(live)) {
      ch <- tryCatch(chol(P[[k]]), error = function(e) NULL)
      if (is.null(ch)) return(1e100)
      out <- out + 0.5 * psi_den[k] *
        (2 * sum(log(diag(ch))) + sum(diag(chol2inv(ch) %*% S[[k]])))
    }
    out
  }
  grad <- function(v) {
    P <- build(v)
    gd <- matrix(0, K, q); gc <- numeric(nc)
    for (k in which(live)) {
      Pi <- tryCatch(chol2inv(chol(P[[k]])), error = function(e) NULL)
      if (is.null(Pi)) return(rep(0, length(v)))
      G <- 0.5 * psi_den[k] * (Pi %*% S[[k]] %*% Pi - Pi)
      gd[k, ] <- diag(G) * diag(P[[k]])
      gc <- gc + 2 * G[lt]
    }
    -c(as.vector(gd), gc)
  }

  P0 <- lapply(psi_old, function(P) (P + t(P)) / 2)
  v0 <- c(as.vector(t(vapply(P0, function(P) log(pmax(diag(P), 1e-8)),
                             numeric(q)))),
          Reduce(`+`, lapply(seq_len(K), function(k)
            psi_den[k] * P0[[k]][lt])) / sum(psi_den))
  f0 <- negQ(v0)
  res <- tryCatch(stats::optim(v0, negQ, grad, method = "BFGS",
                               control = list(maxit = 100L, reltol = 1e-12)),
                  error = function(e) NULL)
  if (is.null(res) || !is.finite(res$value) || res$value >= f0 ||
      res$value >= 1e100) {
    # The entry point itself may not carry the structure yet (a start, or a
    # state the finish moved outside the cone); fall back to it projected.
    return(if (is.finite(f0) && f0 < 1e100) build(v0) else psi_old)
  }
  build(res$par)
}

# The edge of the cone, and the M-step that crosses it.
#
# The Psi updates above treat each case's growth factors as missing data, and a
# covariance estimated from imputed factors is a sum of outer products: it can
# approach the edge of the positive semi-definite cone but never cross it. The
# likelihood itself only needs every Sigma_k positive definite, and its maximum
# can lie past that edge -- a small class whose growth factors vary less than
# their residuals allow, estimated with a negative variance. EM started in such
# a basin crawls along the edge and stops there, well below the maximum.
#
# At the edge the M-step therefore treats only the class as missing. Given the
# posteriors, the expected complete-data log-likelihood in the mean and
# covariance parameters is
#
#     Q = -1/2 sum_k sum_P [ n_kP log|Sigma_kP| + tr(Sigma_kP^{-1} S_kP) ],
#
# with S_kP the class-weighted cross-product of the residuals of the occasions
# pattern P observes. The growth coefficients maximise it in closed form at the
# current Sigma (generalised least squares); Psi and Theta are then climbed
# together by a short BFGS run, Psi as an ordinary symmetric matrix and Theta on
# the log scale, with a point where some Sigma_k is not positive definite scored
# as infinitely bad. The run is kept only if Q rises. Both this step and the one
# above are generalised EM steps for the same likelihood, so EM stays monotone
# whichever runs; inside the cone the closed-form update is kept because it is
# cheaper and leaves every fit there exactly as it was.
#
# "At the edge" is read in the units the Newton finish packs Psi in
# (.qn_sn_psi_unit()), so that a growth factor's variance is measured against
# the residual variance it competes with rather than against its own size.
.sn_edge_ratio <- 1e-3
.sn_psi_at_edge <- function(ms) {
  if (!length(ms$r_cols)) return(FALSE)
  U <- .qn_sn_psi_unit(ms)
  any(vapply(ms$parameters$psi, function(P)
    min(eigen(((P + t(P)) / 2) / U, symmetric = TRUE, only.values = TRUE)$values),
    numeric(1)) < .sn_edge_ratio)
}

# Solve the growth-coefficient normal equations, with the class-equality
# constraint on Gamma. A class that has emptied out, or whose remaining cases
# cannot identify the trajectory, keeps the coefficients it came in with rather
# than producing NaNs that would poison the E-step for every other class.
# Without covariates that is a per-class decision because the systems are
# separate; with a shared Gamma they are not, so the empty class's block is
# dropped from the joint system instead and the rest is solved without it.
.sn_solve_coef <- function(model_state, A_k, rhs_k, den) {
  K <- model_state$n_components; p <- ncol(model_state$design)
  m <- .sn_m(model_state); d <- m + 1L
  alpha <- model_state$parameters$alpha
  gamma <- model_state$parameters$gamma
  if (m == 0L || !model_state$gamma_equal) {
    for (k in seq_len(K)) {
      sol <- tryCatch(solve(A_k[[k]], rhs_k[[k]]), error = function(e) NULL)
      if (!is.null(sol) && all(is.finite(sol))) {
        Bk <- matrix(sol, p, d)
        alpha[k, ] <- Bk[, 1L]
        if (m > 0L) gamma[[k]] <- Bk[, -1L, drop = FALSE]
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
                      unlist(lapply(which(den <= 1e-8),
                                    function(k) (k - 1L) * p + seq_len(p))))
    sub <- tryCatch(solve(A[active, active, drop = FALSE], rhs[active]),
                    error = function(e) NULL)
    if (!is.null(sub) && all(is.finite(sub))) {
      sol <- rep(NA_real_, n_free)
      sol[active] <- sub
      for (k in seq_len(K)) {
        Bk <- matrix(sol[idx[[k]]], p, d)
        if (!all(is.finite(Bk))) next
        alpha[k, ] <- Bk[, 1L]
        gamma[[k]] <- Bk[, -1L, drop = FALSE]
      }
    }
  }
  list(alpha = alpha, gamma = gamma)
}

# The M-step at the edge (see above). `resp` already carries any case weights.
.sn_cov_step <- function(model_state, X, resp, pats) {
  K  <- model_state$n_components
  L  <- model_state$design
  Tn <- nrow(L); p <- ncol(L)
  Lr <- L[, model_state$r_cols, drop = FALSE]; q <- ncol(Lr)
  d  <- .sn_m(model_state) + 1L
  Z  <- .sn_z(model_state, nrow(X))

  # ---- Growth coefficients by generalised least squares at the current Sigma
  A_k   <- lapply(seq_len(K), function(k) matrix(0, p * d, p * d))
  rhs_k <- lapply(seq_len(K), function(k) numeric(p * d))
  den   <- numeric(K)
  for (k in seq_len(K)) {
    Sig <- .sn_sigma(model_state, k)
    for (pat in pats) {
      cols <- pat$cols; rows <- pat$rows
      if (!length(cols)) next
      ch <- .sn_chol(Sig[cols, cols, drop = FALSE])
      if (is.null(ch)) return(NULL)
      wk  <- resp[rows, k]
      LSi <- crossprod(L[cols, , drop = FALSE], chol2inv(ch))   # p x |cols|
      Zp  <- Z[rows, , drop = FALSE]
      A_k[[k]]   <- A_k[[k]] + kronecker(crossprod(Zp, wk * Zp),
                                         LSi %*% L[cols, , drop = FALSE])
      rhs_k[[k]] <- rhs_k[[k]] +
        as.vector(LSi %*% crossprod(X[rows, cols, drop = FALSE], wk * Zp))
      den[k] <- den[k] + sum(wk)
    }
  }
  co <- .sn_solve_coef(model_state, A_k, rhs_k, den)
  model_state$parameters$alpha <- co$alpha
  model_state$parameters$gamma <- co$gamma

  # ---- Residual cross-products at the new coefficients
  live <- which(den > 1e-8)
  SS <- lapply(seq_len(K), function(k) {
    mu <- .sn_mu(model_state, k)
    lapply(pats, function(pat) {
      if (!length(pat$cols) || !(k %in% live)) return(NULL)
      wk <- resp[pat$rows, k]
      R  <- .sn_resid(X, mu, pat$rows, pat$cols)
      list(cols = pat$cols, n = sum(wk), S = crossprod(R, wk * R))
    })
  })

  # ---- Where each class's Psi and Theta read their values in the vector
  ptype <- .sn_psi_type(model_state)
  li  <- which(lower.tri(diag(q), diag = TRUE), arr.ind = TRUE)
  nl  <- nrow(li); dg <- li[, 1L] == li[, 2L]
  npsi <- switch(ptype, equal = nl, free = K * nl, variances = K * q + sum(!dg))
  pid <- function(k) switch(ptype,
    equal = seq_len(nl),
    free  = (k - 1L) * nl + seq_len(nl),
    variances = { v <- integer(nl)
                  v[dg]  <- (k - 1L) * q + li[dg, 1L]
                  v[!dg] <- K * q + seq_len(sum(!dg)); v })
  te <- isTRUE(model_state$theta_equal)
  so <- isTRUE(model_state$theta_shared_occasions)
  tid <- function(k)
    if (te && so) rep(1L, Tn) else if (te) seq_len(Tn) else
      if (so) rep(k, Tn) else (k - 1L) * Tn + seq_len(Tn)
  U <- .qn_sn_psi_unit(model_state)
  Ul <- U[li]

  build <- function(v) {
    psi <- lapply(seq_len(K), function(k) {
      P <- matrix(0, q, q)
      P[li] <- v[pid(k)] * Ul
      P + t(P) - diag(diag(P), q)
    })
    th <- t(vapply(seq_len(K), function(k) exp(v[npsi + tid(k)]), numeric(Tn)))
    list(psi = psi, theta = matrix(th, K, Tn))
  }
  # Value and the derivative of Q in each class's Sigma, in one pass.
  evalQ <- function(v, want_grad) {
    b <- build(v)
    if (any(b$theta < .sn_theta_floor)) return(NULL)
    out <- 0; H <- vector("list", K)
    for (k in live) {
      Sig <- Lr %*% b$psi[[k]] %*% t(Lr) + diag(b$theta[k, ], Tn)
      if (is.null(tryCatch(chol(Sig), error = function(e) NULL))) return(NULL)
      Hk <- matrix(0, Tn, Tn)
      for (s in SS[[k]]) {
        if (is.null(s)) next
        ch <- tryCatch(chol(Sig[s$cols, s$cols, drop = FALSE]),
                       error = function(e) NULL)
        if (is.null(ch)) return(NULL)
        Si  <- chol2inv(ch)
        out <- out + 0.5 * (s$n * 2 * sum(log(diag(ch))) + sum(Si * s$S))
        if (want_grad)
          Hk[s$cols, s$cols] <- Hk[s$cols, s$cols] +
            0.5 * (Si %*% s$S %*% Si - s$n * Si)
      }
      H[[k]] <- Hk
    }
    list(value = out, H = H, b = b)
  }
  negQ <- function(v) {
    e <- evalQ(v, FALSE)
    if (is.null(e)) 1e100 else e$value
  }
  grad <- function(v) {
    e <- evalQ(v, TRUE)
    g <- numeric(length(v))
    if (is.null(e)) return(g)
    for (k in live) {
      dP <- crossprod(Lr, e$H[[k]] %*% Lr)
      gp <- dP[li] * ifelse(dg, 1, 2) * Ul
      ip <- pid(k)
      for (j in seq_along(ip)) g[ip[j]] <- g[ip[j]] + gp[j]
      it <- npsi + tid(k)
      gt <- diag(e$H[[k]]) * e$b$theta[k, ]
      for (t in seq_len(Tn)) g[it[t]] <- g[it[t]] + gt[t]
    }
    -g
  }

  # Start from the current values, averaged where a constraint shares them.
  v0  <- numeric(npsi + max(vapply(seq_len(K), function(k) max(tid(k)),
                                   numeric(1))))
  cnt <- numeric(length(v0))
  th  <- model_state$parameters$theta
  for (k in seq_len(K)) {
    P  <- model_state$parameters$psi[[k]]
    ip <- pid(k); it <- npsi + tid(k)
    v0[ip] <- v0[ip] + ((P + t(P)) / 2)[li] / Ul; cnt[ip] <- cnt[ip] + 1
    for (t in seq_len(Tn)) {
      v0[it[t]] <- v0[it[t]] + log(max(th[k, t], .sn_theta_floor))
      cnt[it[t]] <- cnt[it[t]] + 1
    }
  }
  v0 <- v0 / pmax(cnt, 1)
  f0 <- negQ(v0)
  if (f0 >= 1e100) return(NULL)
  res <- tryCatch(stats::optim(v0, negQ, grad, method = "BFGS",
                               control = list(maxit = 100L, reltol = 1e-12)),
                  error = function(e) NULL)
  b <- build(if (!is.null(res) && is.finite(res$value) && res$value < f0)
    res$par else v0)
  model_state$parameters$psi   <- b$psi
  model_state$parameters$theta <- b$theta
  model_state
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

    # A Psi outside the positive semi-definite cone is a legitimate estimate so
    # long as Sigma_k stays positive definite; once it does not, there is no
    # likelihood at all. Only the Newton finish can put Psi there (EM's update
    # is a sum of outer products), and NaN is what makes the finish's line
    # search refuse the step, where .sn_chol()'s ridge would quietly score a
    # point that does not exist.
    if (length(model_state$r_cols) &&
        is.null(tryCatch(chol(Sig), error = function(e) NULL)) &&
        min(eigen(model_state$parameters$psi[[k]], symmetric = TRUE,
                  only.values = TRUE)$values) < 0) {
      log_eps[, k] <- NaN
      next
    }

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

  # At or past the edge of the cone the update below cannot move Psi outward.
  if (q && .sn_psi_at_edge(model_state)) {
    out <- .sn_cov_step(model_state, X, resp, pats)
    if (!is.null(out)) return(out)
  }

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
  co <- .sn_solve_coef(model_state, A_k, rhs_k, psi_den)
  alpha_new <- co$alpha
  gamma_new <- co$gamma

  # ---- Psi, with the class-equality constraint imposed by pooling ------------
  if (q) {
    if (model_state$psi_equal) {
      num <- Reduce(`+`, psi_num)
      den <- sum(psi_den)
      shared <- .sn_psi_pd(if (den > 1e-8) num / den else psi[[1L]])
      psi <- rep(list(shared), K)
    } else if (identical(.sn_psi_type(model_state), "variances") && q > 1L) {
      psi <- .sn_psi_variances(psi_num, psi_den, psi)
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

  n_psi   <- switch(.sn_psi_type(model_state),
    equal     = q * (q + 1L) / 2L,
    variances = K * q + q * (q - 1L) / 2L,
    q * (q + 1L) / 2L * K)
  n_theta <- (if (model_state$theta_shared_occasions) 1L else Tn) *
    (if (model_state$theta_equal) 1L else K)
  n_gamma <- p * m * (if (model_state$gamma_equal) 1L else K)

  K * p + n_gamma + n_psi + n_theta
}
