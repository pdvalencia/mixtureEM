# ------------------------------------------------------------------------------
# The Vuong-Lo-Mendell-Rubin test of K against K+1 classes.
#
# Two implementations are computed. Both take the eigenvalues of a matrix built
# from the score cross-products and the information of the two models; they
# differ only in which covariance of the parameters enters it — the ordinary one
# or the sandwich. The statistic 2(ll1 - ll0) is then referred to the weighted
# sum of one-degree-of-freedom chi-squares those eigenvalues define, whose tail
# is evaluated by Imhof's (1961) inversion integral.
#
# Vuong, Q. H. (1989). Likelihood ratio tests for model selection and
#   non-nested hypotheses. Econometrica, 57(2), 307-333. Theorem 3.3, eq. (3.6).
# Lo, Y., Mendell, N. R., & Rubin, D. B. (2001). Testing the number of
#   components in a normal mixture. Biometrika, 88(3), 767-778.
# Imhof, J. P. (1961). Computing the distribution of quadratic forms in normal
#   variables. Biometrika, 48(3/4), 419-426.
# Jeffries, N. O. (2003). A note on "Testing the number of components in a
#   normal mixture". Biometrika, 90(4), 991-994.
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# 1. Imhof's inversion integral
# ------------------------------------------------------------------------------

# Upper-tail probability P(Q > q) for Q = sum_j lam_j * chisq_1, the weights
# `lam` real and possibly negative.
#
# This is Imhof (1961) eq. (3.2), p. 422, specialised to one degree of freedom
# and no non-centrality in every term — the same expression Lo et al. (2001)
# write as their eq. (13), p. 773, oriented to the upper tail:
#
#   P(Q > q) = 1/2 + (1/pi) * integral_0^Inf sin(theta(u)) / (u * rho(u)) du
#
# The integrand is 0/0 at u = 0; Imhof (3.3) supplies the limit, and it is
# written in below rather than left to a quadrature rule that happens not to
# evaluate there.
.imhof <- function(q, lam, eps = 1e-10) {
  lam <- lam[lam != 0]
  if (!length(lam)) return(if (q < 0) 1 else 0)

  integrand <- function(u) {
    out <- numeric(length(u))
    at0 <- u == 0
    if (any(at0)) out[at0] <- 0.5 * sum(lam) - 0.5 * q
    if (any(!at0)) {
      uu    <- u[!at0]
      theta <- 0.5 * colSums(atan(outer(lam, uu))) - 0.5 * q * uu
      rho   <- exp(0.25 * colSums(log1p(outer(lam^2, uu^2))))
      out[!at0] <- sin(theta) / (uu * rho)
    }
    out
  }

  min(1, max(0, 0.5 + .imhof_integral(integrand, q, eps) / pi))
}

# The integral itself, accumulated over a finite window at a time rather than
# handed to `integrate(upper = Inf)`.
#
# Imhof (§3) truncates deliberately, and the reason shows up immediately here:
# the integrand decays only as u^(-(1 + m/2)) while oscillating, so with `m`
# small it is still oscillating perceptibly where a general-purpose rule has
# decided it has converged. Passing `upper = Inf` returns with message "OK" and
# an answer wrong in the third decimal for a single unit weight — the exact
# case the unit tests check against a chi-square on one degree of freedom.
# Application-sized weight vectors are nowhere near this badly behaved.
#
# For large u the arctangents saturate and theta'(u) -> -q/2, so the period
# settles at 4*pi/|q|. Integrating a fixed number of periods at a time keeps
# every sub-integral something a quadrature rule can see the shape of, and the
# window's contribution decays monotonically once the cancellation between
# half-periods takes over. Three consecutive negligible windows end it, rather
# than one, so that a window straddling a near-zero crossing cannot stop it
# early.
.imhof_integral <- function(integrand, q, eps) {
  step  <- if (abs(q) > 1e-8) 8 * 4 * pi / abs(q) else 1e3
  step  <- min(max(step, 10), 1e5)
  total <- 0
  a     <- 0
  quiet <- 0L
  for (i in seq_len(20000L)) {
    piece <- stats::integrate(integrand, a, a + step, subdivisions = 500L,
                              rel.tol = 1e-12, stop.on.error = FALSE)$value
    total <- total + piece
    if (abs(piece) < eps) {
      quiet <- quiet + 1L
      if (quiet >= 3L) break
    } else {
      quiet <- 0L
    }
    a <- a + step
  }
  total
}

# ------------------------------------------------------------------------------
# 2. The two W matrices
# ------------------------------------------------------------------------------

# Per-case scores and the Hessian of the summed log-likelihood, both by central
# differences on the step-one parameter vector, which is the unpenalised
# likelihood Vuong's theory is stated for.
#
# Returns NULL when the measurement model has no packable parameter vector.
.vlmr_pieces <- function(fit) {
  par <- .step1_pack(fit)
  if (is.null(par) || !length(par)) return(NULL)
  X <- fit$data
  if (is.null(X)) return(NULL)
  n <- nrow(X)
  w <- fit$sample_weights %||% rep(1, n)

  p <- length(par)
  h <- .step1_fd_step * pmax(1, abs(par))
  s <- vapply(seq_len(p), function(m) {
    a <- par; a[m] <- a[m] + h[m]
    b <- par; b[m] <- b[m] - h[m]
    w * (.step1_ll_case(fit, X, a) - .step1_ll_case(fit, X, b)) / (2 * h[m])
  }, numeric(n))
  dim(s) <- c(n, p)

  ll <- function(v) sum(w * .step1_ll_case(fit, X, v))
  list(s = s, A = .step1_fd_hessian(ll, par), ll = ll(par), p = p, n = n)
}

# The eigenvalues of W for one of the two implementations.
#
# Vuong (1989, Theorem 3.3, eq. 3.6, p. 313) writes W with the information
# matrices post-multiplied. What is formed here is the pre-multiplied
# arrangement D^-1 M where his is M D^-1, with D block-diagonal; the two are
# similar, D^-1 M = D^-1 (M D^-1) D, so they carry identical eigenvalues, which
# is all the test uses. Vermunt (2024) prints the post-multiplied form with the
# two off-diagonal inverses transposed relative to Vuong; that reading is not
# conformable, and Lo et al. (2001, Theorem 1, p. 772) restate Vuong verbatim,
# so Vuong's is followed here.
#
# `B` is a sum of cross-products across cases, not a mean.
#
# The pseudo-inverse is not defensive coding. The larger model reproduces the
# smaller only by emptying or duplicating a class, and at every such point some
# parameters drop out of the likelihood, so A_H1 is singular by construction
# (Jeffries, 2003). That singularity is also why the eigenvalues are not quite
# the ones the theorem describes.
.vlmr_eigen <- function(a0, a1, version) {
  # "standard" is Vuong's own form, on the ordinary covariance; "robust"
  # substitutes the sandwich, a modification another program makes and for which
  # Vermunt (2024) reports having found no theoretical justification in the
  # literature. It is a documented difference between implementations, not a bug
  # in either.
  B0  <- crossprod(a0$s)
  B1  <- crossprod(a1$s)
  B10 <- crossprod(a1$s, a0$s)
  Ai0 <- pinv(a0$A)
  Ai1 <- pinv(a1$A)

  W <- if (version == "standard") {
    rbind(cbind(-Ai1 %*% B1,      -Ai1 %*% B10),
          cbind( Ai0 %*% t(B10),   Ai0 %*% B0))
  } else {
    # "V^-1" is Vermunt's name for the sandwich itself, not an instruction to
    # invert anything further, and the minus signs of the Vuong form are
    # swallowed by the substitution of the robust variance for -A^-1.
    V1 <- Ai1 %*% B1 %*% Ai1
    V0 <- Ai0 %*% B0 %*% Ai0
    rbind(cbind( V1 %*% B1,      V1 %*% B10),
          cbind(-V0 %*% t(B10), -V0 %*% B0))
  }

  lam <- Re(eigen(W, only.values = TRUE)$values)
  lam <- lam[abs(lam) > 1e-8 * max(abs(lam))]
  # An orientation check, not an accuracy check. At a regular maximum B = -A and
  # sum(lam) would be p1 - p0; a mixture is not regular there, so only the sign
  # is trusted. Vermunt's own examples miss the count by several units.
  if (sum(lam) < 0) lam <- -lam
  lam
}

# ------------------------------------------------------------------------------
# 2b. The Satorra-Bentler / Asparouhov scaling factor for a nested comparison
#     under sampling weights or a complex survey design.
#
# c_m = tr(A_m^-1 B_m) / p_m (Satorra & Bentler, 1988; Asparouhov, 2005,
# "Sampling Weights in Latent Variable Modeling", Structural Equation
# Modeling 12(3), 411-434), built from the same step-one score/Hessian pieces
# the VLMR test above uses. `A` is the information matrix (the negated
# Hessian of the packed step-one log-likelihood); `B` is the crossproduct
# meat under simple weighting, or the PSU-aggregated survey meat when a
# design is declared.
#
# Returns NULL -- the caller's cue to refuse rather than compute -- when the
# fit carries a structural model (`.step1_pack()` covers the measurement
# model only, not a covariate or group regression on class membership: using
# it anyway would silently drop the structural model from the packed
# log-likelihood) or when the measurement model has no unconstrained packing
# at all (`.step1_pack_mm()` returns NULL for growth and other structured
# families). Restricted to `mixture_model` fits for the same reason
# `.vlmr_pair()` is: an `lta_model`'s state does not share `.step1_pack()`'s
# field names (`n_components`, `weights`), so calling it in would go straight
# to malformed output at .step1_pack(), never a clean NULL.
.scaling_pieces <- function(info) {
  if (!inherits(info$fit, "mixture_model")) return(NULL)
  if (.supplies_class_probs(info$fit$sm)) return(NULL)
  pieces <- .vlmr_pieces(info$fit)
  if (is.null(pieces)) return(NULL)
  A <- -pieces$A
  B <- if (isTRUE(info$has_survey_design))
    compute_survey_B(pieces$s, info$strata, info$cluster)
  else
    crossprod(pieces$s)
  list(c = sum(diag(pinv(A) %*% B)) / pieces$p, p = pieces$p)
}

# Whether a fit's log-likelihood is a pseudo-likelihood that a plain
# likelihood-ratio-difference test is not asymptotically chi-square on: either
# a complex survey design was declared, or `weights` carries genuine sampling
# weights (not all 1, and not declared as frequencies -- a frequency of 10 is
# ten identical cases and its likelihood is a true likelihood).
.needs_scaling_correction <- function(info) {
  isTRUE(info$has_survey_design) ||
    (identical(info$weight_type, "sampling") && !all(info$weights == 1))
}

# ------------------------------------------------------------------------------
# 3. One K against K+1
# ------------------------------------------------------------------------------

# Returns a one-row list of LR, and per version the mean, sd and p-value of the
# reference distribution; or a `reason` when the test cannot be computed.
.vlmr_pair <- function(fit0, fit1, versions) {
  na <- function(reason) list(reason = reason)

  if (!inherits(fit0, "mixture_model") || !inherits(fit1, "mixture_model"))
    return(na("the test is available for mixture models only"))
  if (is.null(fit0$data) || is.null(fit1$data) ||
      !identical(dim(fit0$data), dim(fit1$data)) ||
      !isTRUE(all.equal(fit0$data, fit1$data, check.attributes = FALSE)))
    return(na("the two models were not fit on the same data"))

  a0 <- .vlmr_pieces(fit0)
  a1 <- .vlmr_pieces(fit1)
  if (is.null(a0) || is.null(a1))
    return(na("this measurement model has no packable parameter vector"))
  if (a0$p + a1$p > .vlmr_max_params)
    return(na(sprintf("the two models have %d parameters between them (limit %d)",
                      a0$p + a1$p, .vlmr_max_params)))

  LR  <- 2 * (a1$ll - a0$ll)
  out <- list(LR = LR)
  for (v in versions) {
    lam <- .vlmr_eigen(a0, a1, v)
    out[[v]] <- list(mean = sum(lam),
                     sd   = sqrt(2 * sum(lam^2)),
                     p    = .imhof(LR, lam))
  }
  out
}

# Above this many parameters between the two models the O(p^2) Hessians stop
# being affordable, and the eigen-decomposition of the (p0 + p1) square W is no
# longer the cheap part either.
.vlmr_max_params <- 200L

# ------------------------------------------------------------------------------
# 4. The columns on a comparison table
# ------------------------------------------------------------------------------

# Appends the VLMR columns to a fit table. Row i tests the model in row i
# against the one in row i + 1, so the last row is NA by construction. Returns
# the table with the columns added and the per-K detail attached as an
# attribute for the comparison object to carry.
.vlmr_augment <- function(tab, models, vlmr) {
  versions <- switch(vlmr,
                     standard   = "standard",
                     robust     = "robust",
                     both       = c("standard", "robust"))
  n <- nrow(tab)
  LR <- p_std <- p_rob <- rep(NA_real_, n)
  detail <- vector("list", n)

  for (i in seq_len(n - 1L)) {
    f0 <- models[[paste0("K", tab$Classes[i])]]
    f1 <- models[[paste0("K", tab$Classes[i + 1L])]]
    res <- .vlmr_pair(f0, f1, versions)
    detail[[i]] <- res
    if (!is.null(res$reason)) {
      message(sprintf("VLMR not computed for %d vs %d classes: %s.",
                      tab$Classes[i], tab$Classes[i + 1L], res$reason))
      next
    }
    LR[i] <- res$LR
    if (!is.null(res$standard)) p_std[i] <- res$standard$p
    if (!is.null(res$robust))   p_rob[i] <- res$robust$p
  }

  tab$VLMR_LR <- LR
  if ("standard" %in% versions) tab$VLMR_p <- p_std
  if ("robust" %in% versions)   tab$VLMR_p_robust <- p_rob
  attr(tab, "vlmr_detail") <- detail
  tab
}
