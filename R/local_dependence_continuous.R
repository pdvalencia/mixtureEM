# ==============================================================================
# Local dependence for continuous indicators
# ==============================================================================
#
# bivariate_residuals()'s categorical statistic is a Pearson chi-square over a
# two-way table, which a continuous item does not have. What it is checking
# for, though, is the same thing: whether the fitted classes reproduce the
# association between a pair of items, or leave some of it as a residual
# dependence the model cannot see.
#
# For a Gaussian measurement model that question has a direct answer: fix a
# residual covariance between items a and b, in one class only, at zero (as
# the fitted conditional-independence model does), and ask how much the
# observed-data log-likelihood wants to move it. That is a score test —
# equivalently, in Sorbom's (1989) terminology, a modification index — and
# unlike a bivariate residual referred to chi-square it keeps its nominal
# size and power (Oberski, van Kollenburg & Vermunt, 2013).
#
# Class-specific, because a residual dependence can run in opposite
# directions in different classes (see the worked example in the package's
# development notes): a class-invariant version is blind by construction to
# any pair whose association reverses sign across classes.
#
# References
#   Sorbom, D. (1989). Model modification. Psychometrika, 54, 371-384.
#   Oberski, D. L., van Kollenburg, G. H., & Vermunt, J. K. (2013). A Monte
#     Carlo evaluation of three methods to detect local dependence in binary
#     data latent class models. Advances in Data Analysis and Classification,
#     7(3), 267-279. \doi{10.1007/s11634-013-0146-2}

# ------------------------------------------------------------------------------
# The augmented case-level log-likelihood
# ------------------------------------------------------------------------------
#
# `par` is a step-one parameter vector on the .step1_pack() scale; `theta` is
# the residual covariance freed between items `a` and `b`, in class `k` only.
# Every other class, and every other item pair within class k, stays exactly
# as the fitted model has it. Returns one log-likelihood value per case.
.local_dep_ll_case <- function(model_state, X, par, theta, a, b, k) {
  ms <- .step1_unpack(model_state, par)
  ll <- log_likelihood(ms$mm, X)

  mu  <- ms$mm$parameters$means[k, ]
  s2  <- ms$mm$parameters$covariances[k, ]
  xa  <- X[, a] - mu[a]
  xb  <- X[, b] - mu[b]

  # What class k's density currently contributes from items a and b alone,
  # to be replaced by the joint bivariate-normal term below.
  da <- stats::dnorm(X[, a], mu[a], sqrt(s2[a]), log = TRUE)
  db <- stats::dnorm(X[, b], mu[b], sqrt(s2[b]), log = TRUE)

  det_s <- s2[a] * s2[b] - theta^2
  quad  <- (xa^2 * s2[b] + xb^2 * s2[a] - 2 * theta * xa * xb) / det_s
  joint <- -log(2 * pi) - 0.5 * log(det_s) - 0.5 * quad

  ll[, k] <- ll[, k] - da - db + joint
  logsumexp(sweep(ll, 2, log(pmax(ms$weights, 1e-300)), "+"), MARGIN = 1)
}

# ------------------------------------------------------------------------------
# The shared Hessian block, H_psi,psi
# ------------------------------------------------------------------------------
#
# At theta = 0 the augmented model is the fitted model for every (pair,
# class), so this block of the augmented Hessian does not depend on which
# statistic is being computed. It is evaluated once, by the same 4-point
# central-difference loop .step1_variance() uses (R/step3_variance.R), copied
# rather than reused because that function bundles in the outer-product
# fallback and a variance (not raw Hessian) return value.
.local_dep_psi_hessian <- function(model_state, X, par, w) {
  p  <- length(par)
  h  <- .step1_fd_step * pmax(1, abs(par))
  ll <- function(v) sum(w * .step1_ll_case(model_state, X, v))

  H <- matrix(0, p, p)
  for (i in seq_len(p)) {
    for (j in i:p) {
      ei <- numeric(p); ei[i] <- h[i]
      ej <- numeric(p); ej[j] <- h[j]
      H[i, j] <- H[j, i] <-
        (ll(par + ei + ej) - ll(par + ei - ej) -
           ll(par - ei + ej) + ll(par - ei - ej)) / (4 * h[i] * h[j])
    }
  }
  H
}

# ------------------------------------------------------------------------------
# One (pair, class) statistic
# ------------------------------------------------------------------------------
#
# H_psi,psi is shared across every (pair, class) and computed once by the
# caller, which is what keeps this affordable: only the theta row/column of
# the augmented Hessian is computed per statistic, roughly 4p likelihood
# evaluations rather than a fresh (p+1)^2 Hessian per pair per class. A
# cheaper one-sided scheme (reusing the stationarity of psi_hat under the
# *baseline* log-likelihood to skip a psi-perturbed baseline evaluation) was
# tried and measured unstable: the mixed partial it produces is swamped by
# d^2 ll / dpsi_m^2, which is not itself suppressed by h_th, and the result
# was wrong by orders of magnitude. The full central 4-point formula below,
# at the same step used everywhere else in the package, is accurate from
# 1e-5 to 1e-2 relative step on the reference data and is what is used.
.local_dep_statistic <- function(model_state, X, par, w, ll0, ll0_case,
                                 a, b, k, h_th, Hpsipsi_inv, S1) {
  llp_case <- .local_dep_ll_case(model_state, X, par,  h_th, a, b, k)
  llm_case <- .local_dep_ll_case(model_state, X, par, -h_th, a, b, k)
  llp <- sum(w * llp_case)
  llm <- sum(w * llm_case)

  s     <- (llp - llm) / (2 * h_th)
  Hthth <- (llp - 2 * ll0 + llm) / h_th^2

  # The mixed partial cannot reuse llp/llm the way Hthth does: a one-sided
  # scheme built from a shared baseline turns out to be swamped by
  # d^2ll/dpsi_m^2, which is not itself suppressed by h_th (measured: garbage
  # results, off by orders of magnitude and not even stable in sign). The
  # full central 4-point formula, at the same step used everywhere else in
  # the package, is what is actually stable here.
  p1 <- length(par)
  h  <- .step1_fd_step * pmax(1, abs(par))
  ll_at <- function(theta, par_v)
    sum(w * .local_dep_ll_case(model_state, X, par_v, theta, a, b, k))
  Hthpsi <- vapply(seq_len(p1), function(m) {
    parp <- par; parp[m] <- parp[m] + h[m]
    parm <- par; parm[m] <- parm[m] - h[m]
    (ll_at(h_th, parp) - ll_at(h_th, parm) -
       ll_at(-h_th, parp) + ll_at(-h_th, parm)) / (4 * h_th * h[m])
  }, numeric(1))

  # The roadmap's H blocks are minus the Hessian (observed information,
  # positive semi-definite); Hpsipsi_inv passed in is the *raw* Hessian's
  # inverse, negative definite, so both the sign of the diagonal term and the
  # sign inside the quadratic form flip relative to the raw Hessian pieces
  # computed above: D = -Hthth + Hthpsi' Hpsipsi_inv Hthpsi.
  D_hess <- -Hthth + as.numeric(Hthpsi %*% Hpsipsi_inv %*% Hthpsi)

  if (D_hess > 0)
    return(list(mi = s^2 / D_hess, epc = s / D_hess, opg = FALSE))

  # Non-positive Schur complement (Part 7.2): recompute with the
  # outer-product information, positive semi-definite by construction. The
  # per-case theta-score reuses the two evaluations already taken above.
  s_case  <- w * (llp_case - llm_case) / (2 * h_th)
  Bthth   <- sum(s_case^2)
  Bthpsi  <- as.numeric(crossprod(s_case, S1))
  Bpsipsi <- crossprod(S1)
  D_opg   <- Bthth - as.numeric(Bthpsi %*% pinv(Bpsipsi) %*% Bthpsi)
  D_opg   <- max(D_opg, .Machine$double.eps)
  s_opg   <- sum(s_case)

  list(mi = s_opg^2 / D_opg, epc = s_opg / D_opg, opg = TRUE)
}

# ------------------------------------------------------------------------------
# The residual covariance matrix (Part 7.3)
# ------------------------------------------------------------------------------

.local_dep_residual_cov <- function(X, mu, s2, pi_k) {
  mbar <- colSums(pi_k * mu)
  J <- ncol(mu)
  implied <- matrix(0, J, J)
  for (k in seq_along(pi_k))
    implied <- implied + pi_k[k] * outer(mu[k, ] - mbar, mu[k, ] - mbar)
  diag(implied) <- diag(implied) + colSums(pi_k * s2)

  observed <- stats::cov(X, use = "pairwise.complete.obs")
  observed - implied
}

# ------------------------------------------------------------------------------
# The public entry point, called from bivariate_residuals()
# ------------------------------------------------------------------------------

.bivariate_mi_gaussian <- function(object) {
  par <- .step1_pack(object)
  if (is.null(par)) {
    message("Local-dependence statistics are not available for this ",
            "measurement model: it has no unconstrained parameter ",
            "representation (growth and other structured emissions).")
    return(NULL)
  }

  X <- object$data
  w <- object$sample_weights
  mm <- object$mm
  K  <- object$n_components
  J  <- ncol(X)
  mu <- mm$parameters$means
  s2 <- mm$parameters$covariances
  pi_k <- object$weights

  ll0_case <- .step1_ll_case(object, X, par)
  ll0 <- sum(w * ll0_case)

  Hpsipsi <- .local_dep_psi_hessian(object, X, par, w)
  Hpsipsi_inv <- pinv(Hpsipsi)

  p1 <- length(par)
  h1 <- .step1_fd_step * pmax(1, abs(par))
  S1 <- vapply(seq_len(p1), function(m) {
    a <- par; a[m] <- a[m] + h1[m]
    b <- par; b[m] <- b[m] - h1[m]
    w * (.step1_ll_case(object, X, a) - .step1_ll_case(object, X, b)) /
      (2 * h1[m])
  }, numeric(nrow(X)))
  dim(S1) <- c(nrow(X), p1)

  nms <- colnames(X) %||% paste0("Item", seq_len(J))
  mi  <- array(NA_real_, c(J, J, K))
  epc <- array(NA_real_, c(J, J, K))
  opg <- array(FALSE, c(J, J, K))

  # Item-pair variance scale sets the finite-difference step for theta: a
  # residual covariance is bounded in magnitude by the geometric mean of the
  # two items' own variances, so that is what "order one" means for it.
  for (k in seq_len(K)) for (b in seq_len(J)) for (a in seq_len(b - 1L)) {
    h_th <- .step1_fd_step * sqrt(s2[k, a] * s2[k, b])
    st <- .local_dep_statistic(object, X, par, w, ll0, ll0_case,
                               a, b, k, h_th, Hpsipsi_inv, S1)
    mi[b, a, k]  <- st$mi
    epc[b, a, k] <- st$epc
    opg[b, a, k] <- st$opg
  }

  p_value  <- stats::pchisq(mi, df = 1, lower.tail = FALSE)
  obs_cov  <- stats::cov(X, use = "pairwise.complete.obs")
  res_cov  <- .local_dep_residual_cov(X, mu, s2, pi_k)
  implied_cov <- obs_cov - res_cov
  res_cor  <- stats::cov2cor(obs_cov) - stats::cov2cor(implied_cov)

  structure(
    list(mi = mi, epc = epc, p_value = p_value, opg = opg,
         residual_cov = res_cov, residual_cor = res_cor, item_names = nms,
         class_labels = paste("Class", seq_len(K)),
         n_opg = sum(opg, na.rm = TRUE)),
    class = c("bivariate_residuals_gaussian", "bivariate_residuals"))
}

# ------------------------------------------------------------------------------
# Print method
# ------------------------------------------------------------------------------

#' @export
print.bivariate_residuals_gaussian <- function(x, digits = 3, ...) {
  J <- length(x$item_names)
  K <- length(x$class_labels)

  cat("=========================================================\n")
  cat("       LOCAL DEPENDENCE (CONTINUOUS INDICATORS)          \n")
  cat("=========================================================\n")
  cat("Modification index (MI) and expected parameter change (EPC)\n")
  cat("for freeing each pair's within-class residual covariance.\n")
  cat("A large MI, referred to chi-square on 1 df, flags a pair\n")
  cat("whose association the class does not reproduce.\n\n")

  for (k in seq_len(K)) {
    cat(x$class_labels[k], "\n")
    for (b in seq_len(J)) for (a in seq_len(b - 1L)) {
      m <- x$mi[b, a, k]
      if (is.na(m)) next
      cat(sprintf("  %s x %s: MI = %s  EPC = %s  p = %s%s\n",
                  x$item_names[a], x$item_names[b],
                  formatC(m, format = "f", digits = digits),
                  formatC(x$epc[b, a, k], format = "f", digits = digits),
                  formatC(x$p_value[b, a, k], format = "f", digits = 4),
                  if (isTRUE(x$opg[b, a, k])) "  [OPG]" else ""))
    }
  }

  if (x$n_opg > 0)
    cat(sprintf("\n%d statistic(s) used the outer-product information: ",
                x$n_opg),
        "the observed-Hessian Schur complement was not positive there.\n")

  cat("\nResidual covariance (observed - model-implied):\n")
  print(round(x$residual_cov, digits))
  cat("=========================================================\n")
  invisible(x)
}
