# ==============================================================================
# Direct covariate effects on ordinal indicators of a cross-sectional model
# ==============================================================================
#
# The ordinal twin of R/categorical_dif.R. Item j has S_j ordered categories
# and, in class k, a cumulative logit with free thresholds:
#
#   logit P(Y_ij >= s | X = k, z_i) = theta_kjs + z_i' delta_j      (uniform)
#                                   = theta_kjs + z_i' delta_kj     (by class)
#
# for s = 2..S_j: one proportional-odds slope per item and covariate, shared by
# every threshold, so a positive slope moves the item's responses towards its
# higher categories. With no covariate the thresholds are a reparameterisation
# of the plain ordinal emission's free category probabilities, which is why an
# item named in no pair keeps that emission's closed-form update and a model
# with no pairs is the unconditional fit exactly.
#
# Everything else follows R/categorical_dif.R: the covariates enter through a
# few distinct patterns (`mm$dif`, built by .dif_emission()); `parameters$pis`
# holds the category probabilities of a case whose item covariates are all
# zero, so every reader of `pis` -- and the categorical prior, applied there --
# works unchanged; and `parameters$dif_slopes` is K x (J * D), one K x J block
# per covariate, a uniform slope's rows identical and an unnamed pair zero.
#
# The thresholds are carried in the increment form R/ordinal.R uses for the
# random-intercept model (.ordinal_theta_from_pis(), .ordinal_cat_probs()):
# the first threshold, then the log of each gap, which is unconstrained.

# Item j's category probabilities for every class at pattern q's shift, K x S_j.
.ordinal_dif_item_probs <- function(theta, cats, j, shift)
  .ordinal_cat_probs(theta[, .ordinal_theta_cols(cats, j), drop = FALSE],
                     shift, cats[j])

#' @exportS3Method
init_params.ordinal_dif <- function(model_state, X, resp, random_state = NULL, ...) {
  model_state <- init_params.ordinal(model_state, X, resp, random_state)
  model_state$parameters$dif_slopes <-
    matrix(0, model_state$n_components, ncol(X) * ncol(model_state$dif$Z))
  model_state
}

#' @exportS3Method
log_likelihood.ordinal_dif <- function(model_state, X, ...) {
  dif <- model_state$dif
  if (nrow(X) != length(dif$pat))
    stop("A model with `predictors_items` can only be evaluated on the cases ",
         "it was fitted to, since each case's item probabilities depend on its ",
         "own covariate values.", call. = FALSE)
  cats  <- model_state$cats
  theta <- .ordinal_theta_from_pis(model_state$parameters$pis, cats)
  X_oh  <- one_hot_ragged(X, cats)
  out   <- matrix(0, nrow(X), model_state$n_components)
  for (q in seq_len(dif$P)) {
    r  <- dif$rows[[q]]
    sh <- .dif_shift(model_state, q)
    P  <- do.call(cbind, lapply(seq_along(cats), function(j)
      .ordinal_dif_item_probs(theta, cats, j, sh[, j])))
    out[r, ] <- X_oh[r, , drop = FALSE] %*% t(log(pmax(P, 1e-300)))
  }
  out
}

#' @exportS3Method
n_parameters.ordinal_dif <- function(model_state, ...) {
  K   <- model_state$n_components
  per <- ifelse(model_state$dif$by_class, K, 1L)
  K * sum(model_state$cats - 1L) + sum(rowSums(model_state$dif$free) * per)
}

# The M-step. Items in no named pair: m_step.ordinal()'s closed form, prior
# included. An item in one: the expected category counts aggregated by (class,
# covariate pattern), plus -- when the categorical prior is on -- each class's
# prior pseudo-counts on a row at z = 0, fitted as a cumulative logit with free
# per-class thresholds and the item's slopes. The fit maximises this item's
# part of the expected complete-data objective from the current values, so
# the ECM ascent property holds.
#' @exportS3Method
m_step.ordinal_dif <- function(model_state, X, resp, weights = NULL,
                               alpha = NULL, ...) {
  alpha <- alpha %||% .bayes_alpha(model_state, "categorical")
  dif   <- model_state$dif
  cats  <- model_state$cats
  K     <- model_state$n_components
  J     <- length(cats)
  D     <- ncol(dif$Z)
  S     <- model_state$parameters$dif_slopes
  if (is.null(S)) S <- matrix(0, K, J * D)
  theta_old <- if (!is.null(model_state$parameters$pis))
    .ordinal_theta_from_pis(model_state$parameters$pis, cats)

  plain <- m_step.ordinal(model_state, X, resp, weights = weights, alpha = alpha)
  pis   <- plain$parameters$pis
  if (!is.null(weights)) resp <- resp * weights
  prior_obs <- alpha / K
  X_oh <- one_hot_ragged(X, cats)

  for (j in which(rowSums(dif$free) > 0)) {
    cols  <- .ordinal_item_cols(model_state, j)
    valid <- !is.na(X[, j])
    if (!any(valid)) next
    Sj <- cats[j]
    # N[k, q, s]: expected count of category s among pattern q's cases in k.
    N <- array(0, c(K, dif$P, Sj))
    for (q in seq_len(dif$P)) {
      r <- dif$rows[[q]]
      r <- r[valid[r]]
      if (length(r))
        N[, q, ] <- crossprod(resp[r, , drop = FALSE], X_oh[r, cols, drop = FALSE])
    }
    zrows <- dif$Zu
    if (prior_obs > 0) {
      marg <- if (is.null(weights)) colMeans(X_oh[valid, cols, drop = FALSE])
              else colSums(X_oh[valid, cols, drop = FALSE] * weights[valid]) /
                sum(weights[valid])
      N2 <- array(0, c(K, dif$P + 1L, Sj))
      N2[, seq_len(dif$P), ] <- N
      N2[, dif$P + 1L, ] <- matrix(prior_obs * marg, K, Sj, byrow = TRUE)
      N <- N2
      zrows <- rbind(zrows, 0)
    }
    ds  <- which(dif$free[j, ])
    th0 <- if (is.null(theta_old))
      .ordinal_theta_from_pis(pis, cats)[, .ordinal_theta_cols(cats, j), drop = FALSE]
    else theta_old[, .ordinal_theta_cols(cats, j), drop = FALSE]
    fit <- .ordinal_dif_item_fit(N, zrows, ds, dif$by_class[j], th0,
                                 .dif_item_slopes(S, K, J, j, ds, dif$by_class[j]))
    pis[, cols] <- .ordinal_cat_probs(fit$theta, 0, Sj)
    for (i in seq_along(ds)) {
      col <- .dif_col(J, j, ds[i])
      S[, col] <- if (dif$by_class[j]) fit$b[(i - 1L) * K + seq_len(K)] else fit$b[i]
    }
  }
  model_state$parameters$pis        <- pis
  model_state$parameters$dif_slopes <- S
  model_state
}

# One item's cumulative logit on aggregated counts N (K x Q x S): thresholds
# th0 (K x (S - 1), increment form) and slopes b0 in .dif_item_slopes()'s
# order. The same +/-15 box as .dif_item_newton() on the first threshold and
# the slopes, for the same reason; a log gap is held in [-20, 4], which leaves
# an empty middle category free to close up and an implausible gap of e^4
# logits out of reach.
.ordinal_dif_item_fit <- function(N, zrows, ds, by_class, th0, b0) {
  K  <- dim(N)[1]; Q <- dim(N)[2]; Sj <- dim(N)[3]
  m  <- Sj - 1L
  nb <- length(b0)
  kk <- rep(seq_len(K), Q)                      # rows: (k, q), class fastest
  qq <- rep(seq_len(Q), each = K)
  Nm <- matrix(N, K * Q, Sj)
  keep <- rowSums(Nm) > 0
  # Shift design: row (k, q) -> the slopes it loads on.
  A <- matrix(0, K * Q, nb)
  at <- 0L
  for (d in ds) {
    z <- zrows[qq, d]
    if (by_class) { A[cbind(seq_len(K * Q), at + kk)] <- z; at <- at + K }
    else { A[, at + 1L] <- z; at <- at + 1L }
  }
  Nm <- Nm[keep, , drop = FALSE]; A <- A[keep, , drop = FALSE]; kk <- kk[keep]
  onehot_k <- outer(kk, seq_len(K), "==") * 1
  unpack <- function(v) list(t = matrix(v[seq_len(K * m)], K, m), b = v[K * m + seq_len(nb)])
  fg <- function(v) {
    u  <- unpack(v)
    th <- u$t
    if (m > 1L) for (s in 2:m) th[, s] <- th[, s - 1L] - exp(u$t[, s])
    eta <- th[kk, , drop = FALSE] + drop(A %*% u$b)
    Fm  <- plogis(eta)
    P   <- cbind(1, Fm) - cbind(Fm, 0)
    P   <- pmax(P, 1e-300)
    val <- -sum(Nm * log(P))
    R   <- Nm / P
    G   <- Fm * (1 - Fm) * (R[, -1L, drop = FALSE] - R[, -Sj, drop = FALSE])
    Gk  <- crossprod(onehot_k, G)                # K x m, d/d threshold s
    gt  <- matrix(0, K, m)
    gt[, 1] <- rowSums(Gk)
    if (m > 1L) for (s in 2:m)
      gt[, s] <- -exp(u$t[, s]) * rowSums(Gk[, s:m, drop = FALSE])
    gb <- drop(crossprod(A, rowSums(G)))
    list(val = val, grad = -c(as.vector(gt), gb))
  }
  # The exact Hessian of fg()$val. Per row the objective depends on the
  # boundaries eta_c = theta_c(k) + A b only; its curvature in the boundary
  # heights F_c is tridiagonal (each P_s is a difference of neighbours), and
  # the chain to `v` adds each boundary's own bend f'(eta) and the bend of
  # the increment form (theta_c = t_1 - sum_{u <= c} exp(t_u)).
  nv  <- K * m + nb
  hfn <- function(v) {
    u  <- unpack(v)
    th <- u$t
    if (m > 1L) for (s in 2:m) th[, s] <- th[, s - 1L] - exp(u$t[, s])
    eta <- th[kk, , drop = FALSE] + drop(A %*% u$b)
    Fm  <- plogis(eta)
    f   <- Fm * (1 - Fm)
    P   <- pmax(cbind(1, Fm) - cbind(Fm, 0), 1e-300)
    R1  <- Nm / P
    Q2  <- Nm / P^2
    dF  <- R1[, -Sj, drop = FALSE] - R1[, -1L, drop = FALSE]   # d nll / d F_c
    ge  <- dF * f                                              # d nll / d eta_c
    nr  <- nrow(Nm)
    # d eta_c / d v, one rows x nv matrix per boundary.
    J <- lapply(seq_len(m), function(c) {
      Jc <- matrix(0, nr, nv)
      Jc[cbind(seq_len(nr), kk)] <- 1
      if (c > 1L) for (uu in 2:c)
        Jc[cbind(seq_len(nr), (uu - 1L) * K + kk)] <- -exp(u$t[kk, uu])
      if (nb) Jc[, K * m + seq_len(nb)] <- A
      Jc
    })
    H <- matrix(0, nv, nv)
    for (c in seq_len(m)) {
      dcc <- f[, c]^2 * (Q2[, c] + Q2[, c + 1L]) + dF[, c] * f[, c] * (1 - 2 * Fm[, c])
      H <- H + crossprod(J[[c]], J[[c]] * dcc)
      if (c < m) {
        dcd <- -f[, c] * f[, c + 1L] * Q2[, c + 1L]
        X   <- crossprod(J[[c]], J[[c + 1L]] * dcd)
        H   <- H + X + t(X)
      }
    }
    # The increment form's own bend: d2 theta_c / d t_u^2 = -exp(t_u), u <= c.
    if (m > 1L) for (uu in 2:m) {
      w <- -exp(u$t[kk, uu]) * rowSums(ge[, uu:m, drop = FALSE])
      idx <- (uu - 1L) * K + seq_len(K)
      H[cbind(idx, idx)] <- H[cbind(idx, idx)] + as.vector(crossprod(onehot_k, w))
    }
    H
  }

  lo <- c(rep(-15, K), rep(-20, K * (m - 1L)), rep(-15, nb))
  hi <- c(rep(15, K), rep(4, K * (m - 1L)), rep(15, nb))
  v0 <- pmin(pmax(c(as.vector(th0), b0), lo), hi)
  v  <- .ordinal_dif_newton(v0, fg, hfn, lo, hi)
  if (is.null(v)) {
    res <- stats::nlminb(v0, function(v) fg(v)$val, function(v) fg(v)$grad,
                         lower = lo, upper = hi,
                         control = list(eval.max = 400, iter.max = 300,
                                        rel.tol = 1e-12, x.tol = 1e-10))
    v <- if (is.finite(res$objective) && res$objective <= fg(v0)$val) res$par else v0
  }
  u <- unpack(v)
  list(theta = u$t, b = u$b)
}

# Newton-Raphson for one item's cumulative logit inside the M-step, from the
# previous iteration's values: a full step along -H^-1 g, halved until the
# objective does not rise and the point stays inside the box, until the
# Newton decrement is below 1e-12 of the objective's scale. NULL -- and the
# caller runs its bounded nlminb() fit, the fit it always ran -- when the
# Hessian is not positive definite or cannot be solved, no halving helps (a
# maximum on the box's edge ends here), or 50 steps pass. The nlminb() fit took most of an ordinal DIF model's time: dozens of
# evaluations per item and M-step, where Newton from a warm start needs two
# or three.
.ordinal_dif_newton <- function(v0, fg, hfn, lo, hi, max_iter = 50L) {
  v   <- v0
  cur <- fg(v)
  if (!is.finite(cur$val)) return(NULL)
  for (it in seq_len(max_iter)) {
    H    <- tryCatch(hfn(v), error = function(e) NULL)
    if (is.null(H) || !all(is.finite(H))) return(NULL)
    step <- tryCatch(solve(H, cur$grad), error = function(e) NULL)
    if (is.null(step) || !all(is.finite(step))) return(NULL)
    dec <- sum(cur$grad * step)
    if (!is.finite(dec) || dec < 0) return(NULL)
    if (dec < 1e-12 * (1 + abs(cur$val))) return(v)
    t <- 1
    repeat {
      nv <- v - t * step
      if (all(nv >= lo) && all(nv <= hi)) {
        new <- fg(nv)
        if (is.finite(new$val) && new$val <= cur$val) break
      }
      t <- t / 2
      if (t < 1e-10) return(NULL)
    }
    v <- nv; cur <- new
  }
  NULL
}

# Packing for the joint numerical Hessian (R/step3_variance.R) and the
# Newton-type finish (R/qn_pack.R): every threshold in increment form, item
# by item, then the slopes in .dif_pack()'s order.
.ordinal_dif_pack <- function(mm) {
  K <- mm$n_components
  J <- length(mm$cats)
  S <- mm$parameters$dif_slopes
  sl <- unlist(lapply(seq_len(J), function(j)
    .dif_item_slopes(S, K, J, j, which(mm$dif$free[j, ]), mm$dif$by_class[j])),
    use.names = FALSE)
  c(as.vector(.ordinal_theta_from_pis(mm$parameters$pis, mm$cats)), sl)
}

.ordinal_dif_unpack <- function(mm, par) {
  K    <- mm$n_components
  cats <- mm$cats
  J    <- length(cats)
  nt   <- K * sum(cats - 1L)
  theta <- matrix(par[seq_len(nt)], K, sum(cats - 1L))
  for (j in seq_len(J))
    mm$parameters$pis[, .ordinal_item_cols(mm, j)] <-
      .ordinal_dif_item_probs(theta, cats, j, 0)
  at <- nt
  for (j in seq_len(J)) {
    for (d in which(mm$dif$free[j, ])) {
      col <- .dif_col(J, j, d)
      if (mm$dif$by_class[j]) {
        mm$parameters$dif_slopes[, col] <- par[at + seq_len(K)]; at <- at + K
      } else {
        mm$parameters$dif_slopes[, col] <- par[at + 1L]; at <- at + 1L
      }
    }
  }
  mm
}
