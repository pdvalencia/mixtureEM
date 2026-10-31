# ==============================================================================
# Per-case scores for the fit_mixture() finish
# ==============================================================================
#
# The finish (R/qn_finish.R) needs, at every iteration, each case's score: the
# derivative of its log-likelihood in every coordinate of the packed vector
# (R/qn_pack.R). Central differences give it at 2p passes over the data per
# iteration, which is what made the finish cost more than the EM it follows on
# a large model: five classes, twelve items and four groups is p = 256, and
# `group_effects = "both"` on it went from nine minutes to over half an hour
# (RECORDS.md, Part 56).
#
# One E-step gives it instead. A case's log-likelihood is
# l_i = log sum_k pi_ik f_k(x_i), and by Fisher's identity its derivative in
# any parameter is sum_k r_ik d/dtheta [log pi_ik + log f_k(x_i)], r_ik the
# posterior. Every family here is an exponential family in the coordinates it
# is packed on, so d/dtheta log f_k is "observed minus expected" times the
# design: x - p for a binary item's logit, 1[x = c] - p_c for a category's
# log-ratio, x - lambda for a log rate, (x - mu) / sigma^2 for a mean and
# (x - mu)^2 / sigma^2 - 1 for a log standard deviation, r_ik - pi_k for a
# class-size logit, and the same forms through a design matrix for the
# class-membership and outcome regressions. A missing cell contributes nothing.
#
# The functions below walk each model in exactly the order its packing does,
# so column m of the score matrix is coordinate m of the packed vector. A
# parameter the packing holds once for several places (an item or a parameter
# equal across groups, a variance equal across classes) sums its places'
# terms. Anything without scores here -- latent class growth (R/lcga.R),
# which carries few parameters and a cheap likelihood -- returns NULL, and the
# finish differences that fit as before; no fit is ever part analytic, part
# differenced. Every family is checked
# against central differences to 1e-5 (tests/testthat/test-qn-mixture-scores.R).

# Observed-cell mask and the data with missing cells at zero.
.qn_obs <- function(X) {
  X <- as.matrix(X)
  M <- (!is.na(X)) * 1
  X[is.na(X)] <- 0
  list(X = X, M = M)
}

# One flat measurement model's scores for all of its items, in the order
# .qn_pack_flat(emis) packs them, with each column's item and the part of the
# packing it belongs to (`mean` or `var` for a Gaussian, `p` otherwise), so a
# composite model can pick out or sum an item's columns.
.qn_flat_scores <- function(emis, X, R) {
  fam <- .step1_family(emis)
  # The binary growth mixture (R/gmm_binary.R): integrated over its nodes,
  # in its own packed order; it is never a block of a composite model.
  if (fam == "binary_growth")
    return(list(S = .qn_bg_scores(emis, X, R), item = NA_integer_, part = "p"))
  if (fam == "structured_normal")
    return(list(S = .qn_sn_scores(emis, X, R), item = NA_integer_, part = "p"))
  K   <- emis$n_components
  n   <- nrow(X)
  p   <- emis$parameters
  if (fam %in% c("bernoulli", "poisson", "gaussian_unit", "gaussian_diag")) {
    o  <- .qn_obs(X)
    E  <- switch(fam, bernoulli = p$pis, poisson = p$rates, p$means)
    J  <- ncol(E)
    ix <- list(kk = rep(seq_len(K), J), jj = rep(seq_len(J), each = K))
    Rk <- R[, ix$kk, drop = FALSE]
    D  <- o$X[, ix$jj, drop = FALSE] -
      sweep(o$M[, ix$jj, drop = FALSE], 2, E[cbind(ix$kk, ix$jj)], "*")
    if (fam != "gaussian_diag")
      return(list(S = Rk * D, item = ix$jj, part = rep("p", K * J)))
    v    <- p$covariances[cbind(ix$kk, ix$jj)]
    Smu  <- Rk * sweep(D, 2, v, "/")
    Ssd  <- Rk * (sweep(D * D, 2, v, "/") - o$M[, ix$jj, drop = FALSE])
    if (isTRUE(emis$variances_equal)) {
      Ssd <- t(rowsum(t(Ssd), ix$jj, reorder = TRUE))
      vj  <- seq_len(J)
    } else vj <- ix$jj
    return(list(S = cbind(Smu, Ssd), item = c(ix$jj, vj),
                part = c(rep("mean", K * J), rep("var", length(vj)))))
  }
  if (fam %in% c("multinoulli", "ordinal")) {
    ordinal <- fam == "ordinal"
    Xo <- if (ordinal) one_hot_ragged(as.matrix(X), emis$cats) else
      one_hot(as.matrix(X), emis$max_val, emis$cats)
    J  <- if (ordinal) length(emis$cats) else .step1_n_items(emis)
    P  <- p$pis
    cols <- list(); item <- integer(0)
    for (j in seq_len(J)) {
      cj  <- .qn_cat_cols(emis, j)
      obs <- rowSums(Xo[, cj, drop = FALSE])
      Pj  <- P[, cj, drop = FALSE]
      ref <- .qn_ref(pmax(Pj, 1e-12))
      # .qn_lr_pack()'s K x (M - 1) matrix, column-major: category slot m of
      # class k is that class's m-th category other than its own anchor.
      for (m in seq_len(length(cj) - 1L)) for (k in seq_len(K)) {
        cat_m <- seq_along(cj)[-ref[k]][m]
        cols[[length(cols) + 1L]] <-
          R[, k] * (Xo[, cj[cat_m]] - obs * Pj[k, cat_m])
      }
      item <- c(item, rep(j, K * (length(cj) - 1L)))
    }
    S <- if (length(cols)) do.call(cbind, cols) else matrix(0, n, 0)
    return(list(S = S, item = item, part = rep("p", length(item))))
  }
  NULL
}

# The growth mixture (R/structured_normal.R), in .qn_pack_sn()'s order. Within
# a class and a pattern of observed occasions O a case is normal with mean
# mu_O and covariance Sigma_O, so with a = Sigma_O^-1 (y_O - mu_O) its
# log-density moves by a in the mean and by (a a' - Sigma_O^-1) / 2 in
# Sigma_O. The mean is Lambda (alpha_k + Gamma_k x_i), so a growth mean takes
# Lambda_O' a and a regression coefficient that times x_i. Sigma is
# Lambda_r Psi Lambda_r' + Theta: a Psi element (a, b) takes
# [Lambda_r' (a a' - Sigma^-1) Lambda_r]_ab / 2, twice over off the diagonal
# for the two places a symmetric element sits, times the unit the packing
# divides it by (.qn_sn_psi_unit()); a residual variance, packed as its log
# standard deviation, takes theta_t (a_t^2 - [Sigma^-1]_tt). The covariance is
# factored as log_likelihood.structured_normal() factors it (.sn_chol()), so
# the scores are those of the log-likelihood the finish evaluates.
.qn_sn_scores <- function(ms, X, R) {
  X  <- as.matrix(X)
  n  <- nrow(X); K <- ms$n_components
  L  <- ms$design; p <- ncol(L)
  q  <- length(ms$r_cols); Lr <- L[, ms$r_cols, drop = FALSE]
  m  <- .sn_m(ms); Z <- .sn_z(ms, n)
  lo <- .qn_sn_layout(ms)
  th <- lo$th; nr <- lo$nr; np <- lo$np; off_th <- lo$off_th
  off_psi <- lo$off_psi; li <- lo$li; fl <- lo$fl
  S  <- matrix(0, n, np)

  for (k in seq_len(K)) {
    Sig <- .sn_sigma(ms, k); mu <- .sn_mu(ms, k)
    pc  <- if (q) lo$psi_col(k)
    tr  <- if (nr == 1L) 1L else k
    gb  <- K * p + if (isTRUE(ms$gamma_equal)) 0L else (k - 1L) * p * m
    for (pat in .sn_patterns(X)) {
      cols <- pat$cols; rows <- pat$rows
      if (!length(cols)) next
      ch <- .sn_chol(Sig[cols, cols, drop = FALSE])
      if (is.null(ch)) next
      Si <- chol2inv(ch)
      A  <- .sn_resid(X, mu, rows, cols, pat$Y) %*% Si
      w  <- R[rows, k]
      G  <- (w * A) %*% L[cols, , drop = FALSE]
      S[rows, (seq_len(p) - 1L) * K + k] <- G
      for (c in seq_len(m)) {
        ci <- gb + (c - 1L) * p + seq_len(p)
        S[rows, ci] <- S[rows, ci] + G * Z[rows, c + 1L]
      }
      if (q) {
        Lo <- Lr[cols, , drop = FALSE]
        C  <- A %*% Lo
        M  <- crossprod(Lo, Si %*% Lo)
        for (j in which(!is.na(pc))) {
          ci <- off_psi + pc[j]
          S[rows, ci] <- S[rows, ci] +
            w * fl[j] * (C[, li[j, 1L]] * C[, li[j, 2L]] - M[li[j, 1L], li[j, 2L]])
        }
      }
      Tt <- w * sweep(A * A, 2, diag(Si), "-") *
        rep(ms$parameters$theta[k, cols], each = length(rows))
      if (length(th$cols) == 1L) {
        S[rows, off_th + tr] <- S[rows, off_th + tr] + rowSums(Tt)
      } else {
        ci <- off_th + (cols - 1L) * nr + tr
        S[rows, ci] <- S[rows, ci] + Tt
      }
    }
  }
  S
}

# Where .qn_pack_sn() puts each piece: the offsets of the Psi and Theta
# blocks, and, for the Psi block, each lower-triangle slot's coordinate class
# by class (NA where the element is fixed) with its multiplier `fl`. A slot
# (a, b) moves Sigma by fl (l_a l_b' + l_b l_a'), l_a column a of Lambda_r,
# which is half its unit on the diagonal and the whole unit off it.
.qn_sn_layout <- function(ms) {
  K  <- ms$n_components
  p  <- ncol(ms$design); q <- length(ms$r_cols); m <- .sn_m(ms)
  th <- .qn_sn_theta_dims(ms); nr <- length(th$rows)
  np <- length(.qn_pack_sn(ms)$par)
  out <- list(th = th, nr = nr, np = np, off_th = np - nr * length(th$cols),
              off_psi = K * p + if (m > 0L)
                (if (isTRUE(ms$gamma_equal)) 1L else K) * p * m else 0L)
  if (!q) return(out)
  li <- which(lower.tri(diag(q), diag = TRUE), arr.ind = TRUE)
  nl <- nrow(li); dg <- li[, 1L] == li[, 2L]
  type <- .sn_psi_type(ms)
  sl <- if (!is.null(ms$psi_fixed)) .sn_psi_slots(ms)
  c(out, list(
    li = li,
    fl = ifelse(dg, 0.5, 1) * (ms$qn_psi_unit %||% .qn_sn_psi_unit(ms))[li],
    psi_col = function(k) {
      if (!is.null(sl)) return(sl$id[k, ])
      if (identical(type, "variances") && q > 1L) {
        out <- integer(nl)
        out[dg]  <- (li[dg, 1L] - 1L) * K + k
        out[!dg] <- K * q + seq_len(sum(!dg))
        return(out)
      }
      if (identical(type, "equal")) seq_len(nl) else (k - 1L) * nl + seq_len(nl)
    }))
}

# The growth mixture's Hessian in closed form, for the finish's Newton steps
# (R/qn_finish.R), over .qn_pack_sn()'s coordinates and then the class-size
# logits: what differencing the scores 2p times would give, from one pass
# over the data instead. Case i's log-likelihood is
# l_i = log sum_k pi_k f_k(y_i), and with l_ik = log pi_k + log f_k(y_i), g_ik
# its gradient, r_ik the posterior and s_i = sum_k r_ik g_ik the score row,
#
#   d2 l_i = sum_k r_ik (d2 l_ik + g_ik g_ik') - s_i s_i'
#
# (Louis, 1982). The second part takes each class's own score rows. The first
# is the normal's Hessian. Within a class and a pattern of observed
# occasions, with a = Sigma^-1 (y - mu), a mean coordinate moving mu by m_a
# and a covariance coordinate moving Sigma by D_b (and, for a log standard
# deviation, by D_bb in its own second derivative):
#
#   mean, mean  -m_a' Sigma^-1 m_b
#   mean, cov   -m_a' Sigma^-1 D_b a
#   cov, cov     tr(Sigma^-1 D_b Sigma^-1 D_c) / 2 - a' D_b Sigma^-1 D_c a
#                + [b = c] (a' D_bb a - tr(Sigma^-1 D_bb)) / 2
#
# Summed over the pattern's cases with weights r_ik w_i these need only
# N = sum r w, abar = sum r w a and Q = sum r w a a': T x T arithmetic per
# class and pattern. A class's score rows touch only its own coordinates and the shared ones, so they are built
# here on those alone, and the score matrix, sum_k r_ik g_ik, is the same
# rows summed, so it is returned too (`S`, in .qn_mixture_score_matrix()'s
# order) and the finish needs no second pass for it. Without growth
# covariates only; NULL otherwise. `w` is the case weights.
.qn_sn_hessian <- function(ms, X, R, w, weights) {
  if (.sn_m(ms) > 0L) return(NULL)
  X  <- as.matrix(X)
  K  <- ms$n_components
  L  <- ms$design; p <- ncol(L)
  q  <- length(ms$r_cols); Lr <- L[, ms$r_cols, drop = FALSE]
  lo <- .qn_sn_layout(ms)
  np <- lo$np; nt <- np + K - 1L
  one_col <- length(lo$th$cols) == 1L
  iw <- np + seq_len(K - 1L)
  pats <- .sn_patterns(X)
  H  <- matrix(0, nt, nt)
  S  <- matrix(0, nrow(X), nt)
  for (k in seq_len(K)) {
    Sig <- .sn_sigma(ms, k); mu <- .sn_mu(ms, k)
    pc  <- if (q) lo$psi_col(k)
    tr  <- if (lo$nr == 1L) 1L else k
    th  <- ms$parameters$theta[k, ]
    rk  <- R[, k] * w
    im  <- (seq_len(p) - 1L) * K + k
    ek  <- -weights[-K]; if (k < K) ek[k] <- ek[k] + 1
    for (pat in pats) {
      cols <- pat$cols; rows <- pat$rows
      if (!length(cols)) next
      ch <- .sn_chol(Sig[cols, cols, drop = FALSE])
      if (is.null(ch)) next
      Si <- chol2inv(ch)
      To <- length(cols)
      rw <- rk[rows]
      N  <- sum(rw)
      A  <- .sn_resid(X, mu, rows, cols, pat$Y) %*% Si
      ab <- colSums(rw * A); Q <- crossprod(A, rw * A)
      Lo <- L[cols, , drop = FALSE]
      # The covariance coordinates this class and pattern move: their place
      # in the packing, D_b, and D_bb (zero for a Psi element).
      id <- integer(0); D <- list(); D2 <- list()
      if (q) {
        Ll <- Lr[cols, , drop = FALSE]
        for (j in which(!is.na(pc))) {
          a <- lo$li[j, 1L]; b <- lo$li[j, 2L]
          id <- c(id, lo$off_psi + pc[j])
          D[[length(D) + 1L]]  <- lo$fl[j] * (tcrossprod(Ll[, a], Ll[, b]) +
                                              tcrossprod(Ll[, b], Ll[, a]))
          D2[[length(D)]] <- matrix(0, To, To)
        }
      }
      if (one_col) {
        id <- c(id, lo$off_th + tr)
        D[[length(D) + 1L]] <- diag(2 * th[cols], To)
        D2[[length(D)]]     <- diag(4 * th[cols], To)
      } else for (t in seq_along(cols)) {
        id <- c(id, lo$off_th + (cols[t] - 1L) * lo$nr + tr)
        E  <- matrix(0, To, To); E[t, t] <- 1
        D[[length(D) + 1L]] <- 2 * th[cols[t]] * E
        D2[[length(D)]]     <- 4 * th[cols[t]] * E
      }
      nc  <- length(D)
      SD  <- lapply(D, function(Db) Si %*% Db)
      Hcc <- matrix(0, nc, nc)
      for (b in seq_len(nc)) for (c in b:nc) {
        v <- N / 2 * sum(SD[[b]] * t(SD[[c]])) - sum((D[[b]] %*% SD[[c]]) * Q)
        if (b == c) v <- v + (sum(D2[[b]] * Q) - N * sum(Si * D2[[b]])) / 2
        Hcc[b, c] <- Hcc[c, b] <- v
      }
      Hmc <- -crossprod(Lo, vapply(SD, function(M) as.vector(M %*% ab),
                                   numeric(To)))
      Hmm <- -N * crossprod(Lo, Si %*% Lo)
      H[im, im] <- H[im, im] + Hmm
      H[im, id] <- H[im, id] + Hmc
      H[id, im] <- H[id, im] + t(Hmc)
      H[id, id] <- H[id, id] + Hcc
      # The class's own score rows on its coordinates -- a' m_a for a mean,
      # (a' D_b a - tr(Sigma^-1 D_b)) / 2 for a covariance, every quadratic
      # form one product of the cases' a a' with the stacked D_b, and
      # 1[k] - pi, the same for every case, for the class-size logits -- and
      # their outer products weighted by r_ik w_i.
      # Each D_b is symmetric, so the quadratic forms need the upper triangle
      # of a a' only, an off-diagonal cell counted twice.
      ut <- which(upper.tri(diag(To), diag = TRUE), arr.ind = TRUE)
      Dm <- matrix(vapply(D, function(Db) Db[ut] * (2 - (ut[, 1L] == ut[, 2L])),
                          numeric(nrow(ut))), nrow(ut))
      qf <- (A[, ut[, 1L], drop = FALSE] * A[, ut[, 2L], drop = FALSE]) %*% Dm
      G  <- cbind(A %*% Lo, (qf - rep(vapply(D, function(Db) sum(Si * Db),
                                              numeric(1)), each = nrow(qf))) / 2)
      ix <- c(im, id)
      cw <- colSums(rw * G)
      H[ix, ix] <- H[ix, ix] + crossprod(G, rw * G)
      H[ix, iw] <- H[ix, iw] + outer(cw, ek)
      H[iw, ix] <- H[iw, ix] + outer(ek, cw)
      H[iw, iw] <- H[iw, iw] + N * tcrossprod(ek)
      S[rows, ix] <- S[rows, ix] + R[rows, k] * G
    }
  }
  S[, iw] <- sweep(R[, -K, drop = FALSE], 2, weights[-K], "-")
  # log pi_k in the logits: -(diag(pi) - pi pi'), the same for every class.
  pw <- weights[-K]
  H[iw, iw] <- H[iw, iw] - sum(w) * (diag(pw, K - 1L) - tcrossprod(pw))
  list(H = H - crossprod(S, w * S), S = S)
}

# A sub-model of a composite model, flat or mixed (`nested`), with its items
# numbered as the composite numbers them.
.qn_items_scores <- function(sub, X, R) {
  if (!inherits(sub, "nested")) return(.qn_flat_scores(sub, X, R))
  out <- list(S = matrix(0, nrow(X), 0), item = integer(0), part = character(0))
  at  <- 0L
  for (nm in names(sub$models)) {
    nj <- sub$columns_per_model[[nm]]
    s  <- .qn_flat_scores(sub$models[[nm]], X[, at + seq_len(nj), drop = FALSE], R)
    if (is.null(s)) return(NULL)
    out <- list(S = cbind(out$S, s$S), item = c(out$item, s$item + at),
                part = c(out$part, s$part))
    at <- at + nj
  }
  out
}

# `blocks` (multiple-group or repeated-measures): each block's items on its
# own columns. An item held equal across blocks is packed by block one alone,
# so its coordinates collect every block's terms; the same for a Gaussian
# parameter held equal (`invariant_params`), which packs means and variances
# as separate parts.
.qn_blocks_scores <- function(mm, X, R) {
  J  <- mm$n_items
  Sb <- lapply(seq_len(mm$n_blocks), function(b)
    .qn_items_scores(mm$models[[b]], X[, .time_block_cols(b, J), drop = FALSE], R))
  if (any(vapply(Sb, is.null, logical(1)))) return(NULL)
  ip <- mm$invariant_params
  if (length(ip)) {
    shared <- c(mean = "means" %in% ip, var = "covariances" %in% ip)
    one <- Sb[[1L]]
    for (b in seq_along(Sb)[-1L]) for (pt in names(shared)[shared]) {
      sel <- one$part == pt
      one$S[, sel] <- one$S[, sel] + Sb[[b]]$S[, Sb[[b]]$part == pt]
    }
    parts <- lapply(seq_along(Sb), function(b) {
      s <- if (b == 1L) one else Sb[[b]]
      keep <- c(mean = b == 1L || !shared[["mean"]],
                var  = b == 1L || !shared[["var"]])
      do.call(cbind, lapply(c("mean", "var"), function(pt)
        if (keep[[pt]]) s$S[, s$part == pt, drop = FALSE]))
    })
    return(do.call(cbind, parts))
  }
  inv <- mm$invariant_items
  one <- Sb[[1L]]
  if (length(inv)) {
    sel <- one$item %in% inv
    for (b in seq_along(Sb)[-1L])
      one$S[, sel] <- one$S[, sel] + Sb[[b]]$S[, Sb[[b]]$item %in% inv]
  }
  do.call(cbind, c(list(one$S), lapply(Sb[-1L], function(s)
    s$S[, !s$item %in% inv, drop = FALSE])))
}

# Binary items with covariate effects (R/categorical_dif.R): each case's logit
# is theta_kj plus its covariate pattern's shift, so x - p uses the case's own
# p; a slope's score is that times the covariate, per class or summed over
# classes for a uniform slope, in .dif_pack()'s order.
.qn_dif_scores <- function(mm, X, R) {
  K   <- mm$n_components
  J   <- ncol(mm$parameters$pis)
  dif <- mm$dif
  o   <- .qn_obs(X)
  kk  <- rep(seq_len(K), J); jj <- rep(seq_len(J), each = K)
  th  <- qlogis(pmax(pmin(mm$parameters$pis, 1 - 1e-15), 1e-15))
  G   <- matrix(0, nrow(X), K * J)
  for (q in seq_len(dif$P)) {
    r  <- dif$rows[[q]]
    pq <- plogis(th + .dif_shift(mm, q))
    G[r, ] <- o$X[r, jj, drop = FALSE] -
      sweep(o$M[r, jj, drop = FALSE], 2, pq[cbind(kk, jj)], "*")
  }
  G  <- R[, kk, drop = FALSE] * G
  sl <- list()
  for (j in seq_len(J)) for (d in which(dif$free[j, ])) {
    Gj <- G[, jj == j, drop = FALSE] * dif$Z[, d]
    sl[[length(sl) + 1L]] <- if (dif$by_class[j]) Gj else rowSums(Gj)
  }
  do.call(cbind, c(list(G), sl))
}

# Ordinal items with covariate effects (R/ordinal_dif.R), in
# .ordinal_dif_pack()'s order: every threshold in increment form (class
# fastest, item by item), then the slopes. With F_s = P(x > s) = plogis(eta_s)
# and P_c = F_{c-1} - F_c, a case in category c has d log P_c / d eta_s =
# F_s (1 - F_s) (1[c = s + 1] - 1[c = s]) / P_c. The first threshold moves
# every eta_s, the log gap of threshold s moves eta_s and those above it by
# -exp(gap), and a slope moves every eta_s by the case's covariate value.
.qn_ordinal_dif_scores <- function(mm, X, R) {
  K    <- mm$n_components
  cats <- mm$cats
  J    <- length(cats)
  dif  <- mm$dif
  X    <- as.matrix(X)
  n    <- nrow(X)
  th   <- .ordinal_theta_from_pis(mm$parameters$pis, cats)
  St   <- matrix(0, n, K * ncol(th))
  sl   <- list()
  for (j in seq_len(J)) {
    m  <- cats[j] - 1L
    tc <- .ordinal_theta_cols(cats, j)
    u  <- th[, tc, drop = FALSE]
    x  <- X[, j]
    G  <- array(0, c(n, K, m))
    for (q in seq_len(dif$P)) {
      r <- dif$rows[[q]]
      r <- r[!is.na(x[r])]
      if (!length(r)) next
      eta <- u
      if (m > 1L) for (s in 2:m) eta[, s] <- eta[, s - 1L] - exp(u[, s])
      eta <- eta + .dif_shift(mm, q)[, j]
      Fm  <- plogis(eta)
      P   <- pmax(cbind(1, Fm) - cbind(Fm, 0), 1e-300)
      xc  <- x[r]
      iP  <- 1 / t(P[, xc, drop = FALSE])
      for (s in seq_len(m))
        G[r, , s] <- sweep(iP, 2, Fm[, s] * (1 - Fm[, s]), "*") *
          ((xc == s + 1L) - (xc == s))
    }
    # cum[[s]]: the case's derivative summed over eta_s and every eta above.
    cum <- lapply(seq_len(m), function(s) matrix(G[, , s], n, K))
    if (m > 1L) for (s in (m - 1L):1L) cum[[s]] <- cum[[s]] + cum[[s + 1L]]
    tot <- cum[[1L]]
    for (s in seq_len(m)) {
      Ts <- if (s == 1L) tot else -sweep(cum[[s]], 2, exp(u[, s]), "*")
      St[, (tc[s] - 1L) * K + seq_len(K)] <- R * Ts
    }
    for (d in which(dif$free[j, ])) {
      Gj <- R * tot * dif$Z[, d]
      sl[[length(sl) + 1L]] <- if (dif$by_class[j]) Gj else rowSums(Gj)
    }
  }
  do.call(cbind, c(list(St), sl))
}

# Any measurement model, in .qn_pack_mm()'s order.
.qn_mm_scores <- function(mm, X, R) {
  if (inherits(mm, "bernoulli_dif")) return(.qn_dif_scores(mm, X, R))
  if (inherits(mm, "ordinal_dif")) return(.qn_ordinal_dif_scores(mm, X, R))
  if (inherits(mm, "blocks")) return(.qn_blocks_scores(mm, X, R))
  if (inherits(mm, "nested")) {
    out <- .qn_items_scores(mm, X, R)
    return(if (is.null(out)) NULL else out$S)
  }
  .qn_flat_scores(mm, X, R)$S
}

# ------------------------------------------------------------------------------
# Structural models, in .qn_pack_sm()'s order
# ------------------------------------------------------------------------------

# Multinomial outcome scores for one class: 1[y = m + 1] - P_m, times the
# design, for the M - 1 non-reference categories (distal_forward()'s
# reference is the first).
.qn_mnl_resid <- function(P, y, valid) {
  M <- ncol(P)
  Yo <- matrix(0, length(y), M)
  if (any(valid)) Yo[cbind(which(valid), y[valid])] <- 1
  (Yo - P * valid)[, -1L, drop = FALSE]
}

.qn_sm_scores <- function(sm, Y, R) {
  p <- sm$parameters
  K <- sm$n_components
  n <- nrow(R)
  if (inherits(sm, "covariate")) {
    Z  <- complete_covariates(as.matrix(Y))
    Zm <- if (isTRUE(sm$intercept)) cbind(1, Z) else Z
    D  <- ncol(Zm)
    E  <- R - softmax_rows(Zm %*% t(p$beta))
    return(E[, rep(seq_len(K - 1L), each = D), drop = FALSE] *
             Zm[, rep(seq_len(D), K - 1L), drop = FALSE])
  }
  if (inherits(sm, "group_prevalence")) {
    g  <- as.integer(as.matrix(Y)[, 1L])
    gm <- p$gamma
    G  <- nrow(gm)
    S  <- sm$frozen
    if (length(S) == K)
      return(sweep(R[, -K, drop = FALSE], 2, gm[1L, -K], "-"))
    per_group <- function(E, m) {
      out <- matrix(0, n, G * m)
      for (h in seq_len(G)) {
        r <- g == h
        out[r, (h - 1L) * m + seq_len(m)] <- E[[h]](r)
      }
      out
    }
    if (!length(S))
      return(per_group(lapply(seq_len(G), function(h) function(r)
        sweep(R[r, -K, drop = FALSE], 2, gm[h, -K], "-")), K - 1L))
    F    <- setdiff(seq_len(K), S)
    rest <- 1 - sum(gm[1L, S])
    Ssh  <- sweep(R[, S, drop = FALSE], 2, gm[1L, S], "-")
    if (length(F) < 2L) return(Ssh)
    RF   <- rowSums(R[, F, drop = FALSE])
    Ff   <- F[-length(F)]
    Ssp  <- per_group(lapply(seq_len(G), function(h) function(r)
      R[r, Ff, drop = FALSE] - outer(RF[r], gm[h, Ff] / rest)), length(Ff))
    return(cbind(Ssh, Ssp))
  }
  if (inherits(sm, "distal_continuous")) {
    y <- as.numeric(as.matrix(Y)[, 1L]); ok <- !is.na(y); y[!ok] <- 0
    mu <- p$means[, 1L]; v <- p$covariances[, 1L]
    e  <- outer(y, mu, "-") * ok
    Sv <- R * (sweep(e * e, 2, v, "/") - ok)
    if (identical(sm$variances, "equal")) Sv <- rowSums(Sv)
    return(cbind(R * sweep(e, 2, v, "/"), Sv))
  }
  if (inherits(sm, "distal_continuous_regression")) {
    Yf <- as.matrix(Y)
    y  <- as.numeric(Yf[, 1L]); ok <- !is.na(y); y[!ok] <- 0
    Z  <- cbind(1, complete_covariates(Yf[, -1L, drop = FALSE]))
    B  <- p$betas
    v  <- p$covariances[1L]
    e  <- (y - Z %*% t(B)) * ok
    Sb <- (R * e / v)[, rep(seq_len(K), ncol(Z)), drop = FALSE] *
      Z[, rep(seq_len(ncol(Z)), each = K), drop = FALSE]
    return(cbind(Sb, rowSums(R * (e * e / v - ok))))
  }
  if (inherits(sm, "distal_continuous_pooled")) {
    Yf  <- as.matrix(Y)
    y   <- as.numeric(Yf[, 1L]); ok <- !is.na(y); y[!ok] <- 0
    Z   <- complete_covariates(Yf[, -1L, drop = FALSE])
    mod <- sm$moderated %||% integer(0)
    pooled <- setdiff(seq_len(ncol(Z)), mod)
    th  <- as.vector(p$beta_pooled)
    Lc  <- length(th)
    v   <- as.vector(p$covariances)
    Sb  <- matrix(0, n, Lc)
    Ssd <- matrix(0, n, K)
    for (k in seq_len(K)) {
      U <- matrix(0, n, Lc)
      U[, k] <- 1
      if (length(pooled)) U[, K + seq_along(pooled)] <- Z[, pooled, drop = FALSE]
      for (j in seq_along(mod)) U[, K + length(pooled) + (j - 1L) * K + k] <- Z[, mod[j]]
      e  <- (y - as.vector(U %*% th)) * ok
      Sb <- Sb + (R[, k] * e / v[k]) * U
      Ssd[, k] <- R[, k] * (e * e / v[k] - ok)
    }
    return(cbind(Sb, if (identical(sm$variances, "class_specific")) Ssd else
      rowSums(Ssd)))
  }
  if (inherits(sm, "distal_regression")) {
    B  <- p$betas
    Yf <- as.matrix(Y)
    y  <- .validate_distal_Y(Yf[, 1L], "distal_regression scores")
    if (is.null(y)) return(matrix(0, n, length(B)))
    ok <- !is.na(y)
    Z  <- cbind(1, complete_covariates(Yf[, -1L, drop = FALSE]))
    M1 <- dim(B)[2]; D <- dim(B)[3]
    S  <- matrix(0, n, length(B))
    for (k in seq_len(K)) {
      Pk <- pmax(pmin(distal_forward(Z, matrix(B[k, , ], M1, D)), 1 - 1e-15), 1e-15)
      Ek <- R[, k] * .qn_mnl_resid(Pk, y, ok)
      for (d in seq_len(D)) for (m in seq_len(M1))
        S[, ((d - 1L) * M1 + (m - 1L)) * K + k] <- Ek[, m] * Z[, d]
    }
    return(S)
  }
  if (inherits(sm, "distal_pooled")) {
    B  <- p$beta_pooled
    Yf <- as.matrix(Y)
    y  <- .validate_pooled_Y(Yf[, 1L], "distal_pooled scores")
    if (is.null(y)) return(matrix(0, n, length(B)))
    ok <- !is.na(y)
    Z  <- complete_covariates(Yf[, -1L, drop = FALSE])
    M1 <- nrow(B); Lc <- ncol(B)
    S  <- matrix(0, n, length(B))
    for (k in seq_len(K)) {
      U <- matrix(0, n, Lc)
      U[, k] <- 1
      if (Lc > K) U[, (K + 1L):Lc] <- Z
      Pk <- pmax(pmin(distal_forward(U, B), 1 - 1e-15), 1e-15)
      Ek <- R[, k] * .qn_mnl_resid(Pk, y, ok)
      S  <- S + Ek[, rep(seq_len(M1), Lc), drop = FALSE] *
        U[, rep(seq_len(Lc), each = M1), drop = FALSE]
    }
    return(S)
  }
  NULL
}

# ------------------------------------------------------------------------------
# The whole fit
# ------------------------------------------------------------------------------

# Each case's joint log-density with each class, as .qn_mixture_problem()'s
# case_ll() sums it: the measurement model, the structural model, and the
# pooled class sizes unless the structural model supplies them.
.qn_mixture_logdens <- function(st, lay, X, Yp) {
  L <- log_likelihood(st$mm, X)
  if (lay$has_sm) L <- L + log_likelihood(st$sm, Yp)
  if (!lay$sm_probs) L <- sweep(L, 2, log(pmax(st$weights, 1e-300)), "+")
  L
}

# The per-case log-likelihood and posteriors at `st`, from one pass.
.qn_mixture_estep <- function(st, lay, X, Yp) {
  L  <- .qn_mixture_logdens(st, lay, X, Yp)
  ll <- logsumexp(L, MARGIN = 1)
  list(ll = ll, R = exp(L - ll))
}

# The n x p score matrix in `lay`'s packed order, or NULL where a block has no
# scores. A block the finish holds fixed (`frozen`, the two-step estimator's
# measurement model) is never read, so it is filled with zeros rather than
# computed -- which also lets a two-step fit on a growth model finish
# analytically.
.qn_mixture_score_matrix <- function(st, lay, X, Yp, R) {
  frozen <- st$frozen
  zero <- function(nm) matrix(0, nrow(R), length(lay$idx[[nm]]))
  Smm <- if ("mm" %in% frozen) zero("mm") else .qn_mm_scores(st$mm, X, R)
  if (is.null(Smm)) return(NULL)
  parts <- list(Smm)
  if (!is.null(lay$idx$weights)) {
    K <- st$n_components
    parts <- c(parts, list(if ("weights" %in% frozen) zero("weights") else
      sweep(R[, -K, drop = FALSE], 2, st$weights[-K], "-")))
  }
  if (lay$has_sm) {
    Ssm <- .qn_sm_scores(st$sm, Yp, R)
    if (is.null(Ssm)) return(NULL)
    parts <- c(parts, list(Ssm))
  }
  S <- do.call(cbind, parts)
  if (ncol(S) != length(lay$par)) return(NULL)
  dimnames(S) <- NULL
  S
}

# Both at once, for callers that want the scores at one point.
.qn_mixture_scores <- function(st, lay, X, Yp) {
  e <- .qn_mixture_estep(st, lay, X, Yp)
  S <- .qn_mixture_score_matrix(st, lay, X, Yp, e$R)
  if (is.null(S)) NULL else list(ll = e$ll, S = S)
}

# ------------------------------------------------------------------------------
# The prior's gradient and curvature
# ------------------------------------------------------------------------------
#
# The finish climbs the penalised objective, so each iteration also needs the
# log-prior's gradient and, for the curvature, its diagonal (R/qn_finish.R).
# Differencing it costs 2p + 1 evaluations of the prior per iteration, and
# each evaluation recomputes the data marginals the prior is centred on: on a
# four-group model with p = 256 that was 93% of the finish once the scores
# were analytic (RECORDS.md, Part 56). Every prior here is a set of
# pseudo-observations at those marginals (.em_log_prior(), .qn_sm_prior()),
# so its gradient has the same "observed minus expected" form as the scores,
# with the marginal as the observation and alpha / K as the weight, and the
# marginals are computed once, when the finish starts. Same walk, same order,
# same NULL for anything not covered as the scores above; checked against
# central differences in tests/testthat/test-qn-mixture-scores.R.

# A flat sub-model's prior gradient and diagonal curvature over all its
# items, in .qn_pack_flat()'s order, from its marginals `mg` (`m` or `s2`, one
# entry per parameter column) and the stacked-update multiplier `scale`.
.qn_flat_prior_gh <- function(sub, mg, scale = 1) {
  fam <- .step1_family(sub)
  K   <- sub$n_components
  p   <- sub$parameters
  # The growth mixture has no measurement prior (.em_family_log_prior()).
  if (fam == "structured_normal") {
    np <- length(.qn_pack_sn(sub)$par)
    return(list(g = numeric(np), h = numeric(np), item = rep(NA_integer_, np),
                part = rep("p", np)))
  }
  if (fam %in% c("bernoulli", "bernoulli_dif", "poisson", "gaussian_unit",
                 "gaussian_diag")) {
    E  <- switch(fam, poisson = p$rates, gaussian_unit = , gaussian_diag = p$means,
                 p$pis)
    J  <- ncol(E)
    kk <- rep(seq_len(K), J); jj <- rep(seq_len(J), each = K)
    ev <- E[cbind(kk, jj)]
    out <- switch(fam,
      bernoulli = , bernoulli_dif = {
        A <- scale * .bayes_alpha(sub, "categorical") / K
        list(g = A * (mg$m[jj] - ev), h = -A * ev * (1 - ev))
      },
      poisson = {
        # (A)[log lambda - lambda / m] on the log-rate scale.
        A <- scale * .bayes_alpha(sub, "poisson") / K
        w <- ifelse(mg$m[jj] > 0, A / pmax(mg$m[jj], 1e-300), 0)
        list(g = w * (mg$m[jj] - ev), h = -w * ev)
      },
      list(g = numeric(K * J), h = numeric(K * J)))
    out$item <- jj
    out$part <- rep(if (fam == "gaussian_diag") "mean" else "p", K * J)
    if (fam == "bernoulli_dif") {
      n_sl <- length(.dif_pack(sub)) - K * J
      out <- list(g = c(out$g, numeric(n_sl)), h = c(out$h, numeric(n_sl)),
                  item = c(jj, rep(NA_integer_, n_sl)),
                  part = c(out$part, rep("slope", n_sl)))
    }
    if (fam != "gaussian_diag") return(out)
    A  <- scale * .bayes_alpha(sub, "variances") / K
    r  <- mg$s2[jj] / p$covariances[cbind(kk, jj)]
    gs <- A * (r - 1); hs <- -2 * A * r
    vj <- jj
    if (isTRUE(sub$variances_equal)) {
      gs <- as.vector(rowsum(gs, jj)); hs <- as.vector(rowsum(hs, jj))
      vj <- seq_len(J)
    }
    return(list(g = c(out$g, gs), h = c(out$h, hs), item = c(jj, vj),
                part = c(out$part, rep("var", length(vj)))))
  }
  if (fam %in% c("multinoulli", "ordinal")) {
    A <- scale * .bayes_alpha(sub, "categorical") / K
    J <- if (fam == "ordinal") length(sub$cats) else .step1_n_items(sub)
    g <- h <- numeric(0); item <- integer(0)
    for (j in seq_len(J)) {
      cj  <- .qn_cat_cols(sub, j)
      mj  <- mg$m[cj]
      Mj  <- sum(mj)
      Pj  <- p$pis[, cj, drop = FALSE]
      ref <- .qn_ref(pmax(Pj, 1e-12))
      for (m in seq_len(length(cj) - 1L)) for (k in seq_len(K)) {
        cat_m <- seq_along(cj)[-ref[k]][m]
        pk <- Pj[k, cat_m]
        g  <- c(g, A * (mj[cat_m] - pk * Mj))
        h  <- c(h, -A * Mj * pk * (1 - pk))
      }
      item <- c(item, rep(j, K * (length(cj) - 1L)))
    }
    return(list(g = g, h = h, item = item, part = rep("p", length(item))))
  }
  NULL
}

# The marginals a block's prior term is centred on, computed as
# .em_flat_family_log_prior_items() computes them (every item at once: each
# marginal is a function of its own column only).
.qn_prior_marg_items <- function(sub, Xv, wt) {
  switch(class(sub)[1],
    bernoulli = list(m = .em_col_marginal(Xv, wt)),
    bernoulli_nan = , bernoulli_dif = list(m = .em_col_marginal_valid(Xv, wt)),
    multinoulli = , multinoulli_nan = list(m = .em_onehot_marginal(
      Xv, one_hot(Xv, sub$max_val), wt,
      function(j) ((j - 1L) * sub$max_val + 1L):(j * sub$max_val))),
    ordinal = , ordinal_nan = list(m = .em_onehot_marginal(
      Xv, one_hot_ragged(Xv, sub$cats), wt,
      function(j) .ordinal_item_cols(sub, j))),
    poisson = , poisson_nan = {
      m <- .em_col_marginal_valid(Xv, wt)
      m[!is.finite(m)] <- 0
      list(m = m)
    },
    gaussian_diag = , gaussian_diag_nan = list(s2 = .marginal_var(Xv, wt)),
    gaussian_unit = , gaussian_unit_nan = list(),
    NULL)
}

.qn_gh_cols <- function(x, sel)
  list(g = x$g[sel], h = x$h[sel], item = x$item[sel], part = x$part[sel])
.qn_gh_join <- function(xs)
  list(g = unlist(lapply(xs, `[[`, "g")), h = unlist(lapply(xs, `[[`, "h")))

# The measurement prior's builder: marginals now, a function of the model
# later. Mirrors .em_family_log_prior(): per sub-model for `nested`; for
# `blocks`, a free item's term per block on that block's data and an item held
# equal across blocks one pooled term on block one's parameters, or, under
# `invariant_params`, every block's full term with the shared parts summed.
.qn_mm_prior_builder <- function(mm, X, wt, marginals) {
  if (inherits(mm, "blocks")) {
    J  <- mm$n_items
    Bn <- mm$n_blocks
    if (any(vapply(mm$models, inherits, logical(1), "nested"))) return(NULL)
    ip     <- mm$invariant_params
    inv    <- mm$invariant_items
    pooled <- !length(ip) && .em_group_pooled_prior(mm)
    if (pooled || (!length(ip) && length(inv)))
      mg_pool <- .qn_prior_marg_items(mm$models[[1L]],
                                      .strip_block_prefix(.stack_blocks(X, J, Bn)),
                                      if (is.null(wt)) NULL else rep(wt, Bn))
    mg <- if (pooled) rep(list(mg_pool), Bn) else
      lapply(seq_len(Bn), function(b) .qn_prior_marg_items(mm$models[[b]],
        .strip_block_prefix(X[, .time_block_cols(b, J), drop = FALSE]), wt))
    if (any(vapply(mg, is.null, logical(1)))) return(NULL)
    if (!length(ip) && length(inv))
      scale <- if (inherits(mm, "time_blocks")) Bn else 1L
    # A free item on group blocks: 1 / Bn of its prior in each group, centred
    # on the pooled marginal (.em_blocks_family_log_prior(), R/em_core.R).
    free_scale <- if (pooled) 1 / Bn else 1
    return(function(mm) {
      gb <- lapply(seq_len(Bn), function(b)
        .qn_flat_prior_gh(mm$models[[b]], mg[[b]], free_scale))
      if (any(vapply(gb, is.null, logical(1)))) return(NULL)
      if (length(ip)) {
        shared <- c(mean = "means" %in% ip, var = "covariances" %in% ip, p = FALSE)
        one <- gb[[1L]]
        for (b in seq_len(Bn)[-1L]) for (pt in c("mean", "var")[shared[1:2]]) {
          sel <- one$part == pt
          one$g[sel] <- one$g[sel] + gb[[b]]$g[gb[[b]]$part == pt]
          one$h[sel] <- one$h[sel] + gb[[b]]$h[gb[[b]]$part == pt]
        }
        gb[[1L]] <- one
        return(.qn_gh_join(lapply(seq_len(Bn), function(b)
          .qn_gh_cols(gb[[b]], b == 1L | !shared[gb[[b]]$part]))))
      }
      if (length(inv)) {
        pool <- .qn_flat_prior_gh(mm$models[[1L]], mg_pool, scale)
        sel  <- gb[[1L]]$item %in% inv
        gb[[1L]]$g[sel] <- pool$g[sel]
        gb[[1L]]$h[sel] <- pool$h[sel]
      }
      .qn_gh_join(c(list(gb[[1L]]), lapply(gb[-1L], function(s)
        .qn_gh_cols(s, !s$item %in% inv))))
    })
  }
  if (inherits(mm, "nested")) {
    at <- 0L
    mg <- list()
    for (nm in names(mm$models)) {
      nj <- mm$columns_per_model[[nm]]
      mg[[nm]] <- .em_prior_marginals(list(mm = mm$models[[nm]], sample_weights = wt),
                                      X[, at + seq_len(nj), drop = FALSE])
      at <- at + nj
    }
    return(function(mm) {
      parts <- lapply(names(mm$models), function(nm)
        .qn_flat_prior_gh(mm$models[[nm]], mg[[nm]]))
      if (any(vapply(parts, is.null, logical(1)))) return(NULL)
      .qn_gh_join(parts)
    })
  }
  marginals <- marginals %||% .qn_prior_marg_items(mm, X, wt)
  if (is.null(marginals)) return(NULL)
  function(mm) {
    out <- .qn_flat_prior_gh(mm, marginals)
    if (is.null(out)) NULL else out[c("g", "h")]
  }
}

# The structural prior's gradient and curvature: .qn_sm_prior()'s three
# terms. Every other structural model has no prior.
.qn_sm_prior_builder <- function(sm, Y, w = NULL) {
  K <- sm$n_components
  zeros <- function(sm) {
    n <- length(.qn_pack_sm(sm)$par)
    list(g = numeric(n), h = numeric(n))
  }
  if (inherits(sm, "covariate")) {
    a <- .bayes_alpha(sm, "latent")
    if (!(a > 0)) return(zeros)
    Z  <- complete_covariates(as.matrix(Y))
    # Cases listwise deleted for a missing predictor (weight zero) place no
    # prior row, as in .fit_mnl().
    if (!is.null(w) && any(w <= 0)) Z <- Z[w > 0, , drop = FALSE]
    Zm <- unique(if (isTRUE(sm$intercept)) cbind(1, Z) else Z)
    wt <- a / (K * nrow(Zm))
    return(function(sm) {
      P <- softmax_rows(Zm %*% t(sm$parameters$beta))
      list(g = wt * as.vector(crossprod(Zm, 1 - K * P)[, -K, drop = FALSE]),
           h = -wt * K * as.vector(crossprod(Zm^2, P * (1 - P))[, -K, drop = FALSE]))
    })
  }
  if (inherits(sm, "group_prevalence")) {
    a <- .bayes_alpha(sm, "latent")
    if (!(a > 0)) return(zeros)
    G  <- sm$n_groups
    wt <- a / (K * G)
    return(function(sm) {
      gm <- sm$parameters$gamma
      S  <- sm$frozen
      if (length(S) == K) {
        q <- gm[1L, -K]
        return(list(g = wt * G * (1 - K * q), h = -wt * G * K * q * (1 - q)))
      }
      if (!length(S)) {
        q <- as.vector(t(gm[, -K, drop = FALSE]))
        return(list(g = wt * (1 - K * q), h = -wt * K * q * (1 - q)))
      }
      F <- setdiff(seq_len(K), S)
      s <- gm[1L, S]
      g <- wt * G * (1 - K * s); h <- -wt * G * K * s * (1 - s)
      if (length(F) > 1L) {
        q <- as.vector(t(gm[, F[-length(F)], drop = FALSE] / (1 - sum(s))))
        g <- c(g, wt * (1 - length(F) * q))
        h <- c(h, -wt * length(F) * q * (1 - q))
      }
      list(g = g, h = h)
    })
  }
  if (inherits(sm, "distal_categorical") && ncol(as.matrix(Y)) == 1L) {
    a <- .bayes_alpha(sm, "categorical")
    y <- .validate_pooled_Y(as.matrix(Y)[, 1], "distal_categorical prior")
    if (!(a > 0) || is.null(y)) return(zeros)
    ok <- !is.na(y)
    wc <- if (is.null(w)) rep(1, length(y)) else w
    q  <- colSums(distal_one_hot(y[ok], sm$max_val) * wc[ok]) / sum(wc[ok])
    return(function(sm) {
      B <- sm$parameters$beta_pooled
      if (ncol(B) != K) return(NULL)
      P <- distal_forward(diag(K), B)[, -1L, drop = FALSE]   # K x (M - 1)
      list(g = (a / K) * as.vector(t(sweep(-P * sum(q), 2, q[-1L], "+"))),
           h = -(a / K) * sum(q) * as.vector(t(P * (1 - P))))
    })
  }
  zeros
}

# The whole objective's prior as a function of the model: its gradient and
# diagonal curvature in `lay`'s packed order, or NULL where some block is not
# covered, and the finish then differences the prior as before.
.qn_mixture_prior_grad <- function(model_state, lay, X, Y, wt) {
  w    <- model_state$sample_weights
  marg <- tryCatch(.em_prior_marginals(model_state, X), error = function(e) NULL)
  mmf  <- .qn_mm_prior_builder(model_state$mm, X, w, marg)
  if (is.null(mmf)) return(NULL)
  smf  <- if (lay$has_sm) .qn_sm_prior_builder(model_state$sm, Y, wt)
  function(st) {
    parts <- list(mmf(st$mm))
    if (is.null(parts[[1L]])) return(NULL)
    if (!is.null(lay$idx$weights)) {
      K <- st$n_components
      a <- .bayes_alpha(st$mm, "latent")
      q <- st$weights[-K]
      parts <- c(parts, list(if (a > 0) list(g = a / K - a * q, h = -a * q * (1 - q))
                             else list(g = numeric(K - 1L), h = numeric(K - 1L))))
    }
    if (lay$has_sm) {
      s <- smf(st$sm)
      if (is.null(s)) return(NULL)
      parts <- c(parts, list(s))
    }
    out <- .qn_gh_join(parts)
    if (length(out$g) != length(lay$par)) return(NULL)
    out
  }
}
