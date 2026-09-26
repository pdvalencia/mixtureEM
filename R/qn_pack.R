# ==============================================================================
# Packings for the Newton-type finish
# ==============================================================================
#
# The finish (R/qn_finish.R) works on one unconstrained vector per fit and needs
# only three things from a model: that vector, the model rebuilt from it, and
# what each coordinate is (the `kind` tags documented at .qn_wall_hit()). Every
# measurement and structural model the EM engine fits describes itself here, so
# the finish reaches a new family by adding its two functions to this file and
# nothing else. A family with no packing is refused by .qn_pack_mm() /
# .qn_pack_sm() returning NULL, and the coverage test
# (tests/testthat/test-qn-finish-coverage.R) fails on any class the package can
# build that has neither a packing nor an entry in .qn_refused.
#
# The families R/step3_variance.R already packs for the standard errors are
# packed by those same functions, so the finish moves exactly the parameters
# the standard errors describe; the rest are packed here and used by the finish
# alone. Two rules hold throughout:
#
#   * One coordinate per free parameter. A parameter held equal across classes,
#     occasions or groups is packed once and written back to every place that
#     shares it, so the finish can never step off a constraint the M-step
#     imposes.
#   * Probabilities are log-ratios against the last category, variances are
#     log standard deviations (a covariance matrix, the log-diagonal Cholesky
#     factor), rates are logs, and everything else is as stored.

.qn_part <- function(par, kind)
  list(par = as.numeric(par), kind = rep_len(kind, length(par)))

.qn_join <- function(parts) {
  if (any(vapply(parts, is.null, logical(1)))) return(NULL)
  list(par  = unlist(lapply(parts, `[[`, "par"),  use.names = FALSE),
       kind = unlist(lapply(parts, `[[`, "kind"), use.names = FALSE))
}

.qn_len <- function(pk) length(pk$par)

# Rows of probabilities <-> log-ratios against the last column.
.qn_log_ratio <- function(P) {
  P <- pmax(P, 1e-12)
  log(P[, -ncol(P), drop = FALSE] / P[, ncol(P)])
}
.qn_softmax_last <- function(L) softmax_rows(cbind(L, 0))

# Rows of response probabilities <-> log-ratios against each row's own most
# probable category. The anchor has to be a category the class actually uses:
# anchored on a category at 0 -- the last one, say, of an item with fewer
# categories than the code range, or one a class never gives -- every other
# log-ratio heads to +infinity, all of them hit the wall together, and the
# finish is left unable to move that row at all (a ragged polytomous fit,
# RECORDS.md, Part 53). The anchor is read off the probabilities the finish
# starts from and is fixed for the whole finish, since unpacking always starts
# from that same state.
.qn_ref <- function(P) max.col(P, ties.method = "last")

.qn_lr_pack <- function(P) {
  P   <- pmax(P, 1e-12)
  ref <- .qn_ref(P)
  matrix(vapply(seq_len(nrow(P)), function(k)
    log(P[k, -ref[k]] / P[k, ref[k]]), numeric(ncol(P) - 1L)),
    nrow = nrow(P), byrow = TRUE)
}

.qn_lr_unpack <- function(L, P0) {
  ref <- .qn_ref(pmax(P0, 1e-12))
  M   <- ncol(P0)
  out <- matrix(0, nrow(P0), M)
  for (k in seq_len(nrow(P0))) out[k, -ref[k]] <- L[k, ]
  softmax_rows(out)
}

# A covariance matrix <-> its lower Cholesky factor, diagonal on the log scale.
.qn_pack_chol <- function(P) {
  R <- tryCatch(chol(P), error = function(e) NULL)
  if (is.null(R)) return(NULL)
  L  <- t(R)
  lt <- lower.tri(L, diag = TRUE)
  d  <- (row(L) == col(L))[lt]
  v  <- L[lt]
  v[d] <- log(v[d])
  list(par = v, kind = ifelse(d, "scale", "free"))
}
.qn_unpack_chol <- function(v, q) {
  L  <- matrix(0, q, q)
  lt <- lower.tri(L, diag = TRUE)
  d  <- (row(L) == col(L))[lt]
  v[d] <- exp(v[d])
  L[lt] <- v
  tcrossprod(L)
}

# ------------------------------------------------------------------------------
# Flat measurement models
# ------------------------------------------------------------------------------

.qn_step1_kind <- c(bernoulli = "logit", multinoulli = "logit",
                    poisson = "log_pos", gaussian_unit = "free")

# `items` restricts the packing to some of the model's items, which is how a
# `blocks` model packs an item held equal across blocks exactly once. The
# structured emissions have no items in that sense and are never split.
.qn_pack_flat <- function(emis, items = NULL) {
  fam <- .step1_family(emis)
  if (fam == "multinoulli") return(.qn_pack_multinoulli(emis, items))
  if (fam %in% .step1_families) {
    par <- .step1_pack_sub(emis, items)
    if (is.null(par)) return(NULL)
    if (fam == "gaussian_diag") {
      nc <- length(.item_param_cols(emis, items %||% seq_len(.step1_n_items(emis))))
      nm <- emis$n_components * nc
      return(list(par = par,
                  kind = c(rep("free", nm), rep("scale", length(par) - nm))))
    }
    return(.qn_part(par, .qn_step1_kind[[fam]]))
  }
  if (fam == "ordinal") return(.qn_pack_ordinal(emis, items))
  if (!is.null(items)) return(NULL)
  switch(fam,
    structured_normal = .qn_pack_sn(emis),
    lcga              = .qn_pack_lcga(emis),
    NULL)
}

.qn_unpack_flat <- function(emis, par, items = NULL) {
  fam <- .step1_family(emis)
  if (fam == "multinoulli") return(.qn_unpack_multinoulli(emis, par, items))
  if (fam %in% .step1_families) return(.step1_unpack_sub(emis, par, items))
  switch(fam,
    ordinal           = .qn_unpack_ordinal(emis, par, items),
    structured_normal = .qn_unpack_sn(emis, par),
    lcga              = .qn_unpack_lcga(emis, par),
    emis)
}

# Polytomous items: each item's category probabilities, per class, as
# log-ratios against that class's most probable category (.qn_lr_pack()).
# Ordinal items without a random intercept are the same thing with a ragged
# number of categories per item.
#
# A nominal item packs only its own categories. The padding a shorter item
# carries up to the block's widest item is fixed at 0, so it is no coordinate
# of the finish, and unpacking leaves it untouched.
.qn_cat_cols <- function(emis, j)
  if (is.null(emis$max_val)) .ordinal_item_cols(emis, j) else
    .multinoulli_item_cols(emis, j)

.qn_pack_cat <- function(emis, items) {
  P <- emis$parameters$pis
  .qn_part(unlist(lapply(items, function(j)
    as.vector(.qn_lr_pack(P[, .qn_cat_cols(emis, j), drop = FALSE]))),
    use.names = FALSE), "logit")
}

.qn_unpack_cat <- function(emis, par, items) {
  K  <- emis$n_components
  at <- 0L
  for (j in items) {
    cols <- .qn_cat_cols(emis, j)
    n    <- K * (length(cols) - 1L)
    emis$parameters$pis[, cols] <- .qn_lr_unpack(
      matrix(par[at + seq_len(n)], K), emis$parameters$pis[, cols, drop = FALSE])
    at <- at + n
  }
  emis
}

.qn_pack_multinoulli <- function(emis, items = NULL)
  .qn_pack_cat(emis, items %||% seq_len(.step1_n_items(emis)))
.qn_unpack_multinoulli <- function(emis, par, items = NULL)
  .qn_unpack_cat(emis, par, items %||% seq_len(.step1_n_items(emis)))
.qn_pack_ordinal <- function(emis, items = NULL)
  .qn_pack_cat(emis, items %||% seq_along(emis$cats))
.qn_unpack_ordinal <- function(emis, par, items = NULL)
  .qn_unpack_cat(emis, par, items %||% seq_along(emis$cats))

# Growth mixture (R/structured_normal.R): the growth-factor means, then the
# covariate regressions, then Psi, then Theta. Psi goes through its Cholesky
# factor, which keeps it positive definite on every step; Theta is one log
# standard deviation per distinct residual variance, so `residual_equal` and
# `residual = "constant"` are packed as the constraints they are.
.qn_sn_theta_dims <- function(ms)
  list(rows = if (isTRUE(ms$theta_equal)) 1L else seq_len(ms$n_components),
       cols = if (isTRUE(ms$theta_shared_occasions)) 1L else
         seq_len(nrow(ms$design)))

.qn_pack_sn <- function(ms) {
  p <- ms$parameters
  q <- length(ms$r_cols)
  parts <- list(.qn_part(p$alpha, "free"))
  if (.sn_m(ms) > 0L) {
    G <- if (isTRUE(ms$gamma_equal)) p$gamma[1L] else p$gamma
    parts <- c(parts, list(.qn_part(unlist(G), "free")))
  }
  if (q > 0L) {
    Ps <- if (isTRUE(ms$psi_equal)) p$psi[1L] else p$psi
    parts <- c(parts, lapply(Ps, .qn_pack_chol))
  }
  th <- .qn_sn_theta_dims(ms)
  parts <- c(parts, list(.qn_part(
    0.5 * log(pmax(p$theta[th$rows, th$cols, drop = FALSE], 1e-300)), "scale")))
  .qn_join(parts)
}

.qn_unpack_sn <- function(ms, par) {
  K  <- ms$n_components
  pp <- ncol(ms$design)
  Tn <- nrow(ms$design)
  q  <- length(ms$r_cols)
  m  <- .sn_m(ms)
  at <- 0L
  take <- function(n) { v <- par[at + seq_len(n)]; at <<- at + n; v }

  ms$parameters$alpha <- matrix(take(K * pp), K, pp)
  if (m > 0L) {
    if (isTRUE(ms$gamma_equal)) {
      G <- matrix(take(pp * m), pp, m)
      ms$parameters$gamma <- rep(list(G), K)
    } else {
      for (k in seq_len(K)) ms$parameters$gamma[[k]] <- matrix(take(pp * m), pp, m)
    }
  }
  if (q > 0L) {
    nq <- q * (q + 1L) / 2L
    if (isTRUE(ms$psi_equal)) {
      ms$parameters$psi <- rep(list(.qn_unpack_chol(take(nq), q)), K)
    } else {
      for (k in seq_len(K)) ms$parameters$psi[[k]] <- .qn_unpack_chol(take(nq), q)
    }
  }
  th <- .qn_sn_theta_dims(ms)
  v  <- matrix(exp(2 * take(length(th$rows) * length(th$cols))),
               length(th$rows), length(th$cols))
  ms$parameters$theta <- v[if (length(th$rows) == 1L) rep(1L, K) else seq_len(K),
                           if (length(th$cols) == 1L) rep(1L, Tn) else seq_len(Tn),
                           drop = FALSE]
  ms
}

# Latent class growth (R/lcga.R): one GLM per class, with the Gaussian
# family's residual variance on the log-sd scale.
.qn_pack_lcga <- function(ms) {
  p <- ms$parameters
  parts <- list(.qn_part(p$coefs, "free"))
  if (isTRUE(ms$fam$has_dispersion))
    parts <- c(parts, list(.qn_part(0.5 * log(pmax(p$dispersion, 1e-300)),
                                    "scale")))
  .qn_join(parts)
}

.qn_unpack_lcga <- function(ms, par) {
  C <- ms$parameters$coefs
  n <- length(C)
  ms$parameters$coefs[] <- par[seq_len(n)]
  if (isTRUE(ms$fam$has_dispersion))
    ms$parameters$dispersion <- exp(2 * par[n + seq_along(ms$parameters$dispersion)])
  ms
}

# ------------------------------------------------------------------------------
# Composite measurement models
# ------------------------------------------------------------------------------

# One block's share of a `blocks` model: its own free items, reaching into a
# mixed-measurement block through the item-to-sub-model map.
.qn_pack_items <- function(sub, items) {
  if (!inherits(sub, "nested")) return(.qn_pack_flat(sub, items))
  map <- .item_submodel_map(sub, items)
  .qn_join(lapply(names(sub$models), function(nm)
    .qn_pack_flat(sub$models[[nm]], map[[nm]] %||% integer(0))))
}

.qn_unpack_items <- function(sub, par, items) {
  if (!inherits(sub, "nested")) return(.qn_unpack_flat(sub, par, items))
  map <- .item_submodel_map(sub, items)
  at <- 0L
  for (nm in names(sub$models)) {
    it <- map[[nm]] %||% integer(0)
    n  <- .qn_len(.qn_pack_flat(sub$models[[nm]], it))
    sub$models[[nm]] <- .qn_unpack_flat(sub$models[[nm]], par[at + seq_len(n)], it)
    at <- at + n
  }
  sub
}

# .copy_item_params() (R/time_blocks.R) with the two parameter tables it does
# not know about: Poisson rates and an ordinal item's ragged columns.
.qn_copy_items <- function(dst, src, items) {
  if (!length(items)) return(dst)
  if (inherits(src, "nested")) {
    map <- .item_submodel_map(src, items)
    for (nm in names(map))
      dst$models[[nm]] <- .qn_copy_items(dst$models[[nm]], src$models[[nm]],
                                         map[[nm]])
    return(dst)
  }
  cols <- .em_item_cols(src, items)
  for (nm in c("pis", "means", "covariances", "rates"))
    if (!is.null(src$parameters[[nm]]))
      dst$parameters[[nm]][, cols] <- src$parameters[[nm]][, cols, drop = FALSE]
  dst
}

# `invariant_params` (R/blocks_constraints.R): Gaussian blocks whose means
# and/or variances are shared by every block. Block one carries a shared
# matrix; later blocks carry only what they do not share.
.qn_blocks_param_parts <- function(mm, b) {
  ip  <- mm$invariant_params
  sub <- mm$models[[b]]
  p   <- sub$parameters
  out <- list()
  if (b == 1L || !"means" %in% ip) out$means <- .qn_part(p$means, "free")
  if (.step1_family(sub) == "gaussian_diag" &&
      (b == 1L || !"covariances" %in% ip)) {
    v <- if (isTRUE(sub$variances_equal)) p$covariances[1L, ] else p$covariances
    out$covariances <- .qn_part(0.5 * log(pmax(v, 1e-300)), "scale")
  }
  out
}

.qn_pack_blocks <- function(mm) {
  if (length(mm$invariant_params)) {
    if (!all(vapply(mm$models, .step1_family, character(1)) %in%
             c("gaussian_diag", "gaussian_unit"))) return(NULL)
    return(.qn_join(unlist(lapply(seq_along(mm$models), function(b)
      .qn_blocks_param_parts(mm, b)), recursive = FALSE)))
  }
  free <- .blocks_free_items(mm)
  .qn_join(lapply(seq_along(mm$models), function(b)
    .qn_pack_items(mm$models[[b]], free[[b]])))
}

.qn_unpack_blocks <- function(mm, par) {
  at <- 0L
  if (length(mm$invariant_params)) {
    ip <- mm$invariant_params
    for (b in seq_along(mm$models)) {
      parts <- .qn_blocks_param_parts(mm, b)
      sub   <- mm$models[[b]]
      K     <- sub$n_components
      if (!is.null(parts$means)) {
        n <- .qn_len(parts$means)
        sub$parameters$means[] <- par[at + seq_len(n)]
        at <- at + n
      } else sub$parameters$means <- mm$models[[1L]]$parameters$means
      if (!is.null(parts$covariances)) {
        n <- .qn_len(parts$covariances)
        v <- exp(2 * par[at + seq_len(n)])
        sub$parameters$covariances[] <- if (isTRUE(sub$variances_equal))
          matrix(v, K, length(v), byrow = TRUE) else v
        at <- at + n
      } else if (.step1_family(sub) == "gaussian_diag")
        sub$parameters$covariances <- mm$models[[1L]]$parameters$covariances
      mm$models[[b]] <- sub
    }
    return(mm)
  }
  free <- .blocks_free_items(mm)
  for (b in seq_along(mm$models)) {
    n <- .qn_len(.qn_pack_items(mm$models[[b]], free[[b]]))
    mm$models[[b]] <- .qn_unpack_items(mm$models[[b]], par[at + seq_len(n)],
                                       free[[b]])
    at <- at + n
  }
  # Items held equal across blocks are carried by block one alone.
  if (length(mm$invariant_items) && mm$n_blocks > 1L)
    for (b in 2:mm$n_blocks)
      mm$models[[b]] <- .qn_copy_items(mm$models[[b]], mm$models[[1L]],
                                       mm$invariant_items)
  mm
}

# ------------------------------------------------------------------------------
# Any measurement model
# ------------------------------------------------------------------------------

.qn_pack_mm <- function(mm) {
  if (inherits(mm, "bernoulli_dif")) {
    v <- .dif_pack(mm)
    return(if (is.null(v)) NULL else .qn_part(v, "logit"))
  }
  if (inherits(mm, "blocks")) return(.qn_pack_blocks(mm))
  if (inherits(mm, "nested")) return(.qn_join(lapply(mm$models, .qn_pack_mm)))
  .qn_pack_flat(mm)
}

.qn_unpack_mm <- function(mm, par) {
  if (inherits(mm, "bernoulli_dif")) return(.dif_unpack(mm, par))
  if (inherits(mm, "blocks")) return(.qn_unpack_blocks(mm, par))
  if (inherits(mm, "nested")) {
    at <- 0L
    for (nm in names(mm$models)) {
      n <- .qn_len(.qn_pack_mm(mm$models[[nm]]))
      mm$models[[nm]] <- .qn_unpack_mm(mm$models[[nm]], par[at + seq_len(n)])
      at <- at + n
    }
    return(mm)
  }
  .qn_unpack_flat(mm, par)
}

# ------------------------------------------------------------------------------
# Structural models
# ------------------------------------------------------------------------------

# Class probabilities per group with the classes in `frozen` held to one share
# across groups (R/group_prevalence.R). With nothing frozen, one simplex per
# group; with everything frozen, one simplex; otherwise the frozen shares and
# the remainder as one simplex, then each group's split of the remainder over
# the free classes.
.qn_pack_gp <- function(sm) {
  g  <- sm$parameters$gamma
  S  <- sm$frozen
  K  <- ncol(g)
  if (!length(S)) return(.qn_part(t(.qn_log_ratio(g)), "logit"))
  if (length(S) == K) return(.qn_part(.qn_log_ratio(g[1L, , drop = FALSE]), "logit"))
  F  <- setdiff(seq_len(K), S)
  s  <- g[1L, S]
  rest <- 1 - sum(s)
  shared <- log(pmax(s, 1e-12) / max(rest, 1e-12))
  split  <- if (length(F) > 1L)
    t(.qn_log_ratio(g[, F, drop = FALSE] / rest)) else numeric(0)
  .qn_part(c(shared, split), "logit")
}

.qn_unpack_gp <- function(sm, par) {
  g <- sm$parameters$gamma
  S <- sm$frozen
  G <- nrow(g); K <- ncol(g)
  if (!length(S)) {
    sm$parameters$gamma <- .qn_softmax_last(matrix(par, G, K - 1L, byrow = TRUE))
    return(sm)
  }
  if (length(S) == K) {
    sm$parameters$gamma <- matrix(.qn_softmax_last(matrix(par, 1L)), G, K,
                                  byrow = TRUE)
    return(sm)
  }
  F  <- setdiff(seq_len(K), S)
  sh <- .qn_softmax_last(matrix(par[seq_along(S)], 1L))
  s  <- sh[seq_along(S)]
  rest <- sh[length(S) + 1L]
  split <- if (length(F) > 1L)
    .qn_softmax_last(matrix(par[-seq_along(S)], G, length(F) - 1L, byrow = TRUE))
  else matrix(1, G, 1L)
  gm <- matrix(0, G, K)
  gm[, S] <- matrix(s, G, length(S), byrow = TRUE)
  gm[, F] <- rest * split
  sm$parameters$gamma <- gm
  sm
}

.qn_pack_sm <- function(sm) {
  p <- sm$parameters
  if (inherits(sm, "covariate")) {
    v <- .step1_pack_sm(sm)
    return(if (is.null(v)) NULL else .qn_part(v, "logit"))
  }
  if (inherits(sm, "group_prevalence")) return(.qn_pack_gp(sm))
  if (inherits(sm, "distal_continuous"))
    return(.qn_join(list(.qn_part(p$means, "free"),
                         .qn_part(0.5 * log(p$covariances), "scale"))))
  if (inherits(sm, "distal_continuous_regression"))
    return(.qn_join(list(.qn_part(p$betas, "free"),
                         .qn_part(0.5 * log(p$covariances[1L]), "scale"))))
  if (inherits(sm, "distal_continuous_pooled")) {
    v <- if (identical(sm$variances, "class_specific")) p$covariances else
      p$covariances[1L]
    return(.qn_join(list(.qn_part(p$beta_pooled, "free"),
                         .qn_part(0.5 * log(v), "scale"))))
  }
  # Categorical outcomes: the class intercepts are logits and are walled like
  # any response probability; covariate slopes are left free.
  if (inherits(sm, "distal_regression")) {
    B <- p$betas
    if (!length(B)) return(.qn_part(numeric(0), "free"))
    return(list(par = as.vector(B),
                kind = ifelse(as.vector(slice.index(B, 3L)) == 1L, "logit", "free")))
  }
  if (inherits(sm, "distal_pooled")) {
    B <- p$beta_pooled
    if (!length(B)) return(.qn_part(numeric(0), "free"))
    return(list(par = as.vector(B),
                kind = ifelse(as.vector(col(B)) <= sm$n_components, "logit", "free")))
  }
  NULL
}

.qn_unpack_sm <- function(sm, par) {
  p <- sm$parameters
  if (inherits(sm, "covariate")) return(.step1_unpack_sm(sm, par))
  if (inherits(sm, "group_prevalence")) return(.qn_unpack_gp(sm, par))
  K <- sm$n_components
  if (inherits(sm, "distal_continuous")) {
    sm$parameters$means[] <- par[seq_len(K)]
    sm$parameters$covariances[] <- exp(2 * par[K + seq_len(K)])
    return(sm)
  }
  if (inherits(sm, "distal_continuous_regression")) {
    n <- length(p$betas)
    sm$parameters$betas[] <- par[seq_len(n)]
    sm$parameters$covariances[] <- exp(2 * par[n + 1L])
    return(sm)
  }
  if (inherits(sm, "distal_continuous_pooled")) {
    n <- length(p$beta_pooled)
    sm$parameters$beta_pooled[] <- par[seq_len(n)]
    sm$parameters$covariances[] <- exp(2 * par[-seq_len(n)])
    return(sm)
  }
  if (inherits(sm, "distal_regression")) {
    sm$parameters$betas[] <- par
    return(sm)
  }
  if (inherits(sm, "distal_pooled")) {
    sm$parameters$beta_pooled[] <- par
    return(sm)
  }
  sm
}

# The log-prior a structural M-step adds, as a function of the structural
# model, for the finisher's objective. .em_log_prior() covers the class weights
# and the measurement model and leaves the structural side out, so each
# structural prior is priced here, in the form its own M-step writes it:
#
#   covariate          Dirichlet pseudo-rows at every distinct covariate
#                      pattern (.fit_mnl(), R/covariate.R);
#   group_prevalence   alpha / (K G) pseudo-cases in every group-by-class cell
#                      (.group_prevalence_counts());
#   distal_categorical alpha / K pseudo-observations per class at the
#                      outcome's observed marginal (m_step.distal_pooled()).
#
# The other outcome models carry no prior.
.qn_sm_prior <- function(sm, Y, w = NULL) {
  K <- sm$n_components
  if (inherits(sm, "covariate")) return(.qn_covariate_prior(sm, Y))
  if (inherits(sm, "group_prevalence")) {
    a <- .bayes_alpha(sm, "latent")
    if (!(a > 0)) return(function(sm) 0)
    wt <- a / (K * sm$n_groups)
    return(function(sm) wt * sum(log(pmax(sm$parameters$gamma, 1e-300))))
  }
  if (inherits(sm, "distal_categorical") && ncol(as.matrix(Y)) == 1L) {
    a <- .bayes_alpha(sm, "categorical")
    y <- .validate_pooled_Y(as.matrix(Y)[, 1], "distal_categorical prior")
    if (!(a > 0) || is.null(y)) return(function(sm) 0)
    ok <- !is.na(y)
    wc <- if (is.null(w)) rep(1, length(y)) else w
    M  <- sm$max_val
    q  <- colSums(distal_one_hot(y[ok], M) * wc[ok]) / sum(wc[ok])
    return(function(sm) {
      P <- distal_forward(diag(K), sm$parameters$beta_pooled)
      (a / K) * sum(sweep(log(pmax(P, 1e-300)), 2, q, "*"))
    })
  }
  function(sm) 0
}
