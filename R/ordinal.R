# ==============================================================================
# S3 Ordinal (Cumulative-Logit) Model
# ==============================================================================
#
# The plain (non-random-intercept) half of the ordinal family: a ragged
# generalisation of multinoulli.R's multinoulli, where item j has its own
# category count cats[j] instead of one shared max_val. This is what lets a
# 3/3/2-category item block live in one measurement model instead of needing a
# mixed (`nested`) specification. `pis` stays the primary parameterisation
# here exactly as it does for multinoulli; a random intercept re-expresses the
# same category probabilities through free per-status thresholds and one
# loading per item (see the roadmap's `### 14.18`), which is not built yet.

#' Constructor for the ordinal emission model
#'
#' @param n_components Integer. Number of latent classes/components.
#' @param n_features Integer. Number of items (raw items, not expanded
#'   columns).
#' @param cats Integer vector of length `n_features`, each entry `>= 2`: the
#'   number of ordered categories for that item.
#' @param type Character. `"ordinal"` or `"ordinal_nan"`; both use the same
#'   methods (see the `_nan` aliasing below).
#' @param ... Unused; present for `build_emission()`'s uniform call shape.
#'
#' @return A list object of class `c(type, "emission")`.
ordinal_model <- function(n_components, n_features = NULL, cats, type = "ordinal", ...) {
  cats <- as.integer(cats)
  if (any(cats < 2L))
    stop("`cats` must be at least 2 for every item.", call. = FALSE)
  state <- list(
    n_components = n_components,
    parameters   = list(),
    cats         = cats,
    offsets      = cumsum(c(0L, cats))
  )
  class(state) <- c(type, "emission")
  state
}

# Columns of `parameters$pis` belonging to item j (1-based items).
.ordinal_item_cols <- function(model_state, j) {
  off <- model_state$offsets
  (off[j] + 1L):(off[j + 1L])
}

# Ragged twin of one_hot() (R/categorical.R): item j has its own category
# count cats[j] instead of one shared max_val. Same two guards -- non-integer
# codes, and codes outside 1..cats[j] for item j's own range -- checked
# per-item since the range is no longer uniform. NA becomes an all-zero row,
# which is what makes the emission FIML by construction (log_likelihood and
# m_step below never special-case a missing cell; it simply contributes zero
# to every category indicator).
one_hot_ragged <- function(X, cats) {
  n <- nrow(X)
  D <- length(cats)
  offsets <- cumsum(c(0L, cats))
  out <- matrix(0, nrow = n, ncol = offsets[D + 1L])

  for (j in seq_len(D)) {
    col <- X[, j]
    valid <- !is.na(col)
    if (!any(valid)) next
    vals <- col[valid]
    if (!all(vals == as.integer(vals)))
      stop(paste(
        "one_hot_ragged: X contains non-integer values in a categorical",
        "column. Categorical items must be integer-valued (e.g. 1L, 2L)."))
    if (min(vals) < 1 || max(vals) > cats[j])
      stop(sprintf(paste(
        "one_hot_ragged: item %d's codes must run 1..%d, but values from",
        "%g to %g were found. Categories are 1-based here."),
        j, cats[j], min(vals), max(vals)))
    rows <- which(valid)
    out[cbind(rows, offsets[j] + col[rows])] <- 1
  }
  out
}

#' @exportS3Method
init_params.ordinal <- function(model_state, X, resp, random_state = NULL, ...) {
  if (!is.null(random_state)) set.seed(random_state)
  if (!is.null(colnames(X))) model_state$item_names <- colnames(X)

  n_cols <- model_state$offsets[length(model_state$offsets)]
  pis <- matrix(runif(model_state$n_components * n_cols),
               nrow = model_state$n_components)
  for (j in seq_along(model_state$cats)) {
    cols <- .ordinal_item_cols(model_state, j)
    pis[, cols] <- sweep(pis[, cols, drop = FALSE], 1,
                        rowSums(pis[, cols, drop = FALSE]), "/")
  }
  model_state$parameters$pis <- pis
  model_state
}

#' @exportS3Method
m_step.ordinal <- function(model_state, X, resp, weights = NULL, alpha = NULL,
                           prior_scale = 1, ...) {
  alpha <- alpha %||% .bayes_alpha(model_state, "categorical")
  if (!is.null(weights)) resp <- sweep(resp, 1, weights, "*")

  if (is.null(model_state$item_names) && !is.null(colnames(X)))
    model_state$item_names <- colnames(X)

  cats  <- model_state$cats
  X_oh  <- one_hot_ragged(X, cats)
  K     <- model_state$n_components
  prior_obs <- prior_scale * alpha / K

  marginal_prob <- numeric(ncol(X_oh))
  for (j in seq_along(cats)) {
    valid <- !is.na(X[, j])
    cols  <- .ordinal_item_cols(model_state, j)
    if (!any(valid)) next
    if (!is.null(weights)) {
      marginal_prob[cols] <- colSums(X_oh[valid, cols, drop = FALSE] * weights[valid]) /
        sum(weights[valid])
    } else {
      marginal_prob[cols] <- colMeans(X_oh[valid, cols, drop = FALSE])
    }
  }

  pis <- t(resp) %*% X_oh
  pis <- sweep(pis, 2, prior_obs * marginal_prob, "+")

  sum_resp <- colSums(resp) + prior_obs
  pis <- sweep(pis, 1, sum_resp, "/")

  for (j in seq_along(cats)) {
    cols <- .ordinal_item_cols(model_state, j)
    pis[, cols] <- sweep(pis[, cols, drop = FALSE], 1,
                        rowSums(pis[, cols, drop = FALSE]), "/")
  }

  model_state$parameters$pis <- pis
  model_state
}

#' @exportS3Method
log_likelihood.ordinal <- function(model_state, X, ...) {
  p   <- pmax(pmin(model_state$parameters$pis, 1 - 1e-15), 1e-15)
  out <- one_hot_ragged(X, model_state$cats) %*% t(log(p))
  dimnames(out) <- NULL
  out
}

#' @exportS3Method
n_parameters.ordinal <- function(model_state, ...) {
  model_state$n_components * sum(model_state$cats - 1L)
}

# ------------------------------------------------------------------------------
# Ordinal NaN: one family, always NA-tolerant. one_hot_ragged() already zeros
# a missing cell and the marginal above is computed only over valid rows, so
# there is nothing a "_nan" variant would do differently -- the same aliasing
# multinoulli_nan uses, for the same reason (R/categorical.R).
# ------------------------------------------------------------------------------
#' @exportS3Method init_params ordinal_nan
init_params.ordinal_nan <- init_params.ordinal
#' @exportS3Method m_step ordinal_nan
m_step.ordinal_nan <- m_step.ordinal
#' @exportS3Method log_likelihood ordinal_nan
log_likelihood.ordinal_nan <- log_likelihood.ordinal
#' @exportS3Method n_parameters ordinal_nan
n_parameters.ordinal_nan <- n_parameters.ordinal

# ==============================================================================
# Random-intercept ordinal helpers (roadmap ### 14.18, W4)
# ==============================================================================
#
# Under a random intercept the primary parameterisation switches from `pis`
# (category probabilities) to free per-status thresholds `theta` plus one
# loading per item `lambda` -- see R/lta_ri.R, which owns the M-step and calls
# these. `theta` is stored ragged, item-major, like `pis`, but with `cats - 1`
# columns per item instead of `cats`: item j owns columns
# `.ordinal_theta_cols(cats, j)`, holding `(theta_2, eta_3, ..., eta_S)`, the
# increment parameterisation from ### 14.18.2 (`theta_s = theta_2 -
# sum_{u=3..s} exp(eta_u)`), which is unconstrained and keeps every existing
# per-coordinate L-BFGS box and finite-difference Hessian machinery unchanged
# once ordinal is added to their whitelists (W9/W10, not this session).

.lta_is_ordinal_model <- function(m) inherits(m, c("ordinal", "ordinal_nan"))

.ordinal_theta_offsets <- function(cats) cumsum(c(0L, cats - 1L))

.ordinal_theta_cols <- function(cats, j) {
  off <- .ordinal_theta_offsets(cats)
  (off[j] + 1L):(off[j + 1L])
}

# One item's category probabilities at a single node, given that item's
# threshold block (K x (cats_j - 1), increment-parameterised) and the scalar
# per-node shift `lambda_j' d_q` (the same shift for every status -- the
# loading is shared across statuses, categories and occasions, ### 14.18.2).
.ordinal_cat_probs <- function(theta_j, shift, cats_j) {
  K <- nrow(theta_j)
  S <- cats_j
  th <- matrix(0, K, S - 1L)
  th[, 1] <- theta_j[, 1]
  if (S > 2L)
    for (s in 3:S) th[, s - 1L] <- th[, s - 2L] - exp(theta_j[, s - 1L])
  Fmat <- cbind(1, plogis(th + shift), 0)
  Fmat[, seq_len(S), drop = FALSE] - Fmat[, seq_len(S) + 1L, drop = FALSE]
}

# Every item's category probabilities at node q, cbind()ed into one
# K x sum(cats) block -- the ordinal analogue of
# `plogis(ri$A + matrix(ri$L %*% ri$Dnode[q, ], K, R, byrow = TRUE))`.
.ordinal_node_pis <- function(theta, lambda, Dnode, cats, q) {
  R  <- length(cats)
  dq <- Dnode[q, ]
  blocks <- lapply(seq_len(R), function(j) {
    cols  <- .ordinal_theta_cols(cats, j)
    shift <- sum(lambda[j, ] * dq)
    .ordinal_cat_probs(theta[, cols, drop = FALSE], shift, cats[j])
  })
  do.call(cbind, blocks)
}

# Node-mixed (reported) category probabilities -- the ordinal analogue of
# .lta_ri_integrated_pis()'s binary loop; what gets written into
# state$mm$models[[t]]$parameters$pis.
.ordinal_pis_from_theta <- function(theta, lambda, Dnode, mass, cats) {
  Q   <- length(mass)
  K   <- nrow(theta)
  out <- matrix(0, K, sum(cats))
  for (q in seq_len(Q))
    out <- out + mass[q] * .ordinal_node_pis(theta, lambda, Dnode, cats, q)
  out
}

# Inverse: recover the increment-parameterised thresholds from a plain
# (lambda = 0) ordinal fit's `pis`, for the random-intercept warm start. At
# lambda = 0 every node sees the same category probabilities, so there is
# nothing to integrate over -- this is the single-node case of the forward
# map, run backwards.
.ordinal_theta_from_pis <- function(pis, cats) {
  K   <- nrow(pis)
  R   <- length(cats)
  off <- cumsum(c(0L, cats))
  theta <- matrix(0, K, .ordinal_theta_offsets(cats)[R + 1L])
  for (j in seq_len(R)) {
    cols <- (off[j] + 1L):(off[j + 1L])
    p    <- pis[, cols, drop = FALSE]
    S    <- cats[j]
    cum  <- t(apply(p, 1L, function(row) rev(cumsum(rev(row)))))
    cum  <- cum[, 2:S, drop = FALSE]              # P(U >= s), s = 2..S
    cum  <- pmin(pmax(cum, 1e-10), 1 - 1e-10)
    th   <- qlogis(cum)                            # K x (S - 1), decreasing across s
    out_cols <- .ordinal_theta_cols(cats, j)
    theta[, out_cols[1]] <- th[, 1]
    if (S > 2L)
      for (s in 3:S)
        theta[, out_cols[s - 1L]] <- log(pmax(th[, s - 2L] - th[, s - 1L], 1e-10))
  }
  theta
}

# Weighted observed marginal of item j's S_r categories, pooled over every
# occasion -- the ordinal analogue of .lta_ri_item_marginal() (R/lta_ri.R),
# which returns a scalar (P(U = 1)) rather than a length-S_r vector.
.lta_ri_item_marginal_ordinal <- function(X, w, R, Tn, j, Sj) {
  num <- numeric(Sj); den <- 0
  for (t in seq_len(Tn)) {
    xj  <- X[, .time_block_cols(t, R)[j]]
    obs <- !is.na(xj)
    for (s in seq_len(Sj)) {
      in_s   <- obs & (xj == s)
      num[s] <- num[s] + sum(w[in_s])
    }
    den <- den + sum(w[obs])
  }
  if (den == 0) rep(1 / Sj, Sj) else num / den
}

# Negative log-likelihood of one status's thresholds, thresholds free /
# lambda fixed (ECM cycle 1). `counts` is the (Q x Sj) slice of the
# sufficient-statistic array for this status; `shift` is length Q
# (`Dnode %*% lambda_j`).
.ordinal_theta_negloglik <- function(par, counts, shift, Sj) {
  th <- numeric(Sj - 1L)
  th[1] <- par[1]
  if (Sj > 2L)
    for (s in 2:(Sj - 1L)) th[s] <- th[s - 1L] - exp(par[s])
  Q <- length(shift)
  Fmat <- cbind(1, matrix(plogis(outer(shift, th, "+")), Q, Sj - 1L), 0)
  p <- Fmat[, seq_len(Sj), drop = FALSE] - Fmat[, seq_len(Sj) + 1L, drop = FALSE]
  -sum(counts * log(pmax(p, 1e-12)))
}

.ordinal_newton_theta <- function(par0, counts, shift, Sj) {
  fit <- stats::nlminb(par0, .ordinal_theta_negloglik,
                       counts = counts, shift = shift, Sj = Sj)
  fit$par
}

# Negative log-likelihood of one item's loading, thresholds fixed / lambda
# free (ECM cycle 2), pooling every status and node. `n_arr` is the full
# K x Q x Sj sufficient-statistic array for this item.
.ordinal_lambda_negloglik <- function(lam, theta_j, Dnode, n_arr, Sj) {
  K     <- nrow(theta_j)
  Q     <- nrow(Dnode)
  shift <- as.vector(Dnode %*% lam)
  total <- 0
  for (k in seq_len(K)) {
    th <- numeric(Sj - 1L)
    th[1] <- theta_j[k, 1]
    if (Sj > 2L)
      for (s in 2:(Sj - 1L)) th[s] <- th[s - 1L] - exp(theta_j[k, s])
    Fmat <- cbind(1, matrix(plogis(outer(shift, th, "+")), Q, Sj - 1L), 0)
    p <- Fmat[, seq_len(Sj), drop = FALSE] - Fmat[, seq_len(Sj) + 1L, drop = FALSE]
    total <- total - sum(matrix(n_arr[k, , ], Q, Sj) * log(pmax(p, 1e-12)))
  }
  total
}

.ordinal_newton_lambda <- function(lam0, theta_j, Dnode, n_arr, Sj) {
  fit <- stats::nlminb(lam0, .ordinal_lambda_negloglik,
                       theta_j = theta_j, Dnode = Dnode, n_arr = n_arr, Sj = Sj)
  fit$par
}
