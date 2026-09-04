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
