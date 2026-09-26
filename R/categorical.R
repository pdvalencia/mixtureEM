# ==============================================================================
# S3 Categorical Models (Bernoulli and Multinoulli)
# ==============================================================================

#' Constructor for Categorical models
#'
#' @description
#' Sets up the initial state and class structure for categorical emission models
#' (like Bernoulli or Multinoulli) before the EM algorithm runs.
#'
#' @param n_components Integer. The number of latent classes/components to estimate.
#' @param type Character. The specific distribution type, usually "bernoulli".
#' @param max_val Integer or NULL. The maximum category value (used for multinoulli).
#' @param cats Integer vector or NULL. The number of categories of each item
#'   (used for multinoulli). When NULL, each item's count is its own highest
#'   observed code, or \code{max_val} for every item when that is given.
#' @param ... Additional arguments passed to the method.
#'
#' @return A list object of class \code{c(type, "emission")} containing the model state.
#' @export
categorical_model <- function(n_components, type = "bernoulli", max_val = NULL,
                              cats = NULL, ...) {
  if (!is.null(cats)) {
    cats <- as.integer(cats)
    if (anyNA(cats) || any(cats < 1L))
      stop("`cats` must hold a positive category count for every item.",
           call. = FALSE)
    max_val <- max(c(max_val, cats))
  }
  state <- list(
    n_components = n_components,
    parameters = list(),
    max_val = max_val,
    cats = cats
  )
  class(state) <- c(type, "emission")
  return(state)
}

# ------------------------------------------------------------------------------
# 1. Bernoulli (Binary) S3 Methods
# ------------------------------------------------------------------------------

# Where a restart starts from, for a binary indicator.
#
# The draw looks narrow and uninformed, and replacing it with a start centred
# on each item's observed marginal (a perturbed-marginal start) was tried,
# measured and rejected on 2026-09-08 (Part 47's W2, `internal/RECORDS.md`).
# **Do not re-propose it without reading that entry**, which records three
# separate reasons, any one of which is enough:
#
#   1. It recovers the generating parameters WORSE. The marginal start reaches
#      a higher log-likelihood on nearly every seed and lands further from the
#      truth, which is this package's own recorded lesson about a stronger
#      search finding more spurious maxima, arriving from a new direction.
#   2. .lta_random_start() ends at init_params(), so any change here moves every
#      LTA fit too. The marginal start turned a converged binary random-intercept
#      fixture into a diverging one: log-likelihood -456.58 to -1956.39 with a
#      node intercept at 2.8e15.
#   3. A start centred on the item's marginal is not invariant to response-
#      pattern collapsing, and .lta_collapse()'s documented safety argument
#      (R/lta.R) rests on this function reading nothing but ncol(X). On a
#      four-item table the sample's column means were 0.108/0.200/0.796/0.898
#      and the pattern table's were 0.467/0.467/0.533/0.533.
#
# There is a coherent reason a marginal start does not help here, and it is
# about the objective rather than about either search: our M-step already shrinks every
# response probability toward the item's observed marginal through the
# `categorical` Bayes constant, so this surface is marginal-anchored in a way a
# plain-ML one is not. A start that is also marginal-anchored adds little the
# penalty does not already supply, while a wide perturbation finds the extra
# local maxima a shrunk surface has more of.
#' @exportS3Method
init_params.bernoulli <- function(model_state, X, resp, random_state = NULL, ...) {
  if (!is.null(random_state)) set.seed(random_state)
  model_state$parameters$pis <- matrix(
    runif(model_state$n_components * ncol(X), 0.25, 0.75),
    nrow = model_state$n_components
  )
  return(model_state)
}

#' @exportS3Method
m_step.bernoulli <- function(model_state, X, resp, weights = NULL, alpha = NULL,
                             prior_scale = 1, ...) {
  alpha <- alpha %||% .bayes_alpha(model_state, "categorical")
  if (!is.null(weights)) {
    resp <- sweep(resp, 1, weights, "*")
    marginal_prob <- colSums(sweep(X, 1, weights, "*"), na.rm = TRUE) / sum(weights)
  } else {
    marginal_prob <- colMeans(X, na.rm = TRUE)
  }

  K <- model_state$n_components
  # `prior_scale` is the number of equations this one update stands in for. It
  # is 1 everywhere except the stacked update m_step.blocks() runs for an item
  # held invariant across time blocks, which sees that many occasions' data at
  # once and must carry that many occasions' prior with it.
  prior_obs <- prior_scale * alpha / K

  pis <- t(resp) %*% X

  # Add Dirichlet prior to numerator (proportional to marginal probability)
  pis <- sweep(pis, 2, prior_obs * marginal_prob, "+")

  # Add Dirichlet prior to denominator
  sum_resp <- colSums(resp) + prior_obs
  pis <- sweep(pis, 1, sum_resp, "/")

  model_state$parameters$pis <- pis
  return(model_state)
}

#' @exportS3Method
log_likelihood.bernoulli <- function(model_state, X, ...) {
  p   <- pmax(pmin(model_state$parameters$pis, 1 - 1e-15), 1e-15)
  out <- X %*% t(log(p) - log1p(-p)) + rep(rowSums(log1p(-p)), each = nrow(X))
  dimnames(out) <- NULL
  return(out)
}

#' @exportS3Method
n_parameters.bernoulli <- function(model_state, ...) {
  return(length(model_state$parameters$pis))
}

# ------------------------------------------------------------------------------
# 2. Bernoulli NaN (Binary with Missing Data) S3 Methods
# ------------------------------------------------------------------------------
#' @exportS3Method init_params bernoulli_nan
init_params.bernoulli_nan <- init_params.bernoulli
#' @exportS3Method n_parameters bernoulli_nan
n_parameters.bernoulli_nan <- n_parameters.bernoulli

#' @exportS3Method
m_step.bernoulli_nan <- function(model_state, X, resp, weights = NULL, alpha = NULL,
                                 prior_scale = 1, ...) {
  alpha <- alpha %||% .bayes_alpha(model_state, "categorical")
  if (!is.null(weights)) {
    resp <- sweep(resp, 1, weights, "*")
  }

  K <- model_state$n_components
  prior_obs <- prior_scale * alpha / K
  pis <- matrix(0, nrow = K, ncol = ncol(X),
                dimnames = list(NULL, colnames(X)))

  for (j in seq_len(ncol(X))) {
    valid <- !is.na(X[, j])
    if (any(valid)) {
      resp_valid <- resp[valid, , drop = FALSE]

      if (!is.null(weights)) {
        marg_prob <- sum(X[valid, j] * weights[valid]) / sum(weights[valid])
      } else {
        marg_prob <- mean(X[valid, j])
      }

      num <- t(resp_valid) %*% X[valid, j] + (prior_obs * marg_prob)
      den <- colSums(resp_valid) + prior_obs
      pis[, j] <- num / den
    }
  }
  model_state$parameters$pis <- pis
  return(model_state)
}

#' @exportS3Method
log_likelihood.bernoulli_nan <- function(model_state, X, ...) {
  p   <- pmax(pmin(model_state$parameters$pis, 1 - 1e-15), 1e-15)
  X0  <- X
  X0[is.na(X0)] <- 0
  M   <- (!is.na(X)) * 1
  out <- X0 %*% t(log(p) - log1p(-p)) + M %*% t(log1p(-p))
  dimnames(out) <- NULL
  return(out)
}

# ------------------------------------------------------------------------------
# 3. Multinoulli (Categorical) S3 Methods
# ------------------------------------------------------------------------------

# The parameter matrix keeps max_val columns per item whatever the item's own
# category count, so every caller that walks items in steps of max_val keeps
# working. Item j's categories above cats[j] are padding: no case can give
# them, the M-step estimates them at exactly 0, and n_parameters() does not
# count them. `cats` is absent on a state built before items had their own
# counts, and every item then has max_val.
.multinoulli_cats <- function(ms) {
  if (!is.null(ms$cats)) return(ms$cats)
  J <- if (!is.null(ms$parameters$pis)) ncol(ms$parameters$pis) %/% ms$max_val
       else length(ms$item_names)
  rep(as.integer(ms$max_val), J)
}

# The columns of item j that hold one of its own categories.
.multinoulli_item_cols <- function(ms, j)
  (j - 1L) * ms$max_val + seq_len(.multinoulli_cats(ms)[j])

# Every item's own-category columns, item-major.
.multinoulli_real_cols <- function(ms) {
  cats <- .multinoulli_cats(ms)
  as.integer(unlist(lapply(seq_along(cats), function(j)
    (j - 1L) * ms$max_val + seq_len(cats[j]))))
}

# Each item's category count as the data show it: its highest observed code,
# and at least 2, since an item needs two categories to carry a parameter.
.observed_cats <- function(X) {
  vapply(seq_len(ncol(X)), function(j) {
    v <- X[, j]
    v <- v[!is.na(v)]
    if (!length(v)) 2L else max(2L, as.integer(max(v)))
  }, integer(1))
}

one_hot <- function(X, max_val, cats = NULL) {
  n <- nrow(X)
  D <- ncol(X)
  out <- matrix(0, nrow = n, ncol = D * max_val)

  valid <- !is.na(X)

  # Guard against float values: R silently truncates non-integer subscripts,
  # so 1.9 is treated as category 1 and 2.9 as category 2 with no error.
  # Catch this before it produces silent wrong results.
  if (any(valid)) {
    valid_vals <- X[valid]
    if (!all(valid_vals == as.integer(valid_vals)))
      stop(paste(
        "one_hot: X contains non-integer values in a categorical column.",
        "Categorical items must be integer-valued (e.g. 1L, 2L, 3L).",
        "Did you forget to round or convert to integer?"
      ))
    # And against out-of-range codes, which were far worse than non-integers
    # because nothing caught them. Category 0 in item j indexes column
    # (j-1)*max_val, which for j = 1 is column 0 and is dropped by matrix
    # indexing, and for every later item is the *last category of the previous
    # item* -- so a 0-based coding produced no error and a wrong likelihood.
    if (min(valid_vals) < 1 || max(valid_vals) > max_val)
      stop(sprintf(paste(
        "one_hot: categorical codes must run 1..%d, but values from %g to %g",
        "were found. Categories are 1-based here: recode a 0-based item by",
        "adding 1 (X + 1), and make sure the codes are contiguous."),
        max_val, min(valid_vals), max(valid_vals)))
    # An item with fewer categories than the block's widest one has its own,
    # lower ceiling. A code above it would land in a padding column, whose
    # probability is 0, and turn the likelihood into -Inf without saying why.
    if (!is.null(cats)) {
      over <- which(colSums(X > rep(cats, each = n), na.rm = TRUE) > 0)
      if (length(over))
        stop(sprintf(paste(
          "one_hot: item %d has %d categories, but a code above %d was found.",
          "Categorical codes must run 1..(number of categories) for each item."),
          over[1L], cats[over[1L]], cats[over[1L]]))
    }
  }

  if (any(valid)) {
    idx <- which(valid)
    r   <- (idx - 1L) %% n + 1L
    cb  <- ((idx - 1L) %/% n) * max_val
    out[cbind(r, cb + X[idx])] <- 1
  }
  return(out)
}

#' @exportS3Method
init_params.multinoulli <- function(model_state, X, resp, random_state = NULL, ...) {
  if (!is.null(random_state)) set.seed(random_state)

  # Retain the original item (column) names so summaries can label each
  # polytomous indicator; the one-hot expansion used internally discards them.
  if (!is.null(colnames(X))) model_state$item_names <- colnames(X)

  # Each item's own category count. A caller that knows the response space --
  # the declared levels of a factor, or the count pooled over every group and
  # occasion -- passes it as `cats`; otherwise an explicit `max_val` stands for
  # every item, as it always has, and failing both each item's count is its own
  # highest observed code. Reading nothing but each column's maximum keeps the
  # start unchanged when the data are collapsed to their distinct patterns.
  if (is.null(model_state$cats)) {
    model_state$cats <- if (!is.null(model_state$max_val))
      rep(as.integer(model_state$max_val), ncol(X))
    else
      .observed_cats(X)
  }
  if (length(model_state$cats) != ncol(X))
    stop(sprintf("`cats` has %d entries but there are %d categorical items.",
                 length(model_state$cats), ncol(X)), call. = FALSE)
  model_state$max_val <- max(c(model_state$max_val, model_state$cats))

  M <- model_state$max_val
  n_features <- ncol(X) * M
  pis <- matrix(runif(model_state$n_components * n_features), nrow = model_state$n_components)

  # Padding columns start where the M-step would put them, at 0, so the first
  # E-step already sees only each item's own categories. The draw itself is
  # unchanged, so a block whose items all share one count starts exactly where
  # it did.
  pad <- setdiff(seq_len(n_features), .multinoulli_real_cols(model_state))
  if (length(pad)) pis[, pad] <- 0

  for (j in seq_len(ncol(X))) {
    cols <- ((j - 1) * M + 1):(j * M)
    pis[, cols] <- sweep(pis[, cols, drop=FALSE], 1, rowSums(pis[, cols, drop=FALSE]), "/")
  }
  model_state$parameters$pis <- pis
  return(model_state)
}

#' @exportS3Method
m_step.multinoulli <- function(model_state, X, resp, weights = NULL, alpha = NULL,
                               prior_scale = 1, ...) {
  alpha <- alpha %||% .bayes_alpha(model_state, "categorical")
  if (!is.null(weights)) {
    resp <- sweep(resp, 1, weights, "*")
  }

  if (is.null(model_state$item_names) && !is.null(colnames(X)))
    model_state$item_names <- colnames(X)

  X_oh <- one_hot(X, model_state$max_val)
  K <- model_state$n_components
  prior_obs <- prior_scale * alpha / K

  # Calculate marginal probabilities ignoring NAs
  marginal_prob <- numeric(ncol(X_oh))
  for (j in seq_len(ncol(X))) {
    valid <- !is.na(X[, j])
    cols <- ((j - 1) * model_state$max_val + 1):(j * model_state$max_val)
    if (any(valid)) {
      if (!is.null(weights)) {
        marginal_prob[cols] <- colSums(X_oh[valid, cols, drop=FALSE] * weights[valid]) / sum(weights[valid])
      } else {
        marginal_prob[cols] <- colMeans(X_oh[valid, cols, drop=FALSE])
      }
    }
  }

  pis <- t(resp) %*% X_oh
  pis <- sweep(pis, 2, prior_obs * marginal_prob, "+")

  sum_resp <- colSums(resp) + prior_obs
  pis <- sweep(pis, 1, sum_resp, "/")

  # Re-normalize just to prevent floating point drift
  for (j in seq_len(ncol(X))) {
    cols <- ((j - 1) * model_state$max_val + 1):(j * model_state$max_val)
    pis[, cols] <- sweep(pis[, cols, drop=FALSE], 1, rowSums(pis[, cols, drop=FALSE]), "/")
  }

  model_state$parameters$pis <- pis
  return(model_state)
}

#' @exportS3Method
log_likelihood.multinoulli <- function(model_state, X, ...) {
  p   <- pmax(pmin(model_state$parameters$pis, 1 - 1e-15), 1e-15)
  out <- one_hot(X, model_state$max_val, model_state$cats) %*% t(log(p))
  dimnames(out) <- NULL
  return(out)
}

#' @exportS3Method
n_parameters.multinoulli <- function(model_state, ...) {
  # Each item is charged for its own categories, not the block's widest item's:
  # the padding above cats[j] is fixed at 0 and is not a parameter.
  model_state$n_components * sum(.multinoulli_cats(model_state) - 1L)
}

# ------------------------------------------------------------------------------
# 4. Multinoulli NaN S3 Methods
# ------------------------------------------------------------------------------
#' @exportS3Method init_params multinoulli_nan
init_params.multinoulli_nan <- init_params.multinoulli
#' @exportS3Method n_parameters multinoulli_nan
n_parameters.multinoulli_nan <- n_parameters.multinoulli

# The base multinoulli m-step now cleanly handles NAs via the one_hot zeroing
# and the targeted valid-row marginal calculations!
#' @exportS3Method m_step multinoulli_nan
m_step.multinoulli_nan <- m_step.multinoulli
#' @exportS3Method log_likelihood multinoulli_nan
log_likelihood.multinoulli_nan <- log_likelihood.multinoulli
