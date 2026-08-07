# ==============================================================================
# S3 Time-Block Measurement Model (repeated-measures LCA / LPA)
# ==============================================================================
#
# A repeated-measures latent class model is an ordinary latent class model whose
# indicators are the J items crossed with the T occasions: a person holds one
# class throughout, and the classes are patterns of change (Collins & Lanza,
# 2010, sec. 7.2). The likelihood is therefore the ordinary one,
#
#   P(y_i) = Σ_k γ_k Π_t Π_j ρ_{j t k}(y_ijt),
#
# and this class exists solely to impose equality of ρ across occasions,
#
#   ρ_{j 1 k} = ρ_{j 2 k} = ... = ρ_{j T k},
#
# which is what makes a class label mean the same thing at every occasion.
#
# The structure mirrors `nested` (R/nested.R): T sub-models, one per occasion,
# each covering all J items, with column blocks laid out time-major. Invariance
# is imposed in the M-step by estimating the constrained items from the T
# occasions stacked together and writing that estimate into every occasion's
# parameter block. The complete-data log-likelihood is additive over items and
# occasions, so the stacked estimate is the constrained MLE and no new optimiser
# is required. Keeping one full sub-model per occasion, rather than splitting
# items into shared and free groups, preserves the column order, which is why
# the BLRT sampler, the class-alignment code and the profile plot can treat this
# model exactly as they treat `nested`.

# Shared constructor behind both `time_blocks_model()` (blocks = occasions) and
# `group_blocks_model()` (blocks = observed groups, see R/group_blocks.R). Both
# are "J items x n_blocks copies" models whose only job is to optionally hold
# some items' response parameters equal across the blocks; nothing here cares
# whether a block is a time occasion or a group, so the S3 methods below are
# shared and dispatch on the common `"blocks"` class.
.blocks_model <- function(n_components, n_items, n_blocks,
                          sub_model = "bernoulli",
                          invariant_items = integer(0),
                          invariant_params = character(0),
                          variances_equal = FALSE,
                          max_val = NULL, prefix = "B", extra_class = "blocks") {
  n_items  <- as.integer(n_items)
  n_blocks <- as.integer(n_blocks)
  invariant_items <- sort(unique(as.integer(invariant_items)))
  if (length(invariant_items) &&
      (min(invariant_items) < 1L || max(invariant_items) > n_items))
    stop("`invariant_items` must index items in 1:n_items.", call. = FALSE)

  # See R/blocks_constraints.R for why the item-wise and parameter-wise axes are
  # offered separately and not crossed.
  invariant_params <- .check_invariant_params(invariant_params, sub_model)
  if (length(invariant_items) && length(invariant_params))
    stop("Use either `invariant_items` (whole items held equal across groups) ",
         "or `invariant_params` (one kind of parameter held equal, for every ",
         "item), not both.", call. = FALSE)

  state <- list(
    n_components     = n_components,
    n_items          = n_items,
    n_blocks         = n_blocks,
    sub_model        = sub_model,
    invariant_items  = invariant_items,
    invariant_params = invariant_params,
    max_val          = NULL,         # kept NULL: this model is not itself polytomous
    models           = list()
  )
  # max_val is meaningful only for the polytomous family; the Gaussian
  # constructors take no such argument and would reject it. `variances_equal`
  # travels the other way, and .construct_emission() drops it where it does not
  # apply.
  sub_args <- list(descriptor = sub_model, n_components = n_components,
                   variances_equal = isTRUE(variances_equal))
  if (!is.null(max_val)) sub_args$max_val <- max_val

  for (b in seq_len(n_blocks))
    state$models[[paste0(prefix, b)]] <- do.call(build_emission, sub_args)

  class(state) <- unique(c(extra_class, "blocks", "emission"))
  state
}

# Constructor. `sub_model` is any measurement descriptor understood by
# build_emission(); `invariant_items` are the item indices (1..n_items) whose
# response parameters are held equal across time.
time_blocks_model <- function(n_components, n_items, n_times,
                              sub_model = "bernoulli",
                              invariant_items = integer(0),
                              invariant_params = character(0),
                              variances_equal = FALSE,
                              max_val = NULL, ...) {
  state <- .blocks_model(n_components, n_items, n_blocks = n_times,
                         sub_model = sub_model,
                         invariant_items = invariant_items,
                         invariant_params = invariant_params,
                         variances_equal = variances_equal,
                         max_val = max_val, prefix = "T",
                         extra_class = "time_blocks")
  state$n_times <- state$n_blocks
  state
}

# Columns of a sub-model's parameter matrix that belong to the given items.
# One column per item for Bernoulli/Gaussian; max_val columns per item for the
# one-hot layout used by multinoulli.
.item_param_cols <- function(sub, items) {
  M <- sub$max_val
  if (is.null(M)) return(as.integer(items))
  as.integer(unlist(lapply(items, function(j) ((j - 1L) * M + 1L):(j * M))))
}

# Copy the parameter columns of the given items from `src` into `dst`.
.copy_item_params <- function(dst, src, items) {
  if (!length(items)) return(dst)
  cols <- .item_param_cols(src, items)
  for (nm in c("pis", "means", "covariances")) {
    if (!is.null(src$parameters[[nm]]))
      dst$parameters[[nm]][, cols] <- src$parameters[[nm]][, cols, drop = FALSE]
  }
  dst
}

# Stack the n_blocks blocks of X on top of each other, giving an (n*n_blocks) x
# J matrix whose rows repeat the responsibilities and case weights n_blocks
# times. This is the design that yields the across-block constrained estimate
# in one call, whether a block is a time occasion or an observed group.
.stack_blocks <- function(X, n_items, n_blocks) {
  do.call(rbind, lapply(seq_len(n_blocks),
                        function(b) X[, .time_block_cols(b, n_items), drop = FALSE]))
}

# Drop the block tag a padded matrix carries on its column names, so a block's
# own parameters come back labelled with the item rather than with the item and
# the block. .pad_group_blocks() and its time-axis counterpart name columns
# "G1.anxiety" / "T2.anxiety"; the block table already has the block in its
# heading, so repeating it in every row is noise. A matrix with no column names,
# or names in any other shape, is returned untouched.
.strip_block_prefix <- function(X) {
  nm <- colnames(X)
  if (is.null(nm)) return(X)
  colnames(X) <- sub("^[GT][0-9]+\\.", "", nm)
  X
}

#' @exportS3Method
init_params.blocks <- function(model_state, X, resp, random_state = NULL, ...) {
  J <- model_state$n_items
  for (b in seq_len(model_state$n_blocks)) {
    X_sub <- .strip_block_prefix(X[, .time_block_cols(b, J), drop = FALSE])
    model_state$models[[b]] <-
      init_params(model_state$models[[b]], X_sub, resp, random_state)
  }
  # Start on the constraint surface so the first E-step already reflects it.
  inv <- model_state$invariant_items
  if (length(inv) && model_state$n_blocks > 1L) {
    for (b in 2:model_state$n_blocks)
      model_state$models[[b]] <-
        .copy_item_params(model_state$models[[b]], model_state$models[[1]], inv)
  }
  .project_invariant_params(model_state)
}

#' @exportS3Method
m_step.blocks <- function(model_state, X, resp, weights = NULL, ...) {
  # A parameter-wise constraint needs its own conditional-maximisation step; the
  # stacked update below is exact only when a shared item shares every one of
  # its parameters. See R/blocks_constraints.R.
  if (length(model_state$invariant_params))
    return(.blocks_gaussian_ecm(model_state, X, resp, weights = weights, ...))

  J   <- model_state$n_items
  Bn  <- model_state$n_blocks
  inv <- model_state$invariant_items
  all_invariant <- length(inv) == J

  # `resp` is either one n x K matrix, used at every block (repeated-measures
  # LCA and multiple-group LCA, where class membership is fixed across the
  # blocks a case can even belong to), or a list of matrices (latent
  # transition analysis, where the posterior over statuses differs by
  # occasion). Both give the same constrained estimator; only the
  # responsibilities attached to each block differ.
  resp_at <- function(b) if (is.list(resp)) resp[[b]] else resp

  # Free items: an ordinary per-block update.
  if (!all_invariant) {
    for (b in seq_len(Bn)) {
      X_sub <- .strip_block_prefix(X[, .time_block_cols(b, J), drop = FALSE])
      model_state$models[[b]] <-
        m_step(model_state$models[[b]], X_sub, resp_at(b), weights = weights, ...)
    }
  }

  # Invariant items: one update on the stacked blocks, copied everywhere.
  if (length(inv)) {
    X_pool    <- .strip_block_prefix(.stack_blocks(X, J, Bn))
    resp_pool <- do.call(rbind, lapply(seq_len(Bn), resp_at))
    w_pool    <- if (is.null(weights)) NULL else rep(weights, Bn)
    pooled    <- m_step(model_state$models[[1]], X_pool, resp_pool,
                        weights = w_pool, ...)
    for (b in seq_len(Bn))
      model_state$models[[b]] <-
        if (all_invariant) pooled
        else .copy_item_params(model_state$models[[b]], pooled, inv)
  }

  model_state
}

#' @exportS3Method
log_likelihood.blocks <- function(model_state, X, ...) {
  J <- model_state$n_items
  log_eps <- matrix(0, nrow = nrow(X), ncol = model_state$n_components)
  for (b in seq_len(model_state$n_blocks)) {
    X_sub   <- X[, .time_block_cols(b, J), drop = FALSE]
    log_eps <- log_eps + log_likelihood(model_state$models[[b]], X_sub)
  }
  log_eps
}

# ------------------------------------------------------------------------------
# Flat view for the L-BFGS refinement in em_core.R
# ------------------------------------------------------------------------------
#
# refine_lbfgs() works on one wide K x (J*n_blocks) parameter matrix. These two
# helpers translate between that view and the per-block sub-models, and build
# the `col_map` that ties columns held equal across blocks to a single free
# parameter. Returns NULL for any sub-model the refinement does not support, in
# which case refine_lbfgs() leaves the fit as EM produced it.
.refine_time_block_view <- function(mm) {
  if (!inherits(mm, "blocks")) return(NULL)
  # `col_map` can tie two stacked columns to one free column, which is exactly
  # an item held equal across blocks. A parameter held equal is a constraint on
  # part of a column and has no expression here, so such a model is left as EM
  # produced it (and EM therefore runs to the tight stopping rule).
  if (length(mm$invariant_params)) return(NULL)

  sub1 <- mm$models[[1]]
  if (!class(sub1)[1] %in% .refine_supported) return(NULL)
  # Polytomous items would need a per-category column map; they are outside the
  # refinement whitelist anyway, so this only guards against future additions.
  if (!is.null(sub1$max_val)) return(NULL)

  J   <- mm$n_items
  Bn  <- mm$n_blocks
  inv <- mm$invariant_items

  col_map   <- integer(J * Bn)
  next_id   <- 0L
  shared_id <- integer(J)
  for (j in inv) { next_id <- next_id + 1L; shared_id[j] <- next_id }
  for (b in seq_len(Bn)) for (j in seq_len(J)) {
    pos <- (b - 1L) * J + j
    if (j %in% inv) col_map[pos] <- shared_id[j]
    else { next_id <- next_id + 1L; col_map[pos] <- next_id }
  }

  flat <- sub1
  flat$parameters <- list()
  for (nm in c("pis", "means", "covariances")) {
    if (!is.null(sub1$parameters[[nm]]))
      flat$parameters[[nm]] <- do.call(
        cbind, lapply(mm$models, function(s) s$parameters[[nm]]))
  }

  list(mm = flat, col_map = col_map)
}

# Write a wide K x (J*n_blocks) refined parameter set back into the sub-models.
.refine_time_block_write <- function(mm, refined) {
  J <- mm$n_items
  for (b in seq_len(mm$n_blocks)) {
    cols <- .time_block_cols(b, J)
    for (nm in names(refined))
      mm$models[[b]]$parameters[[nm]] <- refined[[nm]][, cols, drop = FALSE]
  }
  mm
}

#' @exportS3Method
n_parameters.blocks <- function(model_state, ...) {
  J  <- model_state$n_items
  Bn <- model_state$n_blocks

  # Parameter-wise invariance: count each parameter matrix once if it is shared
  # across the blocks and B times if it is not. The per-matrix sizes come from
  # the sub-model, so an across-class variance constraint inside a block is
  # already reflected in the `covariances` count.
  inv_p <- model_state$invariant_params
  if (length(inv_p)) {
    sub <- model_state$models[[1]]
    K   <- model_state$n_components
    sizes <- c(means = length(sub$parameters$means),
               covariances = if (isTRUE(sub$variances_equal))
                 length(sub$parameters$covariances) / K
               else length(sub$parameters$covariances))
    sizes <- sizes[!is.na(sizes) & sizes > 0]
    return(sum(vapply(names(sizes), function(nm)
      sizes[[nm]] * (if (nm %in% inv_p) 1L else Bn), numeric(1))))
  }

  n_inv    <- length(model_state$invariant_items)
  per_item <- n_parameters(model_state$models[[1]]) / J
  per_item * (n_inv + (J - n_inv) * model_state$n_blocks)
}
