# ==============================================================================
# Parameter-wise invariance across blocks (continuous indicators)
# ==============================================================================
#
# `invariant_items` (R/time_blocks.R) holds an item equal across the blocks:
# every parameter of that item — for a continuous indicator, both its class
# means and its class variances — is shared. That is the right constraint for
# the questions Collins & Lanza ask of categorical indicators, where an item has
# only one kind of parameter, and it is the wrong shape for continuous ones,
# where the applied literature routinely frees one kind and holds the other.
#
# Olivera-Aguilar & Rikoon (2018) is the case in point. Their "unconstrained"
# multiple-group latent profile model — which they note is also what most
# software fits by default — frees the class *means* across groups while holding
# the indicator *variances* invariant, and it is that model, not the fully
# heterogeneous one, that the invariance literature compares against. An
# item-wise constraint cannot express it: an item is either wholly free or
# wholly shared. `invariant_params` adds the other axis, naming the parameter
# matrices held equal across blocks for *all* items.
#
# The two axes are deliberately not combinable. Their intersection would be a
# per-item, per-parameter mask, which is a different (and much larger) interface
# than either, and nothing in the literature this implements asks for one.
#
# ------------------------------------------------------------------------------
# Why this needs its own M-step rather than the stacked update
# ------------------------------------------------------------------------------
#
# An item held equal across blocks is estimated by `m_step.blocks()` in one call
# on the blocks stacked on top of each other, which works because *every*
# parameter of that item is shared: the stacked data are then exactly the data
# that parameter sees. Sharing only some parameters breaks that. If the means
# are free per block and the variance shared, the stacked update would compute
# the sum of squares about the *stacked* mean, which is the wrong quantity — it
# carries the between-block spread of the means into the variance.
#
# The constrained maximiser is instead a conditional-maximisation (ECM) step,
# and for a diagonal Gaussian it is available in closed form on each half:
#
#   means shared across blocks, variances free:
#       mu_kj = sum_b (S1_kjb / v_kjb) / sum_b (n_kjb / v_kjb)
#     — a precision-weighted average, because a block that measures the item
#       more sharply says more about the shared mean. (With the variances also
#       shared this reduces to the ordinary pooled mean, as it must.)
#
#   variances shared across blocks, means free:
#       v_kj = (sum_b SS_kjb(mu_kjb) + prior) / (sum_b n_kjb + prior weight)
#     — each block's sum of squares taken about *its own* mean, which is the
#       part the stacked update gets wrong.
#
# Conditional on the other half each is the exact constrained maximiser, so the
# pair increases the penalised complete-data objective at every iteration and EM
# stays monotone; it is Meng & Rubin's ECM, not an approximation.
#
# The prior follows the parameters. Each block's emission contributes its own
# alpha/K pseudo-observations at its own observed marginal variance, so when B
# blocks share one variance the constrained optimum sums B of those terms — the
# same accounting `.pool_variances_over_classes()` uses across classes, applied
# on the other index. This is what keeps a constrained fit and its unconstrained
# counterpart on one objective, which is the only reason their log-likelihoods
# may be differenced.

# The parameter matrices `invariant_params` may name, and the emissions that
# understand them. Deliberately short: this is a constraint on continuous
# indicators, and a categorical block's `pis` is already covered item-wise by
# `invariant_items`, which for a Bernoulli item is the same constraint.
.blocks_invariant_param_names <- c("means", "covariances")
# Both spellings, since a block's `sub_model` may arrive as either the internal
# emission name or the user-facing alias `.longitudinal_measurement_spec()` uses.
.blocks_invariant_supported   <- c("gaussian_diag", "gaussian_diag_nan",
                                   "gaussian_unit", "gaussian_unit_nan",
                                   "continuous", "continuous_nan",
                                   "gaussian", "gaussian_nan")

.check_invariant_params <- function(invariant_params, sub_model) {
  invariant_params <- unique(as.character(invariant_params))
  if (!length(invariant_params)) return(character(0))

  bad <- setdiff(invariant_params, .blocks_invariant_param_names)
  if (length(bad))
    stop("`invariant_params` must name parameter matrices held equal across ",
         "groups, from: ", paste(sQuote(.blocks_invariant_param_names), collapse = ", "),
         ". Got: ", paste(sQuote(bad), collapse = ", "), ".", call. = FALSE)

  if (is.character(sub_model) && length(sub_model) == 1L &&
      !sub_model %in% .blocks_invariant_supported)
    stop("`invariant_params` is only available for continuous indicators ",
         "(the model here is ", sQuote(sub_model), "). For categorical ",
         "indicators an item has a single kind of parameter, so holding a ",
         "parameter equal across groups and holding the item equal across ",
         "groups are the same constraint: use `group_invariant_items`.",
         call. = FALSE)

  if ("covariances" %in% invariant_params &&
      is.character(sub_model) && grepl("^gaussian(_unit)?(_nan)?$", sub_model))
    stop("`invariant_params = \"covariances\"` has nothing to constrain: ",
         sQuote(sub_model), " fixes every variance at 1.", call. = FALSE)

  intersect(.blocks_invariant_param_names, invariant_params)
}

# Weighted per-(class, item) sufficient statistics for one block: the weight
# behind each cell and the first two moments, over that column's observed cells
# only, so FIML and the structurally-missing other-group blocks are handled by
# the same masking.
.block_gaussian_stats <- function(Xb, R) {
  K <- ncol(R); J <- ncol(Xb)
  n <- S1 <- S2 <- matrix(0, nrow = K, ncol = J)
  for (j in seq_len(J)) {
    ok <- !is.na(Xb[, j])
    if (!any(ok)) next
    Rj <- R[ok, , drop = FALSE]
    xj <- Xb[ok, j]
    n[, j]  <- colSums(Rj)
    S1[, j] <- as.vector(crossprod(Rj, xj))
    S2[, j] <- as.vector(crossprod(Rj, xj^2))
  }
  list(n = n, S1 = S1, S2 = S2)
}

# The constrained M-step. Called by m_step.blocks() whenever
# `invariant_params` is non-empty; the unconstrained parameters are updated here
# too, so this replaces the per-block loop rather than following it.
.blocks_gaussian_ecm <- function(model_state, X, resp, weights = NULL, ...) {
  J   <- model_state$n_items
  Bn  <- model_state$n_blocks
  K   <- model_state$n_components
  inv <- model_state$invariant_params
  shared_mean <- "means" %in% inv
  shared_var  <- "covariances" %in% inv
  var_equal   <- isTRUE(model_state$models[[1]]$variances_equal)
  has_var     <- !is.null(model_state$models[[1]]$parameters$covariances)

  resp_at <- function(b) if (is.list(resp)) resp[[b]] else resp

  stats <- vector("list", Bn)
  s2    <- matrix(1, nrow = Bn, ncol = J)      # per-block marginal variances
  for (b in seq_len(Bn)) {
    cols <- .time_block_cols(b, J)
    Xb <- X[, cols, drop = FALSE]
    R  <- resp_at(b)
    if (!is.null(weights)) R <- sweep(R, 1, weights, "*")
    stats[[b]] <- .block_gaussian_stats(Xb, R)
    s2[b, ]    <- .marginal_var(Xb, weights)
  }

  # --- CM-step 1: means, given the current variances ---------------------------
  old_mu <- lapply(model_state$models, function(m) m$parameters$means)
  var_at <- function(b) {
    v <- model_state$models[[b]]$parameters$covariances
    if (is.null(v)) matrix(1, nrow = K, ncol = J) else v
  }

  if (shared_mean) {
    num <- den <- matrix(0, nrow = K, ncol = J)
    for (b in seq_len(Bn)) {
      v <- var_at(b)
      num <- num + stats[[b]]$S1 / v
      den <- den + stats[[b]]$n  / v
    }
    mu <- num / den
    # A class-item cell no block observes leaves 0/0; the data say nothing about
    # it, so it keeps the value it had.
    mu[!is.finite(mu)] <- old_mu[[1]][!is.finite(mu)]
    mu_at <- rep(list(mu), Bn)
  } else {
    mu_at <- lapply(seq_len(Bn), function(b) {
      m <- stats[[b]]$S1 / stats[[b]]$n
      m[!is.finite(m)] <- old_mu[[b]][!is.finite(m)]
      m
    })
  }

  # --- CM-step 2: variances, given the new means ------------------------------
  if (has_var) {
    alpha     <- .bayes_alpha(model_state$models[[1]], "variances")
    prior_obs <- alpha / K
    ss <- lapply(seq_len(Bn), function(b) {
      st <- stats[[b]]
      st$S2 - 2 * mu_at[[b]] * st$S1 + mu_at[[b]]^2 * st$n
    })

    if (shared_var) {
      # One variance per (class, item) for the whole model: sum the sums of
      # squares over blocks, and with them the blocks' prior terms.
      SS <- Reduce(`+`, ss)
      NN <- Reduce(`+`, lapply(stats, `[[`, "n"))
      s2_tot <- colSums(s2)                       # sum_b s2_jb
      if (var_equal) {
        pooled <- (colSums(SS) + K * prior_obs * s2_tot) /
                  (colSums(NN) + K * Bn * prior_obs)
        V <- matrix(pooled, nrow = K, ncol = J, byrow = TRUE)
      } else {
        V <- (SS + prior_obs * rep(s2_tot, each = K)) / (NN + Bn * prior_obs)
      }
      V_at <- rep(list(V), Bn)
    } else {
      V_at <- lapply(seq_len(Bn), function(b)
        .pool_variances_over_classes(ss[[b]], stats[[b]]$n, prior_obs, s2[b, ],
                                     var_equal))
    }
  }

  for (b in seq_len(Bn)) {
    model_state$models[[b]]$parameters$means <- mu_at[[b]]
    if (has_var)
      model_state$models[[b]]$parameters$covariances <- pmax(V_at[[b]], 1e-12)
  }
  model_state
}

# Put the named parameter matrices of every block onto the constraint surface by
# copying block 1's, so that the first E-step already reflects the constraint.
# Used from init_params.blocks() and from the group warm start, whose per-group
# sub-fits are by construction unconstrained.
.project_invariant_params <- function(model_state, src = NULL) {
  inv <- model_state$invariant_params
  if (!length(inv)) return(model_state)
  ref <- src %||% model_state$models[[1]]
  for (b in seq_len(model_state$n_blocks))
    for (nm in inv)
      if (!is.null(ref$parameters[[nm]]))
        model_state$models[[b]]$parameters[[nm]] <- ref$parameters[[nm]]
  model_state
}
