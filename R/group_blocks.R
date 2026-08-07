# ==============================================================================
# S3 Group-Block Measurement Model (multiple-group LCA)
# ==============================================================================
#
# Collins & Lanza's multiple-group LCA (2010, sec. 5.7-5.12) fits the same
# latent class model to data split into G observed groups, and asks whether
# the item-response probabilities are equal across groups (measurement
# invariance, their sec. 5.8) and/or whether the class prevalences are (their
# sec. 5.11-5.12).
#
# The prevalence question is already answered by `predictors=`/`group=` with
# `group_effects="prevalence"`: entering the group as a covariate on class
# membership via the existing `covariate` structural model (R/covariate.R)
# gives each group its own free prevalences while the measurement model stays
# pooled/invariant (fit_rmlca()'s own docs note the same equivalence for the
# longitudinal case, Collins & Lanza sec. 6.10.2).
#
# The measurement question needs a model in which item-response probabilities
# can differ by group. This reuses exactly the trick `time_blocks_model()`
# already uses across occasions (R/time_blocks.R): stack G copies of the J
# items, one block per group, and hold some or all items equal across blocks
# via `invariant_items`. The one thing that differs from a time block is that
# a case belongs to exactly one group, so its OTHER (G-1) blocks are
# structurally missing rather than merely occasion-unobserved. FIML already
# treats a missing item as contributing nothing to that item's likelihood
# (R/categorical.R, R/gaussian.R), and does so correctly even when an entire
# block is missing for a case, so no change to the EM core is needed: a
# group-blocks case's own block estimates its own group's parameters, exactly
# as intended, and the M-step of `.blocks_model()` (shared with time blocks)
# already restricts each block's update to the rows with real data in it.

# Pad an n x J indicator matrix into an n x (J*G) matrix, one J-column block
# per group, with a case's own group's block populated and every other
# block NA. `group_factor` must already be a factor with one entry per row.
.pad_group_blocks <- function(X, group_factor) {
  n <- nrow(X)
  J <- ncol(X)
  G <- nlevels(group_factor)
  item_names <- colnames(X) %||% paste0("Item", seq_len(J))

  out <- matrix(NA_real_, nrow = n, ncol = J * G)
  g_idx <- as.integer(group_factor)
  Xnum  <- if (is.numeric(X)) X else data.matrix(X)
  for (g in seq_len(G)) {
    rows <- which(g_idx == g)
    if (!length(rows)) next
    out[rows, .time_block_cols(g, J)] <- Xnum[rows, , drop = FALSE]
  }
  colnames(out) <- paste0("G", rep(seq_len(G), each = J), ".", rep(item_names, G))
  out
}

# Constructor. Thin wrapper around the shared `.blocks_model()`; see
# R/time_blocks.R. `invariant_items` are the item indices held equal across
# groups (Collins & Lanza's partial-invariance models, sec. 5.9); leave empty
# for a fully configural (non-invariant) measurement model, i.e. the
# multiple-group equivalent of `measurement_invariance = "none"` in
# `fit_rmlca()`.
group_blocks_model <- function(n_components, n_items, n_groups,
                               sub_model = "bernoulli",
                               invariant_items = integer(0),
                               invariant_params = character(0),
                               variances_equal = FALSE,
                               max_val = NULL, ...) {
  .blocks_model(n_components, n_items, n_blocks = n_groups,
               sub_model = sub_model,
               invariant_items = invariant_items,
               invariant_params = invariant_params,
               variances_equal = variances_equal,
               max_val = max_val, prefix = "G",
               extra_class = "group_blocks")
}

# Restate the fit's missing-data summary in terms of the data the user supplied.
#
# .pad_group_blocks() turns J items into J*G columns of which a case can only
# ever observe its own group's J -- so (G-1)/G of the padded matrix is empty by
# construction, and fit_mixture_internal() builds `missing_data` from that. A
# three-group configural fit therefore reported "67.4% missing, handled via FIML
# (MAR assumption)", which is alarming, wrong, and meaningless: the padding is
# not missing data and MAR has nothing to say about it. Recomputed here from the
# unpadded matrix, with the block structure reported separately as the design
# fact it is.
.fix_group_block_reporting <- function(fit, X_unpadded, group_info) {
  if (!inherits(fit$mm, "blocks")) return(fit)

  keep <- setdiff(seq_len(nrow(X_unpadded)),
                  fit$missing_data$empty_rows %||% integer(0))
  Xu   <- X_unpadded[keep, , drop = FALSE]

  item_missing <- colSums(is.na(Xu))
  md <- fit$missing_data
  md$any_missing      <- anyNA(Xu) || isTRUE(md$y_any_missing)
  md$n_missing        <- sum(is.na(Xu))
  md$n_cells          <- length(Xu)
  md$prop_missing     <- if (length(Xu) > 0) mean(is.na(Xu)) else 0
  md$per_item         <- item_missing
  md$n_items_affected <- sum(item_missing > 0)
  md$handled_by       <- if (anyNA(Xu)) "FIML (MAR assumption)" else NA_character_
  # Named so print() and measurement_summary() can say what the blocks are
  # instead of describing them as absent data.
  md$group_blocks     <- list(n_groups = nlevels(group_info$factor),
                              n_items  = ncol(Xu),
                              levels   = levels(group_info$factor))
  fit$missing_data <- md
  fit
}

# Report the fit on the scale a known-class program uses, alongside its own.
#
# mixtureEM maximises the likelihood of the items *given* the group, and it does
# not count the group's own G-1 proportions as parameters. Software that instead
# treats the grouping variable as a latent class variable observed without error
# adds the group's multinomial term to the
# log-likelihood and its proportions to the parameter count. Both differences are
# fixed constants that cancel in any likelihood-ratio test and in any comparison
# between models fitted here, so they change no inference; they are supplied so a
# number read off such a program's output can be compared directly.
.add_knownclass_scale <- function(fit, group_info) {
  g <- group_info$factor
  n <- table(g)
  n <- n[n > 0]
  if (!length(n) || is.null(fit$metrics$ll)) return(fit)
  fit$metrics$ll_knownclass       <- fit$metrics$ll + sum(n * log(n / sum(n)))
  fit$metrics$n_params_knownclass <- fit$metrics$n_params + (length(n) - 1L)
  fit
}

# ------------------------------------------------------------------------------
# Warm-starting the configural (group-varying measurement) search
# ------------------------------------------------------------------------------
#
# A group-blocks model holds one measurement model per group, tied to the others
# only through the class labels, which are shared. Random starting values give
# each block its own arbitrary labelling, so class 1 in group 1 and class 1 in
# group 2 begin as unrelated things and EM has no mechanism to bring them into
# correspondence: it climbs from wherever it starts, and the aligned basin is a
# small target among G! * ... labellings. More restarts do not help, because
# every restart draws from the same badly-aligned prior. This is why the
# configural fit could land below its own nested restriction, which is
# impossible at the true optimum and makes the invariance LRT a lower bound
# rather than a test.
#
# The difficulty is structural, and Clogg and Goodman (1985, p. 85) name it:
# with no across-group restrictions the model — their *heterogeneous,
# unrestricted* T-class model — "is actually equivalent to applying the
# unrestricted T-class model separately to each of the S groups." Nothing ties
# the groups together except the labels, and the labels are exactly what random
# starts assign arbitrarily. (Restrict any parameter across groups, as
# `invariant_items` does, and "the model cannot be broken up into S independent
# single-group models"; the difficulty softens in proportion.)
#
# The fix is to build the starting values instead of drawing them, in two steps.
#
#   1. Fit the pooled (group-invariant) measurement model and replicate it into
#      every block. That state is aligned by construction — the blocks are
#      identical, so a class means the same thing in each — and it sits at the
#      restricted model's likelihood, so EM's monotonicity makes the configural
#      fit score at least as well as the model it is compared against and the
#      invariance LRT non-negative whatever else happens.
#   2. Fit each group on its own data and put *that* in the group's block, with
#      its classes permuted to match the pooled solution (align_classes(),
#      R/bootstrap.R). This is the stronger start of the two, because with the
#      measurement model and the prevalences both free by group the configural
#      model is very nearly separable: each block's own group is the only data
#      that informs it, so a per-group optimum is close to a fixed point of the
#      joint EM, and the pooled solution is only its starting point for the
#      alignment. Measured on the multi-country validation model (K = 4, three
#      groups, 177 parameters, 50 restarts): random starts alone reach
#      -5213.65, the replicated pooled start -5206.6, and the aligned per-group
#      start -5198.4.
#
# Step 2 needs step 1 anyway — without a common reference there is nothing to
# align the groups *to* — and it degrades to step 1 if a group's own fit fails,
# so both are built together and the second is optional.
#
# **Why this is one restart among many and never the only one.** Shireman,
# Steinley and Brusco (2017) compared initialisation strategies across ~33,000
# fitted mixtures and found the data-informed deterministic ones markedly worse
# at finding the best solution than plain random starts (proportion of data sets
# on which the strategy found the best solution: random .40, k-means .11,
# hierarchical clustering .05, sum scores .02), precisely because an informed
# start "constrains the search space to an even narrower region". Their informed
# starts *replace* the random pool; this one is added to it and competes on
# log-likelihood like any other restart, so a warm start that lands in a bad
# basin costs one restart and changes nothing. That distinction is the whole
# safety argument, and it is why `warm_start` must never become a way to skip
# the random restarts.
#
# Two things this is NOT. It is not multiple-group factor analysis alignment
# (Asparouhov & Muthén, 2014): that resolves a *continuous* indeterminacy — the
# unidentified factor mean and variance within each group of a configural CFA —
# by optimising a simplicity function. The indeterminacy here is *discrete*, a
# permutation of class labels, and the relevant literature is label switching
# (Grün & Leisch, 2009). Nor is it a substantive method: it produces a starting
# value, never an estimate.

# Copy every parameter matrix of `src` into `dst`, recursing through sub-model
# lists. Returns NULL — meaning "these two emissions are not the same shape, do
# not warm-start from this" — rather than erroring, since a caller that cannot
# build a warm start should simply not use one.
.copy_emission_parameters <- function(dst, src) {
  if (is.null(dst) || is.null(src)) return(NULL)
  if (!setequal(names(src$parameters), names(dst$parameters))) return(NULL)

  for (nm in names(src$parameters)) {
    p <- src$parameters[[nm]]
    q <- dst$parameters[[nm]]
    if (!identical(dim(as.matrix(p)), dim(as.matrix(q)))) return(NULL)
    dst$parameters[[nm]] <- p
  }

  # `nested` and mixed measurement models carry their own sub-models; a pooled
  # fit and one block of a group-blocks fit have the same ones in the same
  # order, differing only in the *_nan resolution, which does not reach here.
  if (!is.null(src$models) || !is.null(dst$models)) {
    if (length(src$models) != length(dst$models)) return(NULL)
    for (i in seq_along(src$models)) {
      sub <- .copy_emission_parameters(dst$models[[i]], src$models[[i]])
      if (is.null(sub)) return(NULL)
      dst$models[[i]] <- sub
    }
  }
  dst
}

# Permute the classes of an emission, recursing through sub-model lists. `perm`
# is what align_classes() returns: perm[i] is the class of this emission that
# plays the role of class i in the reference.
.permute_emission_classes <- function(mm, perm) {
  for (nm in names(mm$parameters)) {
    p <- mm$parameters[[nm]]
    if (is.matrix(p) && nrow(p) == length(perm))
      mm$parameters[[nm]] <- p[perm, , drop = FALSE]
  }
  if (!is.null(mm$models))
    mm$models <- lapply(mm$models, .permute_emission_classes, perm = perm)
  mm
}

# Align a group's own fitted measurement model to the pooled one, so that class
# k means the same thing in every block before the joint EM starts. Returns NULL
# if the two cannot be compared, which the caller reads as "fall back to the
# pooled parameters for this block".
.align_to_pooled <- function(mm, pooled_mm) {
  aligned <- try({
    perm <- align_classes(get_mm_alignment_matrix(pooled_mm),
                          get_mm_alignment_matrix(mm))
    .permute_emission_classes(mm, perm)
  }, silent = TRUE)
  if (inherits(aligned, "try-error")) NULL else aligned
}

# Build the `warm_start` function fit_em() expects.
#
# `pooled` is the fitted group-invariant model; `group_mms` is an optional list
# with one fitted measurement model per group, in group-level order, each of
# which may be NULL. Whatever is missing falls back to the pooled parameters,
# so a failure anywhere costs sharpness rather than correctness.
.group_blocks_warm_start <- function(pooled, group_mms = NULL) {
  force(pooled); force(group_mms)
  function(model_state, X, Y) {
    mm <- model_state$mm
    if (!inherits(mm, "blocks") || is.null(pooled$mm)) return(NULL)

    # A freshly built emission carries no parameters at all, so it is
    # initialised here for its shape; every value is then overwritten below.
    mm <- init_params(mm, X, NULL)

    for (b in seq_len(mm$n_blocks)) {
      src <- (if (is.null(group_mms)) NULL else group_mms[[b]]) %||% pooled$mm
      blk <- .copy_emission_parameters(mm$models[[b]], src)
      # A per-group model that does not match the block's shape is not a reason
      # to abandon the warm start — the pooled one always matches, since the
      # blocks were built from the same measurement descriptor.
      if (is.null(blk)) blk <- .copy_emission_parameters(mm$models[[b]], pooled$mm)
      if (is.null(blk)) return(NULL)
      mm$models[[b]] <- blk
    }

    # An `invariant_items` constraint is satisfied by construction when every
    # block holds the pooled parameters, but not when the blocks hold their own
    # groups' fits, so the constrained items are put back on the constraint
    # surface here — taking them from the pooled model, which is where the
    # group-invariant estimate of exactly those items already lives. This is
    # what init_params.blocks() does for a random start, with a better source
    # than its block 1.
    if (length(mm$invariant_items) && !is.null(group_mms))
      for (b in seq_len(mm$n_blocks))
        mm$models[[b]] <- .copy_item_params(mm$models[[b]], pooled$mm,
                                            mm$invariant_items)

    # Same argument for a parameter-wise constraint, and it binds whether or not
    # the per-group fits were used: the pooled fit's own parameters are the
    # group-invariant estimate of exactly the matrices being shared.
    mm <- .project_invariant_params(mm, src = pooled$mm)

    list(weights = pooled$weights, mm = mm)
  }
}

