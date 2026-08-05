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
                               max_val = NULL, ...) {
  .blocks_model(n_components, n_items, n_blocks = n_groups,
               sub_model = sub_model,
               invariant_items = invariant_items,
               max_val = max_val, prefix = "G",
               extra_class = "group_blocks")
}
