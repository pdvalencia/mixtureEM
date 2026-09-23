# ------------------------------------------------------------------------------
# The bias-adjusted three-step estimator for latent transition models:
# steps one and two.
# ------------------------------------------------------------------------------
#
# Step 1 estimates the measurement model with the items equated across
# occasions and NO transition structure, each occasion's status prevalences
# free. Step 2 classifies at each occasion and forms that occasion's
# classification-error matrix. Step 3 -- estimating the initial status and the
# transitions with those errors held fixed -- is a separate piece of work and
# is not here yet.
#
# WHY STEP 1 DROPS THE TRANSITIONS, which is the one decision in this file that
# is not obvious. Step 3's likelihood treats the assigned status W_t as a lone
# indicator of the true status S_t whose only defect is described by D_t. The
# posteriors of a model WITH transitions are smoothed along the chain: occasion
# 2's assignment borrows from occasions 1 and 3, so W_2 would carry information
# about S_1 and S_3 other than through S_2, and the conditional independence
# step 3 is built on fails. Dropping the transitions is what makes each
# occasion's posterior a function of that occasion's items and that occasion's
# base rates alone.
#
# WHY EACH OCCASION KEEPS ITS OWN ERROR MATRIX rather than one pooled over the
# occasions. The classification error is a function of the prior status
# distribution -- the base rate is what decides where an ambiguous response
# pattern is assigned -- and in a transition model that distribution moves by
# construction, because the movement is the estimand. On the reading panel one
# status holds 1.8% of the sample at occasion 1 and 81.2% at occasion 4, and
# P(assigned = that status | true = that status) is .598 at occasion 1 against
# .979 pooled: a pooled matrix believes occasion 1 is classified almost
# perfectly when two of every five true members are assigned elsewhere. One
# matrix per occasion recovers occasion 1's prevalence to 2e-4 where a pooled
# one lands at 39% of it. The pooled variant remains reachable as a constraint
# -- it is this model with the occasions' prevalences forced equal -- rather
# than as an argument.

# Steps one and two.
#
# `cl` is the caller's `match.call()` and `env` its frame. Step 1 is that same
# call with the structural predictors dropped and the transition structure
# removed, evaluated recursively, exactly as `n_steps = 2` fits its own step 1
# (R/lta.R): one call rather than a second copy of the argument handling, so
# the two models cannot drift apart, and so step 1 sees the identical data,
# weights, strata, groups, restart budget and seed.
#
# Returns the step-1 fit, the per-occasion prevalences and assignments, and the
# per-occasion classification-error matrices with their fixed logits.
.lta_threestep_step12 <- function(cl, env,
                                  assignment = c("modal", "proportional"),
                                  zero_floor = 1e-6) {
  assignment <- match.arg(assignment)

  # Two ways of handing over a starting point; silently preferring one is how a
  # user gets a fit they did not ask for. The same refusal `n_steps = 2` makes.
  if (!is.null(cl$refine_from))
    stop("`n_steps = 3` fits its own step 1, so `refine_from` has nothing to ",
         "hand over. Drop one of them.", call. = FALSE)

  # `tie_initial_status` ties the occasion-1 distribution across the classes
  # of a mover-stayer fit, and the three-step fits one class, so it would tie
  # nothing. Refused here, in the user's terms, rather than by step 1's
  # transition-free model in terms of a model the user never asked for.
  if (isTRUE(eval(cl$tie_initial_status, env)))
    stop("`tie_initial_status` applies only to a fit with several latent ",
         "classes (such as `mover_stayer = TRUE`), and `n_steps = 3` fits one. ",
         "Drop it.", call. = FALSE)

  # `indicators` is pinned to the value already in hand rather than left as the
  # expression the caller wrote: an expression that draws or simulates would
  # otherwise hand step 1 a different sample from the one step 3 will use,
  # silently. Every other argument stays as its expression, which is what makes
  # step 1 see the same weights, strata and groups the caller passed.
  cl1 <- cl
  ind <- eval(cl$indicators, env)
  e   <- new.env(parent = env)
  assign(".lta_step1_indicators", ind, envir = e)
  cl1$indicators <- quote(.lta_step1_indicators)

  cl1$n_steps                     <- 1L
  cl1$predictors_initial          <- NULL
  cl1$predictors_transition       <- NULL
  cl1$predictors_random_intercept <- NULL
  cl1$refine_from                 <- NULL
  cl1$transition_invariance       <- NULL
  cl1$transition_effects          <- NULL
  cl1$forbidden_transitions       <- NULL
  cl1$assignment                  <- NULL
  cl1$correction                  <- NULL
  cl1$zero_floor                  <- NULL
  cl1$.transition_free            <- TRUE

  step1 <- eval(cl1, e)

  K  <- step1$n_statuses
  Tn <- step1$n_times
  w  <- step1$weights_vec

  # Under the transition-free model every row of an occasion's matrix is that
  # occasion's prevalence vector, so occasion 1's prevalences are `delta` and
  # the rest are any row of the matrix that ends at that occasion.
  prevalences <- matrix(NA_real_, Tn, K)
  prevalences[1L, ] <- step1$delta
  for (t in seq_len(Tn - 1L)) prevalences[t + 1L, ] <- step1$tau[[t]][1L, ]
  dimnames(prevalences) <- list(step1$longitudinal$time_labels,
                                paste0("Status ", seq_len(K)))

  # `$gamma` is on the full sample -- .lta_expand() puts the fit back on it
  # before anything per-case is read -- so the posteriors and the case weights
  # are row-aligned here without any further expansion.
  errors <- lapply(seq_len(Tn), function(t)
    .classification_error(step1$gamma[[t]], assignment = assignment,
                          weights = w, zero_floor = zero_floor))

  # The hard classification, reported whichever rule was asked for: under
  # `modal` it is the response variable step 3 fits, and under `proportional`
  # step 3 fits the posteriors themselves, so it is a description rather than
  # an input. `ties.method` is left to get_modal_resp(), which is the rule the
  # error matrices were formed under -- a second tie-break here could disagree
  # with the table it is meant to describe.
  modal <- vapply(seq_len(Tn), function(t)
    max.col(get_modal_resp(step1$gamma[[t]]), ties.method = "first"),
    integer(nrow(step1$gamma[[1L]])))
  dimnames(modal) <- list(NULL, step1$longitudinal$time_labels)

  list(step1       = step1,
       assignment  = assignment,
       zero_floor  = zero_floor,
       prevalences = prevalences,
       modal       = modal,
       D           = lapply(errors, `[[`, "D"),
       logits      = lapply(errors, `[[`, "logits"))
}

# ------------------------------------------------------------------------------
# Step three.
# ------------------------------------------------------------------------------
#
# THE REDUCTION, which is the whole reason this step is small. After step 2 the
# items have done their work and are never looked at again. What is left is one
# assigned-status variable W_t per occasion, taking K values, and a K x K table
# D_t saying how often an assigned value misses the true one. A model for that
# is an ordinary latent transition model with ONE categorical indicator per
# occasion whose response probabilities are already known -- so it is fitted by
# holding the measurement block fixed and letting the initial-status and
# transition blocks run free, which is the machinery the two-step estimator
# already put in place.
#
# THE ORIENTATION, which is the one mistake here that is silent and fatal.
# .classification_error() returns D with rows the ASSIGNED status and columns
# the TRUE one, columns summing to 1, because that is the orientation the fixed
# logits are read off. An emission is the other way round: its parameter matrix
# is rows the latent status, columns the response, rows summing to 1. So the
# emission is t(D_t). Written here, at the single place the two conventions
# meet, and checked by a test rather than trusted to this comment -- planted
# backwards it produces a converged, plausible fit at the wrong answer.

# The reduced data. Returns the response matrix and the weight attached to each
# of its rows; the caller passes them straight to fit_lta().
#
# Under `modal` each case contributes one row -- the statuses it was assigned --
# with its own sampling weight. Under `proportional` a case is not assigned to
# one status but spread over all of them in proportion to its posterior, so it
# contributes to every combination (w_1, ..., w_T) with weight
# `prod_t gamma_t[i, w_t]`. That is the longitudinal form of the expanded data
# set fit_ml() builds cross-sectionally, where a case becomes K records weighted
# by its posterior (R/corrections.R).
#
# The expansion does not have to be carried out case by case. The reduced model
# has exactly one indicator per occasion, so the whole expanded data set has at
# most K^T DISTINCT rows however large the sample is -- 81 on a three-status,
# four-occasion panel -- and the weights are what is summed over cases. The cap
# below is on that grid, which is the only thing here that grows with the model.
#
# WHAT A COVARIATE CHANGES, and it is why `Z` is here. Summing the weights over
# all the cases is legitimate only while the cases are exchangeable in step 3's
# likelihood, and a covariate is exactly what stops them being so: two cases
# with different covariate values have different initial-status and transition
# probabilities, so they cannot share a row. They can still share one if their
# covariate values agree. So the grid is built WITHIN each distinct covariate
# pattern rather than over the sample, and that gives the same likelihood the
# case-by-case expansion would -- exactly, not approximately, because the
# covariates enter step 3 only through the pattern. With `Z` NULL there is one
# pattern and the paragraph above is the special case.
#
# The case-by-case expansion was rejected on cost, not on principle: it is
# n * K^T rows, 289,575 on the reading panel, against 162 for one binary
# covariate.
.lta_threestep_reduce <- function(s12, Z = NULL) {
  K  <- s12$step1$n_statuses
  Tn <- s12$step1$n_times
  w  <- s12$step1$weights_vec

  if (identical(s12$assignment, "modal"))
    return(list(W = s12$modal, weights = w, Z = Z))

  n_cells <- K^Tn
  if (n_cells > 1e4)
    stop(sprintf(paste0(
      "Proportional assignment spreads each case over all %d combinations of ",
      "%d statuses at %d occasions, which is too many to enumerate. Use ",
      "assignment = \"modal\", which assigns each case to one."),
      as.integer(n_cells), K, Tn), call. = FALSE)

  grid <- as.matrix(expand.grid(rep(list(seq_len(K)), Tn)))
  dimnames(grid) <- list(NULL, s12$step1$longitudinal$time_labels)

  # Each case's weight on each cell of the grid, n x K^T: the case's own weight
  # times the product over occasions of its posterior on that cell's status.
  # One multiplication per occasion rather than one pass per cell -- the same
  # arithmetic, and it is what lets the pattern sum below be a rowsum().
  cell_w <- matrix(w, length(w), nrow(grid))
  for (t in seq_len(Tn))
    cell_w <- cell_w * s12$step1$gamma[[t]][, grid[, t], drop = FALSE]

  if (is.null(Z)) {
    idx <- rep(1L, length(w))
    rep_row <- 1L
  } else {
    key     <- do.call(paste, c(as.data.frame(Z), sep = "\r"))
    uniq    <- unique(key)
    idx     <- match(key, uniq)
    rep_row <- match(seq_along(uniq), idx)
  }

  # rowsum() sorts its groups by label and the labels are 1..P, so the rows of
  # `agg` line up with `rep_row`. as.vector() then runs down the patterns within
  # a cell, which is the order `cell_of` and `pat_of` decode.
  agg     <- rowsum(cell_w, idx)
  P       <- nrow(agg)
  wts     <- as.vector(agg)
  cell_of <- rep(seq_len(nrow(grid)), each = P)
  pat_of  <- rep(seq_len(P), times = nrow(grid))

  # A combination no case gives any weight to contributes nothing to the
  # likelihood and is dropped rather than carried as a row of zeros. The
  # weights still sum to the number of cases, since every case's posteriors
  # multiply out to one across the grid.
  keep <- wts > 0
  list(W       = grid[cell_of[keep], , drop = FALSE],
       weights = wts[keep],
       Z       = if (is.null(Z)) NULL else
                   Z[rep_row[pat_of[keep]], , drop = FALSE])
}

# Step three: the structural model, fitted on the reduced data with the
# classification error held fixed.
#
# `cl` and `env` are the caller's `match.call()` and frame, as in
# .lta_threestep_step12(). The reduced fit is a fit_lta() call like any other,
# so the settings that describe a search rather than a model -- the restart
# budget, the seed, the cores, the stopping rule, the prior on the structural
# probabilities -- are carried across from the caller and mean the same thing
# they meant in step 1. What is NOT carried is everything describing the items:
# there are no items left. The structural predictors ARE carried, because
# regressing the initial status and the transitions on covariates is what step 3
# exists to do.
#
# THE PREDICTORS ARE PREPARED HERE rather than passed on as the caller's
# expressions, which is the one thing about this function that is not obvious.
# Two reasons, and the second is the binding one.
#
#   The reduced data under proportional assignment has one row per (covariate
#   pattern, status combination), not one row per case, so the caller's
#   covariate vector is the wrong length for it. The expansion has to know the
#   covariate values to group by them, and the fit has to receive the expanded
#   ones. Both need the design matrix in hand, here.
#
#   Preparing it once, on the case scale, is also the only way the two agree.
#   .lta_design() imputes a missing covariate value at the column mean; run on
#   the deduplicated rows instead it would average over patterns rather than
#   over cases and quietly shift the design. So it is run once, before the
#   expansion, and its output -- a plain numeric matrix, which
#   prepare_covariates() returns untouched -- is what fit_lta() receives, where
#   passing it through .lta_design() a second time is the identity.
.lta_threestep_step3 <- function(s12, cl, env) {
  n <- length(s12$step1$weights_vec)

  # The intercept column .lta_design() prepends is dropped: fit_lta() adds its
  # own back. Position, not name, because a user's covariate may be called
  # anything at all.
  prepared <- function(nm) {
    if (is.null(cl[[nm]])) return(NULL)
    .lta_design(eval(cl[[nm]], env), n, nm)[, -1L, drop = FALSE]
  }
  Z_delta <- prepared("predictors_initial")
  Z_tau   <- prepared("predictors_transition")

  # One key over both blocks: two cases are the same pattern only if they agree
  # on everything step 3's structural model reads.
  red <- .lta_threestep_reduce(
    s12, if (is.null(Z_delta) && is.null(Z_tau)) NULL else cbind(Z_delta, Z_tau))

  nd <- if (is.null(Z_delta)) 0L else ncol(Z_delta)
  nt <- if (is.null(Z_tau))   0L else ncol(Z_tau)

  e <- new.env(parent = env)
  assign(".lta_step3_W", red$W,                    envir = e)
  assign(".lta_step3_w", red$weights,              envir = e)
  assign(".lta_step3_D", lapply(s12$D, t),         envir = e)
  assign(".lta_step3_labels", s12$step1$longitudinal$time_labels, envir = e)
  if (nd) assign(".lta_step3_Zd", red$Z[, seq_len(nd), drop = FALSE], envir = e)
  if (nt) assign(".lta_step3_Zt", red$Z[, nd + seq_len(nt), drop = FALSE],
                 envir = e)

  carried <- c("n_init", "refine", "max_iter", "n_cores", "tol", "smoothing",
               "random_state", "standard_errors", "bayes_constants",
               "transition_invariance", "transition_effects",
               "forbidden_transitions")
  carried <- cl[intersect(carried, names(cl))]

  cl3 <- as.call(c(
    list(quote(fit_lta),
         indicators             = quote(.lta_step3_W),
         n_statuses             = s12$step1$n_statuses,
         times                  = s12$step1$n_times,
         measurement            = "categorical",
         measurement_invariance = "none",
         time_labels            = quote(.lta_step3_labels),
         weights                = quote(.lta_step3_w),
         weight_type            = "frequency",
         .fixed_emission        = quote(.lta_step3_D)),
    if (nd) list(predictors_initial    = quote(.lta_step3_Zd)),
    if (nt) list(predictors_transition = quote(.lta_step3_Zt)),
    as.list(carried)))

  eval(cl3, e)
}
