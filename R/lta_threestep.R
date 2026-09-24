# ------------------------------------------------------------------------------
# The bias-adjusted three-step estimator for latent transition models:
# steps one and two.
# ------------------------------------------------------------------------------
#
# Step 1 estimates the measurement model under the caller's invariance
# constraints and with NO transition structure, each occasion's status
# prevalences free. Step 2 classifies at each occasion and forms that
# occasion's classification-error matrix. Step 3, further down, estimates the
# initial status and the transitions with those errors held fixed.
# .lta_threestep() at the bottom of the file is the entry point fit_lta()
# calls for `n_steps = 3`.
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
#
# `correction = "none"` is the naive classify-analyse baseline: every error
# matrix is the identity, so step 3 believes the assigned statuses. It is what
# the correction exists to improve on, and having it available is what makes
# the improvement visible on the user's own data.
.lta_threestep_step12 <- function(cl, env,
                                  assignment = c("modal", "proportional"),
                                  zero_floor = 1e-6,
                                  correction = c("ML", "none")) {
  assignment <- match.arg(assignment)
  correction <- match.arg(correction)

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
  cl1$distal                      <- NULL
  cl1$zero_floor                  <- NULL
  cl1$.transition_free            <- TRUE

  # With `predictors_items`, step 1 also regresses the status at every occasion
  # on the same covariates, and step 3 fits whatever structural predictors the
  # caller asked for. Without the status regression the item slope has to
  # carry the whole association between the covariate and the status, so a
  # covariate that only moves people between statuses would read as bias in
  # the items. The regressions go in through the transition-free model's
  # origin-free design (.lta_tau_design()), which keeps the occasions
  # independent given the covariate: one prevalence regression per occasion.
  pi_items <- eval(cl$predictors_items, env)
  if (!is.null(pi_items)) {
    Zi <- if (is.list(pi_items) && !is.data.frame(pi_items))
      .lta_dif_spec(pi_items, NULL, names(pi_items), NROW(pi_items[[1L]]))$Z
    else pi_items
    assign(".lta_step1_status_covariates", Zi, envir = e)
    cl1$predictors_initial    <- quote(.lta_step1_status_covariates)
    cl1$predictors_transition <- quote(.lta_step1_status_covariates)
  }

  step1 <- eval(cl1, e)

  # With the item parameters free at every occasion and no transitions, nothing
  # in step 1's likelihood says which status at occasion 2 is "the same" as a
  # status at occasion 1: every per-occasion relabelling fits identically. So
  # the labels are matched here, by the rule a reader would use -- status k at
  # occasion t is the one whose item profile is closest to status k's at
  # occasion 1. `"partial"` needs none of this: the items held equal pin the
  # labels inside the likelihood itself.
  mi <- if (is.null(cl$measurement_invariance)) "full" else
    match.arg(eval(cl$measurement_invariance, env), c("full", "none", "partial"))
  alignment <- NULL
  if (identical(mi, "none")) {
    step1     <- .lta_threestep_align(step1)
    alignment <- step1$alignment
  }

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
  if (correction == "none")
    errors <- lapply(errors, function(x) {
      x$D[] <- diag(K)
      x$logits <- NULL
      x
    })

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
       correction  = correction,
       alignment   = alignment,
       zero_floor  = zero_floor,
       prevalences = prevalences,
       modal       = modal,
       D           = lapply(errors, `[[`, "D"),
       logits      = lapply(errors, `[[`, "logits"))
}

# Match the status labels of a transition-free step 1 whose item parameters are
# free at every occasion. Each occasion t >= 2 is permuted onto occasion 1 by
# the smallest total squared distance between item profiles, which is
# align_classes() over the same profile matrix the bootstrap aligns on. The
# permutation reaches every per-occasion quantity the three-step reads: the
# occasion's measurement block, its posteriors, and the transition-free rows
# that carry its prevalences. Standard errors and pairwise posteriors are
# dropped if anything moved rather than permuted piecemeal; nothing downstream
# reads them. `$alignment` records the permutations, identity included.
.lta_threestep_align <- function(step1) {
  Tn   <- step1$n_times
  ref  <- get_mm_alignment_matrix(step1$mm$models[[1L]])
  perm <- c(list(seq_len(step1$n_statuses)), lapply(seq_len(Tn)[-1L], function(t)
    align_classes(ref, get_mm_alignment_matrix(step1$mm$models[[t]]))))

  if (any(vapply(perm, is.unsorted, logical(1)))) {
    for (t in seq_len(Tn)) {
      p <- perm[[t]]
      step1$mm$models[[t]] <- .permute_emission_classes(step1$mm$models[[t]], p)
      step1$gamma[[t]]     <- step1$gamma[[t]][, p, drop = FALSE]
    }
    for (t in seq_len(Tn - 1L))
      step1$tau[[t]] <- step1$tau[[t]][perm[[t]], perm[[t + 1L]], drop = FALSE]
    # One class above the chain: `tau_c[[1]]` is the same list `tau` reports.
    step1$tau_c[[1L]] <- step1$tau
    step1$prevalences <- .lta_prevalences(step1)
    step1$se <- NULL
    step1$xi <- NULL
  }
  names(perm) <- step1$longitudinal$time_labels
  step1$alignment <- perm
  step1
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
#
# The preparation is its own function because the variance below needs the same
# case-scale design to rebuild the reduced data.
.lta_threestep_Z <- function(s12, cl, env) {
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
  list(delta = Z_delta, tau = Z_tau,
       Z = if (is.null(Z_delta) && is.null(Z_tau)) NULL else cbind(Z_delta, Z_tau))
}

.lta_threestep_step3 <- function(s12, cl, env) {
  Zc      <- .lta_threestep_Z(s12, cl, env)
  Z_delta <- Zc$delta
  Z_tau   <- Zc$tau

  red <- .lta_threestep_reduce(s12, Zc$Z)

  nd <- if (is.null(Z_delta)) 0L else ncol(Z_delta)
  nt <- if (is.null(Z_tau))   0L else ncol(Z_tau)
  extra <- .lta_threestep_pattern_start(s12, red, nd, nt, cl, env)

  e <- new.env(parent = env)
  assign(".lta_step3_W", red$W,                    envir = e)
  assign(".lta_step3_w", red$weights,              envir = e)
  assign(".lta_step3_D", lapply(s12$D, t),         envir = e)
  assign(".lta_step3_labels", s12$step1$longitudinal$time_labels, envir = e)

  # The distal outcome, on the user's rows, is cut to step 1's: step 1 drops a
  # case with no observed indicator at any occasion, and step 3's rows are
  # step 1's (modal assignment, the only one `distal` is allowed with).
  distal <- if (is.null(cl$distal)) NULL else eval(cl$distal, env)
  if (!is.null(distal)) {
    md <- s12$step1$missing_data
    distal <- as.data.frame(distal)
    if (!is.null(md) && length(md$empty_rows) && nrow(distal) == md$n_input_rows)
      distal <- distal[-md$empty_rows, , drop = FALSE]
    assign(".lta_step3_Y", .lta_distal_prepare(distal, nrow(red$W)), envir = e)
  }
  if (nd) assign(".lta_step3_Zd", red$Z[, seq_len(nd), drop = FALSE], envir = e)
  if (nt) assign(".lta_step3_Zt", red$Z[, nd + seq_len(nt), drop = FALSE],
                 envir = e)
  if (!is.null(extra)) assign(".lta_step3_extra", list(extra), envir = e)

  # A survey design reaches step 3 only under modal assignment (fit_lta()
  # refuses the other), where the reduced data is one row per case in step 1's
  # order, so step 1's own per-case vectors are the ones step 3 needs.
  design <- isTRUE(s12$step1$has_survey_design)
  if (design) {
    if (length(s12$step1$strata) != nrow(red$W))
      stop("The survey design does not line up with step 3's rows.",
           call. = FALSE)
    assign(".lta_step3_strata",  s12$step1$strata,  envir = e)
    assign(".lta_step3_cluster", s12$step1$cluster, envir = e)
  }

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
    if (!is.null(extra)) list(.extra_starts = quote(.lta_step3_extra)),
    if (!is.null(distal)) list(.distal = quote(.lta_step3_Y)),
    if (design) list(strata  = quote(.lta_step3_strata),
                     cluster = quote(.lta_step3_cluster)),
    as.list(carried)))

  eval(cl3, e)
}

# The pattern start: one extra starting point for a covariate step 3, built
# from the unconditional model fitted within each covariate pattern.
#
# Step 3's emission is a known constant, so its likelihood is a sum over cases
# that couple only through the structural parameters, and a covariate reaches
# those only through the case's pattern. Within one pattern the step-3 model
# is therefore the unconditional one -- whose search replicates reliably --
# and a model saturated in the patterns (`transition_effects = "by_origin"` on
# one categorical covariate is) is exactly the patterns fitted separately.
# For any other design the separate fits are a start, not a solution.
#
# Why it is needed: under "by_origin" the structural surface is multimodal and
# the random restarts miss the highest mode, while the separate fits find it
# in every restart (RECORDS.md, "W4.5's own numbers").
#
# The probabilities are carried onto the fitted model's own coefficients in
# its own design -- .lta_tau_design() builds the rows -- weighted by the
# expected pattern-level counts, so every design -- "common", "by_origin", "full",
# "slopes" -- is converted by the machinery that defines it. NULL, and so no
# extra start, when there are no covariates, one pattern, more than
# `max_patterns` (a continuous covariate), or a pattern fit fails.
.lta_threestep_pattern_start <- function(s12, red, nd, nt, cl, env,
                                         max_patterns = 16L) {
  Z <- red$Z
  if (is.null(Z) || nd + nt == 0L) return(NULL)
  key <- do.call(paste, c(lapply(seq_len(ncol(Z)), function(j)
    match(Z[, j], unique(Z[, j]))), sep = "\r"))
  grp <- match(key, unique(key))
  P   <- max(grp)
  if (P < 2L || P > max_patterns) return(NULL)

  K  <- s12$step1$n_statuses
  Tn <- s12$step1$n_times
  arg <- function(nm, default) {
    v <- if (is.null(cl[[nm]])) NULL else eval(cl[[nm]], env)
    if (is.null(v)) default else v
  }
  inv       <- arg("transition_invariance", "none")[1L]
  effects   <- arg("transition_effects", "common")[1L]
  shared    <- inv %in% c("full", "slopes")

  search <- c("n_init", "max_iter", "tol", "smoothing", "random_state",
              "bayes_constants", "n_cores")
  search <- lapply(cl[intersect(search, names(cl))], eval, env)
  fits <- lapply(seq_len(P), function(g) {
    i <- which(grp == g)
    # A pattern may hold few cases and fail to converge or warn; this is a
    # start, ranked against the random pool, so neither is the user's concern.
    try(suppressWarnings(do.call(fit_lta, c(list(
      red$W[i, , drop = FALSE], n_statuses = K, times = Tn,
      measurement = "categorical", measurement_invariance = "none",
      weights = red$weights[i], weight_type = "frequency",
      transition_invariance = if (inv == "full") "full" else "none",
      standard_errors = FALSE, .fixed_emission = lapply(s12$D, t)),
      search))), silent = TRUE)
  })
  if (any(vapply(fits, inherits, TRUE, "try-error"))) return(NULL)

  # Coefficients by weighted least squares on the log-ratios against the last
  # category, which is the anchor .fit_mnl() pins: exact when the design is
  # saturated in the patterns, a projection otherwise. Not .fit_mnl() itself --
  # its BFGS stops after optim()'s default 100 iterations, which is enough
  # inside EM but leaves cells near a boundary visibly off, and a start that
  # is off by that much was measured to fall into the wrong mode.
  logit_ls <- function(Zm, Pm, w) {
    L  <- log(pmax(Pm, 1e-12))
    L  <- L - L[, ncol(Pm)]
    ok <- w > 0
    sw <- sqrt(w[ok])
    B  <- vapply(seq_len(ncol(Pm) - 1L), function(d)
      qr.coef(qr(Zm[ok, , drop = FALSE] * sw), L[ok, d] * sw), numeric(ncol(Zm)))
    B  <- matrix(B, ncol(Zm))
    B[is.na(B)] <- 0
    rbind(t(B), 0)
  }

  rep_row <- match(seq_len(P), grp)
  n_g     <- vapply(seq_len(P), function(g) sum(red$weights[grp == g]), 0)
  delta_g <- t(vapply(fits, function(f) f$delta_c[[1]], numeric(K)))
  # Expected occupancy of each origin at each occasion, within each pattern.
  occ_g <- lapply(fits, function(f) {
    o <- list(f$delta_c[[1]])
    for (t in seq_len(Tn - 1L)) o[[t + 1L]] <- as.vector(o[[t]] %*% f$tau[[t]])
    o
  })

  delta_beta <- if (nd) logit_ls(
    cbind(1, Z[rep_row, seq_len(nd), drop = FALSE]), delta_g,
    n_g)

  st <- list(Z_tau = cbind(1, Z[rep_row, nd + seq_len(nt), drop = FALSE]),
             n_statuses = K, n_times = Tn, transition_effects = effects,
             tau_occasion_free_intercepts = inv == "slopes")
  rows <- function(t, k) t(vapply(fits, function(f) f$tau[[t]][k, ], numeric(K)))
  wts  <- function(t, k) n_g * vapply(occ_g, function(o) o[[t]][k], 0)
  tau_beta <- if (nt) lapply(seq_len(if (shared) 1L else Tn - 1L), function(m) {
    ts <- if (shared) seq_len(Tn - 1L) else m
    if (identical(effects, "by_origin")) {
      lapply(seq_len(K), function(k) logit_ls(
        do.call(rbind, rep(list(st$Z_tau), length(ts))),
        do.call(rbind, lapply(ts, rows, k = k)),
        unlist(lapply(ts, wts, k = k))))
    } else {
      tk <- expand.grid(k = seq_len(K), t = ts)
      logit_ls(
        do.call(rbind, Map(function(t, k) .lta_tau_design(st, k, t),
                           tk$t, tk$k)),
        do.call(rbind, Map(rows, tk$t, tk$k)),
        unlist(Map(wts, tk$t, tk$k)))
    }
  })

  big <- fits[[which.max(n_g)]]
  list(delta_c = list(colSums(delta_g * n_g) / sum(n_g)),
       tau_c = big$tau_c, delta_beta = delta_beta, tau_beta = tau_beta,
       mm = big$mm, n_classes = 1L, class_weights = 1)
}

# ------------------------------------------------------------------------------
# Standard errors that know step 1 was estimated.
# ------------------------------------------------------------------------------
#
# Step 3 holds the classification error fixed, so the variance its own
# information gives -- V2 -- treats `D_t` as known constants. They are not: they
# were computed from step 1's estimates, which carry sampling error of their
# own. Step 3 is a pseudo-maximum-likelihood estimator in the sense of Gong and
# Samaniego (1981), and the first-order correction Bakk, Oberski and Vermunt
# (2014) derive for the three-step is
#
#   V = V2 + V2 H21' S11 H21 V2,
#
# with S11 the sampling variance of step 1's estimates and H21 the cross
# derivative of step 3's log-likelihood in its structural coefficients and step
# 1's parameters. That is the two-step's equation with step 1's parameters in
# place of the measurement block, so the cross derivative is taken by the same
# function (.lta_twostep_information(), R/twostep_variance.R) on an objective
# that carries step 1's parameters through the three things step 3 reads off
# them: each occasion's posterior, the error matrix built from it, and -- under
# proportional assignment only -- the weights of the reduced data. Modal
# assignment is held where step 2 put it: the derivative of a modal table is
# the slope of its cell shares, and a case crossing a boundary is a jump, not
# a slope.
#
# Step 1's own vector is read through a copy with the shared-row flag off. The
# transition-free model ties every origin row of an occasion to one vector;
# with the rows equal, the ordinary packing describes it exactly, and the
# reduction below maps one free vector per occasion onto the K rows that share
# it. So S11 is the variance of the model step 1 fitted, not of a model with K
# free rows per occasion.
#
# Returns NULL when there is nothing to propagate (modal assignment with
# `correction = "none"`: every table is the identity and every weight a case
# weight) or when either fit falls outside the packing.
.lta_threestep_vcov <- function(fit, s12, Z = NULL, parts = FALSE) {
  step1 <- s12$step1
  modal <- identical(s12$assignment, "modal")
  none  <- identical(s12$correction, "none")
  if (modal && none) return(NULL)
  if (!isTRUE(fit$mm_fixed) || !is.null(fit$ri) || !is.null(step1$ri) ||
      (fit$n_classes %||% 1L) > 1L)
    return(NULL)

  K  <- step1$n_statuses
  Tn <- step1$n_times
  X1 <- as.matrix(step1$data)
  w1 <- step1$weights_vec

  # Step 1, reduced to its free parameters.
  s1 <- step1
  s1$tau_independent <- FALSE
  lay1  <- .lta_par_layout(s1)
  full1 <- .lta_par_pack(s1, lay1)
  len   <- vapply(lay1, function(b) as.integer(b$len), integer(1))
  beg   <- cumsum(len) - len
  src   <- integer(length(full1))
  n_red <- 0L
  first <- list()
  for (i in seq_along(lay1)) {
    b    <- lay1[[i]]
    cols <- beg[i] + seq_len(len[i])
    key  <- if (identical(b$kind, "tau")) paste0("tau", b$i_mat) else NULL
    if (!is.null(key) && !is.null(first[[key]])) {
      src[cols] <- first[[key]]
    } else {
      src[cols] <- n_red + seq_len(len[i])
      n_red <- n_red + len[i]
      if (!is.null(key)) first[[key]] <- src[cols]
    }
  }
  r0 <- full1[match(seq_len(n_red), src)]
  if (max(abs(r0[src] - full1)) > 1e-8) return(NULL)

  # Every nudge below re-runs a forward-backward pass, so both steps run on
  # their distinct rows with the case weights summed, as the search does
  # (.lta_collapse()). Step 1 has no covariates, so its rows are response
  # patterns, and a posterior -- hence a modal assignment -- is a function of
  # the pattern, so the per-case posteriors are the pattern's, expanded back.
  pat1 <- .pattern_index(X1)
  X1c  <- X1[pat1$rep_row, , drop = FALSE]
  w1c  <- as.vector(rowsum(w1, pat1$idx))
  s1$weights_vec <- w1c

  ll1 <- function(r) sum(w1c * .lta_ll_case(s1, X1c, r[src], lay1))
  S11 <- .psd_pinv(-.step1_fd_hessian(ll1, r0))
  S11 <- (S11 + t(S11)) / 2

  # What step 3 reads off step 1, at an arbitrary point of step 1's vector.
  posteriors <- function(r) {
    st <- .lta_par_unpack(r[src], s1, lay1)
    g  <- .lta_forward_backward(.lta_emission_loglik(st$mm, X1c),
                                .lta_log_delta(st), .lta_log_tau(st), w1c)$gamma
    lapply(g, function(m) m[pat1$idx, , drop = FALSE])
  }
  emission <- function(gamma) lapply(seq_len(Tn), function(t) {
    if (none) return(diag(K))
    t(.classification_error(gamma[[t]], assignment = s12$assignment,
                            weights = w1, zero_floor = s12$zero_floor,
                            assigned = if (modal) s12$modal[, t])$D)
  })

  lay3 <- .lta_par_layout(fit)
  par3 <- .lta_par_pack(fit, lay3)
  W3   <- as.matrix(fit$data)
  w3   <- fit$weights_vec
  if (length(par3) != fit$n_params) return(NULL)

  # Under proportional assignment the reduced data's weights are step 1's
  # posteriors multiplied out, so they move with step 1 too. The rows are
  # rebuilt by the function that built them and must come back in the same
  # order, which is checked once at the fitted point.
  weights_at <- function(gamma) {
    if (modal) return(w3)
    s <- s12
    s$step1$gamma <- gamma
    .lta_threestep_reduce(s, Z)$weights
  }
  if (!modal) {
    red0 <- .lta_threestep_reduce(s12, Z)
    if (!identical(dim(red0$W), dim(W3)) || any(red0$W != W3) ||
        max(abs(red0$weights - w3)) > 1e-8 * max(w3))
      return(NULL)
  }

  # Step 1's point changes only across the cross terms, so each point's tables
  # and weights are computed once and reused over the whole row of them.
  cache <- new.env(parent = emptyenv())
  step1_side <- function(r) {
    key <- paste(format(r, digits = 17), collapse = ",")
    hit <- cache[[key]]
    if (!is.null(hit)) return(hit)
    g   <- posteriors(r)
    out <- list(emission = emission(g), weights = weights_at(g))
    if (length(ls(cache)) > 8L) rm(list = ls(cache), envir = cache)
    assign(key, out, envir = cache)
    out
  }

  # Step 3 likewise, on its distinct (assigned statuses, covariate) rows: at
  # most K^T per covariate pattern. The weights are summed within a row at
  # every point, since under proportional assignment they move with step 1.
  pat3 <- .pattern_index(W3, cbind(fit$Z_delta, fit$Z_tau))
  W3c  <- W3[pat3$rep_row, , drop = FALSE]
  st3c <- fit
  if (!is.null(fit$Z_delta))
    st3c$Z_delta <- fit$Z_delta[pat3$rep_row, , drop = FALSE]
  if (!is.null(fit$Z_tau))
    st3c$Z_tau <- fit$Z_tau[pat3$rep_row, , drop = FALSE]

  p3   <- length(par3)
  idx2 <- seq_len(p3)
  idx1 <- p3 + seq_len(n_red)
  ll <- function(v) {
    s   <- step1_side(v[idx1])
    wc  <- as.vector(rowsum(s$weights, pat3$idx))
    st3 <- st3c
    for (t in seq_len(Tn)) st3$mm$models[[t]]$parameters$pis <- s$emission[[t]]
    st3$weights_vec <- wc
    sum(wc * .lta_ll_case(st3, W3c, v[idx2], lay3))
  }

  info <- .lta_twostep_information(fit, W3, w3, lay3, c(par3, r0), idx1, idx2,
                                   ll = ll)
  V2 <- .psd_pinv(-info$H22)
  V2 <- (V2 + t(V2)) / 2
  V  <- V2 + V2 %*% t(info$H12) %*% S11 %*% info$H12 %*% V2
  V  <- (V + t(V)) / 2
  if (!all(is.finite(V))) return(NULL)
  method <- "Three-step pseudo-ML, first order (Bakk, Oberski and Vermunt, 2014)"

  # Under a survey design the two steps are one stacked set of estimating
  # equations, and the variance is their sandwich. The bread is block lower
  # triangular -- step 1 does not read step 3 -- so the step-3 row of its
  # inverse is [V2 H21' S11, V2], and each case's influence on step 3 is V2
  # times U_i = s3_i + s1_i S11 H21: its own step-3 score plus its step-1 score
  # carried through the same cross derivative as above. Summing U within PSU
  # and comparing across PSUs within stratum (compute_survey_B()) then counts
  # the clustering of both steps and of their covariance. With independent
  # cases its expectation is the V above. Scores are taken on the distinct
  # rows, as everything else here is, and expanded back to the cases.
  B11 <- B33 <- NULL
  if (isTRUE(fit$has_survey_design)) {
    if (!modal || nrow(W3) != length(w1) || length(fit$strata) != length(w1))
      return(NULL)
    row_grad <- function(f, v) {
      h <- .step1_fd_step * pmax(1, abs(v))
      vapply(seq_along(v), function(j) {
        vp <- v; vp[j] <- vp[j] + h[j]
        vm <- v; vm[j] <- vm[j] - h[j]
        (f(vp) - f(vm)) / (2 * h[j])
      }, numeric(length(f(v))))
    }
    G1 <- row_grad(function(r) .lta_ll_case(s1, X1c, r[src], lay1), r0)
    s3 <- step1_side(r0)
    st3 <- st3c
    for (t in seq_len(Tn)) st3$mm$models[[t]]$parameters$pis <- s3$emission[[t]]
    st3$weights_vec <- as.vector(rowsum(w3, pat3$idx))
    G3 <- row_grad(function(v) .lta_ll_case(st3, W3c, v, lay3), par3)
    Sc1 <- w1 * G1[pat1$idx, , drop = FALSE]
    Sc3 <- w3 * G3[pat3$idx, , drop = FALSE]

    U <- Sc3 + Sc1 %*% S11 %*% info$H12
    V <- V2 %*% compute_survey_B(U, fit$strata, fit$cluster) %*% V2
    V <- (V + t(V)) / 2
    if (!all(is.finite(V))) return(NULL)
    B11 <- S11 %*% compute_survey_B(Sc1, fit$strata, fit$cluster) %*% S11
    B33 <- V2 %*% compute_survey_B(Sc3, fit$strata, fit$cluster) %*% V2
    method <- paste(method, "with the survey-linearized sandwich over both steps")
  }

  out <- list(V = V, V2 = V2, S11 = S11, method = method)
  # The pieces a test needs to check the propagation term from outside: the
  # cross derivative, step 1's reduced point and what step 3 reads off it.
  if (parts) out <- c(out, list(H12 = info$H12, r0 = r0, par3 = par3,
                                lay3 = lay3, step1_side = step1_side,
                                design_step1 = B11, design_step3 = B33))
  out
}

# Put the corrected variance on the fit: the structural block of `se$vcov`,
# the probability-scale standard errors read off it, and V2 kept alongside.
.lta_threestep_attach_vcov <- function(fit, s12, Z = NULL) {
  if (is.null(fit$se) || is.null(fit$se$vcov)) return(fit)
  ts <- tryCatch(.lta_threestep_vcov(fit, s12, Z), error = function(e) {
    warning("The three-step standard errors could not be corrected for step ",
            "1's uncertainty and treat the classification error as known: ",
            conditionMessage(e), call. = FALSE)
    NULL
  })
  if (is.null(ts) || !identical(dim(ts$V), dim(fit$se$vcov))) return(fit)

  fit$se$vcov         <- ts$V
  fit$se$prob_se      <- .lta_prob_se(fit$se$blocks, ts$V)
  fit$se$conditional  <- FALSE
  fit$se$threestep    <- TRUE
  fit$se$threestep_V2 <- ts$V2
  fit$se$step1_vcov   <- ts$S11
  fit$se$method       <- ts$method
  fit
}

# ------------------------------------------------------------------------------
# The entry point.
# ------------------------------------------------------------------------------
#
# fit_lta(n_steps = 3) hands its own call over here once its refusals have run.
# Steps one and two, then step three, each from the caller's call, so every
# argument the user wrote reaches the step that reads it and no other. The
# returned object is the step-3 fit -- the model the user asked about -- with
# what it was built from attached: the step-1 fit, and the error matrices step 3
# held fixed, in the rows-assigned, columns-true orientation they were formed in.
.lta_threestep <- function(cl, env, correction, assignment) {
  s12 <- .lta_threestep_step12(cl, env, assignment = assignment,
                               correction = correction)
  fit <- .lta_threestep_step3(s12, cl, env)
  fit <- .lta_threestep_attach_vcov(fit, s12, .lta_threestep_Z(s12, cl, env)$Z)

  fit$n_steps   <- 3L
  fit$step1     <- s12$step1
  fit$threestep <- list(correction           = s12$correction,
                        assignment           = s12$assignment,
                        zero_floor           = s12$zero_floor,
                        prevalences_step1    = s12$prevalences,
                        classification_error = s12$D,
                        modal                = s12$modal,
                        alignment            = s12$alignment)
  if (!is.null(fit$mm$distal)) fit$distal <- .lta_distal_table(fit)
  fit
}
