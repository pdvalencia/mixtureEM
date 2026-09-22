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
