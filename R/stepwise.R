# ==============================================================================
# Shared stepwise-estimation internals
#
# These helpers factor the step-wise portions of fit_mixture_internal() into
# reusable pieces so that add_covariates() / add_outcome() can run steps 2-3 on
# an already-fitted measurement model without repeating step 1. Both entry
# points must produce byte-identical model states for the same inputs, so any
# change here affects the fit_mixture() path too.
# ==============================================================================

# Measurement-only fit metrics, computed after step 1 and before any structural
# model touches the state. The pooled class weights contribute K - 1 free
# parameters here even when a covariate model will later replace them, because
# at step 1 they are still estimated.
.step1_metrics <- function(model_state) {
  n_params_s1  <- n_parameters(model_state$mm) + (model_state$n_components - 1)
  ll_s1        <- sum(model_state$sample_weights * model_state$lower_bound)
  resp_s1      <- exp(model_state$log_resp)
  abs_ent_s1   <- sum(model_state$sample_weights *
                        (-resp_s1 * log(resp_s1 + 1e-15)))
  # Use relative_entropy() to handle the K=1 edge case cleanly.
  rel_ent_s1   <- relative_entropy(abs_ent_s1,
                                   sum(model_state$sample_weights),
                                   model_state$n_components)
  n_eff <- model_state$n_eff

  # The sample-size-adjusted BIC uses Sclove's (1987) effective sample size,
  # (n + 2)/24, in place of n.
  list(
    ll       = ll_s1,
    n_params = n_params_s1,
    aic      = -2 * ll_s1 + 2 * n_params_s1,
    bic      = -2 * ll_s1 + log(n_eff) * n_params_s1,
    sabic    = -2 * ll_s1 + log((n_eff + 2) / 24) * n_params_s1,
    entropy  = rel_ent_s1
  )
}

# Mirror the design onto the structural sub-model so that variance code
# running inside m_step methods (which only receive the sub-model) can
# reach the strata and cluster vectors. These are kept row-aligned with the
# data the sub-model is fit on by any caller that subsets rows.
.mirror_design_onto_sm <- function(model_state) {
  if (!is.null(model_state$sm)) {
    model_state$sm$strata            <- model_state$strata
    model_state$sm$cluster           <- model_state$cluster
    model_state$sm$has_survey_design <- model_state$has_survey_design
  }
  model_state
}

# Fit the structural model on top of a completed step-1 state. The measurement
# model is frozen at this point; only $sm (and, for the ML correction, the
# joint posteriors) change. Callers are responsible for having run fit_em()
# with Y = NULL first so $log_resp holds measurement-only posteriors.
.apply_structural_steps <- function(model_state, X, Y, n_steps, correction,
                                    max_iter, se,
                                    assignment = "proportional") {
  if (is.null(Y) || is.null(model_state$sm)) return(model_state)

  if (n_steps == 2) {
    resp <- exp(model_state$log_resp)
    model_state$sm <- init_params(model_state$sm, Y, resp)
    model_state$sm <- m_step(model_state$sm, Y, resp)
    # Two-step estimation is an unadjusted third step: the posteriors act as
    # K weighted records per case, so the variance needs the same treatment
    # as the ML-adjusted path, with no classification table to correct for.
    model_state <- .attach_step3_covariate_vcov(
      model_state, X, Y, resp, NULL, model_state$sample_weights, se = se)

  } else if (n_steps == 3) {
    if (correction == "ML") {
      model_state <- fit_ml(model_state, X, Y, max_iter = max_iter, se = se,
                            assignment = assignment)
    } else if (correction == "BCH") {
      model_state <- fit_bch(model_state, X, Y, assignment = assignment)
    } else {
      # correction = "none": plain 2-step update on the structural model.
      # The measurement model is already frozen at this point; the SM is fit on
      # the posterior responsibilities from step 1 without any bias correction.
      resp <- exp(model_state$log_resp)
      model_state$sm <- init_params(model_state$sm, Y, resp)
      model_state$sm <- m_step(model_state$sm, Y, resp)
      model_state <- .attach_step3_covariate_vcov(
        model_state, X, Y, resp, NULL, model_state$sample_weights, se = se)
    }
  }

  model_state
}

# Post-fit finishing shared by every estimation path: optional class sorting,
# display names on the fitted parameter blocks, and the combined-model metrics.
.finalize_model_state <- function(model_state, X, Y, order_by_size) {
  if (order_by_size) model_state <- sort_model_classes(model_state)

  # Attach column names.
  # Guard: only assign colnames to pis when dimensions match (Bernoulli: J cols;
  # Multinoulli: J*M cols - colnames(X) has length J so would misfire there).
  if (!is.null(colnames(X)) && !is.null(model_state$mm$parameters$pis) &&
      ncol(model_state$mm$parameters$pis) == ncol(X))
    colnames(model_state$mm$parameters$pis) <- colnames(X)

  # Attach column names to betas for distal_continuous_regression.
  if (!is.null(Y) && !is.null(model_state$sm) &&
      inherits(model_state$sm, "distal_continuous_regression") &&
      !is.null(model_state$sm$parameters$betas)) {
    K_br    <- model_state$n_components
    y_names <- if (!is.null(colnames(Y))) colnames(Y) else
      paste0("V", seq_len(ncol(Y)))
    cov_names_br <- if (ncol(Y) > 1L) y_names[-1L] else character(0L)
    br_names     <- c("Intercept", cov_names_br)
    if (length(br_names) == ncol(model_state$sm$parameters$betas))
      colnames(model_state$sm$parameters$betas) <- br_names
  }

  # Attach column names to beta_pooled for distal_continuous_pooled so that
  # summary() displays real variable names instead of Z1, Z2, etc. A
  # class-specific ("moderated") covariate gets one column per class, named
  # "term:Class{k}", following its block in the beta_pooled layout.
  if (!is.null(Y) && !is.null(model_state$sm) &&
      inherits(model_state$sm, "distal_continuous_pooled") &&
      !is.null(model_state$sm$parameters$beta_pooled)) {
    K_bp    <- model_state$n_components
    y_names <- if (!is.null(colnames(Y))) colnames(Y) else
      paste0("V", seq_len(ncol(Y)))
    # First column of Y is the outcome; remaining are covariates
    cov_names_bp <- if (ncol(Y) > 1L) y_names[-1L] else character(0L)
    mod_bp       <- model_state$sm$moderated %||% integer(0)
    pooled_names_bp <- if (length(mod_bp) > 0L) cov_names_bp[-mod_bp]
                       else cov_names_bp
    mod_names_bp <- character(0L)
    for (term in cov_names_bp[mod_bp])
      mod_names_bp <- c(mod_names_bp, paste0(term, ":Class", seq_len(K_bp)))
    bp_names <- c(paste0("Class_", seq_len(K_bp)), pooled_names_bp, mod_names_bp)
    if (length(bp_names) == ncol(model_state$sm$parameters$beta_pooled))
      colnames(model_state$sm$parameters$beta_pooled) <- bp_names
  }

  # Attach covariate names to beta only when the SM has been initialised.
  # Guards against NULL beta on paths where the SM was not fit.
  if (!is.null(Y) && !is.null(model_state$sm) &&
      inherits(model_state$sm, "covariate") &&
      !is.null(model_state$sm$parameters$beta)) {
    intercept_flag <- isTRUE(model_state$sm$intercept)
    expected_D     <- ncol(model_state$sm$parameters$beta)

    y_col_names <- if (!is.null(colnames(Y))) colnames(Y)
    else paste0("V", seq_len(ncol(Y)))
    cov_names   <- if (intercept_flag) c("Intercept", y_col_names)
    else y_col_names
    if (!is.null(cov_names) && !is.null(expected_D) &&
        length(cov_names) == expected_D)
      colnames(model_state$sm$parameters$beta) <- cov_names

    # Which original variable each beta column belongs to. Stored alongside the
    # coefficients rather than re-derived from their names, because the k-1
    # dummies of a k-level factor are one term and no naming convention makes
    # that recoverable unambiguously. The intercept is its own term and is
    # never a candidate for an omnibus test.
    y_terms   <- .covariate_terms(Y)
    cov_terms <- if (intercept_flag) c("Intercept", y_terms) else y_terms
    if (length(cov_terms) == expected_D)
      model_state$sm$parameters$terms <- cov_terms
  }

  # Final metrics
  # When a covariate structural model is active it replaces the pooled class
  # weights in the likelihood entirely (see `covariate_active` in
  # `e_step()`), and `n_parameters.covariate()` already counts the full
  # (K-1)*D free regression parameters, intercepts included. Adding the
  # pooled-weights `K-1` term on top of that would double-count.
  n_params <- n_parameters(model_state$mm)
  if (!has_covariate(model_state$sm)) n_params <- n_params + (model_state$n_components - 1)
  if (!is.null(model_state$sm)) n_params <- n_params + n_parameters(model_state$sm)
  ll       <- sum(model_state$sample_weights * model_state$lower_bound)
  resp     <- exp(model_state$log_resp)
  abs_ent  <- sum(model_state$sample_weights * (-resp * log(resp + 1e-15)))
  # Use relative_entropy() to handle the K=1 edge case cleanly.
  ent      <- relative_entropy(abs_ent,
                               sum(model_state$sample_weights),
                               model_state$n_components)
  n_eff    <- model_state$n_eff

  model_state$metrics <- list(
    ll       = ll,
    n_params = n_params,
    aic      = -2 * ll + 2 * n_params,
    bic      = -2 * ll + log(n_eff) * n_params,
    sabic    = -2 * ll + log((n_eff + 2) / 24) * n_params,
    entropy  = ent,
    # How many restarts found this solution, out of how many were run. Carried
    # up from fit_em(); absent on a state that did not come through it.
    n_starts     = model_state$n_starts,
    n_replicated = model_state$n_replicated,
    # How many were asked for, which on a staged search is not how many ran to
    # convergence. Both numbers are needed: the first is what the user set and
    # what a reporting checklist asks for, the second is what the count of
    # replications is out of.
    n_requested  = model_state$n_requested
  )

  model_state
}

# The block of fit indices under the header, shared by print.mixture_model()
# and print.lta_model() so the two cannot show different sets.
#
# The six are exactly the columns of compare_mixtures()'s fit_table. That is the
# point of the choice: a user who prints one model and a user who compares a
# range must never see two different sets of numbers for the same fit. Nothing
# else belongs here.
#
# `suffix` labels which set of metrics `m` is (a three-step fit has two).
# `flag_bic` marks the BIC of a fit whose variances collapsed, where the number
# is inflated by the spike and is not comparable with a clean fit's.
.print_fit_indices <- function(m, suffix = "", flag_bic = FALSE,
                               entropy_note = "") {
  if (is.null(m)) return(invisible(NULL))
  labs <- paste0(c("Log-Likelihood", "Parameters", "AIC", "BIC", "SABIC",
                   "Rel. Entropy"), suffix)
  w    <- max(nchar(labs))
  line <- function(i, v) cat(sprintf("  %-*s : %s\n", w, labs[i], v))

  line(1, sprintf("%.2f", m$ll))
  line(2, sprintf("%d", as.integer(m$n_params)))
  line(3, sprintf("%.2f", m$aic))
  line(4, paste0(sprintf("%.2f", m$bic),
                 if (isTRUE(flag_bic))
                   " (inflated by the variance collapse; see note below)"
                 else ""))
  line(5, sprintf("%.2f", m$sabic))
  if (!is.null(m$entropy) && is.finite(m$entropy))
    line(6, paste0(sprintf("%.4f", m$entropy), entropy_note))
  invisible(NULL)
}

# How many restarts reached the reported solution, printed under the
# log-likelihood. A maximum found once is weaker evidence than the same maximum
# found repeatedly, and the difference is the main thing a user can act on: if
# it was found once, raise n_init and see whether anything better turns up.
# Silent for a single-start fit, where there is nothing to replicate.
.print_replication_note <- function(x) {
  n  <- x$metrics$n_starts
  nr <- x$metrics$n_replicated
  if (is.null(n) || is.null(nr) || !is.finite(n) || n < 2L) return(invisible(NULL))
  # Both counts where they differ. Reporting only the converged one understates
  # what was asked for - a checklist item in its own right - and reporting only
  # the requested one overstates what the replication count is out of.
  req    <- x$metrics$n_requested
  detail <- if (!is.null(req) && is.finite(req) && req != n)
    sprintf("%d of %d starts that ran to convergence (of %d requested)", nr, n, req)
  else
    sprintf("%d of %d starts", nr, n)
  cat(sprintf("  Best solution  : found by %s%s\n", detail,
              if (nr == 1L)
                paste0(" - ", .replication_advice(req %||% n, n)) else ""))
  invisible(NULL)
}

# Was the reported maximum found by exactly one start, out of enough starts for
# that to mean anything? Shared by the warning below and by the `Unreplicated`
# column of both comparison tables, so the three cannot drift apart.
#
# The threshold is on the *requested* count, not on the number that ran to
# convergence. On the staged searches only three restarts are carried to
# convergence, so a threshold on `n_starts` alone would go quiet on exactly the
# fits that need it: "1 of 3 survivors, out of 50 requested" is the case this
# exists for. On an unstaged search the two counts are equal.
#
# Below ten requested starts a lone replication carries no information, so the
# answer is FALSE rather than TRUE; NA when the fit carries no counts at all.
.is_unreplicated <- function(metrics) {
  nr <- metrics$n_replicated
  n  <- metrics$n_requested %||% metrics$n_starts
  if (is.null(nr) || is.null(n) || !is.finite(nr) || !is.finite(n)) return(NA)
  nr == 1L && n >= 10L
}

# The advice depends on how many starts were actually asked for, and on how many
# of them were run out to convergence. Telling a user who ran n_init = 200 to
# "refit with n_init = 100" is worse than saying nothing.
#
# The two counts must both be consulted, because the strongest reading - that
# the search is large enough that the fault lies in the specification rather
# than in the search - is a claim about restarts that were actually run out, and
# on a staged search those are a fraction of the ones requested (see
# fit_em(): the survivors are the better of `frac` of them, subject to a floor).
# A rule on the requested count alone would tell a user whose 100 requested
# starts were refined 20 at a time that the maximum "does not replicate at this
# many starts", when it was never tested against a hundred of them. So the
# requested count decides whether the search was large enough to be worth
# discussing (matching .is_unreplicated()'s threshold), and the converged count
# decides which of the two readings is available.
#
# `n_conv` may be NULL from a caller that carries only the requested count, in
# which case the two are equal, which is what they are on an unstaged search.
.replication_advice <- function(n_req, n_conv = NULL) {
  if (is.null(n_req) || !is.finite(n_req) || n_req < 100L)
    return("refit with n_init = 100 before reporting")

  if (is.null(n_conv) || !is.finite(n_conv)) n_conv <- n_req

  if (n_conv < 100L)
    return(sprintf(paste0(
      "the maximum failed to replicate among the %d restart%s run out to ",
      "convergence, which is a thinner test than the %d requested makes it ",
      "sound - a staged search refines only the most promising of them - so ",
      "raise n_init further before reading anything into it"),
      n_conv, if (n_conv == 1L) "" else "s", n_req))

  paste0("a search this large that still does not replicate is more likely to ",
         "be telling you about the specification than about the search. Look ",
         "at how well separated the classes are, at how heavily parameterised ",
         "the within-class structure is, and at whether there are more classes ",
         "than the data support. More starts can still help: models with many ",
         "classes or many free parameters may need several hundred")
}

# The same fact as a warning, because the printed note above is invisible to
# anyone working from summary() or from the coefficients. Raised as a classed
# condition so the comparison functions can muffle it by class rather than by
# matching its text.
.check_replication <- function(fit) {
  m <- fit$metrics
  if (is.null(m) || !isTRUE(.is_unreplicated(m))) return(invisible(NULL))
  # A collapsed variance and a growth-factor boundary both say that raising
  # n_init can make matters worse - more starts means more chances to find the
  # spike - and this warning says to raise it. They must never both fire on one
  # fit, and where they compete the collapse is the more urgent diagnosis.
  if (!is.null(fit$degenerate) || length(fit$growth$boundary))
    return(invisible(NULL))

  n_conv <- m$n_starts
  n_req  <- m$n_requested %||% n_conv
  count  <- if (!is.null(n_conv) && is.finite(n_conv) && n_conv != n_req)
    sprintf("1 of %d starts that ran to convergence, out of %d requested",
            n_conv, n_req)
  else
    sprintf("1 of %d starts", n_req)

  msg <- sprintf(paste0(
    "The reported solution was found by %s. EM climbs the peak it starts ",
    "nearest, so a maximum seen once may be the best of a small sample of the ",
    "likelihood surface rather than the best there is: %s."),
    count, .replication_advice(n_req, n_conv))

  warning(structure(class = c("mixtureEM_replication", "warning", "condition"),
                    list(message = msg, call = NULL)))
  invisible(NULL)
}

# Non-convergence, worded once and shared, so the models with their own EM
# driver say the same thing as the ones that go through fit_mixture(). The
# escalation is a doubling rather than a round number: the iteration budget is
# the thing being explored, and it is explored on a multiplicative scale.
.warn_non_convergence <- function(max_iter) {
  warning(sprintf(
    paste0("EM did not converge within max_iter = %d iterations. The ",
           "estimates are wherever the algorithm had reached, which need ",
           "not be a maximum. Refit with `max_iter = %d`, doubling again if ",
           "that is still not enough; if doubling does not help, the model is ",
           "probably weakly identified at this number of classes."),
    max_iter, 2L * max_iter), call. = FALSE)
  invisible(NULL)
}

# Normalise `slopes`: "pooled" and "class_specific" keep their existing
# all-or-nothing meaning, while a character vector of covariate *term* names
# (or a one-sided formula naming them, e.g. ~ loc1 + loc2) selects a subset to
# get a slope per class, the rest staying pooled. Shared by fit_mixture() and
# add_outcome() in place of match.arg(), which cannot express the third form.
.validate_slopes <- function(slopes) {
  if (inherits(slopes, "formula")) {
    if (length(slopes) != 2L)
      stop("`slopes` must be a one-sided formula (e.g. ~ age + sex) when ",
           "naming covariates that way.", call. = FALSE)
    return(list(mode = "moderated", terms = all.vars(slopes)))
  }
  if (is.character(slopes) && length(slopes) == 1L &&
      slopes %in% c("pooled", "class_specific"))
    return(list(mode = slopes, terms = NULL))
  if (is.character(slopes) && length(slopes) >= 1L)
    return(list(mode = "moderated", terms = slopes))
  stop('`slopes` must be "pooled", "class_specific", a character vector of ',
       "covariate names, or a one-sided formula naming them.", call. = FALSE)
}

# Translate a distal-outcome specification into the structural engine string
# and the Y matrix the fitting engine consumes. Column 1 of Y is always the
# outcome; covariates follow. Shared by fit_mixture() and add_outcome() so the
# two paths cannot drift apart.
.build_outcome_spec <- function(outcome, outcome_covariates, outcome_type,
                                slopes, cov_expr) {
  # One distal outcome per call. A multi-column `outcome` otherwise reaches the
  # coercions below as a list and dies there with a message about doubles or
  # xtfrm, which says nothing about what the caller did wrong. A single column
  # arriving as a one-column data frame or matrix is unambiguous, so it is
  # unwrapped rather than rejected.
  if (!is.null(dim(outcome))) {
    if (ncol(outcome) > 1L)
      stop(sprintf(paste("`outcome` must be a single distal outcome, but %d",
                         "columns were supplied (%s). Fit them one at a time:",
                         "the three-step estimates are identical either way,",
                         "because each outcome is regressed on the same frozen",
                         "measurement model and the outcomes do not enter each",
                         "other's structural equation."),
                   ncol(outcome),
                   paste(utils::head(colnames(outcome) %||%
                                       seq_len(ncol(outcome)), 4L),
                         collapse = ", ")),
           call. = FALSE)
    out_label <- .outcome_label(outcome)
    outcome   <- if (is.data.frame(outcome)) outcome[[1L]] else outcome[, 1L]
  } else {
    out_label <- .outcome_label(outcome)
  }

  otype <- .resolve_outcome_type(outcome, outcome_type)
  if (outcome_type == "auto")
    message(sprintf("Outcome treated as %s (set `outcome_type` to override).",
                    otype))

  has_cov <- !is.null(outcome_covariates)
  sv      <- .validate_slopes(slopes)

  if (sv$mode == "moderated" && otype != "continuous")
    stop('A character `slopes` naming covariates (or a formula) is only ',
         "supported for a continuous outcome; use \"pooled\" or ",
         '"class_specific" for a categorical outcome.', call. = FALSE)

  # "moderated" reuses the pooled engine, generalised to accept a subset of
  # covariates as class-specific -- see distal_continuous_pooled_model()'s
  # `moderated` argument.
  engine  <- if (otype == "continuous") {
    if (!has_cov)                          "continuous_outcome"
    else if (sv$mode == "class_specific")  "continuous_outcome_moderated"
    else                                   "continuous_outcome_adjusted"
  } else {
    if (!has_cov)                  "categorical_outcome"
    else if (sv$mode == "pooled")  "categorical_outcome_adjusted"
    else                           "categorical_outcome_moderated"
  }

  # The outcome is coerced to a plain numeric column (categorical outcomes to
  # 1-indexed integer codes) so it is never dummy-coded as though it were a
  # covariate.
  if (otype == "categorical") {
    out_col <- as.integer(as.factor(outcome))
  } else {
    out_col <- suppressWarnings(as.numeric(outcome))
    if (anyNA(out_col) && !anyNA(outcome))
      stop("A continuous `outcome` must be numeric.", call. = FALSE)
  }
  out_mat <- matrix(out_col, ncol = 1L, dimnames = list(NULL, out_label))

  Y <- if (has_cov) .cbind_covariates(out_mat, prepare_covariates(
    .as_named_covariates(outcome_covariates, cov_expr, "covariate")))
  else out_mat

  # Resolve the named terms to covariate-column indices. A factor expands into
  # several dummy columns under the same term, which substring matching would
  # miss; .covariate_terms() already knows the grouping prepare_covariates()
  # set.
  moderated <- integer(0)
  if (sv$mode == "moderated") {
    cov_terms  <- if (has_cov) .covariate_terms(Y)[-1L] else character(0L)
    unresolved <- setdiff(sv$terms, cov_terms)
    if (length(unresolved)) {
      avail <- if (length(cov_terms)) paste(unique(cov_terms), collapse = ", ")
               else "(none)"
      stop(sprintf(
        "`slopes` names %s not found among the outcome covariates: %s. Available terms: %s.",
        if (length(unresolved) > 1L) "terms" else "a term",
        paste(unresolved, collapse = ", "), avail), call. = FALSE)
    }
    moderated <- which(cov_terms %in% sv$terms)
  }

  list(engine = engine, Y = Y, otype = otype, moderated = moderated)
}
