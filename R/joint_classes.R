# ==============================================================================
# Two latent class variables in one model
# ==============================================================================
#
# A loglinear model with two latent class variables (Hagenaars 1993): `x`,
# measured by its own items, and `y`, measured by its own, with the
# association written as a multinomial logit of `y` on `x`. That regression
# is exactly a two-occasion LTA's by-origin transition regression, so the fit
# rides fit_lta(): occasion 1 carries x's items, occasion 2 y's, each padded
# with all-NA columns to the longer list, measurement free by occasion, a
# smaller `y` written as a padded status whose destination column is
# forbidden, and a smaller `x` as a padded status no case starts in. The
# padding is bookkeeping only and must not reach the user: its
# parameters touch no data, so they are taken out of the count (and every
# criterion built on it), and the report speaks of `x` and `y`, never of
# occasions or transitions.

#' Fit a model with two associated latent class variables
#'
#' Fits two latent class variables at once, each measured by its own items:
#' `x` with `n_x` classes and `y` with `n_y` classes, linked by a
#' multinomial logistic regression of `y` on `x` (Hagenaars, 1993). The
#' association is reported against the last class of each variable: the
#' intercept of `y` class j is its log-odds against the last `y` class among
#' members of the last `x` class, and the coefficient of `x` class k is how
#' much membership in `x` class k changes that log-odds. The joint class
#' probabilities and the probabilities of each `y` class given each `x`
#' class come with it.
#'
#' Covariates can predict either variable. On `y` each covariate acts in
#' common across the `x` classes or separately within each, set per
#' covariate by `y_effects`; the second form is a covariate whose effect on
#' `y` is moderated by `x`.
#'
#' `association = "none"` fits the same model with the association removed:
#' the two variables are independent given the covariates. It is the
#' restricted model of the test of association, so
#' `lr_test(fit_none, fit_free)` asks whether the two class variables are
#' related at all, on \eqn{(K_{x} - 1)(K_{y} - 1)} degrees of freedom; a
#' covariate acting on `y` by `x` class against in common is the next test
#' in the same sequence. Without covariates the model is the two separate
#' latent class analyses.
#'
#' @param x_items,y_items The items measuring each latent variable, one row
#'   per case, the same cases in the same order. Missing values are allowed.
#' @param n_x,n_y Integer. Number of classes of each variable.
#' @param measurement_x,measurement_y `"binary"` (items coded 0/1) or
#'   `"ordinal"` (ordered categories coded from 0 or from 1). A binary and an
#'   ordinal side are fitted together as ordinal items, which for a two-category
#'   item is the same model.
#' @param predictors_x A data frame of covariates predicting `x` class
#'   membership, or `NULL`.
#' @param predictors_y A data frame of covariates predicting `y` class
#'   membership, or `NULL`.
#' @param y_effects How each covariate in `predictors_y` acts: `"common"`
#'   (one slope per `y` class, the default), `"by_x"` (a slope per `x`
#'   class for every covariate), or a named vector such as
#'   `c(female = "by_x")`, the unnamed covariates acting in common.
#' @param association `"free"` (the default) estimates the association;
#'   `"none"` holds the two variables independent given the covariates, so
#'   `y_effects` cannot name `"by_x"`.
#' @param predictors_items Direct effects of covariates on items (uniform
#'   differential item functioning): a named list, item name = a data frame of
#'   the covariates acting on that item, e.g. `list(a3 = data.frame(female))`.
#'   Each slope is one proportional-odds shift of the item's response, the same
#'   in every class; binary items take it on the logit of the second category.
#'   Items are named by the column names of `x_items` and `y_items`, which
#'   must then be unique across the two. Covariate patterns are limited as in
#'   [`fit_lta()`]. A DIF covariate should normally also predict the classes
#'   (`predictors_x`, `predictors_y`); with `n_steps = 2` step 1 does this
#'   itself.
#' @param n_steps `1` (the default) fits everything at once. `2` is the
#'   two-step estimator (Bakk & Kuha, 2018): step 1 fits each variable's
#'   measurement model on its own, as two separate latent class analyses;
#'   step 2 holds those item parameters fixed and estimates the association
#'   and the covariate effects. The classes then mean the same thing in every
#'   structural model compared, and the standard errors carry step 1's
#'   uncertainty. The fit keeps step 1 as `step1`, which is the restricted
#'   model of a test of every structural parameter at once. Two two-step fits
#'   are compared by [`lr_test()`] only when they share step 1, so fit them
#'   with the same `random_state` or hand them the same `step1`. With
#'   `predictors_items`, step 1 also regresses each variable on the DIF
#'   covariates, and step 2
#'   holds the item slopes fixed with the rest of the measurement model.
#'   Step 1 is searched as the two latent class analyses it factorises into,
#'   each with its own `n_init` restarts through [`fit_mixture()`], and then
#'   continued as one model from their solution; arguments passed in `...`
#'   (such as `weights`) search it as one model instead.
#' @param step1 With `n_steps = 2`, a step 1 already fitted, so that a
#'   sequence of structural models shares one measurement model and fits it
#'   once: either a previous two-step fit of the same items with the same
#'   numbers of classes and item slopes, whose step 1 is reused as it stands,
#'   or `list(x = , y = )`, the two latent class analyses as
#'   [`fit_mixture()`] fits (classes on the item-slope covariates, with the
#'   same item slopes), from which step 1 is continued. `NULL` (the default)
#'   fits step 1.
#' @param smoothing,bayes_constants,n_init,random_state,n_cores,standard_errors
#'   As in [`fit_lta()`].
#' @param ... Further arguments passed to [`fit_lta()`], such as `weights`.
#'
#' @return An object of class `joint_classes` (inheriting from `lta_model`, so
#'   [`lr_test()`] and the information criteria apply), with `n_x`, `n_y`,
#'   `association` (from [`joint_coefficients()`]), `joint` (the `n_x x n_y`
#'   joint class probabilities) and `conditional` (P(y | x), rows summing to
#'   one).
#'
#' @references
#' Bakk, Z., & Kuha, J. (2018). Two-step estimation of models between latent
#' classes and external variables. *Psychometrika*, 83(4), 871-892.
#' \doi{10.1007/s11336-017-9592-7}
#'
#' Hagenaars, J. A. (1993). *Loglinear Models with Latent Variables*. Sage.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' n  <- 600
#' cx <- sample(1:2, n, TRUE)
#' cy <- ifelse(runif(n) < ifelse(cx == 1, 0.8, 0.2), 1L, 2L)
#' A  <- sapply(1:4, function(j) rbinom(n, 1, ifelse(cx == 1, 0.85, 0.15)))
#' S  <- sapply(1:4, function(j) rbinom(n, 1, ifelse(cy == 1, 0.85, 0.15)))
#' fit <- fit_joint_classes(A, S, n_x = 2, n_y = 2, n_init = 5,
#'                          random_state = 1)
#' fit
#' }
#' @export
fit_joint_classes <- function(x_items, y_items,
                              n_x = 2, n_y = 2,
                              measurement_x = c("binary", "ordinal"),
                              measurement_y = c("binary", "ordinal"),
                              predictors_x = NULL,
                              predictors_y = NULL,
                              y_effects = "common",
                              association = c("free", "none"),
                              predictors_items = NULL,
                              n_steps = 1,
                              step1 = NULL,
                              smoothing = 1.0,
                              bayes_constants = NULL,
                              n_init = 20,
                              random_state = NULL,
                              n_cores = .default_n_cores(),
                              standard_errors = TRUE,
                              ...) {
  measurement_x <- match.arg(measurement_x)
  measurement_y <- match.arg(measurement_y)
  association <- match.arg(association)
  n_steps <- as.integer(n_steps)[1L]
  if (is.na(n_steps) || !(n_steps %in% 1:2))
    stop("`n_steps` must be 1 or 2.", call. = FALSE)
  A <- as.matrix(x_items); S <- as.matrix(y_items)
  if (nrow(A) != nrow(S))
    stop("`x_items` and `y_items` must have the same number of rows.",
         call. = FALSE)
  n_x <- as.integer(n_x); n_y <- as.integer(n_y)
  if (n_x < 2L || n_y < 2L)
    stop("`n_x` and `n_y` must each be at least 2.", call. = FALSE)

  dif <- !is.null(predictors_items)
  if (dif) {
    if (!is.list(predictors_items) || is.data.frame(predictors_items))
      stop("`predictors_items` must be a named list, item name = the ",
           "covariates acting on that item, e.g. list(a3 = data.frame(female)).",
           call. = FALSE)
    nm <- c(colnames(A), colnames(S))
    if (length(nm) != ncol(A) + ncol(S) || anyDuplicated(nm) || any(!nzchar(nm)))
      stop("`predictors_items` names items, so `x_items` and `y_items` need ",
           "column names, unique across the two.", call. = FALSE)
  }
  # A direct item effect is a proportional-odds shift, so DIF fits every item
  # as ordinal; a binary item is the two-category case of the same model.
  ordinal <- dif || "ordinal" %in% c(measurement_x, measurement_y)
  if (ordinal) {
    # The ordinal emission reads categories from 1; a 0/1 item becomes 1/2,
    # the same two-category model.
    up <- function(M) apply(M, 2, function(x)
      if (isTRUE(min(x, na.rm = TRUE) == 0)) x + 1 else x)
    A <- up(A); S <- up(S)
  }
  Ja <- ncol(A); Js <- ncol(S); J <- max(Ja, Js); K <- max(n_x, n_y)
  # Paired: occasion t carries one variable's items, padded to the longer
  # list, measurement free by occasion. Stacked: each occasion carries every
  # item, the other variable's all missing, and the measurement is held equal
  # across occasions: each item is observed at one occasion only, so the
  # equality restricts nothing, and an item is one column of one item list --
  # which is what an item slope is defined on. Same likelihood either way.
  # Stacked is needed with DIF, and wherever two paired items differ in their
  # number of categories: the ordinal emission gives each column position one
  # category count, so a binary item paired with a three-category one would
  # carry an empty third category, its thresholds counted as parameters and
  # priced by the prior.
  n_cat <- function(M) apply(M, 2, function(x)
    if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE))
  m <- min(Ja, Js)
  stacked <- dif || (ordinal &&
    any(n_cat(A)[seq_len(m)] != n_cat(S)[seq_len(m)], na.rm = TRUE))
  if (stacked) {
    na <- function(M) matrix(NA_real_, nrow(M), ncol(M))
    X <- cbind(A, na(S), na(A), S)
  } else {
    pad <- function(M) cbind(M, matrix(NA_real_, nrow(M), J - ncol(M)))
    X <- cbind(pad(A), pad(S))
  }

  # A smaller `y`: its padded classes are destinations no origin may move
  # into. A smaller `x`: its padded classes are statuses no case starts in
  # (`.initial_statuses`), so the moves out of them are never read; without
  # a transition regression each such row is closed to every destination but
  # one, which leaves it no free parameter. With a regression, or with no
  # association (one row for every origin), the engine leaves them out itself.
  forb <- NULL
  if (n_y < K) { forb <- matrix(FALSE, K, K); forb[, (n_y + 1L):K] <- TRUE }
  if (n_x < K && is.null(predictors_y) && association == "free") {
    forb <- matrix(FALSE, K, K); forb[(n_x + 1L):K, -K] <- TRUE
  }

  effects <- "common"
  if (!is.null(predictors_y)) {
    predictors_y <- as.data.frame(predictors_y)
    effects <- .joint_y_effects(y_effects, names(predictors_y))
  }
  # With no association there is no `x` class for a covariate's effect on
  # `y` to differ by: the model is the transition-free fit, whose design for
  # `y` has no origin in it (.lta_tau_design()).
  none <- association == "none"
  if (none && !identical(effects, "common"))
    stop("`y_effects = \"by_x\"` lets a covariate's effect on `y` differ ",
         "by `x` class, which needs the association; with `association = ",
         "\"none\"` every covariate acts in common.", call. = FALSE)
  # Without any covariate the no-association model is exactly the two-step's
  # step 1, the two separate latent class analyses, and step 2 has nothing to
  # estimate that step 1 did not: the one-step fit is the two-step estimate.
  # It is fitted once, with the search a two-step's step 1 runs, so its
  # measurement block is the one a two-step fit with the same `random_state`
  # holds fixed.
  structural <- !is.null(predictors_x) || !is.null(predictors_y) || dif
  redirect <- none && n_steps == 2L && !structural

  # Step 1, as two LCAs or handed over (.joint_fit_lcas(), .joint_step1_arg()).
  # Arguments for fit_lta() that the LCAs would not see (`weights` and the
  # rest of `...`) keep the joint search.
  if (!is.null(step1) && n_steps != 2L)
    stop("`step1` is the measurement model a two-step fit holds fixed; it ",
         "needs `n_steps = 2`.", call. = FALSE)
  s1_fit <- NULL; writer <- NULL
  Zi <- if (dif) .lta_dif_spec(predictors_items, NULL, names(predictors_items),
                               nrow(A))$Z
  # What step 1 is a fit of: the items, and the item slopes it carries.
  key <- list(items = c(.joint_item_key(A), .joint_item_key(S)),
              dif = if (dif) list(sort(names(predictors_items)),
                                  colSums(Zi, na.rm = TRUE)))
  if (n_steps == 2L) {
    lcas <- NULL
    if (!is.null(step1)) {
      got <- .joint_step1_arg(step1, A, S, n_x, n_y, ordinal, key)
      s1_fit <- got$fit; lcas <- got$lcas
    } else if (!length(list(...))) {
      lcas <- .joint_fit_lcas(A, S, n_x, n_y, ordinal, Zi, predictors_items,
                              smoothing, bayes_constants, n_init, random_state,
                              n_cores)
    }
    if (!is.null(lcas)) writer <- .joint_lca_writer(lcas, Ja, Js, stacked)
    # The redirected fit is step 1 itself.
    if (redirect && !is.null(s1_fit)) {
      donor <- s1_fit
      writer <- function(state, X) .lta_refine_start(state, X, donor)
    }
  }
  # A redirected fit started from the LCAs runs once from them; a two-step
  # carries them to its step 1 (fit_lta(), "the two-step estimator").
  start_w <- if (redirect) writer
  if (!redirect) s1_writer <- writer else s1_writer <- NULL
  n_init_lta <- if (!is.null(start_w)) 1L else n_init

  # The two-step's step 1 is the two separate LCAs: with no association and
  # no covariates the joint likelihood factorises into them, which is the
  # transition-free fit (`.twostep_independent`, R/lta.R).
  fit <- fit_lta(X, n_statuses = K, times = 2,
                 n_steps = if (redirect) 1L else n_steps,
                 .twostep_independent = n_steps == 2L && !redirect,
                 .transition_free = none,
                 measurement = if (ordinal) "ordinal" else "binary",
                 measurement_invariance = if (stacked) "full" else "none",
                 item_names = if (dif) nm,
                 predictors_items = predictors_items,
                 forbidden_transitions = forb,
                 .initial_statuses = if (n_x < K) n_x,
                 .prior_statuses = if (!stacked && (n_x != n_y || Ja != Js))
                   c(n_x, n_y),
                 predictors_initial = predictors_x,
                 predictors_transition = predictors_y,
                 transition_effects = effects,
                 smoothing = smoothing, bayes_constants = bayes_constants,
                 n_init = n_init_lta, random_state = random_state,
                 n_cores = n_cores, standard_errors = standard_errors,
                 .start_writer = start_w, .twostep_step1 = s1_fit,
                 .twostep_writer = s1_writer, ...)

  # The padding's parameters come out of the count, and every criterion built
  # on it, on the fit and on its step 1 alike: step 1 is a model a two-step fit
  # can be tested against.
  # A step 1 handed over has been through this already.
  uncount <- function(f) {
    if (!is.null(f$joint_key)) return(f)
    f$n_x <- n_x; f$n_y <- n_y
    f$n_items_x <- Ja; f$n_items_y <- Js; f$joint_stacked <- stacked
    f$joint_key <- key
    f$n_params <- f$n_params - .joint_padding_n_params(f)
    f$metrics$n_params <- f$n_params
    fresh <- .lta_metrics(f)
    for (nm in intersect(c("aic", "bic", "caic", "awe", "aic3", "icl", "sabic"),
                         names(fresh)))
      f$metrics[[nm]] <- fresh[[nm]]
    f
  }
  fit <- uncount(fit)
  if (inherits(fit$step1, "lta_model")) fit$step1 <- uncount(fit$step1)
  if (redirect) fit$n_steps <- 2L
  fit$association_free <- !none
  # The items' own names, for measurement_summary(): the occasion layout
  # above names each column position once, after the first variable's items.
  fit$item_names_x <- colnames(A) %||% sprintf("ca_item%d", seq_len(Ja))
  fit$item_names_y <- colnames(S) %||% sprintf("cs_item%d", seq_len(Js))

  # P(x, y), averaged over the cases. With covariates on both variables and
  # no association the reported `tau` is the plain average of y's
  # distribution, not weighted by who is in each x class, so the product of
  # the two averages is not the joint table; it is summed case by case.
  P <- if (none && !is.null(fit$delta_beta) && !is.null(fit$tau_beta)) {
    w  <- fit$weights_vec
    Pa <- exp(.lta_log_delta(fit))
    Pb <- exp(.lta_log_tau(fit)[[1]][[1]])
    crossprod(Pa * w, Pb)[seq_len(n_x), seq_len(n_y), drop = FALSE] / sum(w)
  } else (fit$delta * fit$tau[[1]])[seq_len(n_x), seq_len(n_y), drop = FALSE]
  dimnames(P) <- list(paste0("X#", seq_len(n_x)), paste0("Y#", seq_len(n_y)))
  fit$joint <- P
  fit$conditional <- P / rowSums(P)
  class(fit) <- c("joint_classes", class(fit))
  fit$association <- joint_coefficients(fit)
  fit
}

# ------------------------------------------------------------------------------
# The two-step's step 1 from the two separate LCAs
# ------------------------------------------------------------------------------
#
# Step 1 is the model with no association, whose likelihood factorises into
# one latent class analysis per variable. So it is searched as two LCAs, each
# with its own restarts, rather than as one stacked model over pairs of
# starts: two small searches instead of one large one, and the best solution
# of each found independently of the other's. The joint step 1 then runs once
# from their solution (.joint_lca_writer()) under its own objective, so what
# it holds fixed in step 2 is still a maximum of the step-1 model as fitted
# here, prior included; the LCAs only decide where it starts.

# The data a step 1 was fitted on, in a form that survives a recoding: item
# by item, the count and the sum of the observed values.
.joint_item_key <- function(M) {
  M <- as.matrix(M)
  c(ncol(M), colSums(!is.na(M)), colSums(M, na.rm = TRUE))
}

# Fit the two LCAs of step 1. `Z` is the design step 1 regresses each
# variable on (the item-slope covariates, or none), `dif` the joint call's
# `predictors_items`.
.joint_fit_lcas <- function(A, S, n_x, n_y, ordinal, Z, dif,
                            smoothing, bayes_constants, n_init, random_state,
                            n_cores) {
  bc <- bayes_constants %||% list()
  if (is.null(bc$latent)) bc$latent <- smoothing
  one <- function(Y, K) {
    args <- list(Y, n_classes = K,
                 measurement = if (ordinal) "ordinal" else "binary",
                 smoothing = smoothing, bayes_constants = bc,
                 n_init = n_init, random_state = random_state,
                 n_cores = n_cores, standard_errors = FALSE)
    if (ordinal)
      args$cats <- as.integer(apply(Y, 2, max, na.rm = TRUE))
    if (!is.null(Z)) {
      Zd <- as.data.frame(Z)
      args$data <- Zd
      args$predictors <- stats::reformulate(colnames(Zd))
      mine <- intersect(names(dif), colnames(Y))
      if (length(mine))
        args$predictors_items <- stats::setNames(lapply(mine, function(it)
          stats::reformulate(colnames(as.data.frame(dif[[it]])))), mine)
    }
    do.call(fit_mixture, args)
  }
  list(x = one(A, n_x), y = one(S, n_y))
}

# A start writer for fit_lta(`.start_writer`): the two LCAs' solution written
# into the joint step 1's freshly drawn state. `x` fills the statuses at
# occasion 1 and the `x` items, `y` the destinations and the `y` items;
# a padded status keeps its draw, which no case reads. A block the LCAs do
# not carry in the same terms (a regression on other covariates, a missing
# item slope) keeps its draw too, and step 1's own EM fits it.
.joint_lca_writer <- function(lcas, Ja, Js, stacked) {
  force(lcas)
  function(state, X) {
    K <- state$n_statuses
    put_items <- function(state, f, v) {
      Jv <- if (v == 1L) Ja else Js
      Kv <- f$n_components
      cols <- if (stacked && v == 2L) Ja + seq_len(Jv) else seq_len(Jv)
      mods <- if (stacked) seq_along(state$mm$models) else v
      p <- f$mm$parameters$pis
      fc <- f$mm$cats
      if (!is.null(state$ri)) {
        th <- .ordinal_theta_from_pis(p, fc)
        for (j in seq_len(Jv))
          state$ri$theta[seq_len(Kv), .ordinal_theta_cols(state$ri$cats, cols[j])] <-
            th[, .ordinal_theta_cols(fc, j), drop = FALSE]
      } else {
        for (t in mods) {
          m <- state$mm$models[[t]]
          for (j in seq_len(Jv)) {
            if (is.null(m$cats)) {
              m$parameters$pis[seq_len(Kv), cols[j]] <- p[, j]
            } else {
              m$parameters$pis[seq_len(Kv), .ordinal_pis_cols(m$cats, cols[j])] <-
                p[, .ordinal_pis_cols(fc, j), drop = FALSE]
            }
          }
          state$mm$models[[t]] <- m
        }
      }
      sl <- f$mm$parameters$dif_slopes
      if (!is.null(state$dif) && !is.null(sl)) {
        for (d in seq_along(f$mm$dif$covariates)) {
          dl <- match(f$mm$dif$covariates[d], state$dif$names)
          if (is.na(dl)) next
          for (j in seq_len(Jv))
            state$dif$beta[seq_len(Kv), cols[j], dl] <- sl[, .dif_col(Jv, j, d)]
        }
      }
      state
    }
    fa <- lcas$x; fs <- lcas$y
    state <- put_items(state, fa, 1L)
    state <- put_items(state, fs, 2L)
    if (!is.null(state$ri))
      for (t in seq_along(state$mm$models))
        state$mm$models[[t]]$parameters$pis <-
          .lta_ri_integrated_pis(state$ri, K, state$n_items)

    pad <- function(w) c(w, numeric(K - length(w)))
    state$delta_c[[1L]] <- pad(fa$weights)
    for (m in seq_along(state$tau_c[[1L]])) {
      P <- matrix(pad(fs$weights), K, K, byrow = TRUE) *
        state$tau_allowed_c[[1L]][[m]]
      ok <- rowSums(P) > 0
      state$tau_c[[1L]][[m]][ok, ] <- P[ok, , drop = FALSE] / rowSums(P)[ok]
    }
    beta_of <- function(f, cn, rows) {
      B <- f$sm$parameters$beta
      if (is.null(B) || is.null(cn) || !all(cn %in% colnames(B)) ||
          nrow(B) != rows) return(NULL)
      B[, cn, drop = FALSE]
    }
    if (!is.null(state$Z_delta))
      state$delta_beta <- beta_of(fa, colnames(state$Z_delta),
                                  length(.lta_delta_keep(state)))
    if (!is.null(state$Z_tau) && isTRUE(state$tau_independent)) {
      B <- beta_of(fs, colnames(.lta_tau_design(state, 1L)),
                   length(.lta_tau_keep(state)))
      if (!is.null(B)) state$tau_beta <- rep(list(B), length(state$tau_c[[1L]]))
    }
    .lta_seed_regressions(state)
  }
}

# `step1` as fit_joint_classes() takes it: a previous two-step fit (its step
# 1 is reused as it stands), or the two LCAs. Checked against this call's
# items and class counts; returns list(fit = <lta step 1> or NULL, lcas =
# <the pair> or NULL).
.joint_step1_arg <- function(step1, A, S, n_x, n_y, ordinal, key) {
  if (inherits(step1, "joint_classes")) {
    s1 <- if (inherits(step1$step1, "lta_model")) step1$step1 else
      if (identical(step1$n_steps, 2L) && !isTRUE(step1$association_free)) step1
    if (is.null(s1))
      stop("`step1` is a one-step fit, which has no step 1 to hand over; give ",
           "a two-step fit (`n_steps = 2`) or the two LCAs.", call. = FALSE)
    if (!identical(c(step1$n_x, step1$n_y), c(n_x, n_y)) ||
        !isTRUE(all.equal(step1$joint_key, key)))
      stop("`step1` was fitted to other items, other item slopes or other ",
           "numbers of classes than this call's.", call. = FALSE)
    return(list(fit = s1, lcas = NULL))
  }
  if (is.list(step1) && !is.null(step1$x) && !is.null(step1$y) &&
      inherits(step1$x, "mixture_model") && inherits(step1$y, "mixture_model")) {
    # The joint fit reads a 0/1 item from 1 when it is fitted as ordinal, so
    # an LCA of the same items may be one off in either direction.
    same <- function(Y, M) any(vapply(-1:1, function(s)
      isTRUE(all.equal(.joint_item_key(Y), .joint_item_key(M + s))), TRUE))
    ok <- function(f, M, K) identical(as.integer(f$n_components), K) &&
      same(f$data, M) && (!is.null(f$mm$cats)) == ordinal
    if (!ok(step1$x, A, n_x) || !ok(step1$y, S, n_y))
      stop("`step1$x` and `step1$y` must be fit_mixture() fits of `x_items` ",
           "with `n_x` classes and of `y_items` with `n_y` classes, ",
           if (ordinal) "ordinal" else "binary", " as this fit is.",
           call. = FALSE)
    return(list(fit = NULL, lcas = step1[c("x", "y")]))
  }
  stop("`step1` must be a two-step fit_joint_classes() fit or ",
       "list(x = , y = ) of two fit_mixture() fits.", call. = FALSE)
}

# `y_effects` in fit_lta()'s terms: "by_x" for every covariate is the named
# vector naming them all, which is the same model as a separate regression per
# origin (the origin dummies already give each origin its own intercept).
.joint_y_effects <- function(y_effects, covs) {
  if (is.null(names(y_effects))) {
    eff <- match.arg(y_effects, c("common", "by_x"))
    if (eff == "common") return("common")
    return(stats::setNames(rep("by_origin", length(covs)), covs))
  }
  bad <- setdiff(unname(y_effects), c("common", "by_x"))
  if (length(bad))
    stop("`y_effects` values must be \"common\" or \"by_x\".", call. = FALSE)
  unknown <- setdiff(names(y_effects), covs)
  if (length(unknown))
    stop("`y_effects` names covariates not in `predictors_y`: ",
         paste(unknown, collapse = ", "), call. = FALSE)
  stats::setNames(ifelse(y_effects == "by_x", "by_origin", "common"),
                  names(y_effects))
}

# Free item parameters the padding added: the padded item positions on each
# side, and every item of the padded classes (`y` items of a padded `y`
# class, `x` items of a padded `x` class). Each cell costs what its
# column costs per class (one for a binary item, categories - 1 for an ordinal
# one). The padded destination column's transition cells are forbidden, so
# fit_lta() never counted them; with covariates on `y` the regression has no
# category for it (.lta_tau_keep()), so the same holds. The stacked layout has
# no padded item positions, only the padded classes' `y` items.
.joint_padding_n_params <- function(fit) {
  K  <- fit$n_statuses
  Ja <- fit$n_items_x; Js <- fit$n_items_y; J <- max(Ja, Js)
  cats <- fit$mm$models[[1]]$cats
  if (isTRUE(fit$joint_stacked))
    return((K - fit$n_y) * sum(as.integer(cats[Ja + seq_len(Js)]) - 1L) +
             (K - fit$n_x) * sum(as.integer(cats[seq_len(Ja)]) - 1L))
  per_item <- if (is.null(cats)) rep(1L, J) else as.integer(cats) - 1L
  pad_a <- if (Ja < J) K * sum(per_item[(Ja + 1L):J]) else 0L
  pad_s <- if (Js < J) K * sum(per_item[(Js + 1L):J]) else 0L
  pad_k <- (K - fit$n_y) * sum(per_item[seq_len(Js)]) +
    (K - fit$n_x) * sum(per_item[seq_len(Ja)])
  pad_a + pad_s + pad_k
}

#' Association coefficients of a two-latent-variable model
#'
#' The multinomial logit of `y` on `x` from [`fit_joint_classes()`], in the
#' reference-class form: every log-odds is against the last `y` class, and
#' `x` enters as indicators against the last `x` class. Rows are the
#' intercept (`term = "Intercept"`, the log-odds among the last `x` class),
#' the `x` effects (`term = "X#k"`), and each covariate on `y`: one row for
#' a common slope, or, for a covariate acting by `x` class, its slope in the
#' last `x` class (`x = NA`) and the difference from it in each other class.
#'
#' @param fit A `joint_classes` fit.
#' @return A data frame with columns `y`, `term`, `x`, `estimate`, `se`,
#'   `z`, `p`. Each row is also a parameter of [`wald_test()`] and
#'   [`derived_estimates()`], under the label `"Y#j Intercept"`,
#'   `"Y#j ON X#k"`, `"Y#j ON x"` or, for the difference of a covariate's
#'   slope in `x` class k from the last class's, `"Y#j ON x | X#k"`.
#' @export
joint_coefficients <- function(fit) {
  stopifnot(inherits(fit, "joint_classes"))
  K1 <- fit$n_x; K2 <- fit$n_y
  V  <- fit$se$vcov
  np <- if (is.null(V)) 0L else ncol(V)
  eta <- .joint_eta(fit)
  cov <- !is.null(fit$tau_beta)
  cn  <- if (cov) colnames(.lta_tau_design(fit, 1L)) else character()

  # log-odds of y class j against the last, for origin a
  lo <- function(a, j, term) {
    e1 <- eta(a, j, term); e2 <- eta(a, K2, term)
    list(v = e1$v - e2$v, r = e1$r - e2$r)
  }
  minus <- function(x, y) list(v = x$v - y$v, r = x$r - y$r)

  rows <- list(); grad <- list()
  add <- function(j, term, origin, e) {
    grad[[length(grad) + 1L]] <<- e$r
    rows[[length(rows) + 1L]] <<- data.frame(
      y = sprintf("Y#%d", j), term = term, x = origin, estimate = e$v,
      se = if (np) sqrt(max(drop(crossprod(e$r, V %*% e$r)), 0)) else NA_real_)
  }
  covs <- if (cov) cn[!grepl("^(Intercept$|from:)", cn)] else character()
  # With no association the `x` rows are zero by construction, not
  # estimates, so they are left out rather than printed with a zero SE.
  free <- !isFALSE(fit$association_free)
  for (j in seq_len(K2 - 1L)) {
    add(j, "Intercept", NA_character_, lo(K1, j, "Intercept"))
    if (free) for (k in seq_len(K1 - 1L))
      add(j, sprintf("X#%d", k), sprintf("X#%d", k),
          minus(lo(k, j, "Intercept"), lo(K1, j, "Intercept")))
    for (v in covs) {
      add(j, v, NA_character_, lo(K1, j, v))
      if (cov && any(sprintf("from:1 x %s", v) == cn))
        for (k in seq_len(K1 - 1L))
          add(j, v, sprintf("X#%d", k), minus(lo(k, j, v), lo(K1, j, v)))
    }
  }
  out <- do.call(rbind, rows)
  out$z <- out$estimate / out$se
  out$p <- 2 * stats::pnorm(-abs(out$z))
  # Each row is a linear function of the packed parameters; the gradient rows
  # give the rows' joint covariance (.wald_parameters(), derived_estimates()).
  attr(out, "gradient") <- do.call(rbind, grad)
  out
}

# The labels .wald_parameters() gives joint_coefficients()'s rows.
.joint_labels <- function(a) {
  is_x <- grepl("^X#", a$term)
  ifelse(a$term == "Intercept", paste(a$y, "Intercept"),
         ifelse(is_x, paste(a$y, "ON", a$term),
                paste0(a$y, " ON ", a$term,
                       ifelse(is.na(a$x), "", paste0(" | ", a$x)))))
}

# eta(a, b, term): the linear predictor of destination b for origin a, as a
# value and a row over the packed parameter vector, so every reported
# coefficient is a difference of etas and its SE one quadratic form. Without
# covariates each origin's row is a softmax anchored at its last allowed
# destination; the block's `ref` lists the free destinations in packing order.
# With covariates the regression's rows are its own destinations
# (.lta_tau_keep(): the `y` classes, the padded ones dropped), the last the
# reference.
.joint_eta <- function(fit) {
  K  <- fit$n_statuses
  V  <- fit$se$vcov
  np <- if (is.null(V)) 0L else ncol(V)
  if (!is.null(fit$tau_beta)) {
    B  <- fit$tau_beta[[1]]; cn <- colnames(.lta_tau_design(fit, 1L))
    blk <- Find(function(b) grepl("^tau_beta", b$name), fit$se$blocks)
    keep <- .lta_tau_keep(fit)
    return(function(a, b, term) {
      nm <- if (term == "Intercept") c("Intercept", if (a < fit$n_x) sprintf("from:%d", a))
            else c(term, intersect(sprintf("from:%d x %s", a, term), cn))
      r <- numeric(np)
      d <- match(b, keep)
      if (d == length(keep)) return(list(v = 0, r = r))
      # No association: the design has no origin columns, so every origin
      # reads the same coefficients.
      idx <- match(intersect(nm, cn), cn)
      if (!is.null(blk)) r[blk$cols[(d - 1L) * length(cn) + idx]] <- 1
      list(v = sum(B[d, idx]), r = r)
    })
  }
  lt <- log(pmax(fit$tau[[1]], 1e-300))
  function(a, b, term) {
    r <- numeric(np)
    blk <- Find(function(x) x$name == sprintf("tau(1)[from %d]", a),
                fit$se$blocks)
    if (!is.null(blk) && b %in% blk$ref) r[blk$cols[match(b, blk$ref)]] <- 1
    list(v = lt[a, b], r = r)
  }
}

#' @export
print.joint_classes <- function(x, digits = 3, ...) {
  cat(sprintf("Two latent class variables: x (%d classes, %d items), y (%d classes, %d items)\n",
              x$n_x, x$n_items_x, x$n_y, x$n_items_y))
  cat(sprintf("Log-likelihood %.3f, %d free parameters, BIC %.3f\n\n",
              x$metrics$ll, as.integer(x$n_params), x$metrics$bic))
  if (identical(x$n_steps, 2L))
    cat("Two-step: item parameters fixed at the separate latent class analyses\n\n")
  if (isFALSE(x$association_free))
    cat("No association: x and y independent given the covariates\n")
  cat("Association: log-odds of y against its last class\n")
  a <- x$association
  a$x[is.na(a$x)] <- ""
  print(format(a, digits = digits), row.names = FALSE)
  cat("\nJoint class probabilities P(x, y)\n")
  print(round(x$joint, digits))
  cat("\nConditional probabilities P(y | x)\n")
  print(round(x$conditional, digits))
  invisible(x)
}

# The inherited LTA summary speaks of occasions and transitions; this model has
# neither, so its summary is its own print.
#' @export
summary.joint_classes <- function(object, digits = 3, ...) {
  print.joint_classes(object, digits = digits, ...)
}

# The item response probabilities of each variable's classes, long: one row
# per item, category and class, with the observed share of the category beside
# it. Read off the occasion layout fit_joint_classes() builds: stacked, one
# block holds every item, `x`'s first; paired, block 1 holds `x`'s items and
# block 2 `y`'s, each padded to the longer list. In both, row k of an item's
# probabilities is class k of the variable the item measures, and the rows
# past that variable's own classes are padding. A two-category item keeps its
# second category only, the probability of a 1 for a 0/1 item.
.joint_item_probabilities <- function(fit) {
  mods <- fit$mm$models
  Ja <- fit$n_items_x; Js <- fit$n_items_y
  stacked <- isTRUE(fit$joint_stacked)
  J <- if (stacked) Ja + Js else max(Ja, Js)
  X <- as.matrix(fit$data)
  w <- fit$weights_vec
  if (length(w) != nrow(X)) w <- rep(1, nrow(X))
  names_of <- function(v) {
    nm <- fit[[paste0("item_names_", v)]]
    n  <- if (v == "x") Ja else Js
    if (length(nm) == n) nm else sprintf("%s_item%d", v, seq_len(n))
  }
  one <- function(v, i, name) {
    K   <- if (v == "x") fit$n_x else fit$n_y
    m   <- if (v == "x" || stacked) mods[[1]] else mods[[2]]
    pos <- if (v == "y" && stacked) Ja + i else i
    col <- if (v == "x") i else if (stacked) J + Ja + i else J + i
    P   <- m$parameters$pis
    if (is.null(m$cats)) {
      pr <- cbind(1 - P[seq_len(K), pos], P[seq_len(K), pos]); val <- 0:1
    } else {
      M  <- as.integer(m$cats[pos])
      pr <- P[seq_len(K), m$offsets[pos] + seq_len(M), drop = FALSE]
      val <- seq_len(M)
    }
    x  <- X[, col]; ok <- !is.na(x)
    overall <- vapply(val, function(u) sum(w[ok] * (x[ok] == u)) / sum(w[ok]),
                      numeric(1))
    show <- if (length(val) == 2L) 2L else seq_along(val)
    data.frame(block = v, parameter = "probability", item = name,
               category = rep(show, each = K), class = rep(seq_len(K), length(show)),
               estimate = as.vector(pr[, show, drop = FALSE]),
               overall = rep(overall[show], each = K),
               label = rep(if (length(val) == 2L) name
                           else sprintf("%s (Cat %d)", name, show), each = K),
               stringsAsFactors = FALSE)
  }
  rows <- lapply(c("x", "y"), function(v) {
    nm <- names_of(v)
    do.call(rbind, lapply(seq_along(nm), function(i) one(v, i, nm[i])))
  })
  do.call(rbind, rows)
}

#' @rdname measurement_summary
#' @export
measurement_summary.joint_classes <- function(object, ...) {
  tab <- .joint_item_probabilities(object)
  cat("=========================================================\n")
  cat("             MEASUREMENT MODEL PARAMETERS                \n")
  cat("=========================================================\n")
  for (v in c("x", "y")) {
    d <- tab[tab$block == v, , drop = FALSE]
    K <- if (v == "x") object$n_x else object$n_y
    cat(sprintf("\n%s: item response probabilities by %s class\n", v, v))
    lab <- unique(d$label)
    w <- max(20L, nchar(lab))
    cat(sprintf("%-*s | Overall", w, "Indicator"))
    for (k in seq_len(K)) cat(sprintf(" | %7s", sprintf("%s#%d", toupper(v), k)))
    cat("\n", paste(rep("-", w + (K + 1L) * 10L), collapse = ""), "\n", sep = "")
    for (l in lab) {
      r <- d[d$label == l, , drop = FALSE]
      cat(sprintf("%-*s | %7.3f", w, l, r$overall[1]))
      for (k in seq_len(K)) cat(sprintf(" | %7.3f", r$estimate[r$class == k]))
      cat("\n")
    }
  }
  cat("\nA two-category item shows its second category: for an item coded 0/1,",
      "the probability of a 1.\n")
  .print_boundary_note(tab)
  cat("=========================================================\n")
  tab$label <- NULL
  invisible(tab)
}
