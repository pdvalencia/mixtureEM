# ==============================================================================
# Latent Transition Analysis - user-facing fitting function
# ==============================================================================

#' Latent Transition Analysis
#'
#' @description
#' Fits a latent transition model: each person occupies a latent *status* at each
#' occasion, measured by the same indicators every time, and may move between
#' statuses from one occasion to the next. Where [`fit_rmlca()`] assigns a person
#' one class for the whole study, `fit_lta()` estimates both the prevalence of
#' each status and the *incidence of change* between them (Collins & Lanza, 2010,
#' ch. 7).
#'
#' Three things are estimated:
#' \describe{
#'   \item{status prevalences}{how common each status is at the first occasion;}
#'   \item{transition probabilities}{the chance of moving from each status to
#'     each other status, as a square table read from row (earlier occasion) to
#'     column (later occasion). There is one table per pair of adjacent
#'     occasions unless `transition_invariance = "full"`;}
#'   \item{item parameters}{what people in each status tend to answer, which is
#'     what gives the statuses their meaning.}
#' }
#'
#' Cases with individual items or whole occasions missing are kept in the
#' analysis rather than dropped.
#'
#' **Why measurement invariance defaults to `"full"`.** If a status is not
#' defined identically at every occasion, an apparent "transition" mixes real
#' change in the person with a change in what the status means, and the two
#' cannot be told apart. Holding the item parameters equal across occasions
#' removes that ambiguity, and is the usual practice (Collins & Lanza,
#' sec. 7.11). It is a testable restriction: fit the model both ways and compare
#' them with [`longitudinal_lrt()`].
#'
#' @param indicators The repeated indicators, in any format accepted by
#'   [`fit_rmlca()`]: a wide matrix, a three-dimensional array, or a long data
#'   frame with `id` and `time`.
#' @param n_statuses Integer. Number of latent statuses.
#' @param times Integer. Number of occasions. Required for wide input.
#' @param measurement Measurement model for one occasion's items: `"binary"`,
#'   `"categorical"`, `"continuous"`, or a named list for a mixed block.
#' @param measurement_invariance Whether the item parameters are held equal
#'   across occasions: `"full"` (the default), `"none"`, or `"partial"` for only
#'   the items named in `invariant_items`. See the note above on why `"full"` is
#'   the sensible starting point.
#' @param invariant_items Items held equal across occasions when
#'   `measurement_invariance = "partial"`, given by name or position.
#' @param transition_invariance Whether the transition probabilities are held
#'   equal across occasions. `"none"` (the default) estimates a separate matrix
#'   for each pair of adjacent occasions; `"full"` shares one matrix throughout.
#'   Whether change happens at a constant rate is usually a substantive question
#'   rather than an assumption, and the two models are nested, so
#'   [`longitudinal_lrt()`] tests it (Collins & Lanza, sec. 7.14).
#' @param forbidden_transitions Transitions that are impossible by design, as in
#'   a stage-sequential process where people cannot move backwards. Give a
#'   \eqn{K \times K} logical or 0/1 matrix, with `TRUE`/`1` marking a forbidden
#'   move, or a list of one matrix per pair of adjacent occasions. Forbidden
#'   cells are fixed at zero and do not count as estimated parameters
#'   (Collins & Lanza, sec. 7.10).
#' @param n_classes Number of latent classes *above* the chain. The default of
#'   1 is the ordinary latent transition model. With more, each class gets its
#'   own initial distribution and its own transition matrices while sharing the
#'   measurement model - a mixture latent Markov model, which asks whether the
#'   population contains several distinct processes rather than one.
#'
#'   Be warned that the unrestricted mixture is demanding of the data. With one
#'   indicator per occasion it is routinely not identified: on a well-known
#'   five-wave life-satisfaction panel, even hundreds of random starts return
#'   a likelihood *worse* than the restricted mover-stayer model nested
#'   inside it. Several indicators per occasion, or a restriction such as
#'   `mover_stayer`, is usually what makes the model findable. Classes that
#'   converge on the same chain are warned about.
#' @param mover_stayer Restrict the **last** latent class to the identity
#'   transition matrix: a "stayer" class with zero probability of change, the
#'   remaining classes being "movers" (Vermunt, 2004). Implies `n_classes = 2`
#'   unless more are asked for,
#'   in which case only the last class is a stayer. The restricted rows cost no
#'   parameters, so the model is nested in the unrestricted mixture and
#'   [`longitudinal_lrt()`] tests it.
#' @param layout,id,time,items,item_names,time_labels Data-shape arguments,
#'   passed through as in [`fit_rmlca()`].
#' @param weights Optional case weights.
#' @param weight_type What the numbers in `weights` mean. `"sampling"` (the
#'   default) treats them as survey or probability weights and rescales them to
#'   sum to the number of cases; `"frequency"` treats them as counts of identical
#'   cases, as in a response-pattern table, and takes their sum as the sample
#'   size behind AIC and BIC.
#' @param strata,cluster Optional complex-survey design variables. When either is
#'   supplied, the standard errors become design-based: the case-level scores are
#'   aggregated to the primary sampling unit within stratum and used in a
#'   linearization sandwich, which protects inference against clustering.
#' @param n_init Number of random starts. Latent transition models have many
#'   local maxima; the default of 20 is a floor, not a recommendation.
#' @param max_iter,tol EM iteration limit and relative convergence tolerance.
#'   A mixture over chains (`n_classes` > 1) converges much more slowly and
#'   defaults to a tighter `tol` of 1e-11 and 5000 iterations, since the shared
#'   rule is a relative one and would otherwise stop the fit mid-climb. The
#'   restarts are then staged - a short first pass ranks them and only the best
#'   three run on to convergence - so the tighter rule does not multiply the
#'   cost of the search. Supplying either argument overrides all of this.
#' @param smoothing How much smoothing to apply to the status prevalences and to
#'   each row of the transition matrices, expressed as a number of pseudo-cases
#'   spread evenly over the possible destinations. Sparse transition tables
#'   otherwise collapse onto probabilities of exactly zero, which are awkward to
#'   interpret and to test. The default of `1` is negligible at any realistic
#'   sample size; set it to `0` for unsmoothed maximum likelihood.
#' @param random_state Optional seed for reproducible starts.
#' @param order_by_size Relabel the statuses from most to least prevalent at the
#'   first occasion. Ignored whenever the labels already carry meaning, which is
#'   the case when `forbidden_transitions` defines an ordering of stages or when
#'   covariates index the statuses through their regression coefficients.
#' @param standard_errors Compute standard errors for \eqn{\delta} and
#'   \eqn{\tau}.
#' @param predictors_initial Optional covariates predicting the latent status at
#'   the first occasion (Collins & Lanza, sec. 8.10.1).
#' @param predictors_transition Optional covariates predicting the transitions
#'   between statuses (sec. 8.10.2).
#' @param transition_effects How covariates act on the transitions.
#'   `"common"` (default) gives each origin status its own intercepts but one
#'   slope per covariate shared across origins, which is the specification in
#'   Wang & Wang eq. 6.28 and by far the better-behaved one. `"by_origin"` fits
#'   a separate regression per origin status, which is saturated and often fails
#'   to converge usefully when transitions are sparse.
#' @param group Optional grouping variable for a multiple-group model
#'   (Collins & Lanza, sec. 8.2-8.3). It is entered as dummy predictors,
#'   saturated over the transition rows, which gives each group its own status
#'   prevalences and its own transition matrices while the measurement model
#'   stays invariant across groups.
#' @param group_effects Which parameters the grouping variable is allowed to
#'   shift: `"both"` (default), `"initial"`, `"transitions"` or `"none"`.
#'   Fitting the same data under two of these and comparing them with
#'   [`longitudinal_lrt()`] gives the group-difference tests of sec. 8.6-8.8.
#' @param ... Ignored.
#'
#' @return An object of class `"lta_model"` with components including `delta`,
#'   `tau` (a list of transition matrices), `prevalences` (status prevalence by
#'   occasion), `gamma` (posterior status probabilities by occasion), `mm`,
#'   `metrics` and `n_params`.
#'
#'   With `n_classes` > 1 those parameters gain a class index: `delta` becomes a
#'   classes-by-statuses matrix, `tau` a list of per-class lists, and
#'   `class_weights`, `class_posterior` and `gamma_by_class` are added.
#'   `prevalences` stays the whole-sample marginal, with the per-class ones in
#'   its `"by_class"` attribute; [`status_prevalences()`] and
#'   [`transition_matrix()`] take a `class` argument to reach them. Standard
#'   errors are not available for a mixture over chains and `se` is `NULL`.
#'
#' @seealso [`transition_matrix()`], [`status_prevalences()`],
#'   [`longitudinal_lrt()`], [`lta_g2()`], [`fit_rmlca()`].
#' @export
fit_lta <- function(indicators,
                    n_statuses = 2,
                    times = NULL,
                    measurement = "binary",
                    measurement_invariance = c("full", "none", "partial"),
                    invariant_items = NULL,
                    transition_invariance = c("none", "full"),
                    forbidden_transitions = NULL,
                    n_classes = 1,
                    mover_stayer = FALSE,
                    layout = c("time_major", "item_major"),
                    id = NULL, time = NULL, items = NULL,
                    item_names = NULL, time_labels = NULL,
                    weights = NULL,
                    weight_type = c("sampling", "frequency"),
                    strata = NULL,
                    cluster = NULL,
                    n_init = 20,
                    max_iter = 1000,
                    tol = 1e-8,
                    smoothing = 1.0,
                    random_state = NULL,
                    order_by_size = TRUE,
                    standard_errors = TRUE,
                    predictors_initial = NULL,
                    predictors_transition = NULL,
                    transition_effects = c("common", "by_origin"),
                    group = NULL,
                    group_effects = c("both", "initial", "transitions", "none"),
                    ...) {

  measurement_invariance <- match.arg(measurement_invariance)
  transition_invariance  <- match.arg(transition_invariance)
  layout                 <- match.arg(layout)
  transition_effects     <- match.arg(transition_effects)
  group_effects          <- match.arg(group_effects)
  weight_type            <- match.arg(weight_type)

  time_invariance <- measurement_invariance
  tau_homogeneous <- transition_invariance == "full"
  tau_zeros       <- forbidden_transitions
  alpha           <- smoothing

  prep <- .prepare_longitudinal(indicators, times = times, items = items,
                                layout = layout, id = id, time = time,
                                item_names = item_names,
                                time_labels = time_labels)
  if (prep$n_times < 2L)
    stop("Latent transition analysis needs at least two occasions.",
         call. = FALSE)

  # --- Cases with no observed indicator at any occasion -----------------------
  # A case observed at no occasion has a flat likelihood at every wave, so its
  # posterior is the status prior and it contributes nothing to the fit, while
  # still inflating n in BIC/SABIC and lowering entropy. Deleted here, before the
  # measurement engine and every row-aligned design matrix are built, so that
  # each of them describes the analysed cases. See .empty_rows().
  empty_rows   <- .empty_rows(prep$X)
  n_input_rows <- nrow(prep$X)

  if (length(empty_rows) > 0L) {
    if (length(empty_rows) == n_input_rows)
      stop("Every case is missing on all indicators at every occasion, so ",
           "there are no data to fit.", call. = FALSE)

    # Length-checked before subsetting so a mis-specified argument still raises
    # its own error rather than being silently truncated.
    if (!is.null(weights) && length(weights) != n_input_rows)
      stop("`weights` must have one entry per case.", call. = FALSE)
    if (!is.null(strata) && length(strata) != n_input_rows)
      stop("`strata` must have one entry per case.", call. = FALSE)
    if (!is.null(cluster) && length(cluster) != n_input_rows)
      stop("`cluster` must have one entry per case.", call. = FALSE)

    keep              <- setdiff(seq_len(n_input_rows), empty_rows)
    prep$X            <- prep$X[keep, , drop = FALSE]
    prep$wave_missing <- prep$wave_missing[keep, , drop = FALSE]
    prep$any_missing  <- anyNA(prep$X)

    weights               <- .subset_cases(weights, keep)
    strata                <- .subset_cases(strata, keep)
    cluster               <- .subset_cases(cluster, keep)
    predictors_initial    <- .subset_cases(predictors_initial, keep)
    predictors_transition <- .subset_cases(predictors_transition, keep)
    group                 <- .subset_cases(group, keep)

    warning(sprintf(
      paste0("%d case%s had no observed value on any indicator at any ",
             "occasion and %s removed before estimation (n = %d analysed). ",
             "Rows: %s."),
      length(empty_rows), if (length(empty_rows) == 1L) "" else "s",
      if (length(empty_rows) == 1L) "was" else "were", length(keep),
      .abbreviate_indices(empty_rows)), call. = FALSE)
  }

  spec   <- .resolve_invariance(time_invariance, invariant_items,
                               prep$item_names, measurement)
  engine <- .longitudinal_measurement_spec(measurement, prep$X, prep$n_items,
                                           prep$n_times)

  X <- prep$X
  n <- nrow(X)

  if (weight_type == "sampling" && .looks_like_frequencies(weights, n))
    message("These weights look like frequency counts (whole numbers summing ",
            "to ", format(sum(weights)), " across ", n, " rows). ",
            "They are being treated as sampling weights; use ",
            "weight_type = \"frequency\" if each row stands for that many cases.")

  wt    <- .resolve_weights(weights, n, weight_type)
  w     <- wt$weights
  n_eff <- wt$n_eff

  if (!is.null(strata) && length(strata) != n)
    stop("`strata` must have one entry per case.", call. = FALSE)
  if (!is.null(cluster) && length(cluster) != n)
    stop("`cluster` must have one entry per case.", call. = FALSE)
  has_design <- !is.null(strata) || !is.null(cluster)

  K  <- as.integer(n_statuses)
  Tn <- prep$n_times

  # --- latent classes above the chain -----------------------------------------
  # `mover_stayer` is a restriction of the mixture, not a model of its own, so
  # it implies at least two classes and simply pins the last one's transitions.
  if (isTRUE(mover_stayer) && n_classes < 2) n_classes <- 2
  C <- as.integer(n_classes)
  if (is.na(C) || C < 1L)
    stop("`n_classes` must be a positive whole number.", call. = FALSE)
  if (C > 1L && Tn < 3L)
    stop("A mixture latent Markov model needs at least three occasions to be ",
         "identified (Vermunt, Mover-Stayer Models); with two, the classes ",
         "cannot be told apart from the transitions themselves.", call. = FALSE)

  # --- covariate and grouping design matrices ---------------------------------
  Z_delta <- .lta_design(predictors_initial,    n, "predictors_initial")
  Z_tau   <- .lta_design(predictors_transition, n, "predictors_transition")
  group_info <- NULL
  if (!is.null(group)) {
    group_info <- .lta_group_design(group, n)
    gd <- group_info$design
    if (group_effects %in% c("both", "initial"))
      Z_delta <- if (is.null(Z_delta)) cbind(Intercept = 1, gd) else
        cbind(Z_delta, gd)
    if (group_effects %in% c("both", "transitions")) {
      Z_tau <- if (is.null(Z_tau)) cbind(Intercept = 1, gd) else cbind(Z_tau, gd)
      # A multiple-group model gives every group its own transition matrix, so
      # the group must enter every transition row separately rather than as one
      # shared shift (Collins & Lanza, sec. 8.3 and the equivalence in 8.14).
      transition_effects <- "by_origin"
    }
  }
  if (!is.null(Z_tau) && !is.null(forbidden_transitions))
    stop("`forbidden_transitions` cannot be combined with covariates on the ",
         "transitions: a logistic regression has no way to hold a probability ",
         "at exactly zero. Model the transitions with covariates, or forbid ",
         "moves, but not both.", call. = FALSE)
  if (C > 1L && (!is.null(Z_delta) || !is.null(Z_tau)))
    stop("Covariates on the initial status or the transitions are not yet ",
         "available for a mixture latent Markov model (`n_classes` > 1 or ",
         "`mover_stayer = TRUE`). Fit the mixture without them, or use one ",
         "class.", call. = FALSE)

  allowed <- .lta_tau_allowed(tau_zeros, K, Tn, C, mover_stayer)

  state <- list(
    n_statuses      = K,
    n_classes       = C,
    mover_stayer    = isTRUE(mover_stayer),
    n_items         = prep$n_items,
    n_times         = Tn,
    class_weights   = rep(1 / C, C),
    delta_c         = rep(list(rep(1 / K, K)), C),
    tau_c           = rep(list(rep(list(matrix(1 / K, K, K)), Tn - 1L)), C),
    tau_allowed_c   = allowed,
    tau_homogeneous = isTRUE(tau_homogeneous),
    weights_vec     = w,
    weight_type       = weight_type,
    n_eff             = n_eff,
    strata            = if (is.null(strata))  rep(1L, n) else strata,
    cluster           = if (is.null(cluster)) seq_len(n) else cluster,
    has_survey_design = has_design,
    Z_delta         = Z_delta,
    Z_tau           = Z_tau,
    transition_effects = transition_effects,
    group_info      = group_info,
    group_effects   = if (is.null(group)) NULL else group_effects,
    mm              = time_blocks_model(K, prep$n_items, Tn,
                                        sub_model       = engine$sub_model,
                                        invariant_items = spec$invariant_items,
                                        max_val         = engine$max_val)
  )

  # A mixture over chains converges much more slowly than a single chain, and
  # the shared stopping rule is a *relative* one, so on a log-likelihood of
  # -15000 it fires at a change of 1.5e-4 and leaves the fit mid-climb: on the
  # five-wave satisfaction panel that costs 0.009 of log-likelihood and 400 of
  # the 800 iterations the model actually needs. The default is tightened for
  # the mixture rather than for everything, and only when the user has not asked
  # for a particular tolerance - the same reasoning that gave the growth
  # emissions their own `em_max_iter` instead of raising the package default.
  #
  # Running *every* restart to that tolerance is what it cannot afford: on
  # ex8.15 that is 20 restarts x 5000 iterations, and most of them are climbing
  # a hill they will lose anyway. The search is therefore staged, as the growth
  # mixture models' `em_stage1` stages theirs - a short first pass ranks the
  # restarts, and only the survivors are run on to convergence, resuming from
  # where they stopped. Three survivors rather than one because a short first
  # pass is a noisy ranking: the winning basin can be the slow one.
  staged <- C > 1L
  if (staged) {
    if (missing(tol))      tol      <- 1e-11
    if (missing(max_iter)) max_iter <- 5000
  }
  n_survivors <- if (staged) min(3L, max(1L, n_init)) else 0L

  best <- NULL
  stage1 <- list()
  for (init in seq_len(max(1L, n_init))) {
    if (!is.null(random_state)) set.seed(random_state + init)
    cand <- try(.lta_em(.lta_random_start(state, X), X,
                        max_iter = if (staged) min(250L, max_iter) else max_iter,
                        tol = if (staged) 1e-7 else tol, alpha = alpha),
                silent = TRUE)
    if (inherits(cand, "try-error")) next
    if (staged) stage1[[length(stage1) + 1L]] <- cand
    else if (is.null(best) || cand$loglik > best$loglik) best <- cand
  }

  if (staged && length(stage1)) {
    ord <- order(vapply(stage1, `[[`, numeric(1), "loglik"), decreasing = TRUE)
    for (i in utils::head(ord, n_survivors)) {
      cand <- try(.lta_em(stage1[[i]], X, max_iter = max_iter, tol = tol,
                          alpha = alpha), silent = TRUE)
      if (inherits(cand, "try-error")) cand <- stage1[[i]]
      if (is.null(best) || cand$loglik > best$loglik) best <- cand
    }
  }
  if (is.null(best))
    stop("Every random start failed; check the data and the model settings.",
         call. = FALSE)

  best$.tau_design_cache <- NULL      # working memory, not part of the fit
  best$n_params <- .lta_n_parameters(best)
  best$data     <- X
  best$longitudinal <- list(
    model           = "lta",
    n_items         = prep$n_items,
    n_times         = Tn,
    item_names      = prep$item_names,
    time_labels     = prep$time_labels,
    time_invariance = time_invariance,
    invariant_items = spec$invariant_items,
    wave_missing    = prep$wave_missing,
    measurement     = measurement
  )
  # Cases deleted for having no observed indicator at any occasion, recorded so
  # the analysed sample size can be reconciled with the input data.
  best$missing_data <- list(
    n_empty_rows = length(empty_rows),
    empty_rows   = empty_rows,
    n_input_rows = n_input_rows
  )
  class(best) <- "lta_model"

  # Reordering by prevalence is only safe when the status labels are arbitrary.
  # They are not when the user has forbidden particular moves, since that
  # declares an ordering of stages, and not when covariates are present, where
  # the statuses are also indexed by the regression coefficients and by the
  # origin dummies in the transition design.
  keep_labels <- !is.null(forbidden_transitions) ||
    !is.null(Z_delta) || !is.null(Z_tau)
  if (order_by_size && !keep_labels)
    best <- .sort_lta_statuses(best)

  best$prevalences <- .lta_prevalences(best)
  best$boundary    <- .lta_boundary_cells(best)
  best$metrics     <- .lta_metrics(best)
  if (isTRUE(standard_errors)) best$se <- .lta_standard_errors(best, X)

  # Recorded on the object as well as warned about: a warning is transient, and
  # someone reading a saved fit months later should still be able to see it.
  collapsed <- .lta_collapsed_classes(best)
  best$collapsed_classes <- collapsed
  if (!is.null(collapsed))
    warning(sprintf(paste0(
      "Latent class%s %s have converged on the same chain, so the model ",
      "describes one process with the parameters of several. Refit with ",
      "fewer classes, or with a restriction such as `mover_stayer = TRUE` ",
      "that tells the classes apart."),
      if (nrow(collapsed) == 1L) "es" else "es",
      paste(apply(collapsed, 1, paste, collapse = " and "), collapse = "; ")),
      call. = FALSE)

  best
}

# ------------------------------------------------------------------------------
# Constraints and parameter counting
# ------------------------------------------------------------------------------

# Expand the user's forbidden-transition specification into one logical
# "allowed" matrix per pair of adjacent occasions, for each latent class. TRUE
# means the move is free.
#
# `mover_stayer` is expressed here rather than anywhere else in the package: a
# stayer is a class whose only admissible move is to stay put, so its mask is
# the identity, and .lta_normalise() over a single admissible cell returns
# exactly 1 without any special case in the M-step. The stayer is the *last*
# class, matching the package's convention that the
# last category is the reference one.
.lta_tau_allowed <- function(forbidden, K, Tn, n_classes = 1L,
                             mover_stayer = FALSE) {
  one <- .lta_tau_allowed_1(forbidden, K, Tn)
  out <- rep(list(one), n_classes)
  if (isTRUE(mover_stayer) && n_classes >= 1L)
    out[[n_classes]] <- rep(list(diag(TRUE, K)), length(one))
  out
}

.lta_tau_allowed_1 <- function(forbidden, K, Tn) {
  n_mat <- max(Tn - 1L, 0L)
  if (is.null(forbidden))
    return(rep(list(matrix(TRUE, K, K)), n_mat))

  as_allowed <- function(m) {
    m <- as.matrix(m)
    if (!identical(dim(m), c(K, K)))
      stop(sprintf(paste("`forbidden_transitions` must be a %d x %d matrix,",
                         "one row and column per latent status."), K, K),
           call. = FALSE)
    allowed <- !(m == 1 | m == TRUE)
    if (any(rowSums(allowed) == 0L))
      stop("`forbidden_transitions` rules out every move from at least one ",
           "status, including staying put. Each row must leave at least one ",
           "destination possible.", call. = FALSE)
    allowed
  }

  if (is.list(forbidden)) {
    if (length(forbidden) != n_mat)
      stop(sprintf(paste("`forbidden_transitions` must be a single matrix or a",
                         "list of %d, one per pair of adjacent occasions."),
                   n_mat), call. = FALSE)
    return(lapply(forbidden, as_allowed))
  }
  rep(list(as_allowed(forbidden)), n_mat)
}

# Free parameters (Collins & Lanza, sec. 7.6):
#   class weights              C - 1
#   latent status prevalences  K - 1, once per class
#   transition probabilities   one (K-1)-vector per row per estimated matrix,
#                              less any cell fixed to zero, once per class - so
#                              a stayer class, whose rows admit one destination
#                              apiece, contributes nothing
#   item-response parameters   counted by the measurement model, which already
#                              accounts for across-time equality constraints
.lta_n_parameters <- function(state) {
  K  <- state$n_statuses
  Tn <- state$n_times
  C  <- state$n_classes %||% 1L

  n_delta <- if (!is.null(state$delta_beta))
    (K - 1L) * ncol(state$delta_beta) else (K - 1L) * C

  n_tau <- 0L
  if (Tn > 1L) {
    if (!is.null(state$tau_beta)) {
      n_tau <- state$tau_n_params
    } else {
      idx <- if (isTRUE(state$tau_homogeneous)) 1L else seq_len(Tn - 1L)
      for (c in seq_len(C)) for (i in idx)
        n_tau <- n_tau +
          sum(pmax(rowSums(state$tau_allowed_c[[c]][[i]]) - 1L, 0L))
    }
  }

  (C - 1L) + n_delta + n_tau + n_parameters(state$mm)
}

# ------------------------------------------------------------------------------
# Derived quantities
# ------------------------------------------------------------------------------

# Marginal status prevalence at each occasion, obtained by propagating delta
# through the transition matrices. This is the model-implied version; compare
# with colMeans of the posteriors for the empirical one.
.lta_prevalences <- function(state) {
  Tn <- state$n_times
  one <- function(delta, tau) {
    P <- matrix(0, Tn, state$n_statuses)
    P[1, ] <- delta
    if (Tn > 1L) for (t in 2:Tn)
      P[t, ] <- as.vector(P[t - 1, ] %*% tau[[t - 1]])
    dimnames(P) <- list(state$longitudinal$time_labels,
                        paste0("Status ", seq_len(state$n_statuses)))
    P
  }
  C <- state$n_classes %||% 1L
  if (C == 1L) return(one(state$delta_c[[1]], state$tau_c[[1]]))

  # Per class, plus the mixture's own marginal - the quantity longitudinal
  # profile plots label "Overall", and the one that should be
  # compared with the observed proportions.
  by_class <- lapply(seq_len(C), function(c)
    one(state$delta_c[[c]], state$tau_c[[c]]))
  names(by_class) <- paste("Class", seq_len(C))
  overall <- Reduce(`+`, Map(function(P, w) P * w,
                             by_class, state$class_weights))
  structure(overall, by_class = by_class)
}

# The whole-sample prevalences without the per-class ones riding along, for
# printing and for arithmetic: a stray attribute under a table is noise, and
# comparisons against a plain matrix fail on it.
.lta_bare_prevalences <- function(P) {
  attr(P, "by_class") <- NULL
  P
}

# Transition cells that have collapsed onto the boundary. Logit-scale software
# reports these as a logit fixed at a large negative value; here the cell is
# simply zero, and it
# is worth flagging because it costs a degree of freedom the count above still
# charges for.
#
# The threshold is 1e-4 rather than something nearer machine zero because EM
# stops when the *likelihood* stops moving, not when a cell reaches zero, and by
# then the cell is merely negligible: a transition a reference run reports as
# 0.000 is 3.2e-06 here, which a 1e-6 threshold missed entirely. Anything below
# 1e-4 is also well under one expected case in any sample this model is
# identified on, so there is nothing there to distinguish from zero.
.lta_boundary_cells <- function(state, tol = 1e-4) {
  if (state$n_times < 2L) return(NULL)
  C <- state$n_classes %||% 1L
  hits <- list()
  for (c in seq_len(C)) for (t in seq_along(state$tau_c[[c]])) {
    idx <- which(state$tau_c[[c]][[t]] < tol & state$tau_allowed_c[[c]][[t]],
                 arr.ind = TRUE)
    if (nrow(idx))
      hits[[length(hits) + 1L]] <- data.frame(
        class = c, occasion = t, from = idx[, "row"], to = idx[, "col"])
  }
  if (!length(hits)) return(NULL)
  out <- do.call(rbind, hits)
  if (C == 1L) out$class <- NULL
  out
}

.lta_metrics <- function(state) {
  ll <- state$loglik
  p  <- state$n_params
  n  <- state$n_eff %||% sum(state$weights_vec)

  # Relative entropy over all occasions: a status assignment is made at every
  # occasion, so the normalising constant counts n * T classifications.
  Tn  <- state$n_times
  abs_ent <- sum(vapply(state$gamma, function(g)
    sum(state$weights_vec * (-g * log(g + 1e-15))), numeric(1)))
  ent <- relative_entropy(abs_ent, n * Tn, state$n_statuses)

  # With classes above the chain there are two latent variables and two
  # separations to report: one
  # classification per case for the class, one per case per occasion for the
  # status. They are not comparable and are not combined.
  class_ent <- NULL
  if ((state$n_classes %||% 1L) > 1L && !is.null(state$class_posterior)) {
    P <- state$class_posterior
    class_ent <- relative_entropy(
      sum(state$weights_vec * rowSums(-P * log(P + 1e-15))),
      n, state$n_classes)
  }

  list(ll = ll, n_params = p,
       aic   = -2 * ll + 2 * p,
       bic   = -2 * ll + log(n) * p,
       sabic = -2 * ll + log((n + 2) / 24) * p,
       entropy = ent,
       class_entropy = class_ent,
       n_eff = n)
}

# Two classes whose chains have converged on each other describe one process
# twice, and the extra parameters are then bought with nothing. This is the
# degenerate optimum a free mixture latent Markov model falls into when the
# data cannot support it - single-indicator designs land there even with
# hundreds of random starts - so it is worth saying out loud rather
# than leaving in the parameter count.
.lta_collapsed_classes <- function(state, tol = 1e-3) {
  C <- state$n_classes %||% 1L
  if (C < 2L) return(NULL)
  flat <- function(c) c(state$delta_c[[c]], unlist(state$tau_c[[c]]))
  pairs <- list()
  for (b in seq_len(C)) for (a in seq_len(b - 1L))
    if (max(abs(flat(a) - flat(b))) < tol)
      pairs[[length(pairs) + 1L]] <- c(a, b)
  if (!length(pairs)) return(NULL)
  do.call(rbind, pairs)
}

# Order statuses by decreasing Time 1 prevalence. A single permutation applies
# to delta, to the rows AND columns of every transition matrix, to the
# measurement parameters, and to the posteriors.
.sort_lta_statuses <- function(state) {
  K <- state$n_statuses
  C <- state$n_classes %||% 1L
  if (K <= 1L) return(state)

  # One permutation for every class: the measurement model is shared, so a
  # status means the same thing in each class and they must be relabelled
  # together. The ordering is by Time 1 prevalence in the mixture as a whole.
  weight <- if (C == 1L) 1 else state$class_weights
  delta1 <- Reduce(`+`, Map(function(d, w) d * w, state$delta_c, weight))
  ord <- order(delta1, decreasing = TRUE)
  if (identical(ord, seq_len(K))) return(state)

  perm_mats <- function(l) lapply(l, function(m) m[ord, ord, drop = FALSE])
  state$delta_c       <- lapply(state$delta_c, function(d) d[ord])
  state$tau_c         <- lapply(state$tau_c, perm_mats)
  state$tau_allowed_c <- lapply(state$tau_allowed_c, perm_mats)
  state$gamma <- lapply(state$gamma, function(g) g[, ord, drop = FALSE])
  state$xi    <- perm_mats(state$xi)
  if (C > 1L) {
    state$gamma_by_class <- lapply(state$gamma_by_class, function(gl)
      lapply(gl, function(g) g[, ord, drop = FALSE]))
    state$xi_by_class <- lapply(state$xi_by_class, perm_mats)
  }

  for (t in seq_along(state$mm$models)) {
    sub <- state$mm$models[[t]]
    for (nm in c("pis", "means", "covariances"))
      if (!is.null(sub$parameters[[nm]]))
        sub$parameters[[nm]] <- sub$parameters[[nm]][ord, , drop = FALSE]
    state$mm$models[[t]] <- sub
  }
  .lta_pack(state)
}

# ------------------------------------------------------------------------------
# Standard errors (outer product of gradients)
# ------------------------------------------------------------------------------
#
# By the Fisher identity the score of the observed-data log-likelihood equals the
# expectation of the complete-data score given the data, so the E-step quantities
# already computed give each case's score directly. Stacking them gives the
# empirical information, and inverting it the covariance matrix. Parameters are
# taken on their multinomial-logit scale (last category anchored), which is where
# the normal approximation behaves; standard errors for the probabilities
# themselves follow by the delta method.
.lta_standard_errors <- function(state, X) {
  K  <- state$n_statuses
  Tn <- state$n_times
  w  <- state$weights_vec
  n  <- nrow(X)
  if (K < 2L || Tn < 2L) return(NULL)
  if (!is.null(state$delta_beta) || !is.null(state$tau_beta)) return(NULL)
  # A mixture over chains adds a class-membership block and makes every other
  # score class-conditional; the Fisher identity still applies but the blocks
  # below are not the right ones. Declined rather than reported wrongly, as for
  # the covariate models on the line above.
  if ((state$n_classes %||% 1L) > 1L) return(NULL)

  logB <- .lta_emission_loglik(state$mm, X)
  fb <- .lta_forward_backward(logB, log(pmax(state$delta, 1e-300)),
                              lapply(state$tau, function(m)
                                log(pmax(m, 1e-300))), w,
                              keep_pairwise = TRUE)

  scores <- list(); blocks <- list(); pos <- 0L
  add_block <- function(s, probs, ref, name) {
    scores[[length(scores) + 1L]] <<- s
    blocks[[length(blocks) + 1L]] <<- list(
      cols = pos + seq_len(ncol(s)), probs = probs, ref = ref, name = name)
    pos <<- pos + ncol(s)
  }

  # delta: the score of the multinomial logit is gamma_1[i, k] - delta_k, with
  # the last status anchored.
  add_block(sweep(state$gamma[[1]][, seq_len(K - 1L), drop = FALSE], 2,
                  state$delta[seq_len(K - 1L)], "-"),
            state$delta, seq_len(K - 1L), "delta")

  # tau: xi_i[t, k, l] - gamma_t[i, k] * tau_t[k, l], anchored on the last
  # admissible destination in the row.
  idx <- if (isTRUE(state$tau_homogeneous)) 1L else seq_len(Tn - 1L)
  for (i_mat in idx) {
    ts <- if (isTRUE(state$tau_homogeneous)) seq_len(Tn - 1L) else i_mat
    for (k in seq_len(K)) {
      allowed <- which(state$tau_allowed[[i_mat]][k, ])
      if (length(allowed) < 2L) next
      free <- allowed[-length(allowed)]
      s <- matrix(0, n, length(free))
      for (tt in ts) {
        pk <- fb$pairwise[[tt]][[k]]          # n x K: P(S_t = k, S_{t+1} = .)
        s <- s + pk[, free, drop = FALSE] -
          outer(state$gamma[[tt]][, k], state$tau[[i_mat]][k, free])
      }
      add_block(s, state$tau[[i_mat]][k, ], free,
                sprintf("tau%s[from %d]",
                        if (isTRUE(state$tau_homogeneous)) "" else
                          paste0("(", i_mat, ")"), k))
    }
  }

  # Measurement parameters. Including them is what makes these standard errors
  # unconditional; when the emission family is not one of the two handled here
  # the measurement model is treated as known, which understates uncertainty,
  # and the flag below says so.
  J <- state$n_items
  conditional <- TRUE
  fam <- class(state$mm$models[[1]])[1]
  if (fam %in% c("bernoulli", "bernoulli_nan")) {
    conditional <- FALSE
    inv <- state$longitudinal$invariant_items
    for (j in seq_len(J)) {
      ts_groups <- if (j %in% inv) list(seq_len(Tn)) else
        lapply(seq_len(Tn), identity)
      for (grp in ts_groups) {
        s <- matrix(0, n, K)
        for (tt in grp) {
          xj  <- X[, .time_block_cols(tt, J)[j]]
          rho <- state$mm$models[[tt]]$parameters$pis[, j]
          obs <- !is.na(xj); xj[!obs] <- 0
          s <- s + state$gamma[[tt]] * (xj - matrix(rho, n, K, byrow = TRUE)) *
            obs
        }
        add_block(s, NULL, NULL, sprintf("rho[item %d]", j))
      }
    }
  } else if (fam %in% c("gaussian_diag", "gaussian_diag_nan",
                        "gaussian_unit", "gaussian_unit_nan")) {
    # Class means enter the information matrix; the residual variances of
    # gaussian_diag do not, which is harmless for the blocks reported here
    # because means and variances are orthogonal in a normal model.
    conditional <- fam %in% c("gaussian_diag", "gaussian_diag_nan")
    inv <- state$longitudinal$invariant_items
    for (j in seq_len(J)) {
      ts_groups <- if (j %in% inv) list(seq_len(Tn)) else
        lapply(seq_len(Tn), identity)
      for (grp in ts_groups) {
        s <- matrix(0, n, K)
        for (tt in grp) {
          xj  <- X[, .time_block_cols(tt, J)[j]]
          sub <- state$mm$models[[tt]]
          mu  <- sub$parameters$means[, j]
          v   <- if (!is.null(sub$parameters$covariances))
            sub$parameters$covariances[, j] else rep(1, K)
          obs <- !is.na(xj); xj[!obs] <- 0
          s <- s + state$gamma[[tt]] *
            sweep(matrix(xj, n, K) - matrix(mu, n, K, byrow = TRUE), 2, v, "/") *
            obs
        }
        add_block(s, NULL, NULL, sprintf("mu[item %d]", j))
      }
    }
  }

  S    <- do.call(cbind, scores)
  info <- t(S) %*% sweep(S, 1, w, "*")
  V    <- tryCatch(pinv(info), error = function(e) NULL)
  if (is.null(V)) return(NULL)

  # With a complex survey design the empirical information is replaced by the
  # linearization sandwich, aggregating the same case-level scores to the
  # primary sampling unit within stratum. This is the machinery already behind
  # the covariate model's design-based variance (compute_survey_B in utils.R).
  design_based <- FALSE
  if (isTRUE(state$has_survey_design)) {
    meat <- tryCatch(
      compute_survey_B(sweep(S, 1, w, "*"), state$strata, state$cluster),
      error = function(e) NULL)
    if (!is.null(meat)) {
      V <- V %*% meat %*% V
      design_based <- TRUE
    }
  }

  # Delta method onto the probability scale. For a multinomial logit with the
  # reference category anchored, d p_l / d eta_m = p_l (1{l = m} - p_m), so the
  # covariance of the whole probability vector - reference category included -
  # is G V_block G'.
  prob_se <- list()
  for (b in blocks) {
    if (is.null(b$probs)) next
    p  <- b$probs
    Kp <- length(p)
    G  <- matrix(0, Kp, length(b$ref))
    for (m in seq_along(b$ref))
      G[, m] <- p * ((seq_len(Kp) == b$ref[m]) - p[b$ref[m]])
    Vb <- V[b$cols, b$cols, drop = FALSE]
    prob_se[[b$name]] <- sqrt(pmax(diag(G %*% Vb %*% t(G)), 0))
  }

  list(vcov = V, blocks = blocks, prob_se = prob_se,
       conditional = conditional, design_based = design_based)
}
