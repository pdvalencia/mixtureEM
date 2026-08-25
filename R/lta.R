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
#' them with [`lr_test()`].
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
#'   [`lr_test()`] tests it (Collins & Lanza, sec. 7.14).
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
#'   [`lr_test()`] tests it.
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
#'   local maxima; the default of 20 is a floor, not a recommendation. The fit
#'   reports how many starts reached the solution it kept, and warns when that
#'   count is 1; the answer there is `n_init = 100`, and a maximum that still
#'   does not replicate at 100 starts points at the specification rather than at
#'   the search. See `vignette("estimation")`.
#' @param n_cores Positive integer. Number of processes to spread the random
#'   starts over. Default `1` (sequential). These are the slowest fits in the
#'   package and `n_init` is high by necessity, so this is where the argument
#'   earns the most. Starting values are drawn in this session before any
#'   fitting begins, so the fit is identical at every `n_cores`.
#' @param max_iter,tol EM iteration limit and relative convergence tolerance.
#'   A mixture over chains (`n_classes` > 1) converges much more slowly and
#'   defaults to a tighter `tol` of 1e-11 and 5000 iterations, since the shared
#'   rule is a relative one and would otherwise stop the fit mid-climb. The
#'   restarts are then staged - a short first pass ranks them and only the best
#'   three run on to convergence - so the tighter rule does not multiply the
#'   cost of the search. Supplying either argument overrides all of this. (Unlike [`fit_mixture()`],
#'   whose EM tolerance is fixed and not user-adjustable, `tol` here is a
#'   real, respected argument, because a chain mixture converges slowly
#'   enough that the fixed rule would not do.)
#' @param smoothing How much smoothing to apply to the status prevalences and to
#'   each row of the transition matrices, expressed as a number of pseudo-cases
#'   spread evenly over the possible destinations. Sparse transition tables
#'   otherwise collapse onto probabilities of exactly zero, which are awkward to
#'   interpret and to test. The default of `1` is negligible at any realistic
#'   sample size; set it to `0` for unsmoothed maximum likelihood. It governs
#'   the status prevalences, the transition matrices and - with `n_classes` > 1 -
#'   the class weights, and it does **not** govern the measurement model, whose
#'   prior is `bayes_constants`.
#'
#'   It also does not reach the transitions once anything predicts them. With
#'   `predictors_transition`, or with a `group` whose `group_effects` include
#'   the transitions, each transition matrix is the fitted value of a
#'   multinomial logit rather than a smoothed table of counts, and `smoothing`
#'   has no effect on it. The initial status prevalences behave the same way
#'   under `predictors_initial`.
#'
#'   The mass is one pseudo-case per *origin row*, which is the prior Chung,
#'   Lanza and Loken (2008) use for this model, and it is spread evenly rather
#'   than in proportion to how often each destination is occupied: a rare origin
#'   row shrunk toward the destination marginal would be asserting that everyone
#'   moves to the prevalent status, which is a confident claim to make about a
#'   row the sample says little about, whereas an even spread is uninformative.
#'
#'   The cost falls on exactly those rows. On a row with few expected cases the
#'   prior carries a visible share of the estimate - at most
#'   \eqn{[\alpha / (m + \alpha)](1 - 1/K_a)} of it, for a row with \eqn{m}
#'   expected cases and \eqn{K_a} reachable destinations - and the fit says so
#'   when that share exceeds five percentage points, naming the worst row. The
#'   remedy worth reaching for first is `transition_invariance = "full"`, which
#'   puts every occasion's cases behind one pseudo-case; `smoothing = 0.5`
#'   simply halves the pull. `smoothing = 0` is not a good answer, since it
#'   removes the protection against transition probabilities of exactly zero
#'   that the prior is there to give.
#' @param bayes_constants Optional named list of prior strengths for the
#'   *measurement* model (`categorical`, `poisson`, `variances`); see
#'   [`fit_mixture()`]. The status and transition probabilities are governed by
#'   `smoothing` instead, so `latent` is not read here.
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
#'   [`lr_test()`] gives the group-difference tests of sec. 8.6-8.8.
#' @param ... Ignored.
#'
#' @return An object of class `"lta_model"` with components including `delta`,
#'   `tau` (a list of transition matrices), `prevalences` (status prevalence by
#'   occasion), `gamma` (posterior status probabilities by occasion), `mm`,
#'   `metrics` and `n_params`.
#'
#'   Two diagnostics come with it. `boundary` lists the transition cells the data
#'   have driven to zero, and `smoothing_influence` gives, for each origin row,
#'   the cases expected in it and the share of the estimate the `smoothing`
#'   prior is carrying. `smoothing_influence` is `NULL` when `smoothing` is 0 or
#'   when covariates predict the transitions, since the prior does not reach
#'   them then.
#'
#'   With `n_classes` > 1 those parameters gain a class index: `delta` becomes a
#'   classes-by-statuses matrix, `tau` a list of per-class lists, and
#'   `class_weights`, `class_posterior` and `gamma_by_class` are added.
#'   `prevalences` stays the whole-sample marginal, with the per-class ones in
#'   its `"by_class"` attribute; [`status_prevalences()`] and
#'   [`transition_matrix()`] take a `class` argument to reach them. Standard
#'   errors are not available for a mixture over chains and `se` is `NULL`.
#'
#' @references
#' Collins, L. M., & Lanza, S. T. (2010). \emph{Latent Class and Latent
#' Transition Analysis: With Applications in the Social, Behavioral, and Health
#' Sciences}. Wiley (chapters 7-8).
#'
#' Nylund-Gibson, K., Grimm, R., Quirk, M., & Furlong, M. (2014). A latent
#' transition mixture model using the three-step specification.
#' \emph{Structural Equation Modeling}, \emph{21}(3), 439-454.
#' \doi{10.1080/10705511.2014.915375}
#'
#' Chung, H., Lanza, S. T., & Loken, E. (2008). Latent transition analysis:
#' inference and estimation. \emph{Statistics in Medicine}, \emph{27}(11),
#' 1834-1854. \doi{10.1002/sim.3130}
#'
#' Fienberg, S. E., & Holland, P. W. (1973). Simultaneous estimation of
#' multinomial cell probabilities. \emph{Journal of the American Statistical
#' Association}, \emph{68}(343), 683-691.
#' \doi{10.1080/01621459.1973.10481405}
#'
#' @seealso [`transition_matrix()`], [`status_prevalences()`],
#'   [`lr_test()`], [`lta_g2()`], [`fit_rmlca()`].
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
                    n_cores = 1L,
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
                    bayes_constants = NULL,
                    ...) {

  measurement_invariance <- match.arg(measurement_invariance)
  transition_invariance  <- match.arg(transition_invariance)
  layout                 <- match.arg(layout)
  transition_effects     <- match.arg(transition_effects)
  group_effects          <- match.arg(group_effects)
  weight_type            <- match.arg(weight_type)

  # `latent` is the one bayes_constants name fit_lta() does not read: the
  # status and transition priors are `smoothing`'s job. Resolving it silently
  # would let a four-way prior specification lose a quarter of itself without
  # a word, so the value that will be ignored is named where it is passed.
  if (is.list(bayes_constants) && "latent" %in% names(bayes_constants))
    warning(sprintf(
      paste0("`bayes_constants$latent` is not read by fit_lta(); the prior on ",
             "the initial-status and transition probabilities is the ",
             "`smoothing` argument. Did you mean smoothing = %s?"),
      format(bayes_constants[["latent"]])), call. = FALSE)

  bayes_constants        <- .resolve_bayes_constants(bayes_constants)

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
  prep$X <- engine$X

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
    bayes_constants = bayes_constants,
    mm              = time_blocks_model(K, prep$n_items, Tn,
                                        sub_model       = engine$sub_model,
                                        invariant_items = spec$invariant_items,
                                        max_val         = engine$max_val)
  )

  # The measurement M-steps read their prior strengths off the emission. LTA has
  # its own EM driver, so the constants are pushed down here rather than by
  # fit_mixture_internal(). `smoothing` continues to govern the status and
  # transition probabilities and is passed separately as `alpha`; a continuous
  # indicator's variance prior comes from `bayes_constants$variances`, which is
  # why the two are not the same knob.
  state$mm <- .attach_bayes_constants(state$mm, bayes_constants)

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
  # Every restart that ran to convergence, kept so the fit can report how many
  # of them found the reported maximum. On the staged path that is the second
  # loop: the first is a ranking pass stopped at 250 iterations, and its
  # log-likelihoods are not the maxima of anything.
  final_lls <- numeric(0)
  # The random starts are drawn here, in restart order, and fitted afterwards.
  # .lta_random_start() is the only RNG consumer on this path -- .lta_em() draws
  # nothing -- so this takes the same numbers the sequential loop took and makes
  # each restart a deterministic function of the start it is given. The fits can
  # then run on workers without moving a value, at any `n_cores`.
  starts <- lapply(seq_len(max(1L, n_init)), function(i) {
    if (!is.null(random_state)) set.seed(random_state + i)
    .lta_random_start(state, X)
  })
  cands <- .par_lapply(starts, function(s)
    try(.lta_em(s, X,
                max_iter = if (staged) min(250L, max_iter) else max_iter,
                tol = if (staged) 1e-7 else tol, alpha = alpha),
        silent = TRUE),
    n_cores = n_cores)

  for (cand in cands) {
    if (inherits(cand, "try-error")) next
    if (staged) stage1[[length(stage1) + 1L]] <- cand
    else {
      final_lls <- c(final_lls, cand$loglik)
      if (is.null(best) || cand$loglik > best$loglik) best <- cand
    }
  }

  if (staged && length(stage1)) {
    ord <- order(vapply(stage1, `[[`, numeric(1), "loglik"), decreasing = TRUE)
    survivors <- .par_lapply(utils::head(ord, n_survivors), function(i) {
      cand <- try(.lta_em(stage1[[i]], X, max_iter = max_iter, tol = tol,
                          alpha = alpha), silent = TRUE)
      if (inherits(cand, "try-error")) stage1[[i]] else cand
    }, n_cores = n_cores)
    for (cand in survivors) {
      final_lls <- c(final_lls, cand$loglik)
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
  best$smoothing_influence <- .lta_smoothing_influence(best, alpha)
  best$metrics     <- .lta_metrics(best)
  # The multi-start report the mixture models already carry. Two restarts count
  # as the same solution when their log-likelihoods are within 1e-2, the rule
  # fit_em() uses (R/em_core.R): genuinely different optima in these models sit
  # whole units apart, and a tighter rule splits one optimum into several.
  if (length(final_lls)) {
    best$metrics$n_starts     <- length(final_lls)
    best$metrics$n_replicated <- sum(abs(final_lls - max(final_lls)) <= 1e-2)
  }
  best$metrics$n_requested <- max(1L, n_init)
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
      if (nrow(collapsed) == 1L) "" else "es",
      paste(apply(collapsed, 1, paste, collapse = " and "), collapse = "; ")),
      call. = FALSE)

  # A small origin row, not a bad fit: the prior is doing what it is there to do
  # and the only question is how much of the estimate is left over for the data.
  # Classed so compare_longitudinal() can muffle it, as it does the replication
  # warning - a K-range sweep would raise it once per model and drown the table.
  .warn_smoothing_influence(best$smoothing_influence, smoothing)

  # Continuous indicators are checked here too: LTA runs its own EM driver, so
  # it does not pass through fit_mixture_internal() where this normally happens.
  # See R/gaussian_boundary.R.
  best <- .check_gaussian_degeneracy(best, X)

  # Both warnings come last, and in this order, for the same reason: the
  # degeneracy check above sets the flag that .check_replication() reads before
  # deciding whether raising n_init is good advice. Convergence was previously
  # reported only by print(), which a user working from the transition matrix
  # never sees.
  if (isFALSE(best$converged)) .warn_non_convergence(max_iter)
  .check_replication(best)

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

# How much of each transition row is prior rather than data. With `alpha`
# pseudo-cases spread over the Ka admissible destinations of a row whose
# expected count is m, the posterior mean is the shrinkage of Fienberg &
# Holland (1973, eq. 2.6) with weight alpha/(m + alpha), so the largest amount
# the prior can move any cell of that row is
#
#     pull = [alpha / (m + alpha)] * (1 - 1/Ka)
#
# which is exact, not an approximation. It is reported per row and the worst row
# is what the message names. Rows are the thing to look at rather than the
# matrix as a whole because the number of classes divides the sample among them
# - a row's expected count is n * pi_c * P(status k | c) - while the number of
# occasions does not.
#
# The counterpart of .lta_boundary_cells() above: that one reports cells that
# have collapsed *onto* the boundary, this one reports rows the prior has pulled
# *away* from it.
#
# Covariates on the transitions take the whole thing out of scope. That M-step
# is a multinomial logit fitted by .lta_mstep_tau_cov() (R/lta_covariates.R),
# which never calls .lta_normalise(), so `smoothing` does not reach the
# transition matrices at all on that path and there is no pull to report. A
# grouping variable enters the same way, as dummy predictors saturated over the
# transition rows, so `group_effects` of "both" or "transitions" is covered by
# the same test. Without this guard the row counts are real but the prior they
# are compared against is not applied, and the message would send the reader
# after a `smoothing` that is doing nothing.
.lta_smoothing_influence <- function(state, alpha) {
  if (state$n_times < 2L || !isTRUE(alpha > 0) || !is.null(state$Z_tau))
    return(NULL)
  C <- state$n_classes %||% 1L
  K <- state$n_statuses
  # Under `transition_invariance = "full"` the M-step pools the occasions before
  # it smooths, so one pseudo-case is spread over a row with several occasions'
  # counts in it. Pooling here too; per occasion the reported pull would be
  # several times the real one.
  pooled <- isTRUE(state$tau_homogeneous)

  rows <- list()
  for (c in seq_len(C)) {
    Xi      <- if (C > 1L) state$xi_by_class[[c]] else state$xi
    allowed <- state$tau_allowed_c[[c]]
    groups  <- if (pooled) list(seq_along(Xi)) else as.list(seq_along(Xi))
    for (g in groups) {
      counts <- Reduce(`+`, Xi[g])
      mask   <- allowed[[g[1]]]
      for (k in seq_len(K)) {
        a  <- mask[k, ]
        Ka <- sum(a)
        # A row with no admissible destination is not smoothed at all - the
        # normaliser returns a uniform vector - and there is nothing to report.
        if (Ka == 0L) next
        m <- sum(counts[k, a])
        rows[[length(rows) + 1L]] <- data.frame(
          class = c, occasion = if (pooled) NA_integer_ else g[1], from = k,
          n_expected = m, pull = (alpha / (m + alpha)) * (1 - 1 / Ka))
      }
    }
  }
  if (!length(rows)) return(NULL)
  out <- do.call(rbind, rows)
  if (C == 1L) out$class <- NULL
  out
}

# The threshold is a reporting-precision choice: above it the prior can move a
# transition probability by more than five percentage points, which is more than
# the precision such probabilities are reported to. It is a calibrated rule of
# thumb rather than a sharp boundary - it was checked against the risk of the
# prior on sparse rows, where it fires while alpha = 1 costs 8-20% against a
# smaller prior and goes quiet as that penalty falls to a few percent - so it is
# never presented as a verdict on the fit.
.lta_smoothing_tol <- 0.05

# The worst row, as a phrase both the warning and the printer can use.
.lta_worst_smoothing_row <- function(influence) {
  if (is.null(influence) || !nrow(influence)) return(NULL)
  w <- influence[which.max(influence$pull), ]
  where <- sprintf("transitions out of status %d", w$from)
  if (!is.na(w$occasion)) where <- sprintf("%s at occasion %d", where, w$occasion)
  if (!is.null(w$class))  where <- sprintf("%s in class %d", where, w$class)
  list(pull = w$pull, n_expected = w$n_expected, where = where)
}

.warn_smoothing_influence <- function(influence, smoothing) {
  w <- .lta_worst_smoothing_row(influence)
  if (is.null(w) || w$pull <= .lta_smoothing_tol) return(invisible(NULL))

  msg <- sprintf(paste0(
    "The transition prior is carrying more than %d%% of some rows: the largest ",
    "effect is on %s, where %.1f cases are expected and `smoothing = %s` can ",
    "move a transition probability by up to %.2f. Those transitions rest partly ",
    "on the prior rather than on the sample, so read them as indicative. ",
    "Reporting them at face value overstates how much the data say about that ",
    "row; `smoothing = %s` reduces the pull, and pooling the occasions with ",
    "`transition_invariance = \"full\"` removes it where the transitions can be ",
    "assumed constant over time."),
    round(100 * .lta_smoothing_tol), w$where, w$n_expected,
    format(smoothing), w$pull, format(smoothing / 2))

  warning(structure(class = c("mixtureEM_smoothing", "warning", "condition"),
                    list(message = msg, call = NULL)))
  invisible(NULL)
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
