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
#'     occasions unless `transition_invariance` restricts them;}
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
#'   `"categorical"`, `"ordinal"`, `"continuous"`, or a named list for a mixed
#'   block. `"ordinal"` fits one cumulative-logit block per item, with its own
#'   number of ordered categories inferred per item (so a 3/3/2-category block
#'   needs no mixed specification); without a random intercept it is
#'   numerically identical to `"categorical"`.
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
#'
#'   `"slopes"` sits between the two and applies only when something predicts
#'   the transitions. It shares the covariate slopes across occasions while
#'   giving each occasion its own intercepts, so the effect of a covariate on
#'   moving between statuses is held constant over time but the underlying rate
#'   of movement is not. It costs `(n_statuses - 1) * (occasions - 2)`
#'   parameters more than `"full"` and is nested inside `"none"`, so
#'   [`lr_test()`] tests both restrictions. It requires
#'   `predictors_transition` (or a `group` acting on the transitions) and
#'   `transition_effects = "common"`, and with only two occasions it is
#'   identical to `"full"`.
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
#' @param tie_initial_status For a mover-stayer fit, hold the occasion-1
#'   status distribution equal across the latent classes instead of
#'   estimating one per class. Drops the initial-status parameter count from
#'   `(K - 1) * C` to `K - 1`. Default `FALSE`.
#' @param random_intercept Add a random intercept to the measurement model
#'   (Muthen & Asparouhov, 2022): a person-level "how likely to endorse items
#'   in general" trait that regular LTA has no way to represent, and that can
#'   otherwise be mistaken for status separation and for stability over time.
#'   `"none"` (the default) fits ordinary LTA. `"continuous"` integrates over a
#'   normally-distributed factor by Gauss-Hermite quadrature (`n_quadrature`
#'   nodes); `"binary"` instead estimates a small number of discrete intercept
#'   classes (`n_ri` of them, 2 by default, the model's own case). Both require
#'   `measurement_invariance = "full"` and binary indicators. A random
#'   intercept combines with `mover_stayer`, with covariates on the initial
#'   status or the transitions (`predictors_initial`/`predictors_transition`),
#'   and with `group` (implemented as covariates on the same two, so this is
#'   one capability, not two); it still does not support `n_classes` > 1.
#'
#'   `"continuous"`'s restart search (`options(mixtureEM.lta_ri_search =
#'   "wide")`, the default) draws wide random starts -- perturbed item logits
#'   and loadings of either sign -- ranks them with a one-step generalised-EM
#'   M-step, and rescores every ranked candidate on the full quadrature grid
#'   before promoting survivors to the full search. On several published
#'   benchmark data sets it reaches the same or a better optimum, several
#'   times faster, than the search this package used before,
#'   and it escapes the inflated-loading local maximum that search could
#'   settle in. `options(mixtureEM.lta_ri_search =
#'   "narrow")` restores the earlier search exactly, for reproducing a fit
#'   made under it.
#'
#'   **Do not test a random intercept against regular LTA with [`lr_test()`]**:
#'   the continuous variant puts the null (loading = 0) on the boundary of the
#'   parameter space, and the binary variant adds a latent class variable, so
#'   the usual chi-squared reference distribution does not apply either way.
#'   Compare the two by BIC instead.
#'
#'   Whether a random intercept is worth adding is a sample-size question more
#'   than a modelling one. Tseng (2024), for the continuous-indicator analogue
#'   (RI-LPTA), puts the requirement at upwards of 2,000 cases for 80% power
#'   and 90% coverage at a between-profile separation of d = 0.75, with more
#'   items, occasions or separation lowering that bar; the asymmetry that
#'   makes trying it worthwhile anyway is that
#'   omitting a random intercept when one belongs costs a lot (inflated
#'   apparent separation and stability), while including one when it does not
#'   belong costs almost nothing (a handful of parameters, and BIC will say
#'   so, as it does on this package's own benchmark replication of the
#'   article's example).
#'
#'   Not built in this release: continuous indicators, a random *slope*,
#'   correlated residuals across time, or lag-2 dependence - the last of which
#'   the article itself reports
#'   as significant in both of its worked examples, so it is a real
#'   simplification and not a hypothetical one. Standard errors and the
#'   post-EM L-BFGS refinement are available for both binary and ordinal
#'   indicators. A random
#'   intercept crossed with several latent classes or `mover_stayer` is
#'   supported for both measurement families.
#'
#'   `predictors_random_intercept` regresses the factor on covariates. The
#'   item probabilities the fit reports are integrated over the factor at a
#'   zero random-intercept mean - the residual grid - not at each case's own
#'   predicted mean; the per-case means are in
#'   `random_intercept_scores()$predicted_mean`.
#' @param n_quadrature Number of Gauss-Hermite nodes for
#'   `random_intercept = "continuous"`. The default of 15 is a starting point
#'   to check, not a settled answer, the same way `n_init`'s default is a
#'   floor: raise it (15-30 nodes is the usual working range, more when the
#'   loadings are large) and confirm the log-likelihood moves by less than
#'   0.01. Make that check over a wide range
#'   of node counts rather than one step up. When a fit's thresholds are
#'   extreme, the item response is almost a step function of the factor, and
#'   the log-likelihood is then not even monotone in the number of nodes: two
#'   nearby small settings can differ by several log-likelihood units while
#'   both sit far from the converged value. One ordinal five-status fit used
#'   in this package's own validation reads -16047.3 at 15 nodes and -16040.5
#'   at 20, but settles at -16041.18 only from roughly 80 nodes upward.
#'   `n_quadrature = 1` is a valid,
#'   deliberate special case - a single node at 0 with mass 1 - under which
#'   the model reduces exactly to regular LTA; it is not a model worth fitting
#'   on its own, but is how the package's own test suite proves the node
#'   machinery is wired correctly. The node count matters more once
#'   `predictors_random_intercept` is used: the person-specific reweighting is
#'   exact only up to whatever quadrature error the fixed grid already has.
#' @param n_ri Number of discrete intercept classes for
#'   `random_intercept = "binary"`. The default of 2 is Muthen & Asparouhov's
#'   own case.
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
#' @param refine Logical. If `TRUE` (default), each start that runs to
#'   convergence is followed by an L-BFGS climb on the same penalised objective
#'   the EM steps maximise, and the starts are then ranked on the refined
#'   log-likelihoods. EM converges to a fixed point of its own surrogate, which
#'   on a near-flat likelihood ridge can sit measurably short of the maximum;
#'   the climb steps past it, and never returns a fit worse than the one it was
#'   given. This is the same refinement [`fit_mixture()`] has always applied.
#'   It is a no-op - silently, and the fit is unchanged - for models whose free
#'   parameters it cannot differentiate: covariate models and measurement
#'   families other than binary and continuous. Mixtures over chains
#'   (`n_classes` > 1) are refined like any other supported model.
#' @param refine_from A model fitted by [`fit_lta()`] on the same data and with
#'   the same shape, whose solution this fit continues from. No random starts
#'   are run and `n_init` is ignored, because the search that produced the donor
#'   already ran and this only carries its winner further - typically to a
#'   tighter `tol` or a larger `max_iter`. Passing `n_init` alongside is an
#'   error rather than a silent override.
#' @param n_cores Positive integer. Number of processes to spread the random
#'   starts over. Default `1` (sequential). These are the slowest fits in the
#'   package and `n_init` is high by necessity, so this is where the argument
#'   earns the most. Starting values are drawn in this session before any
#'   fitting begins, so the fit is identical at every `n_cores`.
#'   `options(mixtureEM.n_cores = )` sets the default for a whole session; an
#'   argument given here overrides it.
#' @param max_iter,tol EM iteration limit and relative convergence tolerance.
#'   A mixture over chains (`n_classes` > 1) converges much more slowly and
#'   defaults to a tighter `tol` of 1e-11 and 5000 iterations, since the shared
#'   rule is a relative one and would otherwise stop the fit mid-climb. The
#'   restarts are then staged - a short first pass ranks them and only the
#'   top tenth (floor of 3) run on to convergence - so the tighter rule does
#'   not multiply the cost of the search. Supplying either argument overrides all of this. (Unlike [`fit_mixture()`],
#'   whose EM tolerance is fixed and not user-adjustable, `tol` here is a
#'   real, respected argument, because a chain mixture converges slowly
#'   enough that the fixed rule would not do.)
#' @param smoothing How much smoothing to apply to the status prevalences and to
#'   the transition matrices, expressed as a number of pseudo-cases spread
#'   evenly over each conditional table. Sparse transition tables
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
#'   The mass is one pseudo-case per *conditional table*, so with `K` statuses
#'   an origin row of a transition matrix carries a `K`th of it. It is spread
#'   evenly rather than in proportion to how often each destination is
#'   occupied: a rare origin row shrunk toward the destination marginal would
#'   be asserting that everyone moves to the prevalent status, which is a
#'   confident claim to make about a row the sample says little about, whereas
#'   an even spread is uninformative.
#'
#'   The cost falls on exactly those rows. On a row with few expected cases the
#'   prior carries a visible share of the estimate - at most
#'   \eqn{[a / (m + a)](1 - 1/K_a)} of it, where \eqn{a} is that row's share of
#'   the mass, \eqn{m} the cases expected in it and \eqn{K_a} the reachable
#'   destinations - and the fit says so when that share exceeds five percentage
#'   points, naming the worst row. The remedy worth reaching for first is
#'   `transition_invariance = "full"`, which puts every occasion's cases behind
#'   one pseudo-case; `smoothing = 0.5` simply halves the pull. `smoothing = 0`
#'   is not a good answer, since it removes the protection against transition
#'   probabilities of exactly zero that the prior is there to give.
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
#'   \eqn{\tau}. `TRUE` (the default) uses the outer product of the case-level
#'   scores; `FALSE` skips them; `"robust"` returns the sandwich estimator, with
#'   the observed information as its bread. `"robust"` costs a
#'   finite-difference Hessian -- \eqn{2p(p + 1)} likelihood evaluations -- and
#'   so takes minutes rather than seconds on a model of any size. It is
#'   supported for random-intercept fits. Where the model still cannot be
#'   packed on an unconstrained scale (covariates, or a measurement family
#'   whose parameters are not all free) `"robust"` falls back silently to the
#'   default estimator, and `summary()` then does not report robust errors.
#' @param predictors_initial Optional covariates predicting the latent status at
#'   the first occasion (Collins & Lanza, sec. 8.10.1).
#' @param predictors_transition Optional covariates predicting the transitions
#'   between statuses (sec. 8.10.2).
#' @param predictors_items Optional covariates entering each indicator
#'   directly, inside each latent status - a test of measurement invariance
#'   with respect to those covariates. One proportional-odds slope per latent
#'   status per item per covariate, shared across occasions, so the cost is
#'   `n_statuses * n_items * ncol(predictors_items)` parameters. **This is not
#'   `group`**: `group` gives every group its own status prevalences and
#'   transition matrices while the measurement model stays invariant across
#'   groups, which is the assumption this argument exists to relax. Requires
#'   `measurement = "ordinal"` and `measurement_invariance = "full"`.
#'   Restricted to covariates taking few distinct values - see
#'   `options(mixtureEM.dif_max_patterns = )`. The item probabilities the fit
#'   reports are those of a case with every one of these covariates at zero.
#' @param predictors_random_intercept Optional covariates predicting the
#'   continuous random intercept itself (the article's own Step 5). The
#'   factor's residual variance is fixed at 1 and its residual mean at 0, so no
#'   intercept is estimated and a constant column is refused: the coefficients
#'   are the regression of the stable trait on the covariates, and each costs
#'   one parameter. Requires `random_intercept = "continuous"`. The factor's
#'   orientation is arbitrary - the fit pins it by making the largest loading
#'   positive - so the sign of every coefficient is meaningful only relative to
#'   the loadings.
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
#'   `metrics$entropy` (the headline number, printed by `print()`) is the
#'   relative entropy of the joint latent-status path across all occasions at
#'   once: one classification per case over the `K^T` possible paths,
#'   normalised by `n * T * log(K)`. The per-occasion alternative -- the
#'   relative entropy of each occasion's status posterior on its own,
#'   normalised by `log(K)` -- is `metrics$entropy_by_occasion` (printed by
#'   `summary()`). A per-occasion entropy can also be normalised by the
#'   entropy of that occasion's estimated status proportions instead of by
#'   `log(K)`; the two are related by `1 - (1 - r2) * log(K) / H(pi_hat)`,
#'   where `r2` is the package's value, `K` is `n_statuses`, and `H(pi_hat)`
#'   is the entropy of that occasion's estimated status proportions (from
#'   [`status_prevalences()`]).
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
#' Muthen, B., & Asparouhov, T. (2022). Latent transition analysis with random
#' intercepts (RI-LTA). \emph{Psychological Methods}, \emph{27}(1), 1-16.
#' \doi{10.1037/met0000370}
#'
#' Tseng, M.-C. (2024). Latent profile transition analysis with random
#' intercepts (RI-LPTA). \emph{Structural Equation Modeling}, \emph{31}(4),
#' 626-634. \doi{10.1080/10705511.2023.2284671} - the sample-size analysis
#' behind the guidance above.
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
                    transition_invariance = c("none", "full", "slopes"),
                    forbidden_transitions = NULL,
                    n_classes = 1,
                    mover_stayer = FALSE,
                    tie_initial_status = FALSE,
                    random_intercept = c("none", "continuous", "binary"),
                    n_quadrature = 15,
                    n_ri = 2,
                    layout = c("time_major", "item_major"),
                    id = NULL, time = NULL, items = NULL,
                    item_names = NULL, time_labels = NULL,
                    weights = NULL,
                    weight_type = c("sampling", "frequency"),
                    strata = NULL,
                    cluster = NULL,
                    n_init = 20,
                    refine = TRUE,
                    refine_from = NULL,
                    max_iter = 1000,
                    n_cores = .default_n_cores(),
                    tol = 1e-8,
                    smoothing = 1.0,
                    random_state = NULL,
                    order_by_size = TRUE,
                    standard_errors = TRUE,
                    predictors_initial = NULL,
                    predictors_transition = NULL,
                    predictors_items = NULL,
                    predictors_random_intercept = NULL,
                    transition_effects = c("common", "by_origin"),
                    group = NULL,
                    group_effects = c("both", "initial", "transitions", "none"),
                    bayes_constants = NULL,
                    ...) {

  # Read before anything can touch `n_init`: `refine_from` refuses to be given
  # a restart budget, and "the user did not ask for one" is only knowable here.
  n_init_default <- missing(n_init)

  measurement_invariance <- match.arg(measurement_invariance)
  transition_invariance  <- match.arg(transition_invariance)
  random_intercept       <- match.arg(random_intercept)
  layout                 <- match.arg(layout)
  transition_effects     <- match.arg(transition_effects)
  group_effects          <- match.arg(group_effects)
  weight_type            <- match.arg(weight_type)
  if (!(isTRUE(standard_errors) || identical(standard_errors, FALSE) ||
        identical(standard_errors, "robust")))
    stop("`standard_errors` must be TRUE, FALSE or \"robust\".", call. = FALSE)

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
  # "slopes" pools the occasions exactly as "full" does and then gives each
  # occasion its own intercept back through the design, so it is a homogeneous
  # fit everywhere except in .lta_tau_design().
  tau_homogeneous <- transition_invariance %in% c("full", "slopes")
  tau_occasion_free_intercepts <- transition_invariance == "slopes"
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
    predictors_items      <- .subset_cases(predictors_items, keep)
    predictors_random_intercept <- .subset_cases(predictors_random_intercept, keep)
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
  if (tau_occasion_free_intercepts) {
    if (is.null(Z_tau))
      stop("`transition_invariance = \"slopes\"` shares the transition slopes ",
           "across occasions while leaving each occasion its own intercepts, ",
           "which is a restriction on a transition regression. Without ",
           "`predictors_transition` (or a `group` acting on the transitions) ",
           "there are no slopes to share; use \"full\" or \"none\".",
           call. = FALSE)
    if (identical(transition_effects, "by_origin"))
      stop("`transition_invariance = \"slopes\"` needs ",
           "`transition_effects = \"common\"`. Under \"by_origin\" every origin ",
           "status already has its own regression, so sharing slopes across ",
           "occasions restricts a different set of coefficients; specify the ",
           "model you want rather than letting it be guessed.", call. = FALSE)
  }

  if (random_intercept != "none") {
    if (measurement_invariance != "full")
      stop("A random intercept needs `measurement_invariance = \"full\"`: a ",
           "time-varying item intercept under a time-constant loading is not a ",
           "random-intercept model (Muthen & Asparouhov 2022, sec. 3.1).",
           call. = FALSE)
    if (!measurement %in% c("binary", "bernoulli", "ordinal"))
      stop("`random_intercept` currently supports binary or ordinal ",
           "indicators only.", call. = FALSE)
  }

  # --- direct covariate effects on the indicators -----------------------------
  # Every refusal here is a model shape the estimator does not cover, stated at
  # the argument rather than discovered inside an E-step. The ordinal one is
  # not a limitation of the covariate: a proportional-odds slope only means
  # something in the cumulative-logit parameterisation, and that is the ordinal
  # family's own.
  if (!is.null(predictors_items)) {
    if (!measurement %in% c("ordinal"))
      stop("`predictors_items` requires `measurement = \"ordinal\"`: a direct ",
           "covariate effect on an item is a proportional-odds slope, which ",
           "only exists in the cumulative-logit parameterisation. Binary ",
           "indicators can be passed as two-category ordinal data (coded 1 ",
           "and 2) to get the same model.", call. = FALSE)
    if (measurement_invariance != "full")
      stop("`predictors_items` requires `measurement_invariance = \"full\"`. ",
           "The covariate slopes are shared across occasions, which is a ",
           "restriction on a measurement model that is itself held equal ",
           "across them; with the occasions free there is nothing shared to ",
           "restrict.", call. = FALSE)
    if (n_classes > 1 || isTRUE(mover_stayer))
      stop("`predictors_items` is not yet available for a mixture latent ",
           "Markov model (`n_classes` > 1 or `mover_stayer = TRUE`). Fit the ",
           "mixture without it, or use one class.", call. = FALSE)
    if (random_intercept == "binary")
      stop("`predictors_items` is not available for a binary random ",
           "intercept, whose nodes are estimated latent classes rather than ",
           "points on a scale. Use `random_intercept = \"continuous\"` or ",
           "\"none\".", call. = FALSE)
  }

  if (random_intercept == "continuous") {
    n_quadrature <- as.integer(n_quadrature)
    # Q = 1 is not a corner case to reject: it is the node loop's own
    # structural test (a single node at 0 with mass 1 must reduce RI-LTA
    # exactly to regular LTA -- roadmap ### 14.10.3, ### 14.10.8 test 4).
    if (length(n_quadrature) != 1L || is.na(n_quadrature) || n_quadrature < 1L)
      stop("`n_quadrature` must be a single whole number of at least 1.",
           call. = FALSE)
  } else if (random_intercept == "binary") {
    n_ri <- as.integer(n_ri)
    if (length(n_ri) != 1L || is.na(n_ri) || n_ri < 2L)
      stop("`n_ri` must be a single whole number of at least 2.", call. = FALSE)
  }

  # --- covariates on the random intercept itself ------------------------------
  Z_ri <- .lta_ri_design(predictors_random_intercept, n,
                         "predictors_random_intercept")
  if (!is.null(Z_ri)) {
    if (random_intercept != "continuous")
      stop("`predictors_random_intercept` requires ",
           "`random_intercept = \"continuous\"`.", call. = FALSE)
    if (C > 1L || isTRUE(mover_stayer))
      stop("Covariates on the random intercept are not yet available for a ",
           "mixture latent Markov model (`n_classes` > 1 or ",
           "`mover_stayer = TRUE`). Fit the mixture without them, or use one ",
           "class.", call. = FALSE)
  }

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
    tau_occasion_free_intercepts = isTRUE(tau_occasion_free_intercepts),
    tie_initial_status = isTRUE(tie_initial_status),
    weights_vec     = w,
    weight_type       = weight_type,
    n_eff             = n_eff,
    strata            = if (is.null(strata))  rep(1L, n) else strata,
    cluster           = if (is.null(cluster)) seq_len(n) else cluster,
    has_survey_design = has_design,
    Z_delta         = Z_delta,
    Z_tau           = Z_tau,
    Z_ri            = Z_ri,
    ri_beta         = NULL,
    transition_effects = transition_effects,
    group_info      = group_info,
    group_effects   = if (is.null(group)) NULL else group_effects,
    bayes_constants = bayes_constants,
    mm              = time_blocks_model(K, prep$n_items, Tn,
                                        sub_model       = engine$sub_model,
                                        invariant_items = spec$invariant_items,
                                        max_val         = engine$max_val,
                                        cats            = engine$cats),
    ri              = .lta_ri_init(random_intercept, n_quadrature, n_ri),
    dif             = if (is.null(predictors_items)) NULL else
      .lta_dif_init(predictors_items, n, K, prep$n_items, "predictors_items")
  )

  # A DIF slope only means something in the cumulative-logit parameterisation,
  # and in this package that parameterisation lives only on the random-
  # intercept path -- the plain ordinal emission is a ragged multinoulli fitted
  # by closed-form counts. Rather than write a second cumulative-logit path,
  # a DIF fit with no random intercept rides the RI machinery with a DEGENERATE
  # factor: one node at z = 0 carrying all the mass, which is already a
  # documented, tested exact reduction to regular LTA (R/quadrature.R). The
  # loading is fixed at zero and flagged as such, so `.lta_ri_loading_free()`
  # keeps it out of the parameter count, the standard errors, the sign
  # normalisation, the random draw and the staged search -- that flag is the
  # difference between a 76-parameter fit and a 79-parameter one.
  if (!is.null(state$dif) && is.null(state$ri))
    state$ri <- list(kind = "continuous", Dnode = matrix(0, 1L, 1L), mass = 1,
                     A = NULL, L = NULL, loading_free = FALSE)

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
  # one LTA benchmark that is 20 restarts x 5000 iterations, and most of them
  # are climbing a hill they will lose anyway. The search is therefore staged, as the growth
  # mixture models' `em_stage1` stages theirs - a short first pass ranks the
  # restarts, and only the survivors are run on to convergence, resuming from
  # where they stopped. Three survivors rather than one because a short first
  # pass is a noisy ranking: the winning basin can be the slow one.
  # `refine_from` continues a fit that has already searched: it seeds one run
  # from that model's own converged parameters and runs on under this call's
  # stopping rule. It is not `fit_mixture(start_from = )`, which *replaces* a
  # search that has not happened; here the pool ran already and this is its
  # winner, so nothing is being skipped. See the argument's documentation.
  if (!is.null(refine_from)) {
    if (!inherits(refine_from, "lta_model"))
      stop("`refine_from` must be a model fitted by fit_lta().", call. = FALSE)
    if (!identical(as.integer(refine_from$n_statuses), as.integer(K)))
      stop(sprintf(paste0(
        "`refine_from` has %d statuses and this fit asks for %d. The donor ",
        "must be the same model."), refine_from$n_statuses, K), call. = FALSE)
    if (!identical(as.integer(refine_from$n_times), as.integer(Tn)))
      stop(sprintf(paste0(
        "`refine_from` has %d occasions and this fit asks for %d. The donor ",
        "must be the same model."), refine_from$n_times, Tn), call. = FALSE)
    if (!identical(as.integer(refine_from$n_classes), as.integer(C)))
      stop(sprintf(paste0(
        "`refine_from` has %d latent classes and this fit asks for %d. The ",
        "donor must be the same model."), refine_from$n_classes, C),
        call. = FALSE)
    if (!n_init_default)
      stop("`refine_from` continues that solution alone and runs no random ",
           "restarts, so `n_init` has nothing to size. Drop it.",
           call. = FALSE)
  }

  # Every restart below runs on `X_fit`, which is the response-pattern table
  # where the fit is eligible for one and `X` itself where it is not. The
  # search state carries the matching weights; both are put back on the full
  # sample once a winner has been chosen, before anything per-case is read off
  # it. See .lta_collapse().
  coll   <- .lta_collapse(state, X)
  X_fit  <- if (is.null(coll)) X else coll$X
  if (!is.null(coll)) state$weights_vec <- coll$w

  # A random intercept mixes over `Q` nodes each iteration and converges as
  # slowly as a mixture over chains does, so it takes the same tightened
  # defaults - and, for the same reason, the same staging. It has one chain, so
  # the `C > 1` test alone left it running every one of its (at least fifty)
  # restarts to `tol = 1e-11`: the tightening without the staging that was
  # written to pay for it. Nothing else about a single-chain fit changes here;
  # plain LTA keeps the unstaged search every locked benchmark figure was
  # measured on.
  # The predicate, not `!is.null(state$ri)`: the degenerate one-node factor a
  # `predictors_items` fit borrows integrates over nothing and converges at
  # plain-LTA speed, so it must keep plain LTA's unstaged search and its
  # untightened defaults rather than inheriting an RI fit's.
  staged <- C > 1L || (!is.null(state$ri) && .lta_ri_loading_free(state))
  if (staged) {
    if (missing(tol))      tol      <- 1e-11
    if (missing(max_iter)) max_iter <- 5000
  }
  if (!is.null(state$ri) && .lta_ri_loading_free(state) && n_init_default)
    n_init <- max(n_init, 50L)
  # Staging ranks a pool of restarts against each other. There is no pool to
  # rank when the start is handed over, so a refine goes straight to the full
  # stopping rule. The tightened `tol`/`max_iter` defaults just above still
  # apply: they are about how slowly a chain mixture converges, not about how
  # the search is organised.
  if (!is.null(refine_from)) staged <- FALSE
  # A tenth of the pool, rounded up, is a common staged-search convention.
  # The floor of 3 means every fit at
  # `n_init <= 30` promotes exactly as many survivors as before (bit-for-bit
  # unchanged); only `n_init >= 40` moves, and only upward, because taking
  # more of the top-ranked candidates can never lower the winner's score.
  # Deliberately NOT a function of `n_cores`: workers never draw random
  # numbers, and a survivor count that did would make a fit's reported number
  # depend on the hardware it ran on. Override with
  # `options(mixtureEM.lta_survivors = <n>)`.
  n_survivors <- if (staged) {
    min(n_init, max(3L, getOption("mixtureEM.lta_survivors",
                                   ceiling(0.10 * n_init))))
  } else 0L

  # The ranking pass integrates on a coarser grid than the fit reports on.
  # .lta_ri_e_step() runs one whole forward-backward pass per node, so an
  # iteration at `n_quadrature = 20` costs twenty times a plain LTA iteration -
  # and the number of nodes is an accuracy setting for one integral, not part
  # of the model. Ranking only has to identify the right basin, which five
  # nodes do; the survivors are promoted back to the full grid below and run to
  # `tol` there, so the reported fit integrates on exactly the grid the user
  # asked for. Continuous only: the binary variant's nodes are estimated
  # classes and `n_ri` IS part of the model.
  #
  # Drawing the starts from the coarse state takes the same random numbers the
  # full grid would: .lta_random_start() sizes `L` by `ncol(Dnode)`, which is 1
  # for a continuous random intercept at every `Q`.
  # ON by default at 5 nodes. When it was first measured, ranking on five
  # nodes halved the search but appeared to cost 0.62 of log-likelihood - but
  # that measurement was taken while the search was landing in the wrong
  # basin entirely, for reasons since fixed in how the restart pool is built.
  # Re-measured against the fixed search on the LTA-FAQ benchmark (continuous,
  # n_init = 50, random_state = 7): ll -14442.0171 against the no-ladder
  # -14442.0171 (gap -0.00001), 3/3 restarts replicating, in 391.8s against
  # 510-590s - no measurable accuracy cost for 25-30% less wall time, so the
  # option now defaults on. Still an option, not a hard-coded 5, so a user
  # who hits a fixture where 5 nodes cannot resolve the right basin can set
  # `options(mixtureEM.ri_rank_nodes = Inf)` to turn it off.
  n_rank_nodes <- getOption("mixtureEM.ri_rank_nodes", 5L)
  ri_ladder <- staged && !is.null(state$ri) &&
    identical(state$ri$kind, "continuous") &&
    length(state$ri$mass) > n_rank_nodes
  state_rank <- state
  if (ri_ladder) {
    coarse <- .lta_ri_init("continuous", n_rank_nodes, NULL)
    state_rank$ri$Dnode <- coarse$Dnode
    state_rank$ri$mass  <- coarse$mass
  }
  # Put a ranked candidate back on the full grid, carrying the parameters it
  # found. `A` is K x R and `L` is R x 1; neither is indexed by node, so both
  # transfer unchanged and the promoted start is the same model evaluated with
  # a more accurate integral.
  #
  # .lta_ri_sign_normalise() may have negated the ranked candidate's `Dnode`
  # alongside its `L`; overwriting `Dnode` with the unflipped full grid is
  # still the same model, because a Gauss-Hermite grid is symmetric in both its
  # nodes and its weights, so negating it only permutes the terms of a sum.
  # `options(mixtureEM.lta_ri_search = "wide")` (the default), on the
  # continuous-RI staged path only: the wide restart construction
  # (.lta_ri_wide_start()), a one-Newton-step M-step during the ranking
  # stage, and one full-grid E-step per ranked candidate before promotion.
  # Measured on the LTA-FAQ benchmark (RECORDS.md, "R12"):
  # the interior-bound candidates rank 1-3 of 100 on the full grid at 250
  # iterations and 7-12 on the five-node ladder, which mis-ranks
  # boundary-bound candidates upward. `"narrow"` leaves every fit bit-for-bit
  # as it was before this search existed; kept selectable to reproduce a
  # specific already-recorded validation run, not as a public compatibility
  # guarantee.
  ri_wide <- !identical(getOption("mixtureEM.lta_ri_search", "wide"), "narrow") &&
    staged && is.null(refine_from) && !is.null(state$ri) &&
    identical(state$ri$kind, "continuous") && .lta_ri_loading_free(state)
  promote <- function(cand) {
    if (inherits(cand, "try-error")) return(cand)
    if (ri_wide) cand$ri$gem <- NULL
    if (!ri_ladder) return(cand)
    cand$ri$Dnode <- state$ri$Dnode
    cand$ri$mass  <- state$ri$mass
    cand
  }

  best <- NULL
  best_degenerate <- TRUE
  # A boundary solution -- a categorical response probability driven towards
  # 0 or 1 by a continuous random intercept -- can score higher than an
  # interior optimum once the prior is off, by construction (RECORDS.md,
  # "R11"). Ranking on score() alone can therefore return a spurious winner.
  # A non-degenerate candidate always outranks a degenerate one regardless of
  # score; among candidates with the same status, score still decides. Never
  # discards a candidate outright -- if every candidate is degenerate the
  # best-scoring one still wins, flagged via fit$degenerate downstream.
  beats_current <- function(deg, s) {
    if (is.null(best)) return(TRUE)
    if (deg != best_degenerate) return(best_degenerate)
    s > best_score
  }
  stage1 <- list()
  # Every restart that ran to convergence, kept so the fit can report how many
  # of them found the reported maximum. On the staged path that is the second
  # loop: the first is a ranking pass stopped at 250 iterations, and its
  # log-likelihoods are not the maxima of anything.
  final_lls <- numeric(0)
  # Restarts are ranked on the objective EM climbed, not on the plain
  # log-likelihood it reports. Two candidates' penalised values and their plain
  # log-likelihoods are not monotonically related, so ranking on the latter can
  # return a point the search did not converge to as the winner. `loglik` still
  # reports the plain log-likelihood -- every information criterion, lr_test(),
  # lta_g2() and every validation target is defined on it. See
  # `### 41.1 re-opened` in internal/ROADMAP.md.
  rank_marg <- tryCatch(.lta_prior_marginals(state, X_fit),
                        error = function(e) NULL)
  score_of <- function(cand) {
    lp <- .lta_log_prior(cand, X_fit, alpha, rank_marg)
    if (is.na(lp)) cand$loglik else cand$loglik + lp
  }
  # The random starts are drawn here, in restart order, and fitted afterwards.
  # .lta_random_start() is the only RNG consumer on this path -- .lta_em() draws
  # nothing -- so this takes the same numbers the sequential loop took and makes
  # each restart a deterministic function of the start it is given. The fits can
  # then run on workers without moving a value, at any `n_cores`.
  starts <- if (!is.null(refine_from)) {
    # A refine is not settled by its donor alone. An RI fit takes its loadings
    # from a random draw the donor has nothing to say about (.lta_refine_start()
    # explains why they cannot start at zero), and on a multimodal surface that
    # draw decides which optimum the run reaches. Seeded here so `random_state`
    # means the same thing on this path as on the random-restart one. Nothing
    # else the draw produces survives - the donor overwrites delta, tau and the
    # measurement model - so no fit without a random intercept moves.
    if (!is.null(random_state)) set.seed(random_state)
    list(.lta_refine_start(state, X_fit, refine_from))
  } else {
    lapply(seq_len(max(1L, n_init)), function(i) {
      if (!is.null(random_state)) set.seed(random_state + i)
      s <- .lta_random_start(state_rank, X_fit)
      # A continuous random intercept's search benefits from mixing two
      # unrelated random constructions in one pool rather than drawing every
      # restart the same way; see .lta_ri_random_start2() for the second
      # construction, and why the binary variant is deliberately left alone.
      # Nothing without a random intercept moves.
      #
      # Only the second half of the pool gets it. Measured one construction
      # against the other, same pool composition, on two real benchmarks
      # that disagree about which half matters. `i` is unchanged from the
      # plain-restart numbering, so this only replaces the second half's
      # draws -- the first half's fits, and every non-continuous-RI fit, are
      # unaffected.
      # The predicate again: the second construction exists to diversify a
      # search over a real loading, and the degenerate factor has none.
      if (ri_wide) {
        # The third construction and the generalised-EM ranking stage
        # (.lta_ri_wide_start(), .wglm_newton_step()); replaces both halves
        # of the pool above. `ri$gem` is cleared again when a survivor is
        # promoted, so the reported fit is converged by the full M-step.
        s <- .lta_ri_wide_start(s, X_fit, i)
        s$ri$gem <- TRUE
      } else if (!is.null(s$ri) && identical(s$ri$kind, "continuous") &&
          .lta_ri_loading_free(s) && i > n_init %/% 2L)
        s <- .lta_ri_random_start2(s, X_fit)
      s
    })
  }
  # The polish belongs after a run that was allowed to converge, never after
  # the staged ranking pass: that pass stops at 250 iterations and its
  # log-likelihoods are, as the comment above says, not the maxima of anything.
  polish <- function(cand) {
    if (!isTRUE(refine) || inherits(cand, "try-error")) return(cand)
    out <- try(.lta_refine_lbfgs(cand, X_fit, alpha = alpha), silent = TRUE)
    if (inherits(out, "try-error")) cand else out
  }
  cands <- .par_lapply(starts, function(s) {
    cand <- try(.lta_em(s, X_fit,
                        max_iter = if (staged) min(250L, max_iter) else max_iter,
                        tol = if (staged) 1e-7 else tol, alpha = alpha),
                silent = TRUE)
    if (staged) cand else polish(cand)
  }, n_cores = n_cores)

  best_score <- -Inf
  for (cand in cands) {
    if (inherits(cand, "try-error")) next
    if (staged) stage1[[length(stage1) + 1L]] <- cand
    else {
      s <- score_of(cand)
      final_lls <- c(final_lls, s)
      deg <- !is.null(.categorical_boundary(cand, X_fit))
      if (beats_current(deg, s)) {
        best <- cand; best_score <- s; best_degenerate <- deg
      }
    }
  }

  if (staged && length(stage1)) {
    rank_score <- if (ri_wide && ri_ladder) {
      # One E-step on the full grid per candidate: rank on the objective the
      # survivors will be converged on rather than the ladder's.
      function(cand) {
        full <- try(.lta_em(promote(cand), X_fit, max_iter = 0L, tol = tol,
                            alpha = alpha), silent = TRUE)
        if (inherits(full, "try-error")) score_of(cand) else score_of(full)
      }
    } else score_of
    ord <- order(vapply(stage1, rank_score, numeric(1)), decreasing = TRUE)
    survivors <- .par_lapply(utils::head(ord, n_survivors), function(i) {
      cand <- try(.lta_em(promote(stage1[[i]]), X_fit, max_iter = max_iter,
                          tol = tol, alpha = alpha), silent = TRUE)
      if (inherits(cand, "try-error")) promote(stage1[[i]]) else polish(cand)
    }, n_cores = n_cores)
    for (cand in survivors) {
      s <- score_of(cand)
      final_lls <- c(final_lls, s)
      deg <- !is.null(.categorical_boundary(cand, X_fit))
      if (beats_current(deg, s)) {
        best <- cand; best_score <- s; best_degenerate <- deg
      }
    }
  }
  if (is.null(best))
    stop("Every random start failed; check the data and the model settings.",
         call. = FALSE)

  # Back onto the full sample before anything per-case is read off the fit:
  # every posterior, the path entropy, the standard errors and the metrics
  # below all describe cases, and the search saw patterns.
  best <- .lta_expand(best, coll, X, alpha)
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

  # A binary random-intercept fit's item table is rebuilt from the factor
  # (`.lta_ri_integrated_pis()`) every time the factor moves, and that rebuild
  # carries no dimnames, so the table reached print(), measurement_summary()
  # and random_intercept_loadings() as `Item_1`, `Item_2`, ... Name it the way
  # a regular fit's table is named, which is also what lets
  # measurement_summary() match the columns to the data for its Overall column.
  if (!is.null(best$ri) && is.null(best$ri$theta)) {
    nm <- .longitudinal_colnames(prep$item_names, prep$time_labels[1L])
    for (t in seq_len(Tn)) {
      p <- best$mm$models[[t]]$parameters$pis
      if (is.matrix(p) && ncol(p) == length(nm)) {
        colnames(p) <- nm
        best$mm$models[[t]]$parameters$pis <- p
      }
    }
  }

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
  # as the same solution when their scores are within 1e-2, the rule fit_em()
  # uses (R/em_core.R): genuinely different optima in these models sit whole
  # units apart, and a tighter rule splits one optimum into several. The score
  # is the ranked quantity -- the penalised objective where it exists -- so
  # "found the same solution" means the same thing here as "won the ranking".
  # It is on the same scale as a log-likelihood, so the 1e-2 rule carries over.
  if (length(final_lls)) {
    best$metrics$n_starts     <- length(final_lls)
    best$metrics$n_replicated <- sum(abs(final_lls - max(final_lls)) <= 1e-2)
  }
  # One start, and it was not a random one. Recorded rather than left at
  # `n_init` so .check_replication() cannot advise raising a restart budget on
  # a fit that never ran a pool; at n_requested = 1 its own rule (>= 10) already
  # declines to warn, and the flag says on the object where the solution came
  # from for anyone reading the fit later.
  best$metrics$n_requested <- if (is.null(refine_from)) max(1L, n_init) else 1L
  best$refined_from <- !is.null(refine_from)
  if (!identical(standard_errors, FALSE))
    best$se <- .lta_standard_errors(
      best, X, robust = identical(standard_errors, "robust"))

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

  # Raised once, on the fit returned, rather than inside the sign
  # normalisation every ranked candidate passes through: a search discards
  # most of its candidates, and a warning about one of those is noise.
  if (!is.null(best$ri) && .lta_ri_loading_free(best) &&
      any(abs(best$ri$L) > 10))
    warning("A random-intercept loading exceeds 10 in absolute value; the ",
            "model may be weakly identified on these data.", call. = FALSE)

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
# The response-pattern table
# ------------------------------------------------------------------------------

# Can this fit run on one row per distinct response pattern instead of one row
# per case?
#
# The forward-backward recursion costs O(n T K^2) an iteration, and a
# categorical panel repeats itself heavily: the benchmark's 3,092 cases hold
# 718 distinct answer patterns, so three quarters of every iteration is the
# same arithmetic done over again. One row per pattern carrying a frequency
# weight leaves the likelihood exactly as it was -- it is a weighted sum over
# cases either way -- for a quarter of the work.
#
# This is the economy fit_mixture() has had since .collapse_patterns()
# (R/em_core.R); the LTA runs its own driver and never got it. It is also
# precisely what `weights = , weight_type = "frequency"` already lets a user do
# by hand, so nothing new is being claimed about the estimator.
#
# NULL where it would not be the same fit:
#   * covariates on the initial status or the transitions, or a grouping
#     variable -- those designs are per case, and two people with the same
#     answers but different covariates are not one row;
#   * a survey design, where strata and cluster are per case for the same
#     reason;
#   * anything but a plain categorical measurement. Continuous data has no
#     duplicate rows to find, and a Gaussian start samples data rows as its
#     class means, so collapsing would move the search rather than shorten it;
#   * fewer than half the rows duplicated, where the bookkeeping is not repaid.
#
# The starts do not move on the paths this does accept. .lta_random_start()
# draws its initial and transition probabilities from Dirichlets that never see
# the data, and finishes at init_params(), which for `bernoulli` reads only
# ncol(X) and for `multinoulli` only ncol(X) and max(X) -- and every distinct
# value survives collapsing. Asserted in test-lta-collapse.R rather than left
# to this argument.
.lta_collapse <- function(state, X) {
  if (!is.null(state$Z_delta) || !is.null(state$Z_tau)) return(NULL)
  if (!is.null(state$Z_ri))             return(NULL)
  # Two cases with the same responses but different covariate values have
  # different emission tables, so they are not one pattern. (The family
  # whitelist below already excludes every model `predictors_items` accepts;
  # this states the reason rather than relying on that coincidence.)
  if (!is.null(state$dif))              return(NULL)
  if (!is.null(state$group_info))       return(NULL)
  if (isTRUE(state$has_survey_design))  return(NULL)

  sub <- state$mm$models[[1L]]
  if (is.null(sub) || !class(sub)[1] %in%
      c("bernoulli", "bernoulli_nan", "multinoulli", "multinoulli_nan"))
    return(NULL)

  n   <- nrow(X)
  pat <- .pattern_index(X)
  Xc  <- X[pat$rep_row, , drop = FALSE]
  if (nrow(Xc) > 0.5 * n) return(NULL)

  # rowsum() sorts its groups by label and the labels are 1..P, so the counts
  # line up with the rows of Xc.
  list(X = Xc, w = as.vector(rowsum(state$weights_vec, pat$idx)),
       w_orig = state$weights_vec, n_patterns = nrow(Xc))
}

# Put the winning fit back on the full sample after the search ran on the
# pattern table.
#
# Everything per case -- ll_case, gamma, the path entropy, the class posterior
# -- describes the rows the recursion saw, and those were patterns. They are
# rebuilt against the rows the user handed in. `max_iter = 0` is .lta_em()'s
# own way of running nothing but its closing E-step, the idiom
# .lta_refine_lbfgs() uses for the same reason. The weighted totals it also
# recomputes -- loglik, xi -- are unchanged by construction, a sum over
# patterns times their counts being the sum over cases.
.lta_expand <- function(best, coll, X, alpha) {
  if (is.null(coll)) return(best)
  best$weights_vec <- coll$w_orig
  out <- .lta_em(best, X, max_iter = 0L, alpha = alpha)
  # No iteration ran, so .lta_em() hands these back re-initialised rather than
  # as the search left them. See the same restoration in .lta_refine_lbfgs().
  out$converged     <- best$converged
  out$n_iter        <- best$n_iter
  out$refined_lbfgs <- best$refined_lbfgs
  out
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
    (K - 1L) * ncol(state$delta_beta)
  else if (isTRUE(state$tie_initial_status))
    (K - 1L)
  else (K - 1L) * C

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

  # A random intercept's integrated `pis` occupy the same `K x R` slot
  # `n_parameters(state$mm)` already counts; the loadings (and, for the binary
  # variant, the node masses beyond the anchor) are the only new parameters.
  #
  # The degenerate factor a `predictors_items` fit borrows has a loading fixed
  # at zero, which is not a parameter: counting its R zeros is exactly the
  # difference between the 76-parameter model and the 79-parameter one.
  n_ri <- if (is.null(state$ri)) 0L else
    (if (.lta_ri_loading_free(state)) length(state$ri$L) else 0L) +
      if (state$ri$kind == "binary") length(state$ri$mass) - 1L else 0L
  n_ri_beta <- if (is.null(state$Z_ri)) 0L else ncol(state$Z_ri)
  # One proportional-odds slope per (latent status, item, covariate), shared
  # across occasions.
  n_dif <- if (is.null(state$dif)) 0L else length(state$dif$beta)

  (C - 1L) + n_delta + n_tau + n_parameters(state$mm) + n_ri + n_ri_beta + n_dif
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

# Transition cells that have collapsed onto the boundary. On the logit scale
# this is a logit at a large negative value; here the cell is simply zero, and
# it is worth flagging because it costs a degree of freedom the count above
# still charges for.
#
# The threshold is 1e-4 rather than something nearer machine zero because EM
# stops when the *likelihood* stops moving, not when a cell reaches zero, and by
# then the cell is merely negligible: a transition that is structurally
# zero came out as 3.2e-06 on one benchmark, which a 1e-6 threshold missed
# entirely. Anything below
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

# How much of each transition row is prior rather than data. `alpha` is the
# mass on the whole transition table, so a row carries a = alpha / (K * C) of
# it; spread over the Ka admissible destinations of a row whose expected count
# is m, the posterior mean is the shrinkage of Fienberg & Holland (1973,
# eq. 2.6) with weight a/(m + a), so the largest amount the prior can move any
# cell of that row is
#
#     pull = [a / (m + a)] * (1 - 1/Ka),   a = alpha / (K * C)
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
  # `alpha` is the mass on the whole transition table, and a row carries
  # `alpha / (K * C)` of it -- .lta_normalise(patterns = K * C). Reporting the
  # pull as though the row carried all of `alpha` overstates it K-fold, which is
  # what this diagnostic did before 2026-08-31.
  a_row <- alpha / (K * C)

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
          n_expected = m, pull = (a_row / (m + a_row)) * (1 - 1 / Ka))
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
  abs_ent_by_occasion <- vapply(state$gamma, function(g)
    sum(state$weights_vec * (-g * log(g + 1e-15))), numeric(1))
  abs_ent <- sum(abs_ent_by_occasion)
  # Relative entropy of the joint status path: one assignment per case over
  # the K^T possible paths, so the normalising constant is n * log(K^T) =
  # n * T * log K. The numerator is the path entropy accumulated in
  # .lta_em(), not the sum of the per-occasion marginal entropies -- the
  # latter is always larger and is not the entropy of any single
  # classification.
  ent <- relative_entropy(state$abs_ent_path %||% abs_ent, n * Tn,
                          state$n_statuses)
  # Per-occasion relative entropy, normalised by n rather than n * T, on the
  # log K normaliser; `?fit_lta` gives the conversion to the
  # estimated-proportions normaliser.
  entropy_by_occasion <- vapply(abs_ent_by_occasion, relative_entropy,
                                numeric(1), n_samples = n,
                                n_classes = state$n_statuses)

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

  bic <- -2 * ll + log(n) * p
  list(ll = ll, n_params = p,
       aic   = -2 * ll + 2 * p,
       bic   = bic,
       caic  = -2 * ll + (log(n) + 1) * p,
       aic3  = -2 * ll + 3 * p,
       icl   = bic + 2 * abs_ent,
       sabic = -2 * ll + log((n + 2) / 24) * p,
       entropy = ent,
       entropy_by_occasion = entropy_by_occasion,
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
  # A random intercept's own measurement parameters ARE the measurement model;
  # the `pis` permuted below is only their integral over the nodes. Permuting
  # `pis` and leaving them where they were leaves the fit carrying two
  # different models at once, and every reader that recomputes from them
  # disagrees with every reader that does not. Both families are indexed by
  # status in their rows -- `A` for a binary indicator, `theta` for an ordinal
  # one -- and exactly one of the two is set. `L` is indexed by item and
  # `mass`/`Dnode` by node, so neither moves with a status.
  if (!is.null(state$ri)) {
    if (!is.null(state$ri$theta))
      state$ri$theta <- state$ri$theta[ord, , drop = FALSE]
    if (!is.null(state$ri$A))
      state$ri$A <- state$ri$A[ord, , drop = FALSE]
  }
  # `dif$beta` is the second status-indexed quantity outside `pis`, and leaving
  # it unpermuted is the same defect `ri$theta` above was paid for once
  # already: everything a user reads stays right and everything recomputed --
  # the scores, the sandwich, the MLR scaling factor -- is evaluated at a
  # scrambled model. Its other two margins are item and covariate and do not
  # move with a status.
  if (!is.null(state$dif))
    state$dif$beta <- state$dif$beta[ord, , , drop = FALSE]
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
.lta_standard_errors <- function(state, X, robust = FALSE) {
  sc <- .lta_score_matrix(state, X)
  if (is.null(sc)) return(NULL)
  S <- sc$S; blocks <- sc$blocks; conditional <- sc$conditional
  w <- state$weights_vec

  # The small-sample correction N / (N - 1) on the case-score outer product is
  # the usual finite-population fix for a sample covariance matrix estimated
  # from N contributions; applying it here (rather than leaving the raw outer
  # product) keeps this estimator's scale consistent with the same correction
  # applied to the sandwich meat below.
  n_cases <- state$n_eff %||% length(w)
  info <- (n_cases / max(n_cases - 1, 1)) * (t(S) %*% sweep(S, 1, w, "*"))
  V    <- tryCatch(pinv(info), error = function(e) NULL)
  if (is.null(V)) return(NULL)

  # With a complex survey design the empirical information is replaced by the
  # linearization sandwich, aggregating the same case-level scores to the
  # primary sampling unit within stratum. This is the machinery already behind
  # the covariate model's design-based variance (compute_survey_B in utils.R).
  meat <- info
  design_based <- FALSE
  if (isTRUE(state$has_survey_design)) {
    m <- tryCatch(
      compute_survey_B(sweep(S, 1, w, "*"), state$strata, state$cluster),
      error = function(e) NULL)
    if (!is.null(m)) {
      meat <- m
      V <- V %*% m %*% V
      design_based <- TRUE
    }
  }

  # Opt-in sandwich with the observed-information bread, the textbook form of
  # the estimator. The outer-product bread above is the
  # cheaper estimator and stays the default here: the finite-difference Hessian
  # this needs is 2p(p+1) forward-backward passes. The fallback is silent by
  # design -- a model the packing cannot describe gets the existing estimator
  # and `robust = FALSE` rather than an error -- which is why the flag is
  # returned and reported.
  robust_used <- FALSE
  if (isTRUE(robust) && .lta_par_packable(state)) {
    layout <- .lta_par_layout(state)
    par    <- .lta_par_pack(state, layout)
    if (length(par) == ncol(S)) {
      A  <- -.step1_fd_hessian(
        function(v) sum(w * .lta_ll_case(state, X, v, layout)), par)
      Ai <- .psd_pinv(A)
      Vr <- Ai %*% meat %*% Ai
      if (all(is.finite(Vr))) {
        V <- Vr
        robust_used <- TRUE
      }
    }
  }

  # The paper's headline estimate is the RI loading with its standard error
  # (`### 14.5`); surface it as a plain R x M matrix keyed the same way
  # `state$ri$L` already is, rather than making a caller parse block-name
  # strings. Raw numbers only -- formatted reporting is a later slice's job.
  loading_se <- NULL
  if (!is.null(state$ri) && .lta_ri_loading_free(state)) {
    n_items    <- nrow(state$ri$L)
    n_dim      <- ncol(state$ri$Dnode)
    loading_se <- matrix(NA_real_, n_items, n_dim)
    for (b in blocks) {
      m <- regmatches(b$name, regexec("^lambda\\[item (\\d+)\\]$", b$name))[[1]]
      if (length(m) == 2L)
        loading_se[as.integer(m[2]), ] <- sqrt(pmax(diag(V)[b$cols], 0))
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

  list(vcov = V, blocks = blocks, prob_se = prob_se, loading_se = loading_se,
       conditional = conditional, design_based = design_based,
       robust = robust_used)
}

# Where the score blocks below are the right ones. Both `.lta_standard_errors()`
# and the L-BFGS refinement in R/lta_core.R read this, so the parameterisation
# one of them trusts cannot drift away from the parameterisation the other
# trusts.
#
# A mixture over chains adds a class-membership block and makes every other
# score class-conditional; both are below, licensed by the Fisher identity
# (Louis, 1982): the gradient of a mixture's observed-data log-likelihood is
# the posterior-class-weighted sum of the gradients of its complete-data
# (class-conditional) log-likelihoods, so every single-chain score is simply
# scaled by the case's posterior probability of that class. Covariate models
# replace delta and tau with case-level regressions (`delta_beta`/`tau_beta`);
# the score of that multinomial logit is the same (observed - predicted)
# residual, now weighted by the covariate row instead of read off a fixed
# probability vector, which is what the "_beta" blocks below compute. `C > 1`
# and covariates are mutually exclusive in `fit_lta()`, so the two never have
# to compose with each other.
# Which items are held equal across occasions. Read off the measurement model,
# which carries the constraint from the moment fit_lta() builds it, NOT off
# `state$longitudinal`, which fit_lta() attaches only once the search is over.
# The refinement runs *during* the search, and reading the empty slot there
# made it treat every item as free: it then bought log-likelihood by breaking
# the invariance restriction the model is defined by, and the next EM step
# re-imposed the restriction and took it straight back.
.lta_invariant_items <- function(state) {
  state$mm$invariant_items %||% state$longitudinal$invariant_items %||%
    integer(0)
}

.lta_scores_supported <- function(state) {
  if (state$n_statuses < 2L || state$n_times < 2L) return(FALSE)
  # The multi-class branch of `.lta_score_matrix()` now runs
  # `.lta_ri_e_step()` per the RI arm above, populating `ri_G` / `ri_pq` for
  # the measurement blocks, so a multi-class RI fit needs no guard here.
  # A tied initial-status distribution is handled by .lta_par_layout()/
  # .lta_par_pack()/.lta_par_unpack() collapsing the C class-level delta
  # blocks into one shared block, so the Jacobian below describes the model
  # actually fitted rather than treating shared numbers as independent.
  TRUE
}

# The measurement families whose parameters the score blocks actually cover.
# `.lta_standard_errors()` does not need this: where the family is not one of
# these it reports delta and tau with the measurement model treated as known
# and says so through `conditional`. The refinement does need it, because a
# parameter it cannot differentiate is a parameter it must not move.
#
# A random intercept fit is no longer excluded here (`### 14.15` W8): once
# `.lta_penalty()` grew an RI branch and the box rule stopped clamping
# `lambda`, the polish climbs the same penalised objective EM does for every
# RI shape too. This predicate and `.lta_par_packable()` below now differ
# only in the measurement-family whitelist -- `.lta_par_packable()` is
# wider only for families this whitelist has not yet been taught, not for
# any structural reason -- and are kept separate so a future family addition
# has to decide deliberately whether it supports the polish, not just the
# packing.
.lta_scores_full <- function(state) {
  .lta_scores_supported(state) &&
    class(state$mm$models[[1]])[1] %in%
      c("bernoulli", "bernoulli_nan", "gaussian_diag", "gaussian_diag_nan",
        "gaussian_unit", "gaussian_unit_nan", "ordinal", "ordinal_nan")
}

# Whether the packed vector describes this model's full free parameter set,
# and therefore whether a finite-difference Hessian on it means anything.
# Wider than `.lta_scores_full()` by exactly the RI case, which `### 14.15`'s
# W2/W3 taught the packing and the likelihood to handle.
#
# The two predicates are deliberately separate. `.lta_scores_full()` is ALSO
# the L-BFGS polish gate (`.lta_refine_lbfgs()`), and the polish additionally
# needs an RI branch in `.lta_penalty()` and a revised box rule. Until that
# has been done, the sandwich may run on an RI fit and the polish may not.
.lta_par_packable <- function(state) {
  if (!.lta_scores_supported(state)) return(FALSE)
  if (!is.null(state$ri)) return(TRUE)
  class(state$mm$models[[1]])[1] %in%
    c("bernoulli", "bernoulli_nan", "gaussian_diag", "gaussian_diag_nan",
      "gaussian_unit", "gaussian_unit_nan", "ordinal", "ordinal_nan")
}

# The n x p matrix of case-level scores, one column per free parameter, on the
# multinomial-logit scale with the last category anchored. Column order is the
# packing order: delta, then the transition rows, then the measurement blocks.
# `.lta_refine_lbfgs()` sums these rows under the case weights to get the
# gradient of the log-likelihood; `.lta_standard_errors()` takes their outer
# product to get the empirical information. Same matrix, two uses.
#
# The forward-backward pass is run here rather than read off `state`, so the
# scores describe the parameters currently in `state` even when no E-step has
# been taken at them. At a converged fit the two agree by construction.
.lta_score_matrix <- function(state, X) {
  K  <- state$n_statuses
  Tn <- state$n_times
  C  <- state$n_classes %||% 1L
  w  <- state$weights_vec
  n  <- nrow(X)
  if (!.lta_scores_supported(state)) return(NULL)

  logB <- .lta_emission_loglik(state$mm, X)

  scores <- list(); blocks <- list(); pos <- 0L
  add_block <- function(s, probs, ref, name) {
    scores[[length(scores) + 1L]] <<- s
    blocks[[length(blocks) + 1L]] <<- list(
      cols = pos + seq_len(ncol(s)), probs = probs, ref = ref, name = name)
    pos <<- pos + ncol(s)
  }

  # One class runs the single chain unchanged: `fb$gamma` is the status
  # posterior directly and `state$tau_allowed` is already the collapsed
  # (non-list) view .lta_pack() maintains for this case.
  ri_G <- NULL
  ri_pq <- NULL
  if (C == 1L) {
    if (!is.null(state$ri)) {
      # A random intercept is shared across occasions and items within a
      # case, so a single forward-backward pass on the node-marginalised
      # emission (`logB` above) would drop that within-case correlation.
      # .lta_ri_e_step() already runs the right thing -- one forward-backward
      # per node, mixed by the node posterior -- and returns the same shape
      # a single-class .lta_em() e_step() does, plus the joint
      # node-and-status posterior `G` the measurement block below needs.
      e   <- .lta_ri_e_step(state, X, w)
      fb  <- e$es[[1L]]
      gam <- fb$gamma
      ll  <- e$ll
      ri_G <- e$ri$G
      ri_pq <- e$ri$pq
    } else {
      # .lta_log_delta()/.lta_log_tau() (R/lta_covariates.R) already return
      # the case-varying log-probabilities a covariate model implies, and
      # fall back to the plain broadcast delta/tau otherwise -- the same
      # functions the ordinary E-step calls, so the log-likelihood here needs
      # no branch on whether covariates are present.
      fb <- .lta_forward_backward(logB, .lta_log_delta(state), .lta_log_tau(state),
                                  w, keep_pairwise = TRUE)
      gam <- fb$gamma
      ll  <- fb$ll
    }

    if (!is.null(state$delta_beta)) {
      # Ordinary multinomial-logit score, covariate-weighted: the same
      # (observed - predicted) residual as the plain "delta" block below,
      # times the covariate row, which is exactly .fit_mnl()'s own gradient
      # (R/covariate.R) kept at the case level instead of summed. Column
      # order is status k (1..K-1) slow, covariate d fast, matching
      # .lta_par_pack()'s `as.vector(t(B[1:(K-1), ]))`.
      phat <- exp(.lta_log_delta(state))
      s <- do.call(cbind, lapply(seq_len(K - 1L), function(k)
        (gam[[1]][, k] - phat[, k]) * state$Z_delta))
      add_block(s, NULL, NULL, "delta_beta")
    } else {
      # delta: the score of the multinomial logit is gamma_1[i, k] - delta_k,
      # with the last status anchored.
      add_block(sweep(gam[[1]][, seq_len(K - 1L), drop = FALSE], 2,
                      state$delta[seq_len(K - 1L)], "-"),
                state$delta, seq_len(K - 1L), "delta")
    }

    idx <- if (isTRUE(state$tau_homogeneous)) 1L else seq_len(Tn - 1L)
    if (!is.null(state$tau_beta)) {
      # tau_beta: same (observed - predicted) ⊗ covariate-row score as delta
      # above, one block per origin status under "by_origin" (mirroring
      # .lta_mstep_tau_cov()'s by-origin loop, R/lta_covariates.R:139-157) or
      # one block per matrix under "common", summing every origin's
      # contribution into the shared coefficient set (mirroring that
      # function's "common" branch, R/lta_covariates.R:159-183, which stacks
      # every (occasion, origin) row into one .fit_mnl() call for the same
      # reason).
      by_origin <- identical(state$transition_effects, "by_origin")
      lt <- .lta_log_tau(state)
      for (i_mat in idx) {
        ts <- if (isTRUE(state$tau_homogeneous)) seq_len(Tn - 1L) else i_mat
        if (by_origin) {
          for (k in seq_len(K)) {
            Zk <- state$Z_tau
            s  <- matrix(0, n, (K - 1L) * ncol(Zk))
            for (tt in ts) {
              phat_k <- exp(lt[[tt]][[k]])
              pk     <- fb$pairwise[[tt]][[k]]
              resid  <- pk[, seq_len(K - 1L), drop = FALSE] -
                gam[[tt]][, k] * phat_k[, seq_len(K - 1L), drop = FALSE]
              s <- s + do.call(cbind, lapply(seq_len(K - 1L), function(d)
                resid[, d] * Zk))
            }
            add_block(s, NULL, NULL,
                      sprintf("tau_beta%s[from %d]",
                              if (isTRUE(state$tau_homogeneous)) "" else
                                paste0("(", i_mat, ")"), k))
          }
        } else {
          D <- ncol(.lta_tau_design(state, 1L))
          s <- matrix(0, n, (K - 1L) * D)
          for (tt in ts) for (k in seq_len(K)) {
            Zk     <- .lta_tau_design(state, k, tt)
            phat_k <- exp(lt[[tt]][[k]])
            pk     <- fb$pairwise[[tt]][[k]]
            resid  <- pk[, seq_len(K - 1L), drop = FALSE] -
              gam[[tt]][, k] * phat_k[, seq_len(K - 1L), drop = FALSE]
            s <- s + do.call(cbind, lapply(seq_len(K - 1L), function(d)
              resid[, d] * Zk))
          }
          add_block(s, NULL, NULL,
                    sprintf("tau_beta%s",
                            if (isTRUE(state$tau_homogeneous)) "" else
                              paste0("(", i_mat, ")")))
        }
      }
    } else {
      # tau: xi_i[t, k, l] - gamma_t[i, k] * tau_t[k, l], anchored on the last
      # admissible destination in the row.
      for (i_mat in idx) {
        ts <- if (isTRUE(state$tau_homogeneous)) seq_len(Tn - 1L) else i_mat
        for (k in seq_len(K)) {
          allowed <- which(state$tau_allowed[[i_mat]][k, ])
          if (length(allowed) < 2L) next
          free <- allowed[-length(allowed)]
          s <- matrix(0, n, length(free))
          for (tt in ts) {
            pk <- fb$pairwise[[tt]][[k]]        # n x K: P(S_t = k, S_{t+1} = .)
            s <- s + pk[, free, drop = FALSE] -
              outer(gam[[tt]][, k], state$tau[[i_mat]][k, free])
          }
          add_block(s, state$tau[[i_mat]][k, ], free,
                    sprintf("tau%s[from %d]",
                            if (isTRUE(state$tau_homogeneous)) "" else
                              paste0("(", i_mat, ")"), k))
        }
      }
    }
  } else {
    # Several classes: run the single-chain recursion once per class -- the
    # same shape .lta_em()'s e_step() uses -- combine into the mixture
    # log-likelihood and the posterior class membership, then scale every
    # single-chain score by that posterior. `.lta_mixed_gamma()` (R/lta_core.R)
    # is the same class-mixed status posterior the M-step's measurement update
    # already uses, reused here for exactly the reason it exists there: the
    # measurement model is shared across classes, so its score reads the
    # mixed posterior rather than any one class's own.
    if (!is.null(state$ri)) {
      # Same reason the C == 1 branch above calls this: the random intercept
      # is shared across occasions and items within a case, so a
      # forward-backward on the node-marginalised emission would drop that
      # within-case correlation. .lta_ri_e_step() runs the per-class,
      # per-node double loop the E-step itself now runs, and returns `es`,
      # `post` and `ll` in exactly the shape the plain branch below builds by
      # hand -- its `es[[c]]$gamma` and `$pairwise` are already mixed over
      # nodes -- plus `ri$G` / `ri$pq`, which the measurement blocks need.
      e     <- .lta_ri_e_step(state, X, w)
      es    <- e$es
      post  <- e$post
      ll    <- e$ll
      ri_G  <- e$ri$G
      ri_pq <- e$ri$pq
    } else {
      es <- lapply(seq_len(C), function(c) {
        sub <- .lta_class_state(state, c)
        .lta_forward_backward(logB, log(pmax(sub$delta, 1e-300)),
                              lapply(sub$tau, function(m) log(pmax(m, 1e-300))),
                              w, keep_pairwise = TRUE)
      })
      lp   <- vapply(seq_len(C), function(c) es[[c]]$ll, numeric(n))
      lp   <- sweep(lp, 2, log(pmax(state$class_weights, 1e-300)), "+")
      ll   <- logsumexp(lp, MARGIN = 1)
      post <- exp(lp - ll)
    }
    gam  <- .lta_mixed_gamma(list(es = es, post = post), Tn, C)

    # class weights: the score of the mixing-proportion logit is
    # P(class = c | y_i) - pi_c, with the last class anchored -- the same
    # multinomial-logit form as delta below, one level up.
    add_block(sweep(post[, seq_len(C - 1L), drop = FALSE], 2,
                    state$class_weights[seq_len(C - 1L)], "-"),
              state$class_weights, seq_len(C - 1L), "class")

    tied <- isTRUE(state$tie_initial_status)
    if (tied) {
      # Summing the per-class delta scores over c collapses to the single-chain
      # score on the class-mixed posterior, because the classes share one delta
      # and the class posteriors sum to 1. See the roadmap's 14.13 derivation.
      add_block(sweep(gam[[1]][, seq_len(K - 1L), drop = FALSE], 2,
                      state$delta_c[[1]][seq_len(K - 1L)], "-"),
                state$delta_c[[1]], seq_len(K - 1L), "delta")
    }

    idx <- if (isTRUE(state$tau_homogeneous)) 1L else seq_len(Tn - 1L)
    for (c in seq_len(C)) {
      gam_c   <- es[[c]]$gamma
      delta_c <- state$delta_c[[c]]

      if (!tied)
        add_block(post[, c] * sweep(gam_c[[1]][, seq_len(K - 1L), drop = FALSE],
                                    2, delta_c[seq_len(K - 1L)], "-"),
                  delta_c, seq_len(K - 1L), sprintf("delta[class %d]", c))

      for (i_mat in idx) {
        ts <- if (isTRUE(state$tau_homogeneous)) seq_len(Tn - 1L) else i_mat
        for (k in seq_len(K)) {
          allowed <- which(state$tau_allowed_c[[c]][[i_mat]][k, ])
          if (length(allowed) < 2L) next
          free <- allowed[-length(allowed)]
          s <- matrix(0, n, length(free))
          for (tt in ts) {
            pk <- es[[c]]$pairwise[[tt]][[k]]
            s <- s + pk[, free, drop = FALSE] -
              outer(gam_c[[tt]][, k], state$tau_c[[c]][[i_mat]][k, free])
          }
          add_block(post[, c] * s, state$tau_c[[c]][[i_mat]][k, ], free,
                    sprintf("tau%s[from %d, class %d]",
                            if (isTRUE(state$tau_homogeneous)) "" else
                              paste0("(", i_mat, ")"), k, c))
        }
      }
    }
  }

  # Measurement parameters. Including them is what makes these standard errors
  # unconditional; when the emission family is not one of the two handled here
  # the measurement model is treated as known, which understates uncertainty,
  # and the flag below says so.
  J <- state$n_items
  conditional <- TRUE
  if (!is.null(state$ri)) {
    # d/d alpha[k,j] = sum_t sum_q G[[q]][[t]][i,k] (y_ijt - pi(k,j,q))
    # d/d lambda[j,]  = the same, times the node value Dnode[q, ], summed over
    # k as well as q since a loading has no status subscript. Mirrors
    # .lta_ri_mstep()'s own succ/tot aggregation (R/lta_ri.R) at the case
    # level instead of aggregated -- same (item, node, occasion) loop, same
    # missing-data mask, the M-step's score rather than a new derivation.
    conditional <- FALSE
    ri <- state$ri
    Q  <- length(ri$mass)
    M  <- ncol(ri$Dnode)
    ordinal_ri <- !is.null(ri$theta)
    cats <- if (ordinal_ri) ri$cats else NULL
    # No block for a loading that is not estimated; the layout omits it too,
    # and the two descriptions of one vector must agree.
    loading_free <- .lta_ri_loading_free(state)
    for (j in seq_len(J)) {
      if (ordinal_ri) {
        Sj      <- cats[j]
        cols    <- .ordinal_theta_cols(cats, j)
        theta_j <- ri$theta[, cols, drop = FALSE]
        s_theta  <- matrix(0, n, K * (Sj - 1L))
        s_lambda <- matrix(0, n, M)
        # With direct covariate effects the shift depends on the case's
        # covariate pattern as well as the node, so the node loop becomes a
        # node AND pattern loop. Two facts already in the code make the DIF
        # block free: shifting the whole linear predictor by c is the same as
        # shifting the base threshold by c (stated at .lta_ri_sign_normalise()
        # and already used for the loading just below), so `blk`'s first K
        # columns ARE d log p / d shift_k; and rows outside pattern p are
        # already zero in them, so scaling by the pattern's covariate row
        # needs no mask.
        dif <- state$dif
        Dd  <- if (is.null(dif)) 0L else ncol(dif$Zu)
        s_dif <- if (is.null(dif)) NULL else matrix(0, n, K * Dd)
        pats <- if (is.null(dif)) 1L else seq_len(dif$P)
        for (q in seq_len(Q)) {
          node <- sum(ri$L[j, ] * ri$Dnode[q, ])
          for (p in pats) {
            rows  <- if (is.null(dif)) seq_len(n) else dif$rows[[p]]
            shift <- if (is.null(dif)) node else
              node + .lta_dif_shift(dif, p)[, j]
            blk   <- matrix(0, n, K * (Sj - 1L))
            for (tt in seq_len(Tn)) {
              xj <- X[, .time_block_cols(tt, J)[j]]
              blk[rows, ] <- blk[rows, ] +
                .ordinal_theta_score_block(theta_j, shift, xj[rows],
                                           ri_G[[q]][[tt]][rows, , drop = FALSE],
                                           Sj)
            }
            s_theta <- s_theta + blk
            base    <- blk[, seq_len(K), drop = FALSE]   # d log p / d shift_k
            s_lambda <- s_lambda + outer(rowSums(base), ri$Dnode[q, ])
            if (!is.null(dif))
              for (dd in seq_len(Dd))
                s_dif[, (dd - 1L) * K + seq_len(K)] <-
                  s_dif[, (dd - 1L) * K + seq_len(K)] + base * dif$Zu[p, dd]
          }
        }
        add_block(s_theta, NULL, NULL, sprintf("theta[item %d]", j))
        if (loading_free)
          add_block(s_lambda, NULL, NULL, sprintf("lambda[item %d]", j))
        if (!is.null(dif))
          add_block(s_dif, NULL, NULL, sprintf("dif[item %d]", j))
        next
      }
      s_alpha  <- matrix(0, n, K)
      s_lambda <- matrix(0, n, M)
      for (q in seq_len(Q)) {
        eta_q <- ri$A[, j] + as.vector(ri$Dnode[q, , drop = FALSE] %*% ri$L[j, ])
        pi_q  <- plogis(eta_q)
        resid <- matrix(0, n, K)
        for (tt in seq_len(Tn)) {
          xj  <- X[, .time_block_cols(tt, J)[j]]
          obs <- !is.na(xj); xj[!obs] <- 0
          resid <- resid + ri_G[[q]][[tt]] *
            (matrix(xj, n, K) - matrix(pi_q, n, K, byrow = TRUE)) * obs
        }
        s_alpha  <- s_alpha + resid
        s_lambda <- s_lambda + outer(rowSums(resid), ri$Dnode[q, ])
      }
      add_block(s_alpha, NULL, NULL, sprintf("alpha[item %d]", j))
      if (loading_free)
        add_block(s_lambda, NULL, NULL, sprintf("lambda[item %d]", j))
    }
    if (identical(ri$kind, "binary") && Q > 1L) {
      # The node masses are a mixing proportion at the case level, so their
      # multinomial-logit score is the same (posterior - prior) residual as
      # the class-weight block one level up: P(node = q | y_i) - mass_q,
      # anchored on the last node. `pq` IS that posterior -- pooled over
      # class with each case's own class posterior by .lta_ri_e_step() -- and
      # is what .lta_ri_mstep() already sums under the case weights to update
      # `mass`. The continuous variant's Gauss-Hermite weights are FIXED and
      # must never get a block here, the same restriction .lta_ri_mstep()
      # observes.
      add_block(sweep(ri_pq[, seq_len(Q - 1L), drop = FALSE], 2,
                      ri$mass[seq_len(Q - 1L)], "-"),
                ri$mass, seq_len(Q - 1L), "ri_mass")
    }
    if (!is.null(state$Z_ri)) {
      # d/dbeta log L_i = x_i (fbar_i - mu_i). Fisher identity on the node
      # prior: d/dbeta log mass_iq = x_i (z_q - mu_i), weighted by the node
      # posterior.
      fbar <- as.vector(ri_pq %*% ri$Dnode[, 1L])
      mu   <- as.vector(state$Z_ri %*% state$ri_beta)
      add_block(state$Z_ri * (fbar - mu), NULL, NULL, "ri_beta")
    }
    fam <- NULL
  } else {
  fam <- class(state$mm$models[[1]])[1]
  if (fam %in% c("bernoulli", "bernoulli_nan")) {
    conditional <- FALSE
    inv <- .lta_invariant_items(state)
    for (j in seq_len(J)) {
      ts_groups <- if (j %in% inv) list(seq_len(Tn)) else
        lapply(seq_len(Tn), identity)
      for (grp in ts_groups) {
        s <- matrix(0, n, K)
        for (tt in grp) {
          xj  <- X[, .time_block_cols(tt, J)[j]]
          rho <- state$mm$models[[tt]]$parameters$pis[, j]
          obs <- !is.na(xj); xj[!obs] <- 0
          s <- s + gam[[tt]] * (xj - matrix(rho, n, K, byrow = TRUE)) *
            obs
        }
        add_block(s, NULL, NULL, sprintf("rho[item %d]", j))
      }
    }
  } else if (fam %in% c("gaussian_diag", "gaussian_diag_nan",
                        "gaussian_unit", "gaussian_unit_nan")) {
    # Both means and (where free) variances are scored below, so this family
    # is fully unconditional -- the same status the Bernoulli branch above
    # already gets.
    conditional <- FALSE
    inv <- .lta_invariant_items(state)
    for (j in seq_len(J)) {
      ts_groups <- if (j %in% inv) list(seq_len(Tn)) else
        lapply(seq_len(Tn), identity)
      for (grp in ts_groups) {
        has_var <- !is.null(state$mm$models[[grp[1]]]$parameters$covariances)
        s_mu <- matrix(0, n, K)
        s_v  <- if (has_var) matrix(0, n, K) else NULL
        for (tt in grp) {
          xj  <- X[, .time_block_cols(tt, J)[j]]
          sub <- state$mm$models[[tt]]
          mu  <- sub$parameters$means[, j]
          v   <- if (has_var) sub$parameters$covariances[, j] else rep(1, K)
          obs <- !is.na(xj); xj[!obs] <- 0
          resid <- matrix(xj, n, K) - matrix(mu, n, K, byrow = TRUE)
          s_mu <- s_mu + gam[[tt]] * sweep(resid, 2, v, "/") * obs
          if (has_var)
            s_v <- s_v + gam[[tt]] *
              (sweep(resid^2, 2, v, "/") - 1) * obs
        }
        add_block(s_mu, NULL, NULL, sprintf("mu[item %d]", j))
        # Log-sd scale: d(loglik)/d(log sigma) = (x - mu)^2 / sigma^2 - 1,
        # the same reparameterisation R/em_core.R's refine_lbfgs() uses for
        # the plain mixture engine's Gaussian variances.
        if (has_var)
          add_block(s_v, NULL, NULL, sprintf("log_sd[item %d]", j))
      }
    }
  } else if (fam %in% c("ordinal", "ordinal_nan")) {
    # Non-RI ordinal: the single-"node" case of the RI arm above (shift = 0,
    # one status-posterior `gam[[tt]]` in place of a node-and-status joint
    # posterior), so it reuses the same .ordinal_theta_score_block() helper.
    conditional <- FALSE
    inv <- .lta_invariant_items(state)
    for (j in seq_len(J)) {
      ts_groups <- if (j %in% inv) list(seq_len(Tn)) else
        lapply(seq_len(Tn), identity)
      for (grp in ts_groups) {
        m       <- state$mm$models[[grp[1]]]
        cats    <- m$cats
        Sj      <- cats[j]
        cols    <- .ordinal_theta_cols(cats, j)
        theta_j <- .ordinal_theta_from_pis(m$parameters$pis, cats)[, cols,
                                                                    drop = FALSE]
        s <- matrix(0, n, K * (Sj - 1L))
        for (tt in grp) {
          xj <- X[, .time_block_cols(tt, J)[j]]
          s  <- s + .ordinal_theta_score_block(theta_j, 0, xj, gam[[tt]], Sj)
        }
        add_block(s, NULL, NULL, sprintf("theta[item %d]", j))
      }
    }
  }
  }

  list(S = do.call(cbind, scores), blocks = blocks,
       conditional = conditional, ll = ll, gamma = gam)
}
