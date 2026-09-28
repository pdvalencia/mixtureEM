# Latent Transition Analysis

Fits a latent transition model: each person occupies a latent *status*
at each occasion, measured by the same indicators every time, and may
move between statuses from one occasion to the next. Where
[`fit_rmlca()`](https://pdvalencia.github.io/mixtureEM/reference/fit_rmlca.md)
assigns a person one class for the whole study, `fit_lta()` estimates
both the prevalence of each status and the *incidence of change* between
them (Collins & Lanza, 2010, ch. 7).

Three things are estimated:

- status prevalences:

  how common each status is at the first occasion;

- transition probabilities:

  the chance of moving from each status to each other status, as a
  square table read from row (earlier occasion) to column (later
  occasion). There is one table per pair of adjacent occasions unless
  `transition_invariance` restricts them;

- item parameters:

  what people in each status tend to answer, which is what gives the
  statuses their meaning.

Cases with individual items or whole occasions missing are kept in the
analysis rather than dropped.

**Why measurement invariance defaults to `"full"`.** If a status is not
defined identically at every occasion, an apparent "transition" mixes
real change in the person with a change in what the status means, and
the two cannot be told apart. Holding the item parameters equal across
occasions removes that ambiguity, and is the usual practice (Collins &
Lanza, sec. 7.11). It is a testable restriction: fit the model both ways
and compare them with
[`lr_test()`](https://pdvalencia.github.io/mixtureEM/reference/lr_test.md).

## Usage

``` r
fit_lta(
  indicators,
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
  id = NULL,
  time = NULL,
  items = NULL,
  item_names = NULL,
  time_labels = NULL,
  weights = NULL,
  weight_type = c("sampling", "frequency"),
  strata = NULL,
  cluster = NULL,
  n_init = 20,
  refine = TRUE,
  refine_from = NULL,
  max_iter = 1000,
  n_cores = .default_n_cores(),
  tol = 1e-08,
  smoothing = 1,
  random_state = NULL,
  order_by_size = TRUE,
  standard_errors = TRUE,
  predictors_initial = NULL,
  predictors_transition = NULL,
  predictors_items = NULL,
  predictors_items_by_status = NULL,
  predictors_random_intercept = NULL,
  n_steps = 1,
  transition_effects = c("common", "by_origin"),
  group = NULL,
  group_effects = c("both", "initial", "transitions", "none"),
  bayes_constants = NULL,
  correction = c("ML", "BCH", "none"),
  assignment = c("proportional", "modal"),
  distal = NULL,
  ...
)
```

## Arguments

- indicators:

  The repeated indicators, in any format accepted by
  [`fit_rmlca()`](https://pdvalencia.github.io/mixtureEM/reference/fit_rmlca.md):
  a wide matrix, a three-dimensional array, or a long data frame with
  `id` and `time`.

- n_statuses:

  Integer. Number of latent statuses.

- times:

  Integer. Number of occasions. Required for wide input.

- measurement:

  Measurement model for one occasion's items: `"binary"`,
  `"categorical"`, `"ordinal"`, `"continuous"`, or a named list for a
  mixed block. `"ordinal"` fits one cumulative-logit block per item,
  with its own number of ordered categories inferred per item (so a
  3/3/2-category block needs no mixed specification); without a random
  intercept it is numerically identical to `"categorical"`.

- measurement_invariance:

  Whether the item parameters are held equal across occasions: `"full"`
  (the default), `"none"`, or `"partial"` for only the items named in
  `invariant_items`. See the note above on why `"full"` is the sensible
  starting point.

- invariant_items:

  Items held equal across occasions when
  `measurement_invariance = "partial"`, given by name or position.

- transition_invariance:

  Whether the transition probabilities are held equal across occasions.
  `"none"` (the default) estimates a separate matrix for each pair of
  adjacent occasions; `"full"` shares one matrix throughout. Whether
  change happens at a constant rate is usually a substantive question
  rather than an assumption, and the two models are nested, so
  [`lr_test()`](https://pdvalencia.github.io/mixtureEM/reference/lr_test.md)
  tests it (Collins & Lanza, sec. 7.14).

  `"slopes"` sits between the two and applies only when something
  predicts the transitions. It shares the covariate slopes across
  occasions while giving each occasion its own intercepts, so the effect
  of a covariate on moving between statuses is held constant over time
  but the underlying rate of movement is not. It costs
  `(n_statuses - 1) * (occasions - 2)` parameters more than `"full"` and
  is nested inside `"none"`, so
  [`lr_test()`](https://pdvalencia.github.io/mixtureEM/reference/lr_test.md)
  tests both restrictions. It requires `predictors_transition` (or a
  `group` acting on the transitions) and
  `transition_effects = "common"`, and with only two occasions it is
  identical to `"full"`.

- forbidden_transitions:

  Transitions that are impossible by design, as in a stage-sequential
  process where people cannot move backwards. Give a \\K \times K\\
  logical or 0/1 matrix, with `TRUE`/`1` marking a forbidden move, or a
  list of one matrix per pair of adjacent occasions. Forbidden cells are
  fixed at zero and do not count as estimated parameters (Collins &
  Lanza, sec. 7.10).

- n_classes:

  Number of latent classes *above* the chain. The default of 1 is the
  ordinary latent transition model. With more, each class gets its own
  initial distribution and its own transition matrices while sharing the
  measurement model - a mixture latent Markov model, which asks whether
  the population contains several distinct processes rather than one.

  Be warned that the unrestricted mixture is demanding of the data. With
  one indicator per occasion it is routinely not identified: on a
  well-known five-wave life-satisfaction panel, even hundreds of random
  starts return a likelihood *worse* than the restricted mover-stayer
  model nested inside it. Several indicators per occasion, or a
  restriction such as `mover_stayer`, is usually what makes the model
  findable. Classes that converge on the same chain are warned about.

- mover_stayer:

  Restrict the **last** latent class to the identity transition matrix:
  a "stayer" class with zero probability of change, the remaining
  classes being "movers" (Vermunt, 2004). Implies `n_classes = 2` unless
  more are asked for, in which case only the last class is a stayer. The
  restricted rows cost no parameters, so the model is nested in the
  unrestricted mixture and
  [`lr_test()`](https://pdvalencia.github.io/mixtureEM/reference/lr_test.md)
  tests it.

- tie_initial_status:

  For a mover-stayer fit, hold the occasion-1 status distribution equal
  across the latent classes instead of estimating one per class. Drops
  the initial-status parameter count from `(K - 1) * C` to `K - 1`.
  Default `FALSE`. With a single latent class there is nothing to tie,
  so it has no effect and a warning says so.

- random_intercept:

  Add a random intercept to the measurement model (Muthen & Asparouhov,
  2022): a person-level "how likely to endorse items in general" trait
  that regular LTA has no way to represent, and that can otherwise be
  mistaken for status separation and for stability over time. `"none"`
  (the default) fits ordinary LTA. `"continuous"` integrates over a
  normally-distributed factor by Gauss-Hermite quadrature
  (`n_quadrature` nodes); `"binary"` instead estimates a small number of
  discrete intercept classes (`n_ri` of them, 2 by default, the model's
  own case). Both require `measurement_invariance = "full"` and binary
  indicators. A random intercept combines with `mover_stayer`, with
  covariates on the initial status or the transitions
  (`predictors_initial`/`predictors_transition`), and with `group`
  (implemented as covariates on the same two, so this is one capability,
  not two); it still does not support `n_classes` \> 1.

  `"continuous"`'s restart search
  (`options(mixtureEM.lta_ri_search = "wide")`, the default) draws wide
  random starts – perturbed item logits and loadings of either sign –
  ranks them with a one-step generalised-EM M-step, and rescores every
  ranked candidate on the full quadrature grid before promoting
  survivors to the full search. On several published benchmark data sets
  it reaches the same or a better optimum, several times faster, than
  the search this package used before, and it escapes the
  inflated-loading local maximum that search could settle in.
  `options(mixtureEM.lta_ri_search = "narrow")` restores the earlier
  search exactly, for reproducing a fit made under it.

  **Do not test a random intercept against regular LTA with
  [`lr_test()`](https://pdvalencia.github.io/mixtureEM/reference/lr_test.md)**:
  the continuous variant puts the null (loading = 0) on the boundary of
  the parameter space, and the binary variant adds a latent class
  variable, so the usual chi-squared reference distribution does not
  apply either way. Compare the two by BIC instead.

  Whether a random intercept is worth adding is a sample-size question
  more than a modelling one. Tseng (2024), for the continuous-indicator
  analogue (RI-LPTA), puts the requirement at upwards of 2,000 cases for
  80% power and 90% coverage at a between-profile separation of d =
  0.75, with more items, occasions or separation lowering that bar; the
  asymmetry that makes trying it worthwhile anyway is that omitting a
  random intercept when one belongs costs a lot (inflated apparent
  separation and stability), while including one when it does not belong
  costs almost nothing (a handful of parameters, and BIC will say so, as
  it does on this package's own benchmark replication of the article's
  example).

  Not built in this release: continuous indicators, a random *slope*,
  correlated residuals across time, or lag-2 dependence - the last of
  which the article itself reports as significant in both of its worked
  examples, so it is a real simplification and not a hypothetical one.
  Standard errors and the post-EM L-BFGS refinement are available for
  both binary and ordinal indicators. A random intercept crossed with
  several latent classes or `mover_stayer` is supported for both
  measurement families.

  `predictors_random_intercept` regresses the factor on covariates. The
  item probabilities the fit reports are integrated over the factor at a
  zero random-intercept mean - the residual grid - not at each case's
  own predicted mean; the per-case means are in
  `random_intercept_scores()$predicted_mean`.

- n_quadrature:

  Number of Gauss-Hermite nodes for `random_intercept = "continuous"`.
  The default of 15 is a starting point to check, not a settled answer,
  the same way `n_init`'s default is a floor: raise it (15-30 nodes is
  the usual working range, more when the loadings are large) and confirm
  the log-likelihood moves by less than 0.01. Make that check over a
  wide range of node counts rather than one step up. When a fit's
  thresholds are extreme, the item response is almost a step function of
  the factor, and the log-likelihood is then not even monotone in the
  number of nodes: two nearby small settings can differ by several
  log-likelihood units while both sit far from the converged value. One
  ordinal five-status fit used in this package's own validation reads
  -16047.3 at 15 nodes and -16040.5 at 20, but settles at -16041.18 only
  from roughly 80 nodes upward. `n_quadrature = 1` is a valid,
  deliberate special case - a single node at 0 with mass 1 - under which
  the model reduces exactly to regular LTA; it is not a model worth
  fitting on its own, but is how the package's own test suite proves the
  node machinery is wired correctly. The node count matters more once
  `predictors_random_intercept` is used: the person-specific reweighting
  is exact only up to whatever quadrature error the fixed grid already
  has.

- n_ri:

  Number of discrete intercept classes for
  `random_intercept = "binary"`. The default of 2 is Muthen &
  Asparouhov's own case.

- layout, id, time, items, item_names, time_labels:

  Data-shape arguments, passed through as in
  [`fit_rmlca()`](https://pdvalencia.github.io/mixtureEM/reference/fit_rmlca.md).

- weights:

  Optional case weights.

- weight_type:

  What the numbers in `weights` mean. `"sampling"` (the default) treats
  them as survey or probability weights and rescales them to sum to the
  number of cases; `"frequency"` treats them as counts of identical
  cases, as in a response-pattern table, and takes their sum as the
  sample size behind AIC and BIC.

- strata, cluster:

  Optional complex-survey design variables. When either is supplied, the
  standard errors become design-based: the case-level scores are
  aggregated to the primary sampling unit within stratum and used in a
  linearization sandwich, which protects inference against clustering.

- n_init:

  Number of random starts. Latent transition models have many local
  maxima; the default of 20 is a floor, not a recommendation. The fit
  reports how many starts reached the solution it kept, and warns when
  that count is 1; the answer there is `n_init = 100`, and a maximum
  that still does not replicate at 100 starts points at the
  specification rather than at the search. See
  [`vignette("estimation")`](https://pdvalencia.github.io/mixtureEM/articles/estimation.md).

- refine:

  Logical. If `TRUE` (default), each start that runs to convergence is
  followed by an L-BFGS climb on the same penalised objective the EM
  steps maximise, and the starts are then ranked on the refined
  log-likelihoods. EM converges to a fixed point of its own surrogate,
  which on a near-flat likelihood ridge can sit measurably short of the
  maximum; the climb steps past it, and never returns a fit worse than
  the one it was given. This is the same refinement
  [`fit_mixture()`](https://pdvalencia.github.io/mixtureEM/reference/fit_mixture.md)
  has always applied. It is a no-op - silently, and the fit is
  unchanged - for models whose free parameters it cannot differentiate:
  covariate models and measurement families other than binary and
  continuous. Mixtures over chains (`n_classes` \> 1) are refined like
  any other supported model.

- refine_from:

  A model fitted by `fit_lta()` on the same data and with the same
  shape, whose solution this fit continues from. No random starts are
  run and `n_init` is ignored, because the search that produced the
  donor already ran and this only carries its winner further - typically
  to a tighter `tol` or a larger `max_iter`. Passing `n_init` alongside
  is an error rather than a silent override.

- max_iter, tol:

  EM iteration limit and relative convergence tolerance. A mixture over
  chains (`n_classes` \> 1) converges much more slowly and defaults to a
  tighter `tol` of 1e-11 and 5000 iterations, since the shared rule is a
  relative one and would otherwise stop the fit mid-climb. The restarts
  are then staged - a short first pass ranks them and only the top tenth
  (floor of 3) run on to convergence - so the tighter rule does not
  multiply the cost of the search. A random intercept with a free
  loading is staged the same way; on binary items its survivors keep the
  ordinary `tol` of 1e-8 and a Newton-type finish takes them to the
  maximum, which reaches the same optimum as the tighter rule several
  times faster, while on ordinal items they run to 1e-11 as above.
  Supplying either argument overrides all of this. (Unlike
  [`fit_mixture()`](https://pdvalencia.github.io/mixtureEM/reference/fit_mixture.md),
  whose EM tolerance is fixed and not user-adjustable, `tol` here is a
  real, respected argument, because a chain mixture converges slowly
  enough that the fixed rule would not do.)

- n_cores:

  Positive integer. Number of processes to spread the random starts
  over. Default `1` (sequential). These are the slowest fits in the
  package and `n_init` is high by necessity, so this is where the
  argument earns the most. Starting values are drawn in this session
  before any fitting begins, so the fit is identical at every `n_cores`.
  `options(mixtureEM.n_cores = )` sets the default for a whole session;
  an argument given here overrides it.

- smoothing:

  How much smoothing to apply to the status prevalences and to the
  transition matrices, expressed as a number of pseudo-cases spread
  evenly over each conditional table. Sparse transition tables otherwise
  collapse onto probabilities of exactly zero, which are awkward to
  interpret and to test. The default of `1` is negligible at any
  realistic sample size; set it to `0` for unsmoothed maximum
  likelihood. It governs the status prevalences, the transition matrices
  and - with `n_classes` \> 1 - the class weights, and it does **not**
  govern the measurement model, whose prior is `bayes_constants`.

  It also does not reach the transitions once anything predicts them.
  With `predictors_transition`, or with a `group` whose `group_effects`
  include the transitions, each transition matrix is the fitted value of
  a multinomial logit rather than a smoothed table of counts, and
  `smoothing` has no effect on it. The initial status prevalences behave
  the same way under `predictors_initial`.

  The mass is one pseudo-case per *conditional table*, so with `K`
  statuses an origin row of a transition matrix carries a `K`th of it.
  It is spread evenly rather than in proportion to how often each
  destination is occupied: a rare origin row shrunk toward the
  destination marginal would be asserting that everyone moves to the
  prevalent status, which is a confident claim to make about a row the
  sample says little about, whereas an even spread is uninformative.

  The cost falls on exactly those rows. On a row with few expected cases
  the prior carries a visible share of the estimate - at most \\\[a /
  (m + a)\](1 - 1/K_a)\\ of it, where \\a\\ is that row's share of the
  mass, \\m\\ the cases expected in it and \\K_a\\ the reachable
  destinations - and the fit says so when that share exceeds five
  percentage points, naming the worst row. The remedy worth reaching for
  first is `transition_invariance = "full"`, which puts every occasion's
  cases behind one pseudo-case; `smoothing = 0.5` simply halves the
  pull. `smoothing = 0` is not a good answer, since it removes the
  protection against transition probabilities of exactly zero that the
  prior is there to give.

- random_state:

  Optional seed for reproducible starts.

- order_by_size:

  Relabel the statuses from most to least prevalent at the first
  occasion. Ignored whenever the labels already carry meaning, which is
  the case when `forbidden_transitions` defines an ordering of stages or
  when covariates index the statuses through their regression
  coefficients.

- standard_errors:

  Compute standard errors for \\\delta\\ and \\\tau\\. `TRUE` (the
  default) uses the outer product of the case-level scores; `FALSE`
  skips them; `"robust"` returns the sandwich estimator, with the
  observed information as its bread. `"robust"` costs a
  finite-difference Hessian – \\2p(p + 1)\\ likelihood evaluations – and
  so takes minutes rather than seconds on a model of any size. It is
  supported for random-intercept fits. Where the model still cannot be
  packed on an unconstrained scale (covariates, or a measurement family
  whose parameters are not all free) `"robust"` falls back silently to
  the default estimator, and
  [`summary()`](https://rdrr.io/r/base/summary.html) then does not
  report robust errors.

- predictors_initial:

  Optional covariates predicting the latent status at the first occasion
  (Collins & Lanza, sec. 8.10.1).

- predictors_transition:

  Optional covariates predicting the transitions between statuses (sec.
  8.10.2).

- predictors_items:

  Optional covariates entering each indicator directly, inside each
  latent status - a test of measurement invariance with respect to those
  covariates. One proportional-odds slope per latent status per item per
  covariate, shared across occasions, so the cost is
  `n_statuses * n_items * ncol(predictors_items)` parameters. **This is
  not `group`**: `group` gives every group its own status prevalences
  and transition matrices while the measurement model stays invariant
  across groups, which is the assumption this argument exists to relax.
  Requires `measurement = "ordinal"` and
  `measurement_invariance = "full"`. Restricted to covariates taking few
  distinct values - see `options(mixtureEM.dif_max_patterns = )`. The
  item probabilities the fit reports are those of a case with every one
  of these covariates at zero.

  Alternatively a named list, item name = a data frame (or a matrix with
  column names) of the covariates acting on that item, e.g.
  `list(item3 = data.frame(female))`, which frees slopes on the named
  items only. Items are named as the fit names them (occasion 1's column
  names with the occasion marker dropped), and a covariate is matched
  across items by its column name. Each slope is then shared by every
  status (uniform DIF, one parameter) unless the item is named in
  `predictors_items_by_status`.

- predictors_items_by_status:

  Names of items in the list form of `predictors_items` whose slopes
  differ by status (non-uniform DIF, one slope per status). The matrix
  form is already by status.

- predictors_random_intercept:

  Optional covariates predicting the continuous random intercept itself
  (the article's own Step 5). The factor's residual variance is fixed at
  1 and its residual mean at 0, so no intercept is estimated and a
  constant column is refused: the coefficients are the regression of the
  stable trait on the covariates, and each costs one parameter. Requires
  `random_intercept = "continuous"`. The factor's orientation is
  arbitrary - the fit pins it by making the largest loading positive -
  so the sign of every coefficient is meaningful only relative to the
  loadings.

- n_steps:

  How the measurement model and the structural model are estimated
  relative to each other. `1` (default) fits both at once: the
  covariates and the latent statuses are estimated in one likelihood, so
  the covariates help decide what the statuses are. `2` is the two-step
  estimator (Bakk & Kuha, 2018; Bartolucci, Montanari & Pandolfi, 2015):
  the measurement model is fitted first on the indicators alone, held
  fixed, and `predictors_initial` / `predictors_transition` /
  `predictors_random_intercept` are then fitted on the full likelihood
  with the statuses no longer free to move. `3` is the bias-adjusted
  three-step estimator (Vermunt, 2010; Nylund-Gibson et al., 2014): the
  statuses are estimated, assigned, and the transitions are then
  estimated from the assignments with their classification error
  corrected for. See the two sections below.

- transition_effects:

  How covariates act on the transitions. `"common"` (default) gives each
  origin status its own intercepts but one slope per covariate shared
  across origins, which is the specification in Wang & Wang eq. 6.28 and
  by far the better-behaved one. `"by_origin"` fits a separate
  regression per origin status, which is saturated and often fails to
  converge usefully when transitions are sparse.

- group:

  Optional grouping variable for a multiple-group model (Collins &
  Lanza, sec. 8.2-8.3). It is entered as dummy predictors, saturated
  over the transition rows, which gives each group its own status
  prevalences and its own transition matrices while the measurement
  model stays invariant across groups.

- group_effects:

  Which parameters the grouping variable is allowed to shift: `"both"`
  (default), `"initial"`, `"transitions"` or `"none"`. Fitting the same
  data under two of these and comparing them with
  [`lr_test()`](https://pdvalencia.github.io/mixtureEM/reference/lr_test.md)
  gives the group-difference tests of sec. 8.6-8.8.

- bayes_constants:

  Optional named list of prior strengths for the *measurement* model
  (`categorical`, `poisson`, `variances`); see
  [`fit_mixture()`](https://pdvalencia.github.io/mixtureEM/reference/fit_mixture.md).
  The status and transition probabilities are governed by `smoothing`
  instead, so `latent` is not read here.

- correction:

  Three-step only. `"ML"` (default; Vermunt, 2010) holds each occasion's
  classification-error matrix fixed while the transitions are estimated,
  which is what removes the bias. `"none"` sets every matrix to the
  identity, which is the naive classify-analyse estimate: useful as the
  baseline the correction is measured against, and biased towards
  whatever the classification gets wrong. `"BCH"` (Bolck, Croon and
  Hagenaars, 2004) reweights instead: each case is spread over the
  status paths by the product, over occasions, of the rows of the
  inverted error tables for the statuses it was assigned, and step 3
  fits the initial status and the transitions to those weighted paths.
  Some weights are negative and are kept. It uses modal assignment,
  reports the BCH-weighted log-likelihood, and its standard errors are
  the case-clustered sandwich; it does not yet support `distal`,
  `predictors_items`, `strata` or `cluster`.

- assignment:

  Three-step only. How each occasion's posteriors become the assigned
  status step 3 reads. `"proportional"` (default) spreads every case
  over the statuses in proportion to its posterior; `"modal"` assigns
  each case to its most likely status. The same argument, with the same
  default and for the same reason, as in
  [`add_covariates()`](https://pdvalencia.github.io/mixtureEM/reference/add_covariates.md).

- distal:

  Three-step only, with `assignment = "modal"`. A distal outcome
  measured after the last occasion: a vector or a data frame with one
  row per case, in the rows of `indicators`. Each column is modelled by
  the status at the last occasion – a 0/1, logical or two-level factor
  column as binary with a status-specific probability, any other numeric
  column as Gaussian with a status-specific mean and variance – and is
  estimated jointly with the structural model in step 3, so the outcome
  informs the last occasion's posteriors and the transitions (Vermunt,
  2010). A missing value contributes nothing to that case's likelihood.
  The estimates, with standard errors that carry step 1's uncertainty,
  are returned as `distal`;
  [`outcome_contrasts()`](https://pdvalencia.github.io/mixtureEM/reference/outcome_contrasts.md)
  gives the pairwise status differences. `NULL` (the default) fits no
  outcome.

- ...:

  Ignored.

## Value

An object of class `"lta_model"` with components including `delta`,
`tau` (a list of transition matrices), `prevalences` (status prevalence
by occasion), `gamma` (posterior status probabilities by occasion),
`mm`, `metrics` and `n_params`.

Two diagnostics come with it. `boundary` lists the transition cells the
data have driven to zero, and `smoothing_influence` gives, for each
origin row, the cases expected in it and the share of the estimate the
`smoothing` prior is carrying. `smoothing_influence` is `NULL` when
`smoothing` is 0 or when covariates predict the transitions, since the
prior does not reach them then.

With `n_classes` \> 1 those parameters gain a class index: `delta`
becomes a classes-by-statuses matrix, `tau` a list of per-class lists,
and `class_weights`, `class_posterior` and `gamma_by_class` are added.
`prevalences` stays the whole-sample marginal, with the per-class ones
in its `"by_class"` attribute;
[`status_prevalences()`](https://pdvalencia.github.io/mixtureEM/reference/status_prevalences.md)
and
[`transition_matrix()`](https://pdvalencia.github.io/mixtureEM/reference/transition_matrix.md)
take a `class` argument to reach them. Standard errors are not available
for a mixture over chains and `se` is `NULL`.

`metrics$entropy` (the headline number, printed by
[`print()`](https://rdrr.io/r/base/print.html)) is the relative entropy
of the joint latent-status path across all occasions at once: one
classification per case over the `K^T` possible paths, normalised by
`n * T * log(K)`. The per-occasion alternative – the relative entropy of
each occasion's status posterior on its own, normalised by `log(K)` – is
`metrics$entropy_by_occasion` (printed by
[`summary()`](https://rdrr.io/r/base/summary.html)). A per-occasion
entropy can also be normalised by the entropy of that occasion's
estimated status proportions instead of by `log(K)`; the two are related
by `1 - (1 - r2) * log(K) / H(pi_hat)`, where `r2` is the package's
value, `K` is `n_statuses`, and `H(pi_hat)` is the entropy of that
occasion's estimated status proportions (from
[`status_prevalences()`](https://pdvalencia.github.io/mixtureEM/reference/status_prevalences.md)).

## The two-step estimator (`n_steps = 2`)

A latent transition model is two models stacked: a measurement model
saying what the statuses are, and a structural model saying who starts
in which one and who moves. Fitting both at once (`n_steps = 1`) lets
the covariates take part in defining the statuses, so adding or dropping
a covariate can change what the statuses mean. The two-step estimator
removes that: step 1 fits the measurement model on the indicators alone,
step 2 holds it fixed and maximises the same full likelihood over the
structural coefficients only. Each E-step still runs on the joint model,
so this is one-step estimation with one block pinned – there is no class
assignment and no classification table anywhere in it.

What that buys: the statuses are fixed before any covariate is looked
at, so several structural models can be compared on one measurement
model and none of them can redefine it. What it costs: like every
stepwise estimator, the coefficients are biased towards zero when the
statuses are poorly separated or the sample is small (Bakk & Kuha, 2018,
Tables 1 and 3). With well-separated statuses the two estimators agree
closely.

Step 1 is this same call with the structural predictors dropped, so it
takes the same `n_init`, `random_state`, invariance constraints and
priors; it is returned on the fitted object as `$step1`, and the item
parameters of the returned fit are identical to its. Step 2 runs no
restarts of its own – with the measurement block fixed there is nothing
left for a restart to search.

Standard errors for the structural coefficients carry the step-1
uncertainty. Step 2 holds the measurement parameters fixed, but they
were estimated rather than known, and their sampling error propagates
into everything built on top of them. What is reported is therefore the
pseudo-maximum-likelihood variance of Bakk & Kuha's equation 5, \\V =
V_2 + V_1\\: \\V_2\\ is the inverse observed information of the full
likelihood in the structural coefficients, and \\V_1\\ carries the
step-1 variance across through the cross-curvature between the two
blocks. \\V_1\\ is negligible when the statuses are well separated and
is most of the variance when they are not, which is the same condition
that governs the attenuation above. The measurement parameters of a
two-step fit report the step-1 fit's own standard errors, since that is
the fit that estimated them.

Two costs come with it. The information matrices are differenced
numerically, which is a few thousand likelihood evaluations and can take
minutes on a large model; pass `standard_errors = FALSE` to skip it. And
`standard_errors = "robust"` is the sandwich on the case-level scores,
which is a step-2-only estimator and does not include \\V_1\\.
`n_steps = 1` is unaffected by any of this.

## The three-step estimator (`n_steps = 3`)

Step 1 fits the measurement model with no transitions at all: each
occasion gets its own status prevalences and nothing links one
occasion's status to the next. That is deliberate. A model with
transitions smooths each occasion's posterior along the chain, so the
status assigned at occasion 2 would carry information about occasions 1
and 3, and step 3 could no longer treat it as a lone, imperfect
measurement of the status at occasion 2. `measurement_invariance`
applies to step 1 as usual.

Step 2 assigns a status at every occasion (see `assignment`) and
records, for each occasion separately, how often an assigned status
differs from the true one. The matrices differ by occasion because the
error depends on the base rates, and in a transition model the base
rates move: on the
[`ecls_reading`](https://pdvalencia.github.io/mixtureEM/reference/ecls_reading.md)
panel one status holds 1.8% of the children at the first occasion and
81% at the last, and it is misclassified two times in five at the first
occasion against one in fifty averaged over all four. A cell the data
leave empty is held at 1e-6 rather than zero, so no transition is ruled
out by a classification table.

Step 3 is a latent transition model with one indicator per occasion, the
assigned status, whose response probabilities are those matrices held
fixed. It takes `predictors_initial`, `predictors_transition`,
`transition_effects`, `transition_invariance` and
`forbidden_transitions` exactly as a one-step fit does. The returned
object is the step-3 fit; `$step1` is step 1, and `$threestep` holds the
error matrices (rows the assigned status, columns the true one), the
step-1 prevalences by occasion and the modal assignments.

The step-3 fit's `loglik` and `metrics` are the likelihood of the
assigned statuses, not of the items, and are on a different scale from a
one- or two-step fit's. [`print()`](https://rdrr.io/r/base/print.html)
and
[`compare_longitudinal()`](https://pdvalencia.github.io/mixtureEM/reference/compare_longitudinal.md)
therefore report step 1's criteria, labelled as such, and
[`lr_test()`](https://pdvalencia.github.io/mixtureEM/reference/lr_test.md)
refuses a three-step fit beside a one- or two-step one.

With `measurement_invariance = "none"` and no transitions, nothing in
step 1 says which status at occasion 2 is the same as a given status at
occasion 1. The labels are matched: status \\k\\ at every occasion is
the one whose item profile is closest, in summed squared distance, to
status \\k\\'s at occasion 1. `$threestep$alignment` records the match.
With `transition_invariance = "none"` step 3's likelihood does not
depend on the labels, only its reading does; with `"full"` or `"slopes"`
it does, so inspect the step-1 profiles before relying on those.

Step 3's standard errors carry step 1's uncertainty. The error matrices
are computed from step 1's estimates, which have sampling error of their
own, so treating them as known understates every standard error in step
3. The reported variance is the first-order pseudo-maximum-likelihood
one of Bakk, Oberski and Vermunt (2014): step 3's own inverse
information plus the variance of step 1's estimates, carried through the
error matrices (and, under proportional assignment, through the weights
the posteriors give the reduced data). The step-3-only part is kept as
`$se$threestep_V2` and step 1's variance as `$se$step1_vcov`. Under
`correction = "none"` with modal assignment step 1 never reaches step 3,
and the standard errors are step 3's own.

With `strata` or `cluster` and `assignment = "modal"`, the variance is
the design-based sandwich over both steps at once: each case's scores in
step 1 and in step 3 are combined through the same first-order term,
summed within primary sampling units and compared across them within
strata. Clustering therefore widens the step-1 part as well as step 3's
own, and any correlation between the two steps' scores within a cluster
is counted rather than assumed away. Proportional assignment refuses a
design, because its step-3 rows are status combinations, not cases.

Two more things to know when reading the result. Step 3's log-likelihood
is that of the assigned statuses, not of the items, so it must never be
compared with a one- or two-step fit's. And the estimator covers one
chain of statuses on a measurement model every case shares: random
intercepts, `n_classes > 1`, `mover_stayer` and `group` are refused.

`predictors_items` is accepted and acts in step 1 only, where the status
at every occasion is also regressed on the same covariates (without that
regression the item slope would have to carry the whole association
between the covariate and the status). Step 3 fits whatever
`predictors_initial` and `predictors_transition` ask for. Step 1 has no
standard errors in this case, so step 3's treat the classification error
as known rather than adding step 1's uncertainty.

## References

Collins, L. M., & Lanza, S. T. (2010). *Latent Class and Latent
Transition Analysis: With Applications in the Social, Behavioral, and
Health Sciences*. Wiley (chapters 7-8).

Nylund-Gibson, K., Grimm, R., Quirk, M., & Furlong, M. (2014). A latent
transition mixture model using the three-step specification. *Structural
Equation Modeling*, *21*(3), 439-454.
[doi:10.1080/10705511.2014.915375](https://doi.org/10.1080/10705511.2014.915375)

Chung, H., Lanza, S. T., & Loken, E. (2008). Latent transition analysis:
inference and estimation. *Statistics in Medicine*, *27*(11), 1834-1854.
[doi:10.1002/sim.3130](https://doi.org/10.1002/sim.3130)

Fienberg, S. E., & Holland, P. W. (1973). Simultaneous estimation of
multinomial cell probabilities. *Journal of the American Statistical
Association*, *68*(343), 683-691.
[doi:10.1080/01621459.1973.10481405](https://doi.org/10.1080/01621459.1973.10481405)

Muthen, B., & Asparouhov, T. (2022). Latent transition analysis with
random intercepts (RI-LTA). *Psychological Methods*, *27*(1), 1-16.
[doi:10.1037/met0000370](https://doi.org/10.1037/met0000370)

Bakk, Z., & Kuha, J. (2018). Two-step estimation of models between
latent classes and external variables. *Psychometrika*, *83*(4),
871-892.
[doi:10.1007/s11336-017-9592-7](https://doi.org/10.1007/s11336-017-9592-7) -
the estimator behind `n_steps = 2`.

Bartolucci, F., Montanari, G. E., & Pandolfi, S. (2015). Three-step
estimation of latent Markov models with covariates. *Computational
Statistics & Data Analysis*, *83*, 287-301.
[doi:10.1016/j.csda.2014.10.017](https://doi.org/10.1016/j.csda.2014.10.017) -
the longitudinal case. Their step 1 pools the occasions into one
cross-sectional latent class model; the step 1 here is the latent
transition model's own measurement block, which already carries whatever
invariance constraints the fit asks for. Both hold the measurement
parameters fixed and maximise the full likelihood over the structural
ones.

Vermunt, J. K. (2010). Latent class modeling with covariates: Two
improved three-step approaches. *Political Analysis*, *18*(4), 450-469.
[doi:10.1093/pan/mpq025](https://doi.org/10.1093/pan/mpq025) - the
correction behind `n_steps = 3`.

Bakk, Z., Oberski, D. L., & Vermunt, J. K. (2014). Relating latent class
assignments to external variables: Standard errors for correct
inference. *Political Analysis*, *22*(4), 520–540.
[doi:10.1093/pan/mpu003](https://doi.org/10.1093/pan/mpu003)

- the standard errors of `n_steps = 3`.

Tseng, M.-C. (2024). Latent profile transition analysis with random
intercepts (RI-LPTA). *Structural Equation Modeling*, *31*(4), 626-634.
[doi:10.1080/10705511.2023.2284671](https://doi.org/10.1080/10705511.2023.2284671) -
the sample-size analysis behind the guidance above.

## See also

[`transition_matrix()`](https://pdvalencia.github.io/mixtureEM/reference/transition_matrix.md),
[`status_prevalences()`](https://pdvalencia.github.io/mixtureEM/reference/status_prevalences.md),
[`lr_test()`](https://pdvalencia.github.io/mixtureEM/reference/lr_test.md),
[`lta_g2()`](https://pdvalencia.github.io/mixtureEM/reference/lta_g2.md),
[`fit_rmlca()`](https://pdvalencia.github.io/mixtureEM/reference/fit_rmlca.md).
