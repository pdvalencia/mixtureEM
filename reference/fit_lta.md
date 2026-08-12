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
  `transition_invariance = "full"`;

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
  transition_invariance = c("none", "full"),
  forbidden_transitions = NULL,
  n_classes = 1,
  mover_stayer = FALSE,
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
  max_iter = 1000,
  tol = 1e-08,
  smoothing = 1,
  random_state = NULL,
  order_by_size = TRUE,
  standard_errors = TRUE,
  predictors_initial = NULL,
  predictors_transition = NULL,
  transition_effects = c("common", "by_origin"),
  group = NULL,
  group_effects = c("both", "initial", "transitions", "none"),
  bayes_constants = NULL,
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
  `"categorical"`, `"continuous"`, or a named list for a mixed block.

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

- max_iter, tol:

  EM iteration limit and relative convergence tolerance. A mixture over
  chains (`n_classes` \> 1) converges much more slowly and defaults to a
  tighter `tol` of 1e-11 and 5000 iterations, since the shared rule is a
  relative one and would otherwise stop the fit mid-climb. The restarts
  are then staged - a short first pass ranks them and only the best
  three run on to convergence - so the tighter rule does not multiply
  the cost of the search. Supplying either argument overrides all of
  this.

- smoothing:

  How much smoothing to apply to the status prevalences and to each row
  of the transition matrices, expressed as a number of pseudo-cases
  spread evenly over the possible destinations. Sparse transition tables
  otherwise collapse onto probabilities of exactly zero, which are
  awkward to interpret and to test. The default of `1` is negligible at
  any realistic sample size; set it to `0` for unsmoothed maximum
  likelihood.

- random_state:

  Optional seed for reproducible starts.

- order_by_size:

  Relabel the statuses from most to least prevalent at the first
  occasion. Ignored whenever the labels already carry meaning, which is
  the case when `forbidden_transitions` defines an ordering of stages or
  when covariates index the statuses through their regression
  coefficients.

- standard_errors:

  Compute standard errors for \\\delta\\ and \\\tau\\.

- predictors_initial:

  Optional covariates predicting the latent status at the first occasion
  (Collins & Lanza, sec. 8.10.1).

- predictors_transition:

  Optional covariates predicting the transitions between statuses (sec.
  8.10.2).

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

- ...:

  Ignored.

## Value

An object of class `"lta_model"` with components including `delta`,
`tau` (a list of transition matrices), `prevalences` (status prevalence
by occasion), `gamma` (posterior status probabilities by occasion),
`mm`, `metrics` and `n_params`.

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

## References

Collins, L. M., & Lanza, S. T. (2010). *Latent Class and Latent
Transition Analysis: With Applications in the Social, Behavioral, and
Health Sciences*. Wiley (chapters 7-8).

Nylund-Gibson, K., Grimm, R., Quirk, M., & Furlong, M. (2014). A latent
transition mixture model using the three-step specification. *Structural
Equation Modeling*, *21*(3), 439-454.
[doi:10.1080/10705511.2014.915375](https://doi.org/10.1080/10705511.2014.915375)

## See also

[`transition_matrix()`](https://pdvalencia.github.io/mixtureEM/reference/transition_matrix.md),
[`status_prevalences()`](https://pdvalencia.github.io/mixtureEM/reference/status_prevalences.md),
[`lr_test()`](https://pdvalencia.github.io/mixtureEM/reference/lr_test.md),
[`lta_g2()`](https://pdvalencia.github.io/mixtureEM/reference/lta_g2.md),
[`fit_rmlca()`](https://pdvalencia.github.io/mixtureEM/reference/fit_rmlca.md).
