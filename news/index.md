# Changelog

## mixtureEM (development version)

### Fixed: two estimation bugs that cost log-likelihood on every affected fit

- **EM now converges before it stops.** Emissions that L-BFGS refines
  afterwards — the binary and Gaussian families — used a much looser
  stopping rule than the others, on the reasoning that the refinement
  would climb the rest of the way. It does not. On a validated
  four-class binary LCA the loose rule stopped 1.27 log-likelihood units
  short and the refinement recovered 0.14 of that. The rule also stopped
  EM at an absolute change of roughly two units on a log-likelihood in
  the thousands, which is far too coarse to *rank* restarts, so it was
  quietly degrading the multi-start search as well. Every emission now
  uses the tight rule, and the refinement starts from a converged fit
  instead of substituting for one. Fits will change, generally for the
  better, and will take more iterations.

- **`bayes_constants` now reaches the L-BFGS refinement.** The
  refinement read `variances` but had `latent` and `categorical`
  hard-coded at `1`, so setting either to `0` switched the prior off in
  the M-step and left it on in the polish. The two stages then optimised
  different objectives and the polish pulled the fit off the
  maximum-likelihood optimum it had been asked for. The documented
  escape hatch for reproducing an unregularized reference analysis now
  works. The analytical gradient is verified against a finite-difference
  one to 1e-10 for complete data and under FIML, at three prior
  settings.

- **The L-BFGS refinement no longer runs when class membership is
  modelled by a regression.** Its parameterisation packs a single pooled
  vector of class weights and has no slot for `P(class | covariates)`,
  so with `predictors` or a `group` on the prevalences it was maximising
  a different model from the one being fitted and then writing its
  measurement parameters back. The damage was large and had been
  invisible: a multiple-group configural model is *separable*, so its
  log-likelihood must equal the sum of the per-group fits, and the
  polish left it 1.10 units short. Removing it closes that to 0.0005 and
  takes the number of restarts reaching the best solution from 1 of 6 to
  4 of 6 — it had been perturbing every restart away from the optimum,
  not just the winner.

Together these close a 3.5-unit gap against an external reference on a
three-group latent class model — 5.1 units on the configural model —
both now matched to within 0.001.

### New: binary indicators are recoded to 0/1 for you

- **`measurement = "binary"` now accepts any two-valued item.** A
  two-level factor or character, a logical, and any numeric pair — 1/2,
  2/5, whatever the source data used — are converted on the way in, and
  the fit is identical to the one hand-coded 0/1 data would have given.
  Only the mapping is announced, once, and
  [`measurement_summary()`](https://pdvalencia.github.io/mixtureEM/reference/measurement_summary.md)
  then names the level each probability belongs to, so “0.87” is never
  ambiguous. Data already in 0/1 is untouched and silent.

- **Two silent-wrong-answer bugs closed on the way.** The old `{0, 1}`
  check was reachable only through a single-string `measurement`, so a
  1/2-coded item inside a mixed
  [`list()`](https://rdrr.io/r/base/list.html) specification reached the
  Bernoulli likelihood unchecked and returned a wrong log-likelihood
  with no error. And a *categorical* item coded from 0 indexed the
  previous item’s last category rather than its own — no error, wrong
  likelihood. Both are now caught, the second with a message saying to
  add 1.

- An item with three or more values is still refused, but the message
  now points at `"categorical"` or `"continuous"` instead of asking for
  0/1 coding that would not have helped.

### Renamed: `longitudinal_lrt()` is now `lr_test()`

- The test was never specific to longitudinal models. It takes any two
  nested fits and differences their log-likelihoods and parameter
  counts, and most uses of it — including the multiple-group
  measurement-invariance test — are cross-sectional.
  [`longitudinal_lrt()`](https://pdvalencia.github.io/mixtureEM/reference/longitudinal_lrt.md)
  still works and warns.

### Changed: what the collapsed-variance warning tells you to do

- The warning used to lead with `bayes_constants = list(variances = 5)`.
  That number was calibrated on one dataset and does not transfer,
  because the constant is divided among the classes: the same value is a
  different amount of prior at every `n_classes`. The warning now leads
  with `variances_equal = TRUE`, which bounds the likelihood so the
  degeneracy cannot arise, then fewer classes, and only then the prior —
  stated as roughly one artificial observation per class and printed as
  the concrete number for the model in hand.

- It also now says two things it should have said before: raising
  `n_init` is **not** a remedy and can make matters worse, since the
  likelihood really is unbounded in that direction and a longer search
  finds a taller spike; and a flagged fit’s BIC must not be compared
  with a clean fit’s, because it is inflated by the spike. Finally, it
  asks for a substantive check — that the flagged variance is no longer
  far below the others and that its class mean has come off the floor or
  ceiling — rather than just that the warning stopped.

### New: measurement invariance one parameter at a time

- **New `group_invariant_params` argument on
  [`fit_mixture()`](https://pdvalencia.github.io/mixtureEM/reference/fit_mixture.md)**,
  for continuous indicators. `group_invariant_items` holds whole *items*
  equal across groups, which is the natural constraint when an item has
  one kind of parameter, as a categorical indicator does. A continuous
  indicator has two — a class mean and a variance — and the
  latent-profile invariance literature routinely frees one and holds the
  other. That model could not be written down before: an item was either
  wholly free across groups or wholly shared.

  `group_invariant_params = "covariances"` frees the class means across
  groups while holding the indicator variances invariant. This is the
  model Olivera-Aguilar and Rikoon (2018) call *unconstrained* and note
  is the default most software fits, and it is the one their invariance
  test compares against — so it is the comparison an applied analysis
  usually wants, and it is smaller and better identified than the fully
  heterogeneous alternative. `group_invariant_params = "means"` is the
  mirror constraint. The two axes are alternatives: pass either
  `group_invariant_items` or `group_invariant_params`, not both.

- **New `variances_equal` argument on
  [`fit_mixture()`](https://pdvalencia.github.io/mixtureEM/reference/fit_mixture.md)**,
  also for continuous indicators: hold each item’s variance equal across
  the classes, so the classes differ in location only. This is the
  homoscedastic latent profile model, and the parameterisation several
  commercial programs estimate by default. It applies to an ordinary
  single-group fit as well, and composes with `group_invariant_params`
  to give a variance shared by the classes but free across groups.

  The constrained variance is still stored once per class, so profiles,
  plots, class alignment and the degeneracy check are unchanged; the
  constraint lives in the M-step and in the parameter count.

  A model carrying either constraint is fitted by EM alone. The L-BFGS
  refinement packs one parameter per class-item cell and has no way to
  express an equality across cells, so it would step off the constraint
  surface; these models instead run EM to the tighter stopping rule the
  unpolished emissions already use.

### Improved: the search for a group-varying measurement model

- **`group_effects = "both"` and `"measurement"` no longer rely on
  random starting values alone.** A group-varying measurement model
  holds one set of item parameters per group, tied together only through
  the class labels, which are shared. Random starts give each group its
  own arbitrary labelling, so class 1 in one group and class 1 in
  another begin as unrelated things and EM has no way to bring them into
  correspondence — and more restarts do not help, since every restart
  draws from the same badly aligned prior. The visible symptom was a
  model that could score *below* the more restricted model nested inside
  it, which is impossible at the optimum and left
  [`longitudinal_lrt()`](https://pdvalencia.github.io/mixtureEM/reference/longitudinal_lrt.md)
  reporting a lower bound instead of a test.

  These models now fit each group on its own data first, permute its
  classes to match the pooled solution, and run one extra restart from
  there. On the `janousch` example the configural fit improves by 15
  log-likelihood units and the measurement-invariance statistic rises
  from 251 to 282 on 112 degrees of freedom; the substantive conclusion
  is unchanged, but the number no longer depends on the draw. The extra
  fits make these models roughly a third slower.

  Everything about it degrades gracefully: a group with fewer than
  `2 * n_classes` cases is not fitted on its own, and any sub-fit that
  fails falls back to the pooled parameters, so the worst case is the
  search as it was before.

### Fixed: collapsed class variances in continuous-indicator models

- **A class variance that has collapsed towards zero is now detected and
  reported.** A mixture of normals with freely estimated variances has
  an unbounded likelihood: a class whose variance on one item is driven
  to zero, with its mean on a few cases sharing a value, produces a
  likelihood that grows without limit while describing nothing. Such a
  solution previously won the restart competition, was returned as the
  best fit, and gave no warning. It now raises one naming the item, the
  class and the group, and reporting whether the class mean sits at a
  data boundary. The flagged cells are stored on the fitted object as
  `$degenerate` and repeated by
  [`print()`](https://rdrr.io/r/base/print.html) and
  [`summary()`](https://rdrr.io/r/base/summary.html), so a saved fit
  still shows the problem months later.

  **This changes results.** A model that previously reported a collapsed
  solution as converged will now warn, and its estimates will differ,
  because the regularisation it is fitted under has changed (see below).
  The `janousch` vignette’s measurement-invariance section is the worked
  example and its numbers have changed accordingly.

- **The Gaussian M-step now carries a prior on the variances instead of
  a hard floor.** `gaussian_diag` and `gaussian_diag_nan` previously
  added `1e-6` to each variance after maximising. That is a constant in
  the units of the data — invisible on a five-point scale, enormous on
  one measured in thousands — and, being applied after the fact, did
  nothing to stop the L-BFGS refinement from re-optimising straight
  through it. The regularisation is now part of the objective: a
  truncated inverse-Wishart prior centred on the item’s observed
  marginal variance, expressed as pseudo-observations, applied
  identically in the M-step and in the refinement so the two cannot
  disagree about what is being maximised. On healthy fits the difference
  is negligible.

- **New `bayes_constants` argument** on
  [`fit_mixture()`](https://pdvalencia.github.io/mixtureEM/reference/fit_mixture.md),
  [`fit_lta()`](https://pdvalencia.github.io/mixtureEM/reference/fit_lta.md)
  and the wrappers that delegate to them, exposing the strength of each
  of the estimator’s priors — `latent`, `categorical`, `poisson`,
  `variances` — all defaulting to `1`, the value the first three were
  already hard-coded to. Raising `variances` is one of the remedies the
  collapsed-variance warning offers; setting a constant to `0` removes
  that prior and recovers plain maximum likelihood for that block, which
  is an escape hatch for reproducing an unregularized reference analysis
  rather than a recommended setting.

- **[`longitudinal_lrt()`](https://pdvalencia.github.io/mixtureEM/reference/longitudinal_lrt.md)
  refuses to interpret a degenerate fit.** It warns when either input
  carries the flag — the statistic is meaningless in either direction,
  since a degenerate log-likelihood reflects a spike in the likelihood
  rather than fit to the data. Its existing warning for a negative
  statistic now names both things that can cause one: a failed search on
  the full model, or a degenerate restricted model outscoring it.

- **The `janousch` measurement-invariance comparison was testing the
  wrong restriction.** Measurement invariance is the hypothesis that the
  item parameters are equal across groups, which is
  `group_effects = "prevalence"` (measurement pooled, prevalences free)
  against `"both"`. The vignette and test compared `"measurement"`
  against `"both"`, which frees the item parameters and pools the
  prevalences — the opposite restriction. Both have been rewritten.

  This supersedes the note in the 0.2.0-era commit that introduced a
  fixed seed for that comparison: the seed was pinning a degenerate
  configural solution rather than curing a search failure, and it has
  been removed.

## mixtureEM 0.2.0

### New: stepwise analyses on a fitted model

- **[`add_covariates()`](https://pdvalencia.github.io/mixtureEM/reference/add_covariates.md)
  and
  [`add_outcome()`](https://pdvalencia.github.io/mixtureEM/reference/add_outcome.md).**
  Once you have chosen an unconditional model, relating its classes to
  external variables no longer requires re-specifying (and
  re-estimating) the whole measurement model:

  ``` r

  fit <- fit_mixture(items, n_classes = 3)
  summary(add_covariates(fit, data[, c("age", "sex")]))
  summary(add_outcome(fit, data$bmi))
  ```

  Both verbs run only steps 2–3 of the bias-adjusted three-step approach
  (Bakk et al., 2014; Vermunt, 2010) on the stored step-1 solution, so
  the classes you inspected are exactly the classes the structural model
  describes — same solution, same class order, and no risk of a re-run
  landing in a different optimum. `fit_mixture(predictors = )` and
  `fit_mixture(outcome = )` remain available for one-call workflows and
  one-step (simultaneous) estimation.

### New: summary functions return their numbers

- [`measurement_summary()`](https://pdvalencia.github.io/mixtureEM/reference/measurement_summary.md)
  now returns (invisibly) a long-format data frame of the item
  parameters — `item`, `category`, `class`, `estimate` — alongside the
  printed tables, so results can be reused programmatically without
  touching the fitted object’s internals.
- [`summary()`](https://rdrr.io/r/base/summary.html) on a model with a
  structural part returns (invisibly) a list of tidy tables:
  `$coefficients` (class contrasts with SEs, p-values, odds ratios and
  confidence limits), `$omnibus` (per-covariate Wald tests), and
  `$outcome` (predicted probabilities or class means with their tests).
- New
  [`class_sizes()`](https://pdvalencia.github.io/mixtureEM/reference/class_sizes.md):
  one row per class with the model-estimated proportion, the expected
  count, and the modal-assignment count — the three numbers applied
  papers report.

### Output and plotting improvements

- Long variable and category names no longer break table alignment.
  Label columns size themselves to their content, and names past ~28
  characters are abbreviated in the display with a key printed under the
  table (the returned data frames always carry the full names).
- Profile and trajectory plots no longer clip long class labels on the
  right: margins are computed from the actual label widths, the legend
  position scales with the x-range, and overly long labels are
  shortened. The bottom margin now adapts to long rotated indicator
  names as well.
- Factor covariates: levels with no observations (e.g. after subsetting
  a pooled dataset) are dropped automatically instead of producing a
  degenerate `OR = 1.00 [1.00, 1.00]` row, and levels observed fewer
  than 5 times trigger a warning that their coefficients will be
  unstable.

### Breaking changes

- **Standard errors for covariate effects in 2- and 3-step models have
  changed, and every p-value and confidence interval for a class
  predictor moves with them.** They were computed from the Hessian the
  M-step stores, which measures the curvature of the *Q function* rather
  than of the step-3 log-likelihood the estimates actually maximise, and
  which takes no account of the fact that proportional assignment gives
  each case K weighted records. Both errors run the same way, so the
  reported standard errors were too small and the intervals too narrow.

  The new default, `se = "corrected"` on
  [`fit_mixture()`](https://pdvalencia.github.io/mixtureEM/reference/fit_mixture.md),
  is the first-order corrected estimator of Bakk, Oberski and Vermunt
  (2014): the step-3 sandwich plus the variance carried over from step
  1, which is estimated and not known. `se = "robust"` gives the
  sandwich alone; `se = "hessian"` inverts the step-3 observed
  information. See
  [`?covariate_se`](https://pdvalencia.github.io/mixtureEM/reference/covariate_se.md).

  A 250-replication coverage study on the design of Bakk et al.
  (`data-raw/covariate_se_simulation.R`) says how much this mattered.
  With n = 500 and entropy R² = 0.64, the old nominal 95% intervals
  covered between **54% and 88%** of the time, worst for the strongest
  coefficient; the new default covers 92–98%. With entropy R² near 0.90
  the old intervals covered 87–95% and the new ones 94–98%, and there
  the correction term is negligible — the sandwich alone does the work.
  **If you have published covariate p-values from this package with
  small or poorly separated classes, they were anti-conservative and are
  worth re-running.**

  [`summary()`](https://rdrr.io/r/base/summary.html),
  [`confint()`](https://rdrr.io/r/stats/confint.html) and
  [`analytical_wald_test()`](https://pdvalencia.github.io/mixtureEM/reference/analytical_wald_test.md)
  now print the name of the estimator behind the numbers they show.
  `correction = "BCH"`, one-step fits, and a covariate combined with a
  distal outcome in one structural model keep the uncorrected Hessian,
  and say so.

### New features

- **[`summary()`](https://rdrr.io/r/base/summary.html) prints an omnibus
  Wald test for each class predictor.** For a covariate with more than
  two categories, none of the per-dummy rows tests the covariate itself;
  the new block, printed under the coefficients, gives one test per
  covariate on (K−1) × (levels−1) degrees of freedom and inherits
  whichever variance estimator the model was fitted with.
  `prepare_covariates()` now records which dummy columns belong to one
  variable, and
  [`analytical_wald_test()`](https://pdvalencia.github.io/mixtureEM/reference/analytical_wald_test.md)
  matches `term_name` against the *variable* rather than by substring.
  The printed block carries the Hauck–Donner caveat: a covariate strong
  enough to nearly separate a class drives its own Wald statistic back
  towards zero.

- **Mixture latent Markov models, and the mover-stayer restriction.**
  [`fit_lta()`](https://pdvalencia.github.io/mixtureEM/reference/fit_lta.md)
  gains `n_classes`, which puts latent classes above the chain: each
  class gets its own initial distribution and its own transition
  matrices while they share one measurement model. `mover_stayer = TRUE`
  restricts the last class to the identity transition matrix — a group
  with zero probability of change (Vermunt, 2004). It answers a question
  a single chain cannot: is this one process, or a mobile group
  alongside a stable one?

  [`transition_matrix()`](https://pdvalencia.github.io/mixtureEM/reference/transition_matrix.md),
  [`status_prevalences()`](https://pdvalencia.github.io/mixtureEM/reference/status_prevalences.md)
  and [`plot()`](https://rdrr.io/r/graphics/plot.default.html) take a
  `class` argument; [`print()`](https://rdrr.io/r/base/print.html) and
  [`summary()`](https://rdrr.io/r/base/summary.html) report the class
  sizes and each class’s chain. Standard errors are **not** available
  for a mixture over chains, and the *unrestricted* mixture is demanding
  of the data — with a single indicator per occasion it is routinely not
  identified, and classes that collapse onto the same chain raise a
  warning.

- **Absolute and local fit diagnostics for categorical models.**
  [`absolute_fit()`](https://pdvalencia.github.io/mixtureEM/reference/absolute_fit.md)
  compares observed response-pattern frequencies with the model-implied
  ones (likelihood-ratio L², Pearson X², Cressie–Read).
  [`bivariate_residuals()`](https://pdvalencia.github.io/mixtureEM/reference/bivariate_residuals.md)
  is the local counterpart: the Pearson chi-square of each *pair* of
  items against the model-implied two-way table, divided by its degrees
  of freedom — the conditional-independence assumption showing its seams
  (Oberski et al., 2013).
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws them as
  a heat table anchored at 1.
  [`classification_table()`](https://pdvalencia.github.io/mixtureEM/reference/classification_table.md)
  yields the classification error that the bias-adjusted 3-step
  estimators exist to undo, and
  [`classification_diagnostics()`](https://pdvalencia.github.io/mixtureEM/reference/classification_diagnostics.md)
  prints it alongside the AvePP matrix. Neither diagnostic enumerates
  the response-pattern table, so both stay usable as indicators are
  added; the tests re-derive each statistic the slow way on small tables
  and require agreement to machine precision.

- [`classification_diagnostics()`](https://pdvalencia.github.io/mixtureEM/reference/classification_diagnostics.md)
  now uses the case weights, which matters for frequency-weighted
  pattern files where one row stands for many cases.

- **[`fit_gmm()`](https://pdvalencia.github.io/mixtureEM/reference/fit_gmm.md)
  fits growth mixture models**: class-specific growth curves in which
  cases also vary *within* their class through random growth factors
  (Muthén & Shedden, 1999). Because the outcome is continuous the random
  effects integrate out in closed form, so everything the package
  already offers applies unchanged:
  [`compare_mixtures()`](https://pdvalencia.github.io/mixtureEM/reference/compare_mixtures.md)
  and
  [`blrt()`](https://pdvalencia.github.io/mixtureEM/reference/blrt.md),
  predictors and distal outcomes through the three-step machinery,
  survey designs, and FIML for unobserved waves. The covariance
  structure is set by `random_effects`, `psi`, `residual`, and
  `residual_equal`, with conventional defaults. Solutions at the edge of
  the parameter space — a residual variance driven to zero, or a
  singular growth-factor covariance — warn rather than being reported as
  ordinary estimates.

- [`fit_gmm()`](https://pdvalencia.github.io/mixtureEM/reference/fit_gmm.md)
  also takes `growth_predictors`: covariates on the growth factors
  themselves, answering who *within* a trajectory class starts higher or
  grows faster — a different question from `predictors`, which asks who
  is in the class at all. Both can be asked in the same model with
  `n_steps = 1`. Growth-factor covariates must be complete: a case
  missing one has no model-implied mean, and
  [`fit_gmm()`](https://pdvalencia.github.io/mixtureEM/reference/fit_gmm.md)
  says so rather than deleting the case.

- **[`fit_lcga()`](https://pdvalencia.github.io/mixtureEM/reference/fit_lcga.md)
  fits latent class growth analysis**: each class follows its own
  fixed-effect polynomial trajectory (Nagin, 2005), with
  `family = "binomial"`, `"gaussian"`, or `"poisson"` setting the link
  scale. `degree` sets the shape of the curve and `time_scores` the
  spacing of the occasions. `flexmix` is a new `Suggests` dependency,
  used in the tests as independent corroboration on simulated designs.

- [`plot()`](https://rdrr.io/r/graphics/plot.default.html) for
  [`fit_lcga()`](https://pdvalencia.github.io/mixtureEM/reference/fit_lcga.md)
  and
  [`fit_gmm()`](https://pdvalencia.github.io/mixtureEM/reference/fit_gmm.md)
  draws the observed data behind the estimated curves — without it there
  is no way to tell a class that describes a real subgroup from one that
  has split a single cloud in half.

- `measurement = "count"` fits a Poisson measurement model: each
  indicator is Poisson-distributed within class with its own rate.
  Counts work everywhere the other indicator types do, and
  [`measurement_summary()`](https://pdvalencia.github.io/mixtureEM/reference/measurement_summary.md)
  reports a rate (lambda) per item and class.

- An emission can ask the EM driver for a higher iteration ceiling and a
  two-stage multi-start search, in which a short first pass ranks the
  restarts and only the survivors are run to convergence. Only the
  growth mixture emission uses them; every other model runs exactly the
  loop it did before.

### Bug fixes

- [`fit_lta()`](https://pdvalencia.github.io/mixtureEM/reference/fit_lta.md)
  reported a transition as sitting on the zero boundary only when it
  fell below `1e-6`; EM stops when the log-likelihood stops moving, not
  when a cell reaches zero, so practically-zero cells went unflagged.
  The threshold is now `1e-4`.

- **EM now converges for count, polytomous and mixed measurement
  models.** The default stopping rule is loose on purpose for binary and
  continuous indicators, whose fits are polished by an L-BFGS
  refinement; every other measurement model was outside that refinement,
  so where EM stopped was the answer — after five to ten iterations,
  short of the optimum by enough to reverse a BIC comparison between
  class solutions. Affected models now run the iterations they always
  needed (and are accordingly slower), including
  [`blrt()`](https://pdvalencia.github.io/mixtureEM/reference/blrt.md)
  replicates on those models.

- [`fit_mixture()`](https://pdvalencia.github.io/mixtureEM/reference/fit_mixture.md)
  warns when EM stops at `max_iter` without converging, instead of
  reporting it only through
  [`print()`](https://rdrr.io/r/base/print.html).

- Cases missing on *every* indicator are deleted before estimation and
  reported by [`print()`](https://rdrr.io/r/base/print.html) and in
  `$missing_data$n_empty_rows`. Previously they were retained, which
  left the parameter estimates correct but inflated *n* in BIC/SABIC,
  assigned every empty case modally to the largest class, and dragged
  relative entropy down. Two of the bundled datasets (`janousch`,
  `yrbs2005`) contain such rows, so their reported sample sizes and
  entropies change accordingly.

### Datasets and vignettes

- The vignette lineup now covers the package’s main analyses end to end,
  using only the exported accessor functions: LCA with covariates and a
  distal outcome (`ventura_leon`), latent profile analysis with missing
  data and multiple groups (`janousch`), class enumeration,
  repeated-measures LCA, latent transition analysis with a mover-stayer
  class, growth mixture modelling, and complex-survey LCA (`yrbs2005`).
- The `nlsy79_drinking` dataset and its vignette were removed; the
  repeated-measures LCA vignette now uses a bundled simulated dataset
  with a documented generating script.

## mixtureEM 0.1.0

- Initial release: latent class and profile analysis via mixture
  modelling, including binary, categorical, and continuous measurement
  models, structural models (predictors, distal outcomes), multi-step
  estimation with BCH/ML bias corrections, complex survey designs,
  missing data handling, and longitudinal models
  ([`fit_lta()`](https://pdvalencia.github.io/mixtureEM/reference/fit_lta.md),
  [`fit_rmlca()`](https://pdvalencia.github.io/mixtureEM/reference/fit_rmlca.md)).
