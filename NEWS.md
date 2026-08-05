# mixtureEM 0.2.0

## New: stepwise analyses on a fitted model

* **`add_covariates()` and `add_outcome()`.** Once you have chosen an
  unconditional model, relating its classes to external variables no longer
  requires re-specifying (and re-estimating) the whole measurement model:

  ```r
  fit <- fit_mixture(items, n_classes = 3)
  summary(add_covariates(fit, data[, c("age", "sex")]))
  summary(add_outcome(fit, data$bmi))
  ```

  Both verbs run only steps 2–3 of the bias-adjusted three-step approach
  (Bakk et al., 2014; Vermunt, 2010) on the stored step-1
  solution, so the classes you inspected are exactly the classes the
  structural model describes — same solution, same class order, and no risk
  of a re-run landing in a different optimum. `fit_mixture(predictors = )`
  and `fit_mixture(outcome = )` remain available for one-call workflows and
  one-step (simultaneous) estimation.

## New: summary functions return their numbers

* `measurement_summary()` now returns (invisibly) a long-format data frame of
  the item parameters — `item`, `category`, `class`, `estimate` — alongside
  the printed tables, so results can be reused programmatically without
  touching the fitted object's internals.
* `summary()` on a model with a structural part returns (invisibly) a list of
  tidy tables: `$coefficients` (class contrasts with SEs, p-values, odds
  ratios and confidence limits), `$omnibus` (per-covariate Wald tests), and
  `$outcome` (predicted probabilities or class means with their tests).
* New `class_sizes()`: one row per class with the model-estimated proportion,
  the expected count, and the modal-assignment count — the three numbers
  applied papers report.

## Output and plotting improvements

* Long variable and category names no longer break table alignment. Label
  columns size themselves to their content, and names past ~28 characters are
  abbreviated in the display with a key printed under the table (the returned
  data frames always carry the full names).
* Profile and trajectory plots no longer clip long class labels on the right:
  margins are computed from the actual label widths, the legend position
  scales with the x-range, and overly long labels are shortened. The bottom
  margin now adapts to long rotated indicator names as well.
* Factor covariates: levels with no observations (e.g. after subsetting a
  pooled dataset) are dropped automatically instead of producing a degenerate
  `OR = 1.00 [1.00, 1.00]` row, and levels observed fewer than 5 times
  trigger a warning that their coefficients will be unstable.

## Breaking changes

* **Standard errors for covariate effects in 2- and 3-step models have
  changed, and every p-value and confidence interval for a class predictor
  moves with them.** They were computed from the Hessian the M-step stores,
  which measures the curvature of the *Q function* rather than of the step-3
  log-likelihood the estimates actually maximise, and which takes no account
  of the fact that proportional assignment gives each case K weighted
  records. Both errors run the same way, so the reported standard errors were
  too small and the intervals too narrow.

  The new default, `se = "corrected"` on `fit_mixture()`, is the first-order
  corrected estimator of Bakk, Oberski and Vermunt (2014): the step-3
  sandwich plus the variance carried over from step 1, which is estimated and
  not known. `se = "robust"` gives the sandwich alone; `se = "hessian"`
  inverts the step-3 observed information. See `?covariate_se`.

  A 250-replication coverage study on the design of Bakk et al.
  (`data-raw/covariate_se_simulation.R`) says how much this mattered. With
  n = 500 and entropy R² = 0.64, the old nominal 95% intervals covered
  between **54% and 88%** of the time, worst for the strongest coefficient;
  the new default covers 92–98%. With entropy R² near 0.90 the old intervals
  covered 87–95% and the new ones 94–98%, and there the correction term is
  negligible — the sandwich alone does the work. **If you have published
  covariate p-values from this package with small or poorly separated
  classes, they were anti-conservative and are worth re-running.**

  `summary()`, `confint()` and `analytical_wald_test()` now print the name of
  the estimator behind the numbers they show. `correction = "BCH"`, one-step
  fits, and a covariate combined with a distal outcome in one structural
  model keep the uncorrected Hessian, and say so.

## New features

* **`summary()` prints an omnibus Wald test for each class predictor.** For a
  covariate with more than two categories, none of the per-dummy rows tests
  the covariate itself; the new block, printed under the coefficients, gives
  one test per covariate on (K−1) × (levels−1) degrees of freedom and
  inherits whichever variance estimator the model was fitted with.
  `prepare_covariates()` now records which dummy columns belong to one
  variable, and `analytical_wald_test()` matches `term_name` against the
  *variable* rather than by substring. The printed block carries the
  Hauck–Donner caveat: a covariate strong enough to nearly separate a class
  drives its own Wald statistic back towards zero.

* **Mixture latent Markov models, and the mover-stayer restriction.**
  `fit_lta()` gains `n_classes`, which puts latent classes above the chain:
  each class gets its own initial distribution and its own transition
  matrices while they share one measurement model. `mover_stayer = TRUE`
  restricts the last class to the identity transition matrix — a group with
  zero probability of change (Vermunt, 2004). It answers a
  question a single chain cannot: is this one process, or a mobile group
  alongside a stable one?

  `transition_matrix()`, `status_prevalences()` and `plot()` take a `class`
  argument; `print()` and `summary()` report the class sizes and each class's
  chain. Standard errors are **not** available for a mixture over chains, and
  the *unrestricted* mixture is demanding of the data — with a single
  indicator per occasion it is routinely not identified, and classes that
  collapse onto the same chain raise a warning.

* **Absolute and local fit diagnostics for categorical models.**
  `absolute_fit()` compares observed response-pattern frequencies with the
  model-implied ones (likelihood-ratio L², Pearson X², Cressie–Read).
  `bivariate_residuals()` is the local counterpart: the Pearson chi-square of
  each *pair* of items against the model-implied two-way table, divided by
  its degrees of freedom — the conditional-independence assumption showing
  its seams (Oberski et al., 2013). `plot()` draws them
  as a heat table anchored at 1. `classification_table()` yields the
  classification error that the bias-adjusted 3-step estimators exist to
  undo, and `classification_diagnostics()` prints it alongside the AvePP
  matrix. Neither diagnostic enumerates the response-pattern table, so both
  stay usable as indicators are added; the tests re-derive each statistic the
  slow way on small tables and require agreement to machine precision.

* `classification_diagnostics()` now uses the case weights, which matters for
  frequency-weighted pattern files where one row stands for many cases.

* **`fit_gmm()` fits growth mixture models**: class-specific growth curves in
  which cases also vary *within* their class through random growth factors
  (Muthén & Shedden, 1999). Because the outcome is continuous the random
  effects integrate out in closed form, so everything the package already
  offers applies unchanged: `compare_mixtures()` and `blrt()`, predictors and
  distal outcomes through the three-step machinery, survey designs, and FIML
  for unobserved waves. The covariance structure is set by `random_effects`,
  `psi`, `residual`, and `residual_equal`, with conventional defaults.
  Solutions at the edge of the parameter space — a residual variance driven
  to zero, or a singular growth-factor covariance — warn rather than being
  reported as ordinary estimates.

* `fit_gmm()` also takes `growth_predictors`: covariates on the growth
  factors themselves, answering who *within* a trajectory class starts higher
  or grows faster — a different question from `predictors`, which asks who is
  in the class at all. Both can be asked in the same model with
  `n_steps = 1`. Growth-factor covariates must be complete: a case missing
  one has no model-implied mean, and `fit_gmm()` says so rather than deleting
  the case.

* **`fit_lcga()` fits latent class growth analysis**: each class follows its
  own fixed-effect polynomial trajectory (Nagin, 2005), with
  `family = "binomial"`, `"gaussian"`, or `"poisson"` setting the link scale.
  `degree` sets the shape of the curve and `time_scores` the spacing of the
  occasions. `flexmix` is a new `Suggests` dependency, used in the tests as
  independent corroboration on simulated designs.

* `plot()` for `fit_lcga()` and `fit_gmm()` draws the observed data behind
  the estimated curves — without it there is no way to tell a class that
  describes a real subgroup from one that has split a single cloud in half.

* `measurement = "count"` fits a Poisson measurement model: each indicator is
  Poisson-distributed within class with its own rate. Counts work everywhere
  the other indicator types do, and `measurement_summary()` reports a rate
  (lambda) per item and class.

* An emission can ask the EM driver for a higher iteration ceiling and a
  two-stage multi-start search, in which a short first pass ranks the
  restarts and only the survivors are run to convergence. Only the growth
  mixture emission uses them; every other model runs exactly the loop it did
  before.

## Bug fixes

* `fit_lta()` reported a transition as sitting on the zero boundary only when
  it fell below `1e-6`; EM stops when the log-likelihood stops moving, not
  when a cell reaches zero, so practically-zero cells went unflagged. The
  threshold is now `1e-4`.

* **EM now converges for count, polytomous and mixed measurement models.**
  The default stopping rule is loose on purpose for binary and continuous
  indicators, whose fits are polished by an L-BFGS refinement; every other
  measurement model was outside that refinement, so where EM stopped was the
  answer — after five to ten iterations, short of the optimum by enough to
  reverse a BIC comparison between class solutions. Affected models now run
  the iterations they always needed (and are accordingly slower), including
  `blrt()` replicates on those models.

* `fit_mixture()` warns when EM stops at `max_iter` without converging,
  instead of reporting it only through `print()`.

* Cases missing on *every* indicator are deleted before estimation and
  reported by `print()` and in `$missing_data$n_empty_rows`. Previously they
  were retained, which left the parameter estimates correct but inflated *n*
  in BIC/SABIC, assigned every empty case modally to the largest class, and
  dragged relative entropy down. Two of the bundled datasets (`janousch`,
  `yrbs2005`) contain such rows, so their reported sample sizes and entropies
  change accordingly.

## Datasets and vignettes

* The vignette lineup now covers the package's main analyses end to end,
  using only the exported accessor functions: LCA with covariates and a
  distal outcome (`ventura_leon`), latent profile analysis with missing data
  and multiple groups (`janousch`), class enumeration, repeated-measures LCA,
  latent transition analysis with a mover-stayer class, growth mixture
  modelling, and complex-survey LCA (`yrbs2005`).
* The `nlsy79_drinking` dataset and its vignette were removed; the
  repeated-measures LCA vignette now uses a bundled simulated dataset with a
  documented generating script.

# mixtureEM 0.1.0

* Initial release: latent class and profile analysis via mixture modelling,
  including binary, categorical, and continuous measurement models,
  structural models (predictors, distal outcomes), multi-step estimation
  with BCH/ML bias corrections, complex survey designs, missing data
  handling, and longitudinal models (`fit_lta()`, `fit_rmlca()`).
