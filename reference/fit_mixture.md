# Fit a Latent Class or Latent Profile Mixture Model

Fits a finite mixture (latent class / latent profile) model. The latent
classes are defined by a set of measurement `indicators`; optionally, a
structural model relates the classes to external variables — either as
`predictors` of class membership, or as a distal `outcome` caused by the
classes.

## Usage

``` r
fit_mixture(
  indicators = NULL,
  n_classes = 2,
  measurement = "binary",
  predictors = NULL,
  outcome = NULL,
  outcome_covariates = NULL,
  outcome_type = c("auto", "continuous", "categorical"),
  slopes = c("pooled", "class_specific"),
  group = NULL,
  group_effects = c("both", "measurement", "prevalence", "none"),
  group_invariant_items = NULL,
  n_steps = 1,
  correction = "none",
  n_init = 1,
  max_iter = 1000,
  random_state = NULL,
  order_by_size = TRUE,
  weights = NULL,
  weight_type = c("sampling", "frequency"),
  strata = NULL,
  cluster = NULL,
  refine = TRUE,
  se = c("corrected", "robust", "hessian"),
  X = NULL,
  Y = NULL,
  n_components = NULL,
  structural = NULL,
  ...
)
```

## Arguments

- indicators:

  Matrix or data frame of measurement items that define the latent
  classes (rows are observations, columns are items).

- n_classes:

  Number of latent classes/profiles to estimate.

- measurement:

  Either a single type string applied to every indicator (`"binary"`,
  `"categorical"`, `"continuous"`, `"gaussian"`, `"count"`, and `"_nan"`
  missing-data variants), or, for a mixed-type model, a named list
  mapping each type to the indicator columns it governs by name or
  index, e.g. `list(binary = c("q1","q2"), continuous = "score")`.

  `"count"` fits a Poisson measurement model: each item is
  Poisson-distributed within class with its own rate, reported by
  [`measurement_summary`](https://pdvalencia.github.io/mixtureEM/reference/measurement_summary.md)
  as a rate (lambda) per item and class. It requires non-negative
  integer indicators.

- predictors:

  Optional covariates that predict latent class membership. Supplying
  this fits a class-membership regression (the "predict class"
  structural model). Mutually exclusive with `outcome`.

- outcome:

  Optional distal outcome caused by the latent classes. Mutually
  exclusive with `predictors`.

- outcome_covariates:

  Optional covariates that adjust the distal `outcome`.

- outcome_type:

  One of `"auto"`, `"continuous"`, or `"categorical"`. With `"auto"`
  (default) the type is inferred from `outcome`.

- slopes:

  When `outcome_covariates` are supplied, whether their effect is
  `"pooled"` (one slope shared across classes) or `"class_specific"` (a
  separate slope per class).

- group:

  Optional observed grouping variable for a multiple-group model
  (Collins & Lanza, 2010, sec. 5.7-5.12), e.g. grade or gender. Unlike
  `predictors`, which only ever shifts class membership, a group can
  also be allowed to shift the item-response probabilities themselves;
  see `group_effects`.

- group_effects:

  Which parameters `group` is allowed to shift. `"both"` (default) frees
  both the item-response probabilities and the class prevalences across
  groups, i.e. fits each group's own model. `"measurement"` frees only
  the item-response probabilities (prevalences stay pooled).
  `"prevalence"` frees only the class prevalences, by entering `group`
  as a class-membership covariate exactly like `predictors`
  (item-response probabilities stay invariant across groups) — this is
  the same equivalence
  [`fit_rmlca()`](https://pdvalencia.github.io/mixtureEM/reference/fit_rmlca.md)
  documents for its own `predictors`, sec. 6.10.2. `"none"` ignores
  `group` for estimation. Fit the same data under two of these and
  compare them with
  [`longitudinal_lrt()`](https://pdvalencia.github.io/mixtureEM/reference/longitudinal_lrt.md)
  to get Collins & Lanza's measurement-invariance test (sec. 5.8,
  comparing `"both"` against `"prevalence"`) and prevalence-equivalence
  test (sec. 5.11, comparing `"prevalence"` against `"none"`). Because a
  multiple-group model like Collins & Lanza's is fit as one simultaneous
  model rather than through the auxiliary-variable 3-step approximation,
  pass `n_steps = 1` explicitly to reproduce their numbers exactly.

  With `"both"` or `"measurement"`, each group's item-response
  probabilities are estimated from that group's cases alone, tied to the
  other groups only by shared initialization — unlike occasions in
  [`fit_rmlca()`](https://pdvalencia.github.io/mixtureEM/reference/fit_rmlca.md),
  which every case informs simultaneously. A likelihood- ratio
  comparison via
  [`longitudinal_lrt()`](https://pdvalencia.github.io/mixtureEM/reference/longitudinal_lrt.md)
  is unaffected (it only compares total log-likelihoods), but the *class
  labels* a configural fit assigns are not guaranteed to line up across
  groups: "Class 1" in one group's profile need not be the same kind of
  class as "Class 1" in another's. Use a larger `n_init` and a fixed
  `random_state` for configural fits, and when reading per-group
  profiles, match classes by their item-response pattern rather than by
  position.

- group_invariant_items:

  Item indices or names held equal across groups even when
  `group_effects` frees the measurement model (Collins & Lanza's
  partial-invariance models, sec. 5.9). `NULL` (the default) leaves
  every item free, i.e. a fully configural model.

- n_steps:

  Estimation strategy: 1 (simultaneous), 2, or 3 (recommended when a
  structural model is present). Defaults to 3 when `predictors` or
  `outcome` is supplied and left unset, otherwise 1.

- correction:

  Bias correction for 3-step estimation: `"none"`, `"ML"`, or `"BCH"`.
  When left unset for a 3-step structural model, a recommended default
  is chosen (ML for predictors and categorical outcomes, BCH for
  continuous outcomes).

- n_init, max_iter, random_state, order_by_size, refine:

  Estimation controls: number of random starts, maximum EM iterations,
  RNG seed, whether to order classes by size, and whether to run L-BFGS
  refinement.

- weights, strata, cluster:

  Optional survey design: sampling `weights`, and `strata`/`cluster`
  identifiers enabling design-based (linearization) standard errors.

- weight_type:

  What the numbers in `weights` mean. `"sampling"` (the default) treats
  them as survey or probability weights, saying how much of the
  population each case represents; only their relative sizes matter, so
  they are rescaled to sum to the number of cases. `"frequency"` treats
  them as counts of identical cases, as in a response-pattern table
  where one row stands for many respondents; the sample size behind AIC
  and BIC is then the sum of the counts. Getting this wrong changes BIC
  but not the parameter estimates.

- se:

  How standard errors for `predictors` are computed in a 2- or 3-step
  model. `"corrected"` (the default) adds the variance carried over from
  step 1 to the step-3 sandwich, following Bakk et al. (2014);
  `"robust"` reports the sandwich alone; `"hessian"` inverts the step-3
  observed information only. See
  [`covariate_se`](https://pdvalencia.github.io/mixtureEM/reference/covariate_se.md).

- X, Y, n_components, structural:

  Deprecated legacy arguments retained for backward compatibility;
  prefer `indicators`, `n_classes`, `predictors`, and `outcome`.

- ...:

  Passed through to the measurement-model constructors.

## Value

A fitted `mixture_model` object.

## Examples

``` r
set.seed(1)
X <- matrix(rbinom(500, 1, 0.5), nrow = 100)
fit <- fit_mixture(X, n_classes = 2, measurement = "binary")

if (FALSE) { # \dontrun{
# Class membership predicted by a covariate (3-step, ML by default)
fit_mixture(X, n_classes = 2, predictors = age)

# Distal outcome with a class-specific covariate slope
fit_mixture(X, n_classes = 2, outcome = bmi,
            outcome_covariates = age, slopes = "class_specific")

# Mixed-type indicators
fit_mixture(items, n_classes = 3,
            measurement = list(binary = 1:5, continuous = 6:8))
} # }
```
