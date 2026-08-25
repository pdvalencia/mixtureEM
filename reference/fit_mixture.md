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
  measurement,
  predictors = NULL,
  outcome = NULL,
  outcome_covariates = NULL,
  outcome_type = c("auto", "continuous", "categorical"),
  slopes = "pooled",
  group = NULL,
  group_effects = c("both", "measurement", "prevalence", "none"),
  group_invariant_items = NULL,
  group_invariant_params = NULL,
  group_prevalence_equal = NULL,
  variances_equal = NULL,
  n_steps = 1,
  correction = "none",
  assignment = c("proportional", "modal"),
  n_init = 20,
  max_iter = 1000,
  random_state = NULL,
  order_by_size = TRUE,
  weights = NULL,
  weight_type = c("sampling", "frequency"),
  strata = NULL,
  cluster = NULL,
  refine = TRUE,
  bayes_constants = NULL,
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

  What kind of variables the indicators are. Either one type string for
  all of them — `"binary"`, `"categorical"`, `"continuous"`,
  `"gaussian"`, `"count"` — or, for a mixed-type model, a named list
  mapping each type to the columns it governs by name or index, e.g.
  `list(binary = c("q1","q2"), continuous = "score")`.

  **Required; there is no default.** The storage mode of a column does
  not determine its measurement model: a 1-5 column is a legitimate
  `"categorical"`, `"continuous"` or `"count"` indicator, and the class
  solution differs across the three. Inferring the type from the data
  would settle a modelling question by inspecting storage mode, and
  would make a script's meaning depend on the data it happens to be run
  against; a constant default is the same guess with the data-dependence
  removed. Omitting the argument is an error that lists the valid types
  and suggests one based on your columns — as a hint for you to confirm,
  not a choice the package makes.

  **Missing values need no special handling.** Any indicator containing
  `NA` is estimated by full-information maximum likelihood under the
  usual missing-at-random assumption, chosen automatically; cases
  missing on every indicator are dropped, with a warning saying how
  many. There is nothing to switch on. (The `"*_nan"` spellings are
  still accepted, and force the same estimator, but they are a leftover
  and you do not need them.)

  **Binary items need not be coded 0/1.** A two-level factor or
  character, a logical, or any pair of numbers — 1/2 is the commonest —
  is converted for you, with a message saying which level became the 1,
  since that is what the reported probabilities refer to.

  `"categorical"` expects integer codes running 1, 2, 3, ...; add 1 to a
  0-based item. `"count"` fits a Poisson model, one rate per item and
  class, and needs non-negative integers.

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
  `"pooled"` (default; one slope shared across classes),
  `"class_specific"` (a separate slope per class for every covariate),
  or a character vector of covariate names (or a one-sided formula
  naming them, e.g. `~ loc1 + loc2`) giving a separate slope per class
  to just those covariates while the rest stay pooled – letting the
  class moderate some covariates while adjusting for others. The last
  form is continuous-outcome only.

- group:

  Optional observed grouping variable for a multiple-group model
  (Collins & Lanza, 2010, sec. 5.7-5.12), e.g. grade or gender. Unlike
  `predictors`, which only ever shifts class membership, a group can
  also be allowed to shift the item-response probabilities themselves;
  see `group_effects`.

  **A covariate whose effect on class membership differs by group** –
  moderation, not just adjustment – does not need `group` at all: build
  the interaction into `predictors` directly, e.g.
  `predictors = model.matrix(~ grade * factor(year))[, -1]` handed to
  `fit_mixture()` in place of a plain covariate matrix. Leave `group`
  unset for this (or use `group_effects = "measurement"` if the item
  parameters should also be free by group): a `"prevalence"` or `"both"`
  group effect appends the group's own design to `predictors`
  internally, so the group columns would then appear twice. A covariate
  matrix built this way also carries none of the column-to-variable
  bookkeeping a data-frame `predictors` does, so functions like
  [`analytical_wald_test`](https://pdvalencia.github.io/mixtureEM/reference/analytical_wald_test.md)
  need the individual interaction column names rather than a single
  variable name.

- group_effects:

  What you are allowing the groups to differ in. Pick by the question
  you are asking:

  `"prevalence"`

  :   Do the groups differ in *how many* people fall in each class? The
      classes themselves mean the same thing in every group; only their
      sizes move. This is the model most applied analyses want.

  `"both"` (default)

  :   Does each group need its *own* class solution? Both the class
      sizes and what the classes look like are free to differ — the
      configural model.

  `"none"`

  :   Ignore the grouping variable when estimating.

  **The usual workflow** is to fit `"prevalence"` and `"both"` and
  compare them with
  [`lr_test()`](https://pdvalencia.github.io/mixtureEM/reference/lr_test.md).
  That is the test of measurement invariance (Collins & Lanza, 2010,
  sec. 5.8): it asks whether letting the classes mean different things
  in different groups buys a significantly better fit. If it does not —
  the common outcome, and the one you want — keep `"prevalence"` and
  report the group differences in class sizes, which are then comparable
  across groups because the classes are. Comparing `"none"` against
  `"prevalence"` tests whether the sizes differ at all (sec. 5.11). Pass
  `n_steps = 1` for both fits.

  There is a fourth setting, `"measurement"`, which frees the item
  parameters while holding the class sizes pooled. It answers an unusual
  question and is rarely what is wanted; `"both"` is the configural
  model the invariance literature actually compares against.

- group_invariant_items:

  Item indices or names held equal across groups even when
  `group_effects` frees the measurement model (Collins & Lanza's
  partial-invariance models, sec. 5.9). `NULL` (the default) leaves
  every item free, i.e. a fully configural model.

- group_invariant_params:

  For continuous indicators, which *kind* of parameter is held equal
  across groups when `group_effects` frees the measurement model:
  `"means"`, `"covariances"`, or both. `NULL` (the default) frees
  everything. Where `group_invariant_items` shares whole items, this
  shares one parameter matrix for every item — which is the model the
  latent-profile invariance literature actually fits: Olivera-Aguilar
  and Rikoon's (2018) "unconstrained" model frees the class means across
  groups while holding the indicator variances invariant, i.e.
  `group_invariant_params = "covariances"`, and it is that model, not
  the fully heterogeneous one, that their invariance test compares
  against. The two are alternatives, not combinable. Categorical
  indicators have a single kind of parameter, so for them this
  constraint and `group_invariant_items` coincide and only the latter is
  offered.

- group_prevalence_equal:

  Class indices (or `TRUE` for all classes) whose prevalence is held to
  one shared value across every group, while the remaining classes'
  prevalences stay free within each group (Collins & Lanza's restricted
  multiple-group models, sec. 5.11-5.12). Requires `group_effects` to be
  `"prevalence"` or `"both"`. `NULL` (the default) leaves every class's
  prevalence free by group, fit through the ordinary class-membership
  regression on `group` — the two are numerically equivalent when
  unconstrained, but only the regression route has `predict_class`'s
  usual Wald-testable coefficients. This is why: an ordinary regression
  coefficient of zero on the group dummies pins a class to the
  *reference group*'s prevalence, not to one shared across every group,
  so there is no coefficient value that expresses this constraint.
  Standard errors for the frozen prevalences are not produced in this
  first pass; compare nested models with
  [`lr_test()`](https://pdvalencia.github.io/mixtureEM/reference/lr_test.md)
  instead. An out-of-range class index is an error, never silently
  ignored.

- variances_equal:

  Logical, for continuous indicators only: hold each item's variance
  equal across the classes, so the classes differ in location only. This
  is the homoscedastic latent profile model and the default
  parameterisation of several commercial programs, and it combines with
  `group_invariant_params` to give a variance that is free across groups
  but shared by the classes within each. Passed through to the
  measurement model, so it is also available on an ordinary single-group
  fit.

  **This is the default when `measurement` is continuous**, and
  `variances_equal = FALSE` recovers the class-varying parameterisation
  the package used previously. The reason is that the unrestricted
  normal-mixture likelihood is *unbounded*: send a class mean to any
  single data point and that class's variance to zero and the likelihood
  diverges, so no maximum likelihood estimate exists and what the EM
  reports is a local optimum whose properties are not guaranteed.
  Holding the variances equal across classes bounds the likelihood, and
  the constrained estimator is consistent (Day, 1969; Hathaway, 1985).
  Freeing them also invites classes that describe non-normality in a
  single population rather than distinct subgroups (Bauer & Curran,
  2003).

  The restriction is substantive and testable, and you are expected to
  fit both and compare rather than accept either blindly. It is not
  offered as the *safe* choice but as the well-posed one: when the
  homoscedastic model is wrong, it fails visibly — a genuinely
  heteroscedastic class splits into two, and the comparison against the
  free model says so. When the free model is wrong it fails silently, as
  a boundary solution that gets written up as a finding.

- n_steps:

  Estimation strategy: 1 (simultaneous), 2, or 3 (recommended when a
  structural model is present). Defaults to 3 when `predictors` or
  `outcome` is supplied and left unset, otherwise 1.

- correction:

  Bias correction for 3-step estimation: `"none"`, `"ML"`, or `"BCH"`.
  When left unset for a 3-step structural model, a recommended default
  is chosen (ML for predictors and categorical outcomes, BCH for
  continuous outcomes).

- assignment:

  How step 1's posteriors are turned into the assigned-class variable
  whose classification error the correction inverts: `"proportional"`
  (default; each case carries its posterior probability in every class)
  or `"modal"` (each case is assigned to its most likely class
  outright). See
  [`add_covariates`](https://pdvalencia.github.io/mixtureEM/reference/add_covariates.md).

- n_init, max_iter, random_state, order_by_size, refine:

  Estimation controls: number of random starts (default 20), maximum EM
  iterations, RNG seed, whether to order classes by size, and whether to
  run L-BFGS refinement.

  A mixture likelihood usually has several local maxima, so a single
  start is a coin toss rather than an estimate; the fit reports how many
  of the starts reached the solution it kept. If that count is 1, refit
  with `n_init = 100` — a maximum found once may simply be the best of a
  small sample of the surface, and if it does not replicate at 100
  starts the problem is more likely the specification than the search
  (Hipp & Bauer, 2006). Set `random_state` to make the search
  reproducible.

  The default of 20 is a floor for small, well-separated models rather
  than a guarantee: under a correct specification the maximum is
  typically found by 39\\ there is one, but local optima multiply with
  more classes and poorer separation. `max_iter` defaults to 1000,
  following the iteration budget in Biernacki et al. (2003); when a fit
  stops at the cap, double it.
  [`vignette("estimation")`](https://pdvalencia.github.io/mixtureEM/articles/estimation.md)
  gives the figures behind both numbers.

  The cost is roughly linear in `n_init`, so it is easy to budget for a
  larger search: `n_init = 20` took about 16 seconds on a 4-class,
  7-item binary model with 2,587 cases, so `n_init = 200` on the same
  model is a couple of minutes and `n_init = 1000` is closer to ten.

  More starts are not *monotonically* better on a large model. Where
  each EM iteration is expensive the search runs in two stages: a short
  pass ranks the starts and only the best fraction is run to
  convergence. That ranking is taken before the starts have converged,
  so it can discard a slow-climbing basin that would have won, and
  raising `n_init` changes which starts survive rather than simply
  adding to them. A larger `n_init` can therefore land on a slightly
  worse solution than a smaller one. On a big model, compare two or
  three values rather than assuming the largest is best.

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

- bayes_constants:

  Optional list adjusting the strength of the weak priors the estimator
  places on each block of parameters. Named `latent` (class weights and,
  in the transition models, the initial and transition probabilities),
  `categorical` (item-response probabilities), `poisson` (count rates),
  and `variances` (class-specific variances of continuous indicators);
  all default to `1`. Each is a number of pseudo-observations spread
  over the classes, so its influence shrinks as the sample grows.

  The default of one observation, and the fact that it is spread in
  agreement with each item's own observed marginal rather than uniformly
  over the cells, are both taken from Galindo Garre and Vermunt (2006).
  Their simulation compares Jeffreys, normal and two Dirichlet priors
  against maximum likelihood and the parametric bootstrap, at a prior
  strength of a single added case, and the marginal-preserving Dirichlet
  has the lowest root median squared error and the best coverage in
  every condition they study — most clearly where a true parameter is
  near the boundary and small samples send the maximum likelihood
  estimate to infinity. The constant Dirichlet, which spreads the same
  one observation uniformly, degrades as the table grows: with nine
  binary items that observation is spread over 512 cells rather than 32,
  and its estimates shrink too far toward zero. This is why
  `categorical` is not simply an add-one adjustment.

  The `variances` prior is centred on each item's own observed marginal
  variance rather than on a fixed number, which is what makes it mean
  the same thing on a five-point scale and on an income variable:
  rescaling an indicator rescales the prior with it. That invariance is
  the reason the penalised-likelihood literature recommends a
  data-scaled penalty over a constant one (Chen, Tan, & Zhang, 2008,
  sec. 4).

  Two uses. **Reproducing an unregularized fit:** setting a constant to
  `0` removes that prior and gives plain maximum likelihood for that
  block. This is an escape hatch for matching a reference analysis, not
  a recommended setting — the unpenalised mixture likelihood for a
  mixture of normals is unbounded, so a *global* maximum likelihood
  estimate does not exist (Day, 1969; Kiefer & Wolfowitz, 1956) and what
  an unregularized program reports is a local maximum.

  **Rescuing a collapsed fit:** this situation is now rare, because a
  continuous measurement model holds the variances equal across classes
  by default and that bounds the likelihood; it arises when you have
  explicitly passed `variances_equal = FALSE`. If a fit warns that a
  class variance has collapsed, the prior is one of three remedies, and
  not the first to reach for. Constraining the variances to be equal
  across classes (`variances_equal = TRUE`) bounds the likelihood so the
  problem cannot arise at all; fitting fewer classes often removes the
  class that was describing a spike. Where neither is acceptable
  substantively, raise `variances`. A useful starting point is **one
  artificial observation per class**, i.e. `variances = n_classes`,
  increasing it a little at a time if the warning persists. State it
  that way rather than as a bare number: the constant is divided among
  the classes, so `variances = n_classes` is what holds the prior at one
  pseudo-observation per class at *every* number of classes — the way
  the penalised-mixture literature applies it, the same amount to every
  component — while a fixed value drifts as classes are added. Then
  check the result rather than just that the warning stopped — the
  flagged variance should no longer be far below the others in the
  model, and its class mean should have come off the floor or ceiling of
  the response scale. It is worth looking at the distribution of the
  named item as well, for the floor, ceiling or spike the class latched
  onto.

  Raising `n_init` is *not* a remedy here. A collapsed variance is not a
  convergence failure but a property of the likelihood, which really is
  unbounded in that direction, so a longer search can find a taller
  spike and report a better log-likelihood for a worse solution. For the
  same reason a flagged fit's BIC is not comparable with a clean fit's.

  The theory behind the prior settles its form and not its size: within
  the conditions Chen, Tan, & Zhang (2008) impose — which any fixed
  positive constant meets — a penalty that diverges as a variance
  approaches zero and grows more slowly than the sample yields a
  consistent estimator whatever its constant. That result is proved for
  *univariate* normal mixtures; the multivariate case, which is what
  this package fits, they leave open (sec. 5). The choice of constant is
  therefore all the more a finite-sample judgement, which is why the
  default is deliberately weak and the guidance above is a rule with a
  check rather than a magic value.

  This is not a tuning menu. The defaults are the intended settings.

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

## Details

**The EM convergence rule is fixed and not user-adjustable.** Each
random start stops once the log-likelihood improves by less than `1e-4`
in absolute terms or a relative `1e-8`, whichever is looser, and there
is no argument that changes this. A looser rule was tried and measured
to cost real log-likelihood — on a validated 4-class binary example,
stopping early and then polishing with a numerical optimizer recovered
only a ninth of what a tighter EM rule found outright — and, more
importantly, it degrades the multi-start search itself: a loose rule
ranks candidate starts on log-likelihood values too coarse to tell a
good basin from a mediocre one, so `n_init` stops doing its job. If a
fit seems slow, raise `n_init` or trim `max_iter` rather than looking
for a tolerance argument; there is not one to find.

## Multiple-group models

Under `group_effects = "both"` or `"measurement"` each group's item
parameters are estimated from that group's own cases, so the *class
labels* need not line up across groups: "Class 1" in one group's profile
is not guaranteed to be the same kind of class as "Class 1" in
another's. Read the per-group profiles by their item patterns rather
than by position. The likelihood-ratio comparison is unaffected, since
it only differences total log-likelihoods.

The log-likelihood is that of the indicators *given* the group; the
grouping variable's own distribution is not modelled and its proportions
are not counted as parameters. Software that treats a known grouping
variable as a latent class variable observed without error adds both.
For comparison with such output, a `group` fit also carries
`metrics$ll_knownclass` and `metrics$n_params_knownclass`; the
difference is a fixed constant and cancels in
[`lr_test()`](https://pdvalencia.github.io/mixtureEM/reference/lr_test.md).

## References

Biernacki, C., Celeux, G., & Govaert, G. (2003). Choosing starting
values for the EM algorithm for getting the highest likelihood in
multivariate Gaussian mixture models. *Computational Statistics & Data
Analysis*, *41*(3-4), 561-575.
[doi:10.1016/S0167-9473(02)00163-9](https://doi.org/10.1016/S0167-9473%2802%2900163-9)
(the iteration budget behind `max_iter`).

Hipp, J. R., & Bauer, D. J. (2006). Local solutions in the estimation of
growth mixture models. *Psychological Methods*, *11*(1), 36-53.
[doi:10.1037/1082-989X.11.1.36](https://doi.org/10.1037/1082-989X.11.1.36)
(the replication rates behind `n_init`).

Shireman, E., Steinley, D., & Brusco, M. J. (2017). Examining the effect
of initialization strategies on the performance of Gaussian mixture
modeling. *Behavior Research Methods*, *49*(1), 282-293.
[doi:10.3758/s13428-015-0697-6](https://doi.org/10.3758/s13428-015-0697-6)

Collins, L. M., & Lanza, S. T. (2010). *Latent Class and Latent
Transition Analysis: With Applications in the Social, Behavioral, and
Health Sciences*. Wiley.

Masyn, K. E. (2013). Latent class analysis and finite mixture modeling.
In T. D. Little (Ed.), *The Oxford Handbook of Quantitative Methods*
(Vol. 2, pp. 551-611). Oxford University Press.

Vermunt, J. K., & Magidson, J. (2002). Latent class cluster analysis. In
J. A. Hagenaars & A. L. McCutcheon (Eds.), *Applied Latent Class
Analysis* (pp. 89-106). Cambridge University Press.

Galindo Garre, F., & Vermunt, J. K. (2006). Avoiding boundary estimates
in latent class analysis by Bayesian posterior mode estimation.
*Behaviormetrika*, *33*(1), 43-59.
[doi:10.2333/bhmk.33.43](https://doi.org/10.2333/bhmk.33.43) (the priors
behind `bayes_constants`).

Chen, J., Tan, X., & Zhang, R. (2008). Inference for normal mixtures in
mean and variance. *Statistica Sinica*, *18*(2), 443-465 (the penalty
behind `bayes_constants$variances`).

Day, N. E. (1969). Estimating the components of a mixture of normal
distributions. *Biometrika*, *56*(3), 463-474 (the unbounded likelihood
behind the `variances_equal` default).

Hathaway, R. J. (1985). A constrained formulation of maximum-likelihood
estimation for normal mixture distributions. *The Annals of Statistics*,
*13*(2), 795-800 (consistency of the constrained estimator).

Bauer, D. J., & Curran, P. J. (2003). Distributional assumptions of
growth mixture models: Implications for overextraction of latent
trajectory classes. *Psychological Methods*, *8*(3), 338-363.
[doi:10.1037/1082-989X.8.3.338](https://doi.org/10.1037/1082-989X.8.3.338)

Kiefer, J., & Wolfowitz, J. (1956). Consistency of the maximum
likelihood estimator in the presence of infinitely many incidental
parameters. *The Annals of Mathematical Statistics*, *27*(4), 887-906.

McLachlan, G. J., & Peel, D. (2000). *Finite Mixture Models*. Wiley (on
the unbounded likelihood and spurious solutions).

Olivera-Aguilar, M., & Rikoon, S. H. (2018). Assessing measurement
invariance in multiple-group latent profile analysis. *Structural
Equation Modeling*, *25*(3), 439-452.
[doi:10.1080/10705511.2017.1408015](https://doi.org/10.1080/10705511.2017.1408015)

## Examples

``` r
set.seed(1)
X <- matrix(rbinom(500, 1, 0.5), nrow = 100)
fit <- fit_mixture(X, n_classes = 2, measurement = "binary")

if (FALSE) { # \dontrun{
# Class membership predicted by a covariate (3-step, ML by default)
fit_mixture(X, n_classes = 2, measurement = "binary", predictors = age)

# Distal outcome with a class-specific covariate slope
fit_mixture(X, n_classes = 2, measurement = "binary", outcome = bmi,
            outcome_covariates = age, slopes = "class_specific")

# Class moderates level-of-care while age and gender are only adjusted
# for: a mix of "class_specific" and "pooled" in one model.
fit_mixture(X, n_classes = 4, measurement = "binary", outcome = cannabis_days,
            outcome_covariates = data.frame(loc1, loc2, loc3, age, gender),
            slopes = c("loc1", "loc2", "loc3"))

# Mixed-type indicators
fit_mixture(items, n_classes = 3,
            measurement = list(binary = 1:5, continuous = 6:8))
} # }
```
