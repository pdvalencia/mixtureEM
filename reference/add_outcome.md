# Examine a Distal Outcome on a Fitted Model

Takes the latent class model you have already chosen and relates the
classes to a distal outcome with the bias-adjusted three-step approach.
The measurement model is reused exactly as fitted — no re-estimation,
and no risk of landing on a different solution.

## Usage

``` r
add_outcome(
  fit,
  outcome,
  covariates = NULL,
  outcome_type = c("auto", "continuous", "categorical"),
  slopes = "pooled",
  predictors = NULL,
  variances = c("equal", "class_specific"),
  correction = c("auto", "BCH", "ML", "none"),
  steps = c(3, 2),
  se = c("corrected", "robust", "hessian"),
  assignment = c("proportional", "modal"),
  max_iter = 1000,
  data = NULL,
  missing_predictors = c("listwise", "mean"),
  ...
)
```

## Arguments

- fit:

  A fitted unconditional model from
  [`fit_mixture()`](https://pdvalencia.github.io/mixtureEM/reference/fit_mixture.md)
  (or
  [`fit_rmlca()`](https://pdvalencia.github.io/mixtureEM/reference/fit_rmlca.md),
  [`fit_lcga()`](https://pdvalencia.github.io/mixtureEM/reference/fit_lcga.md),
  [`fit_gmm()`](https://pdvalencia.github.io/mixtureEM/reference/fit_gmm.md)).

- outcome:

  The distal outcome: a numeric vector (continuous) or a
  factor/character/integer vector (categorical). Must have one value per
  case of the data the model was fit to.

- covariates:

  Optional covariates that adjust the outcome.

- outcome_type:

  One of `"auto"` (default; inferred from `outcome`), `"continuous"`, or
  `"categorical"`.

- slopes:

  When `covariates` are supplied, whether their effect is `"pooled"`
  (default; one slope shared across classes), `"class_specific"` (every
  covariate gets its own slope per class), or a character vector of
  covariate names (or a one-sided formula naming them, e.g.
  `~ loc1 + loc2`) giving a slope per class to just those covariates
  while the rest stay pooled. The last form – letting the class moderate
  some covariates while adjusting for others – is continuous-outcome
  only. For a continuous outcome, `"class_specific"` is the same model
  as naming every covariate, and is fitted as one: its standard errors
  are the sandwich clustered on the case, and
  [`outcome_contrasts()`](https://pdvalencia.github.io/mixtureEM/reference/outcome_contrasts.md)
  can compare its slopes.

- predictors:

  Optional covariates that predict class membership, in any form
  [`add_covariates()`](https://pdvalencia.github.io/mixtureEM/reference/add_covariates.md)
  accepts. When supplied, the class-membership regression and the
  outcome are estimated together in one step-3 model: each case's
  contribution is \\\sum_k P(k \mid z) f(y \mid k, x)\\ times the
  classification-error term, so the class probabilities the outcome is
  related to are the covariate-specific ones rather than the overall
  class sizes (Vermunt, 2010). `outcome` may then name several distal
  outcomes (a data frame, or a formula naming several columns of
  `data`); each is specified as it would be on its own and all share
  `covariates`. Estimated separately, the outcomes and the predictors
  would each answer a different question from different class
  probabilities; estimated jointly, the paths from the predictors to the
  classes, from the classes to each outcome, and from the covariates to
  each outcome are adjusted for one another. Available with
  `correction = "ML"` (the default here) or `"BCH"`, and pooled slopes.
  Under ML the standard errors are those of the joint step-3
  log-likelihood (`se = "hessian"` for the inverse observed information,
  otherwise the sandwich); under BCH the model is fitted to the
  BCH-weighted log-likelihood and its standard errors are always the
  sandwich clustered on the case. Under `se = "corrected"` both add step
  1's uncertainty; see `se`.

- variances:

  For a continuous outcome: `"equal"` (default; one variance, or
  residual variance when there are `covariates`, shared by the classes)
  or `"class_specific"` (one per class). The shared variance is the more
  stable choice when a class is small, because one estimate pools every
  case; one per class is the better description when the classes really
  do differ in spread, and the two can be compared by their BIC on the
  same fit.

- correction:

  Bias correction for the third step: `"auto"` (default) picks `"BCH"`
  for continuous outcomes (Bakk & Vermunt, 2016) and `"ML"` for
  categorical outcomes; or set `"BCH"`, `"ML"`, `"none"` directly.
  Three-step only; an error with `steps = 2`.

- steps:

  `3` (default) for the bias-adjusted three-step, or `2` for the
  two-step estimator of Bakk and Kuha (2018): `fit`'s measurement model
  and class sizes are held fixed and the outcome model is estimated by
  maximising the full likelihood, every case's class probabilities
  recomputed under the joint model at each iteration. No classification
  step, no correction. See `n_steps` in
  [`fit_mixture()`](https://pdvalencia.github.io/mixtureEM/reference/fit_mixture.md);
  for a distal outcome the two-step's standard errors do not yet carry
  the step-1 uncertainty, and the printed output says so.

- se:

  Standard-error estimator passed on to the third step: `"corrected"`
  (default), `"robust"`, or `"hessian"`. It governs the covariate part
  of the third step. Under `correction = "BCH"` every outcome, and a BCH
  model with `predictors`, reports a sandwich clustered on the case
  whatever this is set to, with `"corrected"` adding step 1's
  uncertainty: the expanded data set carries one weighted record per
  class per case, so a case-clustered sandwich is the only estimator
  that prices the information the correction gives up. With
  `predictors`, under either correction, `"corrected"` adds the step-1
  term of
  [covariate_se](https://pdvalencia.github.io/mixtureEM/reference/covariate_se.md)
  over the whole joint model, so the outcome blocks carry step 1's
  uncertainty too; under BCH it enters through the BCH weights, which
  are rebuilt from the moved posteriors, and is added to the
  case-clustered sandwich. Under BCH a distal outcome fitted without
  `predictors` takes the same term, as the one-block case of the same
  model. The term is first-order: with poorly separated classes (entropy
  R-squared below about .70) and a few hundred cases, the corrected BCH
  standard errors tend to run wide, by up to about 20% in simulation,
  where the sandwich alone runs narrow by about as much; the excess
  shrinks as the sample grows. Under ML a distal outcome fitted without
  `predictors` treats step 1's estimates as known, so with poorly
  separated classes its intervals are somewhat too narrow.

- assignment:

  How step 1's posteriors are turned into the assigned-class variable
  whose classification error the correction inverts. `"proportional"`
  (default) gives every case a weight in every class equal to its
  posterior probability; `"modal"` assigns each case to its most likely
  class outright. The default follows Bakk, Tekle and Vermunt (2013),
  who compared the two rules across 54 simulation conditions and found
  proportional at least as accurate everywhere and clearly better when
  the classes are poorly separated. Use `"modal"` when reproducing an
  analysis whose classes were assigned that way.

- max_iter:

  Maximum iterations for the step-3 estimation.

- data:

  Optional data frame to take the variables from, in which case
  `outcome` may be a one-sided formula naming one column (`~ bmi`), and
  `covariates` a one-sided formula (interactions included, e.g.
  `~ age * sex`) or a vector of column names.

- missing_predictors:

  What to do with a case missing one of `covariates` or `predictors`.
  `"listwise"` (the default) leaves it out of the step-3 model, with a
  message saying how many cases that was; `"mean"` keeps it and
  completes the missing value under the class-invariant Gaussian
  marginal of the covariates, which for a single covariate is its mean.
  A missing *outcome* is never deleted for this reason: that case simply
  contributes no outcome information. See "Missing covariates" below.

- ...:

  Currently unused.

## Value

A `mixture_model` with the distal-outcome model attached. Use
[`summary()`](https://rdrr.io/r/base/summary.html) for class-specific
means or probabilities and their tests, and
[`outcome_contrasts()`](https://pdvalencia.github.io/mixtureEM/reference/outcome_contrasts.md)
for which classes differ from which, rather than whether any of them do.

## Missing covariates

Under `"listwise"` a case missing a covariate of the outcome, or a class
predictor, is still part of step 1 and of the classification-error table
the correction inverts, and is left out of step 3 only; its class
probabilities stay those of the measurement model. For a covariate of
the outcome, deletion leaves the regression unbiased when whether the
covariate is missing depends on the covariate itself or on the class,
and completing it with the mean does not; when missingness is unrelated
to anything the two agree (Little, 1992). For class predictors see
"Missing predictors" in
[`add_covariates()`](https://pdvalencia.github.io/mixtureEM/reference/add_covariates.md).

## References

Vermunt, J. K. (2010). Latent class modeling with covariates: Two
improved three-step approaches. *Political Analysis*, *18*(4), 450–469.
[doi:10.1093/pan/mpq025](https://doi.org/10.1093/pan/mpq025)

Bakk, Z., & Vermunt, J. K. (2016). Robustness of stepwise latent class
modeling with continuous distal outcomes. *Structural Equation
Modeling*, *23*(1), 20–31.
[doi:10.1080/10705511.2014.955104](https://doi.org/10.1080/10705511.2014.955104)

Bolck, A., Croon, M., & Hagenaars, J. (2004). Estimating latent
structure models with categorical variables: One-step versus three-step
estimators. *Political Analysis*, *12*(1), 3–27.
[doi:10.1093/pan/mph001](https://doi.org/10.1093/pan/mph001)

Bakk, Z., Tekle, F. B., & Vermunt, J. K. (2013). Estimating the
association between latent class membership and external variables using
bias-adjusted three-step approaches. *Sociological Methodology*,
*43*(1), 272–311.
[doi:10.1177/0081175012470644](https://doi.org/10.1177/0081175012470644)

Little, R. J. A. (1992). Regression with missing X's: A review. *Journal
of the American Statistical Association*, *87*(420), 1227–1237.
[doi:10.1080/01621459.1992.10476282](https://doi.org/10.1080/01621459.1992.10476282)

## See also

[`outcome_contrasts()`](https://pdvalencia.github.io/mixtureEM/reference/outcome_contrasts.md)
for class-vs-class differences on the outcome;
[`add_covariates()`](https://pdvalencia.github.io/mixtureEM/reference/add_covariates.md)
for predictors of class membership;
[`fit_mixture()`](https://pdvalencia.github.io/mixtureEM/reference/fit_mixture.md)
to fit the unconditional model.

## Examples

``` r
set.seed(1)
items <- matrix(rbinom(300, 1, 0.5), nrow = 100)
bmi   <- rnorm(100, mean = 25)
fit   <- fit_mixture(items, n_classes = 2, measurement = "binary")
fit_out <- add_outcome(fit, bmi)
#> Outcome treated as continuous (set `outcome_type` to override).
#> Using 'BCH' bias correction (set `correction` to override).
summary(fit_out)
#> =========================================================
#>              STRUCTURAL MODEL SUMMARY                    
#> =========================================================
#>   Best solution  : found by 20 of 20 starts
#> 
#> CONTINUOUS DISTAL OUTCOME (MEANS)
#> ---------------------------------------------------------
#> 
#> Omnibus test (class differences): Wald chi^2(1) = 0.04, p   0.851
#> 
#>                  Mean       [95% CI]        SE
#>   Class 1       25.063  [24.601, 25.526]     0.236
#>   Class 2       24.967  [24.328, 25.606]     0.326
#> 
#> Pairwise class differences:
#>                     Difference       [95% CI]        P-Value
#>   Class 2 vs 1        -0.096  [-1.105,  0.912]     0.851
#> =========================================================

# The same outcome named in a formula against its data frame
df <- data.frame(bmi = bmi)
fit_out2 <- add_outcome(fit, ~ bmi, data = df)
#> Outcome treated as continuous (set `outcome_type` to override).
#> Using 'BCH' bias correction (set `correction` to override).

if (FALSE) { # \dontrun{
# Class moderates level-of-care while age and gender are only adjusted for:
# a mix of "class_specific" and "pooled" in one model.
fit_mod <- add_outcome(fit, cannabis_days,
                       covariates = data.frame(loc1, loc2, loc3, age, gender),
                       slopes = c("loc1", "loc2", "loc3"))
} # }
```
