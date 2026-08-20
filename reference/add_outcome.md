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
  correction = c("auto", "BCH", "ML", "none"),
  se = c("corrected", "robust", "hessian"),
  assignment = c("proportional", "modal"),
  max_iter = 1000,
  data = NULL,
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
  only.

- correction:

  Bias correction for the third step: `"auto"` (default) picks `"BCH"`
  for continuous outcomes (Bakk & Vermunt, 2016) and `"ML"` for
  categorical outcomes; or set `"BCH"`, `"ML"`, `"none"` directly.

- se:

  Standard-error estimator passed on to the third step: `"corrected"`
  (default), `"robust"`, or `"hessian"`. It governs the covariate part
  of the third step. A continuous distal outcome under
  `correction = "BCH"` always reports a sandwich clustered on the case,
  whatever this is set to: the expanded data set carries one weighted
  record per class per case, so a case-clustered sandwich is the only
  estimator that prices the information the correction gives up.

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
  `covariates` a one-sided formula or a vector of column names.

- ...:

  Currently unused.

## Value

A `mixture_model` with the distal-outcome model attached. Use
[`summary()`](https://rdrr.io/r/base/summary.html) for class-specific
means or probabilities and their tests, and
[`outcome_contrasts()`](https://pdvalencia.github.io/mixtureEM/reference/outcome_contrasts.md)
for which classes differ from which, rather than whether any of them do.

## References

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
#> 
#> CONTINUOUS DISTAL OUTCOME (MEANS)
#> ---------------------------------------------------------
#> 
#> Omnibus test (class differences): Wald chi^2(1) = 0.05, p   0.832
#> 
#>                  Mean       [95% CI]        SE
#>   Class 1       25.063  [24.627, 25.500]     0.223
#>   Class 2       24.967  [24.422, 25.512]     0.278
#> 
#> Pairwise class differences:
#>                     Difference       [95% CI]        P-Value
#>   Class 2 vs 1        -0.096  [-0.987,  0.794]     0.832
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
