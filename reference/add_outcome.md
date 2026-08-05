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
  slopes = c("pooled", "class_specific"),
  correction = c("auto", "BCH", "ML", "none"),
  se = c("corrected", "robust", "hessian"),
  max_iter = 1000,
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
  (one slope shared across classes) or `"class_specific"`.

- correction:

  Bias correction for the third step: `"auto"` (default) picks `"BCH"`
  for continuous outcomes (Bakk & Vermunt, 2016) and `"ML"` for
  categorical outcomes; or set `"BCH"`, `"ML"`, `"none"` directly.

- se:

  Standard-error estimator passed on to the third step: `"corrected"`
  (default), `"robust"`, or `"hessian"`.

- max_iter:

  Maximum iterations for the step-3 estimation.

- ...:

  Currently unused.

## Value

A `mixture_model` with the distal-outcome model attached. Use
[`summary()`](https://rdrr.io/r/base/summary.html) for class-specific
means or probabilities and their tests.

## References

Bakk, Z., & Vermunt, J. K. (2016). Robustness of stepwise latent class
modeling with continuous distal outcomes. *Structural Equation
Modeling*, *23*(1), 20–31.
[doi:10.1080/10705511.2014.955104](https://doi.org/10.1080/10705511.2014.955104)

## See also

[`add_covariates()`](https://pdvalencia.github.io/mixtureEM/reference/add_covariates.md)
for predictors of class membership;
[`fit_mixture()`](https://pdvalencia.github.io/mixtureEM/reference/fit_mixture.md)
to fit the unconditional model.

## Examples

``` r
set.seed(1)
items <- matrix(rbinom(300, 1, 0.5), nrow = 100)
bmi   <- rnorm(100, mean = 25)
fit   <- fit_mixture(items, n_classes = 2)
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
#> Omnibus test (class differences): Wald chi^2(1) = 0.04, p   0.832
#> 
#>                  Mean       [95% CI]        SE
#>   Class 1       25.063  [24.627, 25.500]     0.223
#>   Class 2       24.967  [24.423, 25.511]     0.278
#> =========================================================
```
