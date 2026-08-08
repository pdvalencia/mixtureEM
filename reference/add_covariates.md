# Examine Predictors of Class Membership on a Fitted Model

Takes the latent class model you have already chosen and relates
covariates to class membership with the bias-adjusted three-step
approach (Vermunt, 2010). The measurement model is reused exactly as
fitted — no re-estimation, and no risk of landing on a different
solution — so this is both faster and conceptually cleaner than
re-specifying the model with `predictors`.

## Usage

``` r
add_covariates(
  fit,
  predictors,
  correction = c("ML", "BCH", "none"),
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

- predictors:

  Covariates that predict class membership: a data frame, matrix, or
  single vector/factor. Factors are dummy-coded with the first level as
  reference. Must have one row per case of the data the model was fit
  to.

- correction:

  Bias correction for the third step: `"ML"` (default; Vermunt, 2010),
  `"BCH"`, or `"none"`.

- se:

  Standard-error estimator passed on to the third step: `"corrected"`
  (default), `"robust"`, or `"hessian"`.

- max_iter:

  Maximum iterations for the step-3 estimation.

- ...:

  Currently unused.

## Value

A `mixture_model` with the class-membership regression attached. Use
[`summary()`](https://rdrr.io/r/base/summary.html) for odds ratios and
omnibus tests; [`coef()`](https://rdrr.io/r/stats/coef.html),
[`confint()`](https://rdrr.io/r/stats/confint.html), and
[`wald_omnibus_test()`](https://pdvalencia.github.io/mixtureEM/reference/wald_omnibus_test.md)
also apply.

## References

Vermunt, J. K. (2010). Latent class modeling with covariates: Two
improved three-step approaches. *Political Analysis*, *18*(4), 450–469.
[doi:10.1093/pan/mpq025](https://doi.org/10.1093/pan/mpq025)

Bakk, Z., Oberski, D. L., & Vermunt, J. K. (2014). Relating latent class
assignments to external variables: Standard errors for correct
inference. *Political Analysis*, *22*(4), 520–540.
[doi:10.1093/pan/mpu003](https://doi.org/10.1093/pan/mpu003)

## See also

[`add_outcome()`](https://pdvalencia.github.io/mixtureEM/reference/add_outcome.md)
for distal outcomes;
[`fit_mixture()`](https://pdvalencia.github.io/mixtureEM/reference/fit_mixture.md)
to fit the unconditional model.

## Examples

``` r
set.seed(1)
items <- matrix(rbinom(300, 1, 0.5), nrow = 100)
age   <- rnorm(100)
fit   <- fit_mixture(items, n_classes = 2)
fit_cov <- add_covariates(fit, age)
#> Using 'ML' bias correction (set `correction` to override).
summary(fit_cov)
#> =========================================================
#>              STRUCTURAL MODEL SUMMARY                    
#> =========================================================
#> 
#> CATEGORICAL LATENT VARIABLE REGRESSION (CLASS PREDICTORS)
#> Reference Class: 1
#> Standard errors: Bakk-Oberski-Vermunt corrected (robust step 3, hessian step 1)
#> ---------------------------------------------------------
#>                               OR         [95% CI]         P-Value
#> 
#> Class 2 ON
#>   Intercept                0.730  [    0.000, 18650.316]     0.952
#>   age                      0.919  [    0.271,     3.123]     0.893
#> =========================================================
```
