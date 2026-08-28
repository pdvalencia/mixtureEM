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

  Optional data frame to take the covariates from, in which case
  `predictors` may be a one-sided formula (`~ age + sex`, or
  `~ age * sex` for an interaction) or a vector of column names instead
  of the columns themselves. A formula's terms – a factor's dummies, an
  interaction's several columns – are recognised as one term by the
  omnibus Wald test in
  [`analytical_wald_test()`](https://pdvalencia.github.io/mixtureEM/reference/analytical_wald_test.md).

- ...:

  Currently unused.

## Value

A `mixture_model` with the class-membership regression attached. Use
[`summary()`](https://rdrr.io/r/base/summary.html) for odds ratios and
omnibus tests; [`coef()`](https://rdrr.io/r/stats/coef.html),
[`confint()`](https://rdrr.io/r/stats/confint.html), and
[`wald_omnibus_test()`](https://pdvalencia.github.io/mixtureEM/reference/wald_omnibus_test.md)
also apply.

## Details

A case missing a predictor is retained, not deleted: the missing value
is completed under the class-invariant Gaussian marginal of the
predictors (Sterba, 2014), so the analysis keeps its full N. An analysis
that listwise deletes them is fitted to fewer cases; check the reported
N before comparing coefficients with a published set.

`se = "corrected"` (the default) is the Bakk, Oberski and Vermunt (2014)
estimator, which propagates the uncertainty in the step-1 estimates as
well as the step-3 sampling variability. `se = "robust"` reports only
the latter; use it when reproducing an analysis whose standard errors
were computed that way.

**Why a three-step function at all, when the measurement model could
just be refit with `predictors` supplied to
[`fit_mixture()`](https://pdvalencia.github.io/mixtureEM/reference/fit_mixture.md)
directly.** Jiang, Elliott, Sammel and Wang (2016) name the failure mode
of that one-step alternative: when a covariate participates in forming
the classes, the joint models "tend to have a high chance of
artificially creating spurious mixture components to enhance predictive
accuracy for the sample that is used to derive the model," and the
apparent gain does not survive an independent validation sample. Their
simulation also found the deviance-based criteria (AIC, BIC) more
reliable for choosing the number of classes than predictive ones, which
over-selected components in pursuit of fit. `add_covariates()`'s
two-stage design – fit the measurement model first, decide on classes,
only then look at what predicts them – keeps that choice from being made
by a covariate that happens to correlate with the outcome.

## References

Vermunt, J. K. (2010). Latent class modeling with covariates: Two
improved three-step approaches. *Political Analysis*, *18*(4), 450–469.
[doi:10.1093/pan/mpq025](https://doi.org/10.1093/pan/mpq025)

Bakk, Z., Oberski, D. L., & Vermunt, J. K. (2014). Relating latent class
assignments to external variables: Standard errors for correct
inference. *Political Analysis*, *22*(4), 520–540.
[doi:10.1093/pan/mpu003](https://doi.org/10.1093/pan/mpu003)

Bolck, A., Croon, M., & Hagenaars, J. (2004). Estimating latent
structure models with categorical variables: One-step versus three-step
estimators. *Political Analysis*, *12*(1), 3–27.
[doi:10.1093/pan/mph001](https://doi.org/10.1093/pan/mph001)

Bakk, Z., Tekle, F. B., & Vermunt, J. K. (2013). Estimating the
association between latent class membership and external variables using
bias-adjusted three-step approaches. *Sociological Methodology*,
*43*(1), 272–311.
[doi:10.1177/0081175012470644](https://doi.org/10.1177/0081175012470644)

Jiang, Y., Elliott, M. R., Sammel, M. D., & Wang, N. (2016). Joint
modeling of cross-sectional health outcomes and longitudinal predictors
via mixtures of latent classes. *Statistics and Its Interface*, *9*,
183–201.

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
fit   <- fit_mixture(items, n_classes = 2, measurement = "binary")
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
#>   Intercept                0.760  [    0.000, 18310.293]     0.958
#>   age                      0.929  [    0.280,     3.077]     0.904
#> =========================================================

# The same covariate named in a formula against its data frame
df <- data.frame(age = age, sex = rbinom(100, 1, 0.5))
fit_cov2 <- add_covariates(fit, ~ age + sex, data = df)
#> Using 'ML' bias correction (set `correction` to override).
```
