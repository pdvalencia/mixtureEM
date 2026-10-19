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

- predictors:

  Covariates that predict class membership: a data frame, matrix, or
  single vector/factor. Factors are dummy-coded with the first level as
  reference. Must have one row per case of the data the model was fit
  to.

- correction:

  Bias correction for the third step: `"ML"` (default; Vermunt, 2010),
  `"BCH"`, or `"none"`. Three-step only; an error with `steps = 2`.

- steps:

  `3` (default) for the bias-adjusted three-step, or `2` for the
  two-step estimator of Bakk and Kuha (2018): `fit`'s measurement model
  is held fixed and the class-membership regression is estimated by
  maximising the full likelihood, with every case's class probabilities
  recomputed under the joint model at each iteration. No classification
  step, no correction. `fit` is exactly the step-1 estimate the two-step
  starts from. See `n_steps` in
  [`fit_mixture()`](https://pdvalencia.github.io/mixtureEM/reference/fit_mixture.md)
  for what the two-step is and is not. Its standard errors carry the
  step-1 uncertainty (Bakk and Kuha, 2018, eq. 5); see
  [covariate_se](https://pdvalencia.github.io/mixtureEM/reference/covariate_se.md).

- se:

  Standard-error estimator for the structural coefficients:
  `"corrected"` (default), `"robust"`, or `"hessian"`. What each is
  under the three-step and under the two-step is in
  [covariate_se](https://pdvalencia.github.io/mixtureEM/reference/covariate_se.md).

- assignment:

  How step 1's posteriors are turned into the assigned-class variable
  whose classification error the correction inverts. `"proportional"`
  (default) gives every case a weight in every class equal to its
  posterior probability; `"modal"` assigns each case to its most likely
  class outright. The default follows Bakk, Tekle and Vermunt (2013),
  who compared the two rules across 54 simulation conditions and found
  proportional at least as accurate everywhere and clearly better when
  the classes are poorly separated. Use `"modal"` when reproducing an
  analysis whose classes were assigned that way. Three-step only.

- max_iter:

  Maximum iterations for the structural estimation.

- data:

  Optional data frame to take the covariates from, in which case
  `predictors` may be a one-sided formula (`~ age + sex`, or
  `~ age * sex` for an interaction) or a vector of column names instead
  of the columns themselves. A formula's terms – a factor's dummies, an
  interaction's several columns – are recognised as one term by the
  omnibus Wald test in
  [`analytical_wald_test()`](https://pdvalencia.github.io/mixtureEM/reference/analytical_wald_test.md).

- missing_predictors:

  What to do with a case missing a predictor. `"listwise"` (the default)
  leaves it out of the class-membership regression, with a message
  saying how many cases that was; `"mean"` keeps it and completes the
  missing value under the class-invariant Gaussian marginal of the
  predictors (Sterba, 2014), which for a single predictor is its mean.
  See "Missing predictors" below.

- ...:

  Currently unused.

## Value

A `mixture_model` with the class-membership regression attached. Use
[`summary()`](https://rdrr.io/r/base/summary.html) for odds ratios and
omnibus tests; [`coef()`](https://rdrr.io/r/stats/coef.html),
[`confint()`](https://rdrr.io/r/stats/confint.html), and
[`wald_omnibus_test()`](https://pdvalencia.github.io/mixtureEM/reference/wald_omnibus_test.md)
also apply.

## Missing predictors

Under `"listwise"` a case missing a predictor is still part of step 1
and of the classification-error table the correction inverts, because
its indicators are observed and tell us as much about that table as
anyone else's do; it is only left out of the regression of class on the
predictors, the one step that needs its predictor values. Its class
probabilities stay those of the measurement model.

Listwise deletion is the default because it does not bias the regression
when whether a predictor is missing depends on the predictor itself, a
common pattern (income left blank more often by high earners), whereas
completing the value with the mean does: it shrinks the slopes toward
zero, and the shrinkage does not go away as the sample grows. When
missingness is unrelated to anything (missing completely at random) both
rules are unbiased and the mean fill, which keeps every case, can be
somewhat more precise. When missingness depends on the class itself
neither rule is unbiased. Little (1992) sets out why deleting incomplete
cases is unbiased in a regression whose missingness depends only on the
predictors.

`se = "corrected"` (the default) is the Bakk, Oberski and Vermunt (2014)
estimator, which propagates the uncertainty in the step-1 estimates as
well as the step-3 sampling variability. `se = "robust"` reports only
the latter; use it when the assigned classes are to be treated as given.

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

**Why the classes are counted before any covariate is added.** Arch and
Nylund-Gibson (2026) went looking for conditions under which counting
classes without covariates goes wrong when items behave differently
across covariate values (differential item functioning) that the model
does not include, and largely did not find them. It was enumeration
*with* the covariates in the model, but without their direct effects on
the items, that over-extracted classes, most of all when those effects
were strong and the sample was large. That is the order this function
assumes: the number of classes is settled on the unconditional fit you
pass in. One caveat travels with the result: even the study's
"realistic" separation condition had an entropy of about .84, so how
well the advice holds for poorly separated classes is untested there.

## References

Arch, D. A. N., & Nylund-Gibson, K. (2026). Unmodeled DIF in latent
class analysis: Consequences for class enumeration. *Methods in
Psychology*, *14*, 100243.
[doi:10.1016/j.metip.2026.100243](https://doi.org/10.1016/j.metip.2026.100243)

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

Bakk, Z., & Kuha, J. (2018). Two-step estimation of models between
latent classes and external variables. *Psychometrika*, *83*(4),
871–892.
[doi:10.1007/s11336-017-9592-7](https://doi.org/10.1007/s11336-017-9592-7)

Jiang, Y., Elliott, M. R., Sammel, M. D., & Wang, N. (2016). Joint
modeling of cross-sectional health outcomes and longitudinal predictors
via mixtures of latent classes. *Statistics and Its Interface*, *9*,
183–201.

Little, R. J. A. (1992). Regression with missing X's: A review. *Journal
of the American Statistical Association*, *87*(420), 1227–1237.
[doi:10.1080/01621459.1992.10476282](https://doi.org/10.1080/01621459.1992.10476282)

Sterba, S. K. (2014). Handling missing covariates in conditional mixture
models under missing at random assumptions. *Multivariate Behavioral
Research*, *49*(6), 614–632.
[doi:10.1080/00273171.2014.950719](https://doi.org/10.1080/00273171.2014.950719)

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
#>   Best solution  : found by 20 of 20 starts
#> 
#> CATEGORICAL LATENT VARIABLE REGRESSION (CLASS PREDICTORS)
#> Reference Class: 1
#> Standard errors: Bakk-Oberski-Vermunt corrected (robust step 3, hessian step 1)
#> ---------------------------------------------------------
#>                               OR         [95% CI]         P-Value
#> 
#> Class 2 ON
#>   Intercept                0.761  [    0.041,    13.984]     0.854
#>   age                      0.929  [    0.387,     2.229]     0.869
#> =========================================================

# The same covariate named in a formula against its data frame
df <- data.frame(age = age, sex = rbinom(100, 1, 0.5))
fit_cov2 <- add_covariates(fit, ~ age + sex, data = df)
#> Using 'ML' bias correction (set `correction` to override).
```
