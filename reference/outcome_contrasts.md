# Class-vs-Class Contrasts on a Distal Outcome

Tests which latent classes differ on a distal outcome attached by
[`add_outcome()`](https://pdvalencia.github.io/mixtureEM/reference/add_outcome.md),
rather than whether any of them do. The omnibus Wald test printed by
[`summary()`](https://rdrr.io/r/base/summary.html) is a single joint
statistic; this is the set of contrasts behind it, each with a standard
error, a confidence interval and a p-value.

Each contrast is reported as `class` minus `reference`, so a positive
estimate means the class scores above the class it is compared with.
With `ref = NULL` every pair of classes is reported once; naming a class
in `ref` holds it fixed as the comparison.

The standard errors account for the covariance between the class
parameters, which is not optional here. Two class means from one fit are
estimated from the same posteriors, and under the BCH correction from
the same inverted classification-error matrix, so they are correlated;
the standard error of their difference is not the root of the sum of
their squared standard errors, and the naive version is conservative
exactly when the correlation is positive. The full sandwich covariance
the fit already carries is used wherever it is available.

## Usage

``` r
outcome_contrasts(
  fit,
  ref = NULL,
  adjust = c("none", "holm", "bonferroni"),
  level = 0.95,
  outcome = NULL,
  term = NULL,
  ...
)
```

## Arguments

- fit:

  A fitted model with a distal outcome attached, as returned by
  [`add_outcome()`](https://pdvalencia.github.io/mixtureEM/reference/add_outcome.md),
  or a three-step
  [`fit_lta()`](https://pdvalencia.github.io/mixtureEM/reference/fit_lta.md)
  fit with `distal`, where the classes are the statuses at the last
  occasion.

- ref:

  Optional reference class, as an integer between 1 and the number of
  classes. `NULL` (the default) reports all pairs.

- adjust:

  Multiplicity adjustment across the reported contrasts, passed to
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html): `"none"`
  (the default), `"holm"` or `"bonferroni"`. With `"none"` the p-values
  are the per-contrast ones, which is what to report when the contrasts
  were named in advance; an adjustment belongs on the all-pairs table
  read as a family.

- level:

  Confidence level for the intervals. Default `0.95`.

- outcome:

  For a fit holding several distal outcomes (see `predictors` in
  [`add_outcome()`](https://pdvalencia.github.io/mixtureEM/reference/add_outcome.md)),
  the outcome to contrast, by name or position. `NULL` (the default) is
  the first; the `outcome` attribute of the result names the one used.

- term:

  For an outcome regressed on covariates with class-specific slopes
  (`slopes` in
  [`add_outcome()`](https://pdvalencia.github.io/mixtureEM/reference/add_outcome.md)
  naming the covariate), the covariate whose slopes to contrast, by its
  column name. `NULL` (the default) contrasts the class intercepts
  instead. See "Slope contrasts" below.

- ...:

  Currently unused.

## Value

A data frame of class `outcome_contrasts`, one row per contrast, with
columns `category` (the outcome category the contrast is on, `NA` for a
continuous outcome), `class`, `reference`, `estimate`, `se`, `lower`,
`upper`, `z` and `p`, plus `p_adj` when `adjust` is not `"none"` and
`OR`, `OR_lower`, `OR_upper` for a categorical outcome, whose estimates
are log odds. The `method` attribute names the covariance the standard
errors came from.

## Slope contrasts

When the class moderates a covariate's effect on the outcome – each
class has its own regression line – the question after "do the slopes
differ at all" (the Wald test printed by
[`summary()`](https://rdrr.io/r/base/summary.html)) is "which classes
differ in slope". A slope contrast answers it: the estimate is how much
more the outcome rises per unit of the covariate in `class` than in
`reference`. This is the second stage of the test sequence for latent
class moderation in Arch, Nylund-Gibson and Ing (2026); the first is the
omnibus test and the third is each class's slope against zero, both
already in [`summary()`](https://rdrr.io/r/base/summary.html).

Fit the outcome with `slopes` naming the covariate, e.g.
`add_outcome(fit, y, covariates = df["x"], slopes = "x")`, and pass
`term = "x"`. That fit estimates every class's line in one regression,
so its covariance includes the covariance between the slopes of
different classes, which the contrast standard error needs. On a
continuous outcome `slopes = "class_specific"` is the same fit, with
every covariate moderated, and can be contrasted the same way.

With `term` left `NULL` on such a fit the contrasts are between the
class intercepts, i.e. the class means at covariate values of zero;
centring the covariate first moves that comparison to its mean.

## References

Arch, D. A. N., Nylund-Gibson, K., & Ing, M. (2026). Moderation with a
latent class variable. *Behavior Research Methods*, 58, 108.

## See also

[`add_outcome()`](https://pdvalencia.github.io/mixtureEM/reference/add_outcome.md)
for fitting the outcome model and
[`summary()`](https://rdrr.io/r/base/summary.html) for the omnibus test;
[`wald_omnibus_test()`](https://pdvalencia.github.io/mixtureEM/reference/wald_omnibus_test.md)
and
[`bootstrap_covariates()`](https://pdvalencia.github.io/mixtureEM/reference/bootstrap_covariates.md)
for the bootstrap route, which is what to use for an outcome whose
parameters are estimated one class at a time.

## Examples

``` r
set.seed(1)
items <- matrix(rbinom(300, 1, 0.5), nrow = 100)
bmi   <- rnorm(100, mean = 25)
fit   <- fit_mixture(items, n_classes = 2, measurement = "binary")
fit_out <- add_outcome(fit, bmi)
#> Outcome treated as continuous (set `outcome_type` to override).
#> Using 'BCH' bias correction (set `correction` to override).
outcome_contrasts(fit_out)
#> Class contrasts on the distal outcome
#> Standard errors: sandwich covariance of the class means
#> 
#>   Contrast       Estimate        SE  [95% CI]            P-Value
#>   2 vs 1           -0.096     0.454  [ -0.987,   0.794]     0.832

# Which classes differ in the slope of the outcome on a covariate
x     <- rnorm(100)
fit_m <- add_outcome(fit, bmi + x, covariates = data.frame(x = x),
                     slopes = "x")
#> Outcome treated as continuous (set `outcome_type` to override).
#> Using 'BCH' bias correction (set `correction` to override).
outcome_contrasts(fit_m, term = "x")
#> Class contrasts on the distal outcome
#> Standard errors: covariance of the class-specific slopes on x
#> 
#>   Contrast       Estimate        SE  [95% CI]            P-Value
#>   2 vs 1           -0.063     0.410  [ -0.867,   0.741]     0.877
```
