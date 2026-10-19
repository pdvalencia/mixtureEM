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
  adjust = c("none", "BH", "holm", "bonferroni"),
  level = 0.95,
  outcome = NULL,
  term = NULL,
  effect_size = FALSE,
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
  (the default), `"BH"`, `"holm"` or `"bonferroni"`. With `"none"` the
  p-values are the per-contrast ones, which is what to report when the
  contrasts were named in advance; an adjustment belongs on the
  all-pairs table read as a family. See "Choosing an adjustment" below.

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

- effect_size:

  If `TRUE`, add standardized effect sizes for a continuous outcome:
  Cohen's d for every contrast and the LTB-\\\omega\\ across all
  classes. Default `FALSE`. See "Effect sizes" below.

- ...:

  Currently unused.

## Value

A data frame of class `outcome_contrasts`, one row per contrast, with
columns `category` (the outcome category the contrast is on, `NA` for a
continuous outcome), `class`, `reference`, `estimate`, `se`, `lower`,
`upper`, `z` and `p`, plus `p_adj` when `adjust` is not `"none"` and
`OR`, `OR_lower`, `OR_upper` for a categorical outcome, whose estimates
are log odds. With `effect_size = TRUE`, also `d`, `d_se`, `d_lower` and
`d_upper`, and the `omega` attribute. The `method` attribute names the
covariance the standard errors came from.

## Effect sizes

A contrast's p-value says whether two classes can be told apart, not how
far apart they are; with enough cases a trivial gap is significant, and
a large one can miss significance in a small class. `effect_size = TRUE`
reports the size of the differences on a scale that does not depend on
the outcome's units.

For each pair, Cohen's d is the difference between the two class means
divided by the square root of the average of the two class variances,
\\d\_{jk} = (\mu_j - \mu_k) / \sqrt{(\sigma^2_j + \sigma^2_k)/2}\\
(Carter, 2026); with a variance shared by the classes it is the gap
measured in within-class standard deviations. The point estimate is the
standardized mean difference of Hedges and Olkin (2014).

Across all classes, the LTB-\\\omega\\ of Lanza, Tan and Bray (2013) is
the spread of the class means around their overall mean, each weighted
by its class size: \\\omega = \sqrt{\sum_k \pi_k (\mu_k - \bar\mu)^2}\\,
\\\bar\mu = \sum_k \pi_k \mu_k\\. It is in the outcome's units. The
printed value is standardized by the outcome's total standard deviation,
\\\omega\_{std} = \omega / \sqrt{\sum_k \pi_k \sigma^2_k + \omega^2}\\,
so that \\\omega\_{std}^2\\ is the share of the outcome's variance that
lies between the classes. Both are in the `omega` attribute of the
result (`$standardized` and `$raw`, each with `estimate`, `se`, `lower`
and `upper`); the raw one is the \\\omega\\ of Carter's reporting
workflow. The effect size answers a different question from the omnibus
test in [`summary()`](https://rdrr.io/r/base/summary.html) and does not
replace it.

Standard errors and intervals are by the delta method from the joint
covariance of the class sizes, means and variances, so the uncertainty
in the sizes and the variances is carried along with that of the means.
That covariance is stored by a one-step
[`fit_gmm()`](https://pdvalencia.github.io/mixtureEM/reference/fit_gmm.md)
or
[`fit_lcga()`](https://pdvalencia.github.io/mixtureEM/reference/fit_lcga.md)
outcome (the robust sandwich of
[`growth_estimates()`](https://pdvalencia.github.io/mixtureEM/reference/growth_estimates.md))
and by `add_outcome(correction = "ML")`; for other outcome models,
including the BCH correction, `effect_size = TRUE` is an error. No
thresholds such as "small" or "large" are attached: what size of
difference matters depends on the outcome and the field.

## Choosing an adjustment

The adjustments answer different questions. `"holm"` and `"bonferroni"`
control the family-wise error rate, the chance of even one false
positive among the contrasts; that is the right guarantee for a narrowly
confirmatory comparison, where a single wrong claim is costly. `"BH"`
(Benjamini and Hochberg, 1995) controls the false discovery rate, the
expected share of false positives among the contrasts declared
significant; that is the right guarantee when the aim is to describe the
pattern of differences between the classes, which is the usual aim of a
distal-outcome analysis.

The pairwise contrasts are positively dependent by construction: each
class mean enters several of them (in a three-class model, class 1's
mean is in both the 1-vs-2 and the 1-vs-3 contrast). Family-wise
procedures do not use that dependence and are conservative here, so they
suppress more of a real pattern than the number of contrasts alone
suggests. Benjamini and Yekutieli (2001) prove that the BH procedure
keeps its false discovery rate under positive regression dependence,
which covers positively correlated one-sided tests; the p-values here
are two-sided, for which that proof does not apply as it stands.

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

Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery
rate: A practical and powerful approach to multiple testing. *Journal of
the Royal Statistical Society: Series B*, *57*(1), 289–300.
[doi:10.1111/j.2517-6161.1995.tb02031.x](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x)

Benjamini, Y., & Yekutieli, D. (2001). The control of the false
discovery rate in multiple testing under dependency. *The Annals of
Statistics*, *29*(4), 1165–1188.
[doi:10.1214/aos/1013699998](https://doi.org/10.1214/aos/1013699998)

Carter, D. B. (2026). Distal outcomes in mixture modeling: A guide for
pairwise comparisons, multiplicity control, and effect size reporting.
*Behavior Research Methods*, *58*, 268.
[doi:10.3758/s13428-026-03144-4](https://doi.org/10.3758/s13428-026-03144-4)

Hedges, L. V., & Olkin, I. (2014). *Statistical methods for
meta-analysis*. Academic Press.

Lanza, S. T., Tan, X., & Bray, B. C. (2013). Latent class analysis with
distal outcomes: A flexible model-based approach. *Structural Equation
Modeling*, *20*(1), 1–26.
[doi:10.1080/10705511.2013.742377](https://doi.org/10.1080/10705511.2013.742377)

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
#>   2 vs 1           -0.096     0.515  [ -1.105,   0.912]     0.851

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
#>   2 vs 1           -0.063     0.517  [ -1.077,   0.950]     0.903
```
