# Measurement Parameter Estimates with Standard Errors

Returns the measurement parameters of a
[`fit_mixture()`](https://pdvalencia.github.io/mixtureEM/reference/fit_mixture.md)
fit – the item probabilities, the class means and variances, the count
rates – and the class-size logits, each with a standard error, a z
statistic and a two-sided p-value. It is the standard-error companion to
[`measurement_summary()`](https://pdvalencia.github.io/mixtureEM/reference/measurement_summary.md),
which prints the same estimates without them.

## Usage

``` r
measurement_estimates(fit, se = c("robust", "hessian", "opg"))
```

## Arguments

- fit:

  A model returned by
  [`fit_mixture()`](https://pdvalencia.github.io/mixtureEM/reference/fit_mixture.md).

- se:

  Which standard-error estimator: `"robust"` (default), `"hessian"` or
  `"opg"`. See Details.

## Value

A data frame with one row per parameter and class: `block` (the
sub-model for a mixed or blocked measurement model, `NA` otherwise),
`parameter` (`"probability"`, `"mean"`, `"variance"`, `"rate"` or
`"class_logit"`), `item`, `category` (polytomous items only), `class`,
`estimate`, `se`, `z`, `p_value` and `boundary`. The covariance matrix
of the reported parameters is attached as the attribute `"vcov"`.

## Details

The standard errors come from the log-likelihood of the fitted model,
differentiated with respect to the parameters on the scale the fit
estimates them on (logits for probabilities, logs for rates and standard
deviations, means as they are) and carried to the reported scale by the
delta method. The three estimators are those of
[`growth_estimates()`](https://pdvalencia.github.io/mixtureEM/reference/growth_estimates.md):

- `"robust"` (the default) is the sandwich \\H^{-1} B H^{-1}\\ (White,
  1982), which stays valid when the within-class model is only
  approximately right – for instance, normal indicators whose real
  within-class distribution is skewed, or variances held equal across
  classes when they are not.

- `"hessian"` is the inverse observed information, the textbook
  maximum-likelihood variance when the model is correct.

- `"opg"` is the inverse outer product of the case scores.

A gap between `"robust"` and `"hessian"` is a sign of misspecification.
Sample `weights` and a `cluster` variable are honoured as in
[`growth_estimates()`](https://pdvalencia.github.io/mixtureEM/reference/growth_estimates.md).

A parameter held equal across classes – a variance under
`variances_equal = TRUE` – is repeated once per class, with the same
standard error, as in
[`measurement_summary()`](https://pdvalencia.github.io/mixtureEM/reference/measurement_summary.md).

The estimates maximise the log-likelihood plus the prior that keeps them
off the edge of the parameter space (`bayes_constants` in
[`fit_mixture()`](https://pdvalencia.github.io/mixtureEM/reference/fit_mixture.md)),
so the information is that sum's: the prior's curvature is added to
\\H\\, and to \\B\\ as well, so that `"robust"` reduces to `"hessian"`
when the model holds. In \\B\\ the prior's curvature is taken on each
parameter's natural scale (logits, log rates, and for a variance the
precision \\1/\sigma^2\\, on which the variance prior's curvature does
not depend on the data). With `bayes_constants = 0` these are the
textbook likelihood estimators. With covariates on class membership
fitted in one step, the class probabilities belong to the covariate
model and no class-size logits are reported (see
[`coef()`](https://rdrr.io/r/stats/coef.html)); under `n_steps = 3` the
standard errors are those of the step-1 model, which is where the
measurement parameters are estimated.

## Probabilities at the boundary

An item probability that has run to 0 or 1 (within 0.001) is flagged in
the `boundary` column, and a message says how many there are. Its
standard error is not interpretable: the likelihood is flat in that
direction, so any number would measure how far the optimiser stopped
from the edge rather than how precisely the data pin the probability
down. A probability whose logit has run past \\\pm 15\\ is therefore
held where it stopped, as a fixed value, and its standard error is `NA`;
the other parameters' standard errors are those of the model with it
fixed.

## References

White, H. (1982). Maximum likelihood estimation of misspecified models.
*Econometrica*, *50*(1), 1-25.
[doi:10.2307/1912526](https://doi.org/10.2307/1912526)

Louis, T. A. (1982). Finding the observed information matrix when using
the EM algorithm. *Journal of the Royal Statistical Society: Series B*,
*44*(2), 226-233.
[doi:10.1111/j.2517-6161.1982.tb01203.x](https://doi.org/10.1111/j.2517-6161.1982.tb01203.x)

## See also

[`measurement_summary()`](https://pdvalencia.github.io/mixtureEM/reference/measurement_summary.md),
[`growth_estimates()`](https://pdvalencia.github.io/mixtureEM/reference/growth_estimates.md)
for growth models, and
[`plot.mixture_model()`](https://pdvalencia.github.io/mixtureEM/reference/plot.mixture_model.md),
whose `error_bars` draw these standard errors on the profile bar chart.

## Examples

``` r
set.seed(1)
X <- rbind(matrix(rnorm(300, 0), 100), matrix(rnorm(300, 2), 100))
fit <- fit_mixture(X, n_classes = 2, measurement = "continuous",
                   n_init = 5)
measurement_estimates(fit)
#>    block   parameter       item category class    estimate         se
#> 1   <NA>        mean     Cont_1       NA     1 0.115005689 0.08965207
#> 2   <NA>        mean     Cont_1       NA     2 2.162196420 0.10559643
#> 3   <NA>        mean     Cont_2       NA     1 0.006602604 0.10315322
#> 4   <NA>        mean     Cont_2       NA     2 2.031920995 0.12217526
#> 5   <NA>        mean     Cont_3       NA     1 0.135073104 0.11388157
#> 6   <NA>        mean     Cont_3       NA     2 1.953767073 0.10576491
#> 7   <NA>    variance     Cont_1       NA     1 0.790528560 0.09592641
#> 8   <NA>    variance     Cont_1       NA     2 0.790528560 0.09592641
#> 9   <NA>    variance     Cont_2       NA     1 1.112948932 0.13183461
#> 10  <NA>    variance     Cont_2       NA     2 1.112948932 0.13183461
#> 11  <NA>    variance     Cont_3       NA     1 1.098648458 0.11896909
#> 12  <NA>    variance     Cont_3       NA     2 1.098648458 0.11896909
#> 13  <NA> class_logit log(p1/p2)       NA     1 0.113577140 0.15331913
#>              z      p_value boundary
#> 1   1.28280014 1.995621e-01    FALSE
#> 2  20.47603757 3.521682e-93    FALSE
#> 3   0.06400773 9.489641e-01    FALSE
#> 4  16.63119908 4.142210e-62    FALSE
#> 5   1.18608394 2.355892e-01    FALSE
#> 6  18.47273433 3.422988e-76    FALSE
#> 7   8.24099013 1.707913e-16    FALSE
#> 8   8.24099013 1.707913e-16    FALSE
#> 9   8.44200907 3.119292e-17    FALSE
#> 10  8.44200907 3.119292e-17    FALSE
#> 11  9.23473860 2.589174e-20    FALSE
#> 12  9.23473860 2.589174e-20    FALSE
#> 13  0.74078909 4.588213e-01    FALSE
```
