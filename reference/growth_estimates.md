# Growth Mixture Parameter Estimates with Standard Errors

Returns every parameter of a
[`fit_gmm()`](https://pdvalencia.github.io/mixtureEM/reference/fit_gmm.md)
fit – the growth-factor means (intercepts, when there are
`growth_predictors`), the growth-factor regressions, the growth-factor
variances and covariances, the residual variances, the fitted trajectory
and the class-size logits – with a standard error, a z statistic and a
two-sided p-value. For a
[`fit_lcga()`](https://pdvalencia.github.io/mixtureEM/reference/fit_lcga.md)
fit it returns the growth coefficients of each class (on the link scale:
logits for `"binomial"`, log rates for `"poisson"`), the residual
variance (`"gaussian"` only), the fitted trajectory on the response
scale and the class-size logits.

## Usage

``` r
growth_estimates(fit, se = c("robust", "hessian", "opg"))
```

## Arguments

- fit:

  A model returned by
  [`fit_gmm()`](https://pdvalencia.github.io/mixtureEM/reference/fit_gmm.md)
  or
  [`fit_lcga()`](https://pdvalencia.github.io/mixtureEM/reference/fit_lcga.md).

- se:

  Which standard-error estimator: `"robust"` (default), `"hessian"` or
  `"opg"`. See Details.

## Value

A data frame with one row per parameter and class: `parameter`
(`"growth_mean"`, `"growth_regression"`, `"growth_variance"`,
`"growth_covariance"`, `"residual_variance"`, `"fitted"`,
`"class_logit"`, or, for a continuous `outcome` fitted in one step,
`"outcome_mean"` and `"outcome_variance"`), `item`, `class`, `estimate`,
`se`, `z` and `p_value`. The covariance matrix of the reported
parameters is attached as the attribute `"vcov"`.

## Details

The standard errors are computed from the log-likelihood of the fitted
model itself, by differentiating it with respect to the parameters the
fit estimated, and are then carried to the reported scale by the delta
method. Three estimators are offered:

- `"robust"` (the default) is the sandwich \\H^{-1} B H^{-1}\\, where
  \\H\\ is the observed information (minus the matrix of second
  derivatives of the log-likelihood) and \\B\\ the sum over cases of the
  outer product of each case's score (White, 1982). It stays a valid
  variance when the within-class normal model is only approximately
  right, which for a growth mixture is the usual situation: the classes
  are partly a device for describing a non-normal distribution of
  trajectories (Bauer & Curran, 2003).

- `"hessian"` is \\H^{-1}\\, the inverse observed information, the
  textbook maximum-likelihood variance when the model is correct.

- `"opg"` is \\B^{-1}\\, the inverse outer product of the case scores,
  which estimates the same quantity as `"hessian"` under a correct model
  from first derivatives alone.

When the model is right the three agree in large samples; a marked gap
between `"robust"` and the other two is a sign of misspecification. With
sample `weights` every case's contribution is weighted, and with a
`cluster` variable the scores are summed within clusters before \\B\\ is
formed.

Constrained parameters are repeated across the classes they are shared
by, with the same standard error, as in
[`measurement_summary()`](https://pdvalencia.github.io/mixtureEM/reference/measurement_summary.md).
The prior that keeps the variances away from zero during estimation is
not part of the variance: these are the standard errors of the
likelihood.

A growth-factor covariance matrix outside the positive semi-definite
cone still gets standard errors (see `psi` in
[`fit_gmm()`](https://pdvalencia.github.io/mixtureEM/reference/fit_gmm.md));
the point is a maximum of the likelihood like any other. Where the
information matrix is singular – a parameter the data do not identify –
the standard errors are `NA` and a warning says so.

The class-size logits are \\\log(\pi_k / \pi_K)\\. When the class sizes
are modelled by `predictors` they are not reported here; see
[`coef()`](https://rdrr.io/r/stats/coef.html). Under `n_steps = 3` the
standard errors here are those of the step-1 model, which is where the
growth parameters are estimated.

## References

White, H. (1982). Maximum likelihood estimation of misspecified models.
*Econometrica*, *50*(1), 1-25.
[doi:10.2307/1912526](https://doi.org/10.2307/1912526)

Louis, T. A. (1982). Finding the observed information matrix when using
the EM algorithm. *Journal of the Royal Statistical Society: Series B*,
*44*(2), 226-233.
[doi:10.1111/j.2517-6161.1982.tb01203.x](https://doi.org/10.1111/j.2517-6161.1982.tb01203.x)

Bauer, D. J., & Curran, P. J. (2003). Distributional assumptions of
growth mixture models: Implications for overextraction of latent
trajectory classes. *Psychological Methods*, *8*(3), 338-363.
[doi:10.1037/1082-989X.8.3.338](https://doi.org/10.1037/1082-989X.8.3.338)

## See also

[`fit_gmm()`](https://pdvalencia.github.io/mixtureEM/reference/fit_gmm.md),
[`fit_lcga()`](https://pdvalencia.github.io/mixtureEM/reference/fit_lcga.md),
[`measurement_summary()`](https://pdvalencia.github.io/mixtureEM/reference/measurement_summary.md).

## Examples

``` r
# \donttest{
set.seed(1)
n <- 400
cls <- rbinom(n, 1, 0.5)
icept <- ifelse(cls == 1, 3, 1) + rnorm(n, 0, 1)
slope <- ifelse(cls == 1, 1.0, 0.2) + rnorm(n, 0, 0.4)
y <- outer(icept, rep(1, 4)) + outer(slope, 0:3) + matrix(rnorm(n * 4, 0, 0.7), n, 4)
fit <- fit_gmm(y, times = 4, n_classes = 2, n_init = 10, random_state = 1)
growth_estimates(fit)
#>            parameter                  item class    estimate         se
#> 1        growth_mean             intercept     1 1.022643213 0.11856188
#> 2        growth_mean             intercept     2 2.942648364 0.14195794
#> 3        growth_mean                linear     1 0.176867998 0.05319109
#> 4        growth_mean                linear     2 1.000575848 0.05366114
#> 5    growth_variance             intercept     1 1.205134599 0.20754515
#> 6    growth_variance             intercept     2 1.205134599 0.20754515
#> 7    growth_variance                linear     1 0.189013287 0.03385954
#> 8    growth_variance                linear     2 0.189013287 0.03385954
#> 9  growth_covariance intercept with linear     1 0.005391572 0.04137730
#> 10 growth_covariance intercept with linear     2 0.005391572 0.04137730
#> 11 residual_variance                    T1     1 0.508533257 0.07885768
#> 12 residual_variance                    T1     2 0.508533257 0.07885768
#> 13 residual_variance                    T2     1 0.536970588 0.05124768
#> 14 residual_variance                    T2     2 0.536970588 0.05124768
#> 15 residual_variance                    T3     1 0.572257386 0.05940420
#> 16 residual_variance                    T3     2 0.572257386 0.05940420
#> 17 residual_variance                    T4     1 0.336185730 0.08540090
#> 18 residual_variance                    T4     2 0.336185730 0.08540090
#> 19            fitted                    T1     1 1.022643213 0.11856188
#> 20            fitted                    T1     2 2.942648364 0.14195794
#> 21            fitted                    T2     1 1.199511210 0.12415150
#> 22            fitted                    T2     2 3.943224212 0.14536513
#> 23            fitted                    T3     1 1.376379208 0.14976267
#> 24            fitted                    T3     2 4.943800059 0.16694018
#> 25            fitted                    T4     1 1.553247206 0.18735715
#> 26            fitted                    T4     2 5.944375907 0.20091307
#> 27       class_logit            log(p1/p2)     1 0.181563706 0.15809949
#>             z       p_value
#> 1   8.6253966  6.387090e-18
#> 2  20.7290157  1.896288e-95
#> 3   3.3251432  8.837313e-04
#> 4  18.6461890  1.356057e-77
#> 5   5.8066141  6.374884e-09
#> 6   5.8066141  6.374884e-09
#> 7   5.5822752  2.373924e-08
#> 8   5.5822752  2.373924e-08
#> 9   0.1303026  8.963270e-01
#> 10  0.1303026  8.963270e-01
#> 11  6.4487470  1.127786e-10
#> 12  6.4487470  1.127786e-10
#> 13 10.4779494  1.090832e-25
#> 14 10.4779494  1.090832e-25
#> 15  9.6332820  5.785158e-22
#> 16  9.6332820  5.785158e-22
#> 17  3.9365594  8.265817e-05
#> 18  3.9365594  8.265817e-05
#> 19  8.6253966  6.387090e-18
#> 20 20.7290157  1.896288e-95
#> 21  9.6616730  4.386412e-22
#> 22 27.1263413 4.815973e-162
#> 23  9.1904026  3.913737e-20
#> 24 29.6142009 9.808932e-193
#> 25  8.2903011  1.129621e-16
#> 26 29.5868049 2.209060e-192
#> 27  1.1484142  2.507976e-01
# }
```
