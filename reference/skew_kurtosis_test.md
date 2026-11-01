# Skewness and Kurtosis Tests of Model Fit

Checks whether a fitted mixture of normal distributions reproduces the
skewness and kurtosis of the data. A mixture of normals can be skewed
and heavy-tailed overall, through the spread of its classes; what it
cannot produce is a shape its classes do not add up to. A data set whose
shape the fitted model does not reproduce is evidence that the number of
classes, or the normality of the outcomes within each class, is wrong
(Muthén, 2003).

Each statistic is computed on the data and on `nsim` data sets simulated
from the fitted model; the simulated values give the statistic's
distribution under the model, and the p-value is two-sided: twice the
smaller of the shares of simulated values at or above, and at or below,
the sample value.

The fitted parameters are estimates, not known values, and that
uncertainty belongs in the reference distribution. Each simulated data
set is therefore generated from its own parameter draw, taken from the
normal approximation to the estimates' sampling distribution (centred on
the estimates, with the inverse of the observed information as its
covariance, on a scale where variances are log standard deviations and
class sizes log-ratios), with classes drawn from the drawn class sizes
and the estimation sample size. Holding the parameters at the estimates
instead understates the spread of the skewness statistics markedly,
because a mixture's skewness moves with the separation and sizes of its
classes, which are the least precisely estimated parts of the model; the
kurtosis statistics are dominated by sampling noise and change little.

## Usage

``` r
skew_kurtosis_test(fit, nsim = 100, seed = NULL)
```

## Arguments

- fit:

  A model with continuous outcomes from
  [`fit_gmm()`](https://pdvalencia.github.io/mixtureEM/reference/fit_gmm.md),
  [`fit_lcga()`](https://pdvalencia.github.io/mixtureEM/reference/fit_lcga.md)
  or
  [`fit_mixture()`](https://pdvalencia.github.io/mixtureEM/reference/fit_mixture.md),
  fitted to complete data.

- nsim:

  Number of simulated data sets.

- seed:

  `NULL`, or an integer for the simulation; a supplied seed leaves the
  caller's random number state as it was.

## Value

A list of class `"skew_kurtosis_test"` with data frames `multivariate`,
`univariate` and `bivariate`, each with columns `statistic`
(`"skewness"` or `"kurtosis"`), the variable(s), `sample`, `mean`, `sd`
and `p_value`.

## Details

The statistics, for complete data:

- univariate: the skewness \\m_3 / s^3\\ and the excess kurtosis \\m_4 /
  s^4 - 3\\ of each outcome, \\m_r\\ its \\r\\-th central sample moment
  and \\s^2\\ its sample variance on \\n - 1\\ degrees of freedom;

- bivariate and multivariate: Mardia's (1970) skewness \\b\_{1,p} =
  n^{-2} \sum\_{i,j} (d_i^\top S^{-1} d_j)^3\\ and kurtosis \\b\_{2,p} =
  n^{-1} \sum_i (d_i^\top S^{-1} d_i)^2\\ of each pair of outcomes and
  of all of them, \\d_i\\ the centred case and \\S\\ the covariance
  matrix on \\n\\ degrees of freedom. Under multivariate normality
  \\b\_{2,p}\\ is about \\p(p + 2)\\.

The p-values carry the simulation's own noise; with the default 100 data
sets they are good to about two decimals.

## References

Mardia, K. V. (1970). Measures of multivariate skewness and kurtosis
with applications. *Biometrika*, *57*(3), 519-530.

Muthén, B. (2003). Statistical and substantive checking in growth
mixture modeling: Comment on Bauer and Curran (2003). *Psychological
Methods*, *8*(3), 369-377.

## Examples

``` r
# \donttest{
set.seed(1)
n <- 300
cls <- rbinom(n, 1, 0.5)
y <- outer(ifelse(cls == 1, 3, 1), rep(1, 4)) +
  outer(ifelse(cls == 1, 1, 0.2), 0:3) + matrix(rnorm(n * 4, 0, 0.7), n, 4)
fit <- fit_gmm(y, n_classes = 2, times = 4, random_effects = "none",
               n_init = 3, random_state = 1)
skew_kurtosis_test(fit, nsim = 20, seed = 1)
#> 
#> Skewness and kurtosis tests of model fit (20 simulated data sets)
#> 
#> Multivariate
#>  statistic variables sample   mean    sd p_value
#>   skewness       all  0.478  0.349 0.090     0.1
#>   kurtosis       all 22.467 22.048 0.715     0.8
#> 
#> Univariate
#>  statistic variables sample   mean    sd p_value
#>   skewness  Item1@T1 -0.018  0.062 0.111     0.6
#>   skewness  Item1@T2  0.217  0.096 0.150     0.5
#>   skewness  Item1@T3  0.164  0.124 0.131     0.7
#>   skewness  Item1@T4  0.110  0.138 0.153     0.8
#>   kurtosis  Item1@T1 -0.946 -0.833 0.151     0.4
#>   kurtosis  Item1@T2 -1.099 -1.128 0.073     0.6
#>   kurtosis  Item1@T3 -1.492 -1.437 0.067     0.4
#>   kurtosis  Item1@T4 -1.596 -1.600 0.050     0.9
#> 
#> Bivariate
#>  statistic             variables sample  mean    sd p_value
#>   skewness Item1@T1 and Item1@T2  0.126 0.079 0.040     0.1
#>   skewness Item1@T1 and Item1@T3  0.057 0.085 0.059     0.9
#>   skewness Item1@T1 and Item1@T4  0.055 0.098 0.059     0.5
#>   skewness Item1@T2 and Item1@T3  0.214 0.093 0.059     0.1
#>   skewness Item1@T2 and Item1@T4  0.071 0.088 0.071     0.9
#>   skewness Item1@T3 and Item1@T4  0.037 0.092 0.047     0.3
#>   kurtosis Item1@T1 and Item1@T2  6.909 6.541 0.351     0.3
#>   kurtosis Item1@T1 and Item1@T3  6.382 6.529 0.326     0.6
#>   kurtosis Item1@T1 and Item1@T4  6.379 6.414 0.258     0.8
#>   kurtosis Item1@T2 and Item1@T3  6.232 6.374 0.261     0.5
#>   kurtosis Item1@T2 and Item1@T4  6.050 6.371 0.320     0.2
#>   kurtosis Item1@T3 and Item1@T4  6.646 6.228 0.325     0.3
#> 
# }
```
