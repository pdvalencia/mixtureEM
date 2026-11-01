# Vuong-Lo-Mendell-Rubin Test of K Against K + 1 Classes

Compares a mixture with one class fewer (`restricted`) against one with
an extra class (`full`), fitted to the same data, by the
likelihood-ratio statistic \\2(\ell_1 - \ell_0)\\ referred to the
weighted sum of chi-squares of Vuong (1989), as Lo, Mendell and Rubin
(2001) proposed for mixtures, and by their finite-sample adjusted
statistic \\LR / (1 + 1 / (d \ln n))\\, with \\d\\ the difference in
free parameters and \\n\\ the number of cases, referred to the same
distribution. It takes any pair of fits
[`compare_mixtures()`](https://pdvalencia.github.io/mixtureEM/reference/compare_mixtures.md)
tests, growth mixtures from
[`fit_gmm()`](https://pdvalencia.github.io/mixtureEM/reference/fit_gmm.md)
included, so a growth model's number of classes can be tested on fits
already in hand. It also takes two latent transition models from
[`fit_lta()`](https://pdvalencia.github.io/mixtureEM/reference/fit_lta.md)
that differ by one status, one status against two included, with the
same items, occasions and options otherwise (`smoothing = 0` and
`bayes_constants` all zero for a published comparison).

Read
[`compare_mixtures()`](https://pdvalencia.github.io/mixtureEM/reference/compare_mixtures.md)'s
section on the test before relying on it: the theorem behind the
reference distribution does not hold at the point where the larger
mixture reduces to the smaller one, and neither version's p-values are
uniform under the null. BIC and
[`blrt()`](https://pdvalencia.github.io/mixtureEM/reference/blrt.md) are
the better guides. Fit both models without priors (`bayes_constants` all
zero) when the p-value is to be compared with a published one; the test
is stated for the plain likelihood.

## Usage

``` r
vlmr_test(restricted, full, version = c("robust", "standard"))
```

## Arguments

- restricted:

  The model with fewer classes (statuses, for an LTA).

- full:

  The model with one class (status) more, fitted to the same data.

- version:

  `"robust"` (the default) builds the reference distribution from the
  sandwich covariance; `"standard"` from the ordinary one, Vuong's own
  form (see
  [`compare_mixtures()`](https://pdvalencia.github.io/mixtureEM/reference/compare_mixtures.md)).

## Value

A list of class `"vlmr_test"`: `statistic` (\\2 \Delta\ell\\), `df`,
`mean` and `sd` of the reference distribution, `p_value`, `adjusted`
(the adjusted statistic) and `p_adjusted`, plus `version`,
`ll_restricted`, `ll_full`, `n` and `classes` (the two numbers of
statuses, for an LTA).

## References

Lo, Y., Mendell, N. R., & Rubin, D. B. (2001). Testing the number of
components in a normal mixture. *Biometrika*, *88*(3), 767-778.

Vuong, Q. H. (1989). Likelihood ratio tests for model selection and
non-nested hypotheses. *Econometrica*, *57*(2), 307-333.

## See also

[`compare_mixtures()`](https://pdvalencia.github.io/mixtureEM/reference/compare_mixtures.md),
[`blrt()`](https://pdvalencia.github.io/mixtureEM/reference/blrt.md).

## Examples

``` r
# \donttest{
set.seed(1)
n <- 400
cls <- rbinom(n, 1, 0.5)
y <- outer(ifelse(cls == 1, 3, 1), rep(1, 4)) +
  outer(ifelse(cls == 1, 1, 0.2), 0:3) + matrix(rnorm(n * 4, 0, 0.7), n, 4)
off <- list(categorical = 0, latent = 0, variances = 0)
f1 <- fit_gmm(y, n_classes = 1, times = 4, random_effects = "none",
              n_init = 5, random_state = 1, bayes_constants = off)
f2 <- fit_gmm(y, n_classes = 2, times = 4, random_effects = "none",
              n_init = 5, random_state = 1, bayes_constants = off)
vlmr_test(f1, f2)
#> 
#> Vuong-Lo-Mendell-Rubin test, 1 vs 2 classes (robust)
#> ---------------------------------------------------------
#>   LL: -3120.806 vs -2046.537   2 x diff = 2148.538   df = 3
#>   Reference distribution: mean -3.304, sd 15.688
#>   p = < 1e-04
#>   Adjusted (Lo, Mendell & Rubin): 2035.304, p = < 1e-04
#> 
# }
```
