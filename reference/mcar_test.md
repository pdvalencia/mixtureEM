# Test Whether Data Are Missing Completely at Random

Tests the missingness mechanism itself, separately from whether the
fitted model fits.
[`absolute_fit()`](https://pdvalencia.github.io/mixtureEM/reference/absolute_fit.md)
already reports a model comparison under the weaker missing-at-random
(MAR) assumption; `mcar_test(fit)` asks the stronger question a reviewer
sometimes wants answered on its own – whether the pattern of missing
values could plausibly be unrelated to the data, rather than depending
on it.

The test compares the observed response frequencies, partitioned by
which items each case answered, against a saturated model fit to that
same partition. It does not depend on the mixture model fitting well:
the saturated model is as flexible as the data allow, so what remains is
a statement about the missingness pattern, not about the number of
classes.

## Usage

``` r
mcar_test(object)
```

## Arguments

- object:

  A model fitted by
  [`fit_mixture()`](https://pdvalencia.github.io/mixtureEM/reference/fit_mixture.md),
  with categorical indicators and at least one missing value.

## Value

A list with `stat`, `df` and `p_value`; or `NULL` (with a message) when
the data are complete, since there is then nothing to test.

## Reading the result

A small `p_value` says the missingness is **not** missing completely at
random – whether a value is missing depends on the data in some way.
That is a common and often unsurprising finding (people who skip a
question about drug use are not a random subset of respondents), and it
does **not** by itself mean the weaker, more common missing-at-random
assumption fails too; nothing here tests that. It also says nothing
about whether the mixture model itself fits – that question belongs to
[`absolute_fit()`](https://pdvalencia.github.io/mixtureEM/reference/absolute_fit.md),
which enters only through the size of the table this test uses, not
through its own log-likelihood.

## See also

[`absolute_fit()`](https://pdvalencia.github.io/mixtureEM/reference/absolute_fit.md),
whose missing-data branch this function reuses.

## Examples

``` r
set.seed(1)
X <- matrix(rbinom(600, 1, 0.5), ncol = 6)
X[sample(length(X), 30)] <- NA
fit <- fit_mixture(X, n_components = 2, measurement = "binary")
#> Note: `X`, `Y`, `n_components`, and `structural` are the legacy interface. The current arguments are `indicators`, `n_classes`, `predictors`, and `outcome` / `outcome_covariates`.
mcar_test(fit)
#> =========================================================
#>                MISSING COMPLETELY AT RANDOM              
#> =========================================================
#> Chi-square: 109.1950   df: 246   p-value: 1.0000
#> =========================================================
```
