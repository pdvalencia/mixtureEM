# Absolute Fit of a Categorical Mixture Model

Compares the observed response-pattern frequencies with those the model
implies, over the contingency table formed by crossing every categorical
indicator. Three members of the power-divergence family are reported:
the likelihood-ratio statistic \\G^2\\ (also written \\L^2\\), the
Pearson \\X^2\\, and the Cressie-Read statistic (\\\lambda = 2/3\\),
each on \\df = W - P - 1\\ degrees of freedom, where \\W\\ is the number
of cells in the table and \\P\\ the number of free parameters.

The statistics require fully categorical indicators. Even then they
should be read with care: the table has \\W\\ cells and is usually
extremely sparse, so the chi-square reference distribution is unreliable
and the value is best used to compare models rather than to test one in
isolation.
[`blrt()`](https://pdvalencia.github.io/mixtureEM/reference/blrt.md)
tests a model against one with fewer classes without relying on that
reference distribution.

## Usage

``` r
absolute_fit(object)
```

## Arguments

- object:

  A model fitted by
  [`fit_mixture()`](https://pdvalencia.github.io/mixtureEM/reference/fit_mixture.md),
  [`fit_rmlca()`](https://pdvalencia.github.io/mixtureEM/reference/fit_rmlca.md)
  or
  [`fit_lta()`](https://pdvalencia.github.io/mixtureEM/reference/fit_lta.md).

## Value

An object of class `absolute_fit` with elements `g2`, `x2`,
`cressie_read`, `df`, the corresponding `p_value`s, `dissimilarity`,
`n_cells` and `n_patterns`; or `NULL` (with a message) when the
statistics do not apply. With missing data, also `g2_mcar`, `x2_mcar`,
`cressie_read_mcar`, `df_mcar` and `p_value_mcar` for the block computed
under MCAR, `ll_sat` for the saturated model's log-likelihood, and
`mar = TRUE`.

## Missing data

With one or more missing values (categorical, plain
[`fit_mixture()`](https://pdvalencia.github.io/mixtureEM/reference/fit_mixture.md)
models only), the statistics are computed under the missing-at-random
(MAR) assumption instead: the model is compared not to the raw response
table, which no longer exists once cases have different items observed,
but to a saturated model fit to the same partition of the data by which
items each case observed. `df` is smaller than in the complete-data case
(\\df = W - 1 - P\\) because the saturated baseline already accounts for
the missingness pattern. A short block giving the model's fit jointly
with the stronger missing-completely-at-random (MCAR) assumption is
printed underneath; use
[`mcar_test()`](https://pdvalencia.github.io/mixtureEM/reference/mcar_test.md)
to test that assumption on its own. This can be slow, or refused
outright, once the number of indicator categories crossed together grows
large – the same \\W\\ that already makes the complete-data table
sparse.

## References

Langeheine, R., Pannekoek, J., & van de Pol, F. (1996). Bootstrapping
goodness-of-fit measures in categorical data analysis. *Sociological
Methods & Research*, *24*(4), 492-516.
[doi:10.1177/0049124196024004004](https://doi.org/10.1177/0049124196024004004)

## See also

[`bivariate_residuals()`](https://pdvalencia.github.io/mixtureEM/reference/bivariate_residuals.md)
for the local counterpart,
[`classification_table()`](https://pdvalencia.github.io/mixtureEM/reference/classification_table.md),
[`blrt()`](https://pdvalencia.github.io/mixtureEM/reference/blrt.md),
[`mcar_test()`](https://pdvalencia.github.io/mixtureEM/reference/mcar_test.md)
for testing the missingness mechanism on its own.

## Examples

``` r
set.seed(1)
X <- matrix(rbinom(600, 1, 0.5), ncol = 6)
fit <- fit_mixture(X, n_components = 2, measurement = "binary")
#> Note: `X`, `Y`, `n_components`, and `structural` are the legacy interface. The current arguments are `indicators`, `n_classes`, `predictors`, and `outcome` / `outcome_covariates`.
absolute_fit(fit)
#> =========================================================
#>                   ABSOLUTE FIT                           
#> =========================================================
#> Table: 64 cells, 48 observed response patterns
#> Free parameters: 13   df: 50
#> 
#> Statistic               Value    p-value
#> ---------------------------------------- 
#> L-squared             62.4171     0.1118
#> X-squared             50.0325     0.4721
#> Cressie-Read          51.3465     0.4207
#> Dissimilarity          0.2819           
#> =========================================================
```
