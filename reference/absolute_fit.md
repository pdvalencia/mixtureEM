# Absolute Fit of a Categorical Mixture Model

Compares the observed response-pattern frequencies with those the model
implies, over the contingency table formed by crossing every categorical
indicator. Three members of the power-divergence family are reported:
the likelihood-ratio statistic \\G^2\\ (also written \\L^2\\), the
Pearson \\X^2\\, and the Cressie-Read statistic (\\\lambda = 2/3\\),
each on \\df = W - P - 1\\ degrees of freedom, where \\W\\ is the number
of cells in the table and \\P\\ the number of free parameters.

The statistics are defined only for fully categorical indicators
observed without missingness. Even then they should be read with care:
the table has \\W\\ cells and is usually extremely sparse, so the
chi-square reference distribution is unreliable and the value is best
used to compare models rather than to test one in isolation.
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
`cressie_read`, `df`, the corresponding `p_value`s, `n_cells` and
`n_patterns`; or `NULL` (with a message) when the statistics do not
apply.

## References

Langeheine, R., Pannekoek, J., & van de Pol, F. (1996). Bootstrapping
goodness-of-fit measures in categorical data analysis. *Sociological
Methods & Research*, *24*(4), 492-516.
[doi:10.1177/0049124196024004004](https://doi.org/10.1177/0049124196024004004)

## See also

[`bivariate_residuals()`](https://pdvalencia.github.io/mixtureEM/reference/bivariate_residuals.md)
for the local counterpart,
[`classification_table()`](https://pdvalencia.github.io/mixtureEM/reference/classification_table.md),
[`blrt()`](https://pdvalencia.github.io/mixtureEM/reference/blrt.md).

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
#> =========================================================
```
