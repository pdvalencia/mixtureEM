# Bivariate Residuals

A local measure of fit: for each pair of categorical indicators, the
Pearson chi-square of the observed two-way table against the table the
model implies, divided by its degrees of freedom \\(R_a - 1)(R_b - 1)\\.
If the model were true, a bivariate residual should not be substantially
larger than 1.

Where
[`absolute_fit()`](https://pdvalencia.github.io/mixtureEM/reference/absolute_fit.md)
says whether the model fits, this says *where* it fails: a large value
flags the specific pair of items whose association the latent classes do
not reproduce, which is the conditional-independence assumption showing
its seams. This is the classic local-dependence diagnostic (Oberski et
al., 2013), and what the usual referee question about local dependence
asks for.

Unlike the absolute-fit statistics, bivariate residuals do not require
the full response-pattern table and so remain usable with many
indicators: the model-implied two-way margin follows in closed form from
conditional independence, \\P(y_a = r, y_b = s) = \sum_k \gamma_k\\
p_a(r\|k)\\ p_b(s\|k)\\.

With missing data each pair is computed on the cases observing both
items, and the expected counts are scaled to that pair's total. This is
a pairwise-complete statistic rather than a full-information one, so
read it as descriptive when missingness is heavy.

## Usage

``` r
bivariate_residuals(object)
```

## Arguments

- object:

  A model fitted by
  [`fit_mixture()`](https://pdvalencia.github.io/mixtureEM/reference/fit_mixture.md)
  or
  [`fit_rmlca()`](https://pdvalencia.github.io/mixtureEM/reference/fit_rmlca.md)
  with categorical indicators.

## Value

An object of class `bivariate_residuals`: a lower-triangular
indicator-by-indicator matrix, `NA` on and above the diagonal; or `NULL`
(with a message) when the statistic does not apply.

## See also

[`absolute_fit()`](https://pdvalencia.github.io/mixtureEM/reference/absolute_fit.md),
[`classification_table()`](https://pdvalencia.github.io/mixtureEM/reference/classification_table.md).

## Examples

``` r
set.seed(1)
X <- matrix(rbinom(600, 1, 0.5), ncol = 6)
fit <- fit_mixture(X, n_components = 2, measurement = "binary")
#> Note: `X`, `Y`, `n_components`, and `structural` are the legacy interface. The current arguments are `indicators`, `n_classes`, `predictors`, and `outcome` / `outcome_covariates`.
bivariate_residuals(fit)
#> =========================================================
#>                BIVARIATE RESIDUALS                       
#> =========================================================
#> Pearson chi-square per item pair, divided by its df.
#> Values well above 1 flag a pair whose association the
#> classes do not reproduce (local dependence).
#> 
#>          Item1    Item2    Item3    Item4    Item5
#> Item2   0.0419
#> Item3   0.0120   0.0002
#> Item4   0.6327   1.1869   0.4264
#> Item5   0.0697   0.3653   0.0726   0.0279
#> Item6   0.0149   0.2691   0.0671   0.5234   0.1122
#> 
#> Largest: Item4 x Item2 = 1.1869
#> =========================================================
```
