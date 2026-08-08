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

## References

Oberski, D. L., van Kollenburg, G. H., & Vermunt, J. K. (2013). A Monte
Carlo evaluation of three methods to detect local dependence in binary
data latent class models. *Advances in Data Analysis and
Classification*, *7*(3), 267-279.
[doi:10.1007/s11634-013-0146-2](https://doi.org/10.1007/s11634-013-0146-2)

bivariate_residuals(fit)

## See also

[`absolute_fit()`](https://pdvalencia.github.io/mixtureEM/reference/absolute_fit.md),
[`classification_table()`](https://pdvalencia.github.io/mixtureEM/reference/classification_table.md).

## Examples

``` r
set.seed(1)
X <- matrix(rbinom(600, 1, 0.5), ncol = 6)
fit <- fit_mixture(X, n_components = 2, measurement = "binary")
#> Note: `X`, `Y`, `n_components`, and `structural` are the legacy interface. The current arguments are `indicators`, `n_classes`, `predictors`, and `outcome` / `outcome_covariates`.
```
