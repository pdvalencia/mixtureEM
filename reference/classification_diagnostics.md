# Print Classification Diagnostics

Prints the two tables that describe how cleanly the model assigns cases.

The **Average Posterior Probability (AvePP)** matrix has one row per set
of observations modally assigned to a class and one column per class,
holding the mean posterior probability. High values on the diagonal, low
values off it, indicate well-separated classes.

The **classification table** cross-classifies the probabilistic
memberships against the modal assignment, and yields the classification
error: the proportion of cases the modal rule is expected to place in
the wrong class. See
[`classification_table()`](https://pdvalencia.github.io/mixtureEM/reference/classification_table.md)
for the details, and
[`absolute_fit()`](https://pdvalencia.github.io/mixtureEM/reference/absolute_fit.md)
and
[`bivariate_residuals()`](https://pdvalencia.github.io/mixtureEM/reference/bivariate_residuals.md)
for the fit of the model itself rather than the quality of its
assignments.

Both tables use the case weights when the model was fitted with any.

## Usage

``` r
classification_diagnostics(object, ...)

# Default S3 method
classification_diagnostics(object, ...)
```

## Arguments

- object:

  A fitted `mixture_model` object returned by
  [`fit_mixture`](https://pdvalencia.github.io/mixtureEM/reference/fit_mixture.md).

- ...:

  Passed to methods.

## Value

Invisibly, a list with `ave_pp` (the K x K matrix), `table` (the
classification table) and `error` (the classification error). All are
also printed to the console.

## References

Celeux, G., & Soromenho, G. (1996). An entropy criterion for assessing
the number of clusters in a mixture model. *Journal of Classification*,
*13*(2), 195-212.
[doi:10.1007/BF01246098](https://doi.org/10.1007/BF01246098)

Nagin, D. S. (2005). *Group-Based Modeling of Development*. Harvard
University Press.

classification_diagnostics(fit)

## Examples

``` r
set.seed(1)
X <- matrix(rbinom(500, 1, 0.5), nrow = 100)
fit <- fit_mixture(X, n_components = 2, measurement = "binary")
#> Note: `X`, `Y`, `n_components`, and `structural` are the legacy interface. The current arguments are `indicators`, `n_classes`, `predictors`, and `outcome` / `outcome_covariates`.
```
