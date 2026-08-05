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

## Examples

``` r
set.seed(1)
X <- matrix(rbinom(500, 1, 0.5), nrow = 100)
fit <- fit_mixture(X, n_components = 2, measurement = "binary")
#> Note: `X`, `Y`, `n_components`, and `structural` are the legacy interface. The current arguments are `indicators`, `n_classes`, `predictors`, and `outcome` / `outcome_covariates`.
classification_diagnostics(fit)
#> =========================================================
#>           AVERAGE POSTERIOR PROBABILITIES (AvePP)        
#> =========================================================
#> Rows: Modal Assignment | Columns: Mean Probability
#> 
#>                  Prob C 1 Prob C 2
#> Assigned Class 1    0.822    0.178
#> Assigned Class 2    0.281    0.719
#> =========================================================
#> 
#> =========================================================
#>                CLASSIFICATION TABLE                      
#> =========================================================
#> Rows: model-expected membership | Columns: modal assignment
#> 
#>         Modal 1 Modal 2    Total
#> Class 1 46.0528 12.3858  58.4385
#> Class 2  9.9472 31.6142  41.5615
#> Total   56.0000 44.0000 100.0000
#> 
#> Classification error: 0.2233 (22.33% of 100 cases)
#> =========================================================
```
