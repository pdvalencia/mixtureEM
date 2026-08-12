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

The diagonal of the AvePP matrix reads as a percentage of the cases
assigned to that class: a value of 0.91 "suggests that 91% of subjects
in the assigned class fit that category, while 9% of the subjects in
that class do not accurately fit that category" (Fanti & Henrich, 2010,
as reported by Lee et al., 2023, p. 653).

This function does not compute entropy, which is the other number
solutions are judged by. That is `fit$metrics$entropy`, printed as
`Rel. Entropy` and carried as the `Entropy` column of
[`compare_mixtures()`](https://pdvalencia.github.io/mixtureEM/reference/compare_mixtures.md)
and
[`compare_longitudinal()`](https://pdvalencia.github.io/mixtureEM/reference/compare_longitudinal.md),
where it is documented.

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

Lee, T. K., Wickrama, K. A. S., & O'Neal, C. W. (2023). An introduction
to growth mixture models (GMM). In *International Encyclopedia of
Education* (4th ed., Vol. 14, pp. 646-655). Elsevier.
[doi:10.1016/B978-0-12-818630-5.10076-4](https://doi.org/10.1016/B978-0-12-818630-5.10076-4)

Nagin, D. S. (2005). *Group-Based Modeling of Development*. Harvard
University Press.

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
#> Assigned Class 1    0.823    0.177
#> Assigned Class 2    0.283    0.717
#> =========================================================
#> 
#> =========================================================
#>                CLASSIFICATION TABLE                      
#> =========================================================
#> Rows: model-expected membership | Columns: modal assignment
#> 
#>         Modal 1 Modal 2   Total
#> Class 1  46.094  12.433  58.527
#> Class 2   9.906  31.567  41.473
#> Total    56.000  44.000 100.000
#> 
#> Classification error: 0.2234 (22.34% of 100 cases)
#> =========================================================
```
