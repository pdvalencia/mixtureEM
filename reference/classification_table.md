# Classification Table and Classification Error

Cross-classifies the probabilistic class memberships against the modal
assignment, which is what quantifies the cost of treating a fitted class
as though it were an observed group. Entry \\(k, m)\\ is
\\\sum\_{i:\\\mathrm{modal}(i) = m} w_i P(k \mid y_i)\\, so the rows sum
to the model-expected class sizes and the columns to the modal counts.
The two sets of totals disagree by exactly the amount modal assignment
distorts the class proportions.

The classification error is \\1 - \mathrm{trace}/N\\: the proportion of
cases the modal rule is expected to place in the wrong class. It is the
quantity that motivates the bias-adjusted 3-step estimators, since it is
the error those corrections exist to undo — see
[`fit_mixture()`](https://pdvalencia.github.io/mixtureEM/reference/fit_mixture.md)'s
`n_steps` and `correction` arguments.

## Usage

``` r
classification_table(object)
```

## Arguments

- object:

  A model fitted by
  [`fit_mixture()`](https://pdvalencia.github.io/mixtureEM/reference/fit_mixture.md),
  or by
  [`fit_lta()`](https://pdvalencia.github.io/mixtureEM/reference/fit_lta.md),
  in which case one table per occasion is returned.

## Value

An object of class `classification_table`: the K x K matrix, with the
classification `error`, `n` and expected/modal class sizes attached as
attributes. For a latent transition model, a list of such tables.

## See also

[`classification_diagnostics()`](https://pdvalencia.github.io/mixtureEM/reference/classification_diagnostics.md),
which prints this alongside the average posterior probabilities.

## Examples

``` r
set.seed(1)
X <- matrix(rbinom(600, 1, 0.5), ncol = 6)
fit <- fit_mixture(X, n_components = 2, measurement = "binary")
#> Note: `X`, `Y`, `n_components`, and `structural` are the legacy interface. The current arguments are `indicators`, `n_classes`, `predictors`, and `outcome` / `outcome_covariates`.
classification_table(fit)
#> =========================================================
#>                CLASSIFICATION TABLE                      
#> =========================================================
#> Rows: model-expected membership | Columns: modal assignment
#> 
#>         Modal 1 Modal 2    Total
#> Class 1 57.8988  7.5058  65.4046
#> Class 2 11.1012 23.4942  34.5954
#> Total   69.0000 31.0000 100.0000
#> 
#> Classification error: 0.1861 (18.61% of 100 cases)
#> =========================================================
```
