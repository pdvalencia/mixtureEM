# Print Measurement Model Parameters

Prints a formatted table of the fitted measurement model parameters:
item-response probabilities for categorical models, or means for
Gaussian models. Results are broken down by latent class. Handles both
flat and nested (mixed) measurement models.

## Usage

``` r
measurement_summary(object, ...)

# Default S3 method
measurement_summary(object, ...)
```

## Arguments

- object:

  A fitted `mixture_model` object returned by
  [`fit_mixture`](https://pdvalencia.github.io/mixtureEM/reference/fit_mixture.md).

- ...:

  Passed to methods.

## Value

Invisibly, a data frame in long format with one row per item, response
category (polytomous items only, `NA` otherwise), and class: columns
`block` (sub-model name for mixed measurement models, `NA` otherwise),
`parameter` (`"probability"`, `"mean"`, or `"rate"`), `item`,
`category`, `class`, and `estimate`. The same numbers are printed as
formatted tables.

## Examples

``` r
set.seed(1)
X <- matrix(rbinom(500, 1, 0.5), nrow = 100)
fit <- fit_mixture(X, n_components = 2, measurement = "binary")
#> Note: `X`, `Y`, `n_components`, and `structural` are the legacy interface. The current arguments are `indicators`, `n_classes`, `predictors`, and `outcome` / `outcome_covariates`.
measurement_summary(fit)
#> =========================================================
#>              MEASUREMENT MODEL PARAMETERS                
#> =========================================================
#> 
#> CATEGORICAL PROBABILITIES
#> Indicator            | Class 1 | Class 2
#> ---------------------------------------- 
#> Item_1               |   0.253 |   0.800
#> Item_2               |   0.574 |   0.492
#> Item_3               |   0.270 |   0.534
#> Item_4               |   0.495 |   0.339
#> Item_5               |   0.516 |   0.405
#> =========================================================
params <- measurement_summary(fit)   # reuse the table programmatically
#> =========================================================
#>              MEASUREMENT MODEL PARAMETERS                
#> =========================================================
#> 
#> CATEGORICAL PROBABILITIES
#> Indicator            | Class 1 | Class 2
#> ---------------------------------------- 
#> Item_1               |   0.253 |   0.800
#> Item_2               |   0.574 |   0.492
#> Item_3               |   0.270 |   0.534
#> Item_4               |   0.495 |   0.339
#> Item_5               |   0.516 |   0.405
#> =========================================================
```
