# Print Measurement Model Parameters

Prints a formatted table of the fitted measurement model parameters:
item-response probabilities for categorical models, or means for
Gaussian models. Results are broken down by latent class. Handles both
flat and nested (mixed) measurement models.

For a growth model —
[`fit_gmm`](https://pdvalencia.github.io/mixtureEM/reference/fit_gmm.md)
or
[`fit_lcga`](https://pdvalencia.github.io/mixtureEM/reference/fit_lcga.md)
— the measurement parameters are the growth-factor means, their
variances and covariances, the residual variances and the fitted
trajectory, and those are what the table holds. A parameter held equal
across classes is repeated once per class rather than reported once, so
the table can be joined to anything else indexed by class; the
constraint is stated in the printed heading.

## Usage

``` r
# S3 method for class 'gmm'
measurement_summary(object, ...)

# S3 method for class 'lcga'
measurement_summary(object, ...)

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
`parameter` (`"probability"`, `"mean"`, or `"rate"`; for a growth model
`"growth_mean"`, `"growth_variance"`, `"growth_covariance"`,
`"growth_regression"`, `"residual_variance"` or `"fitted"`), `item`,
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
