# Compare Mixture Models Across a Range of Class Numbers

Fits a sequence of measurement-only mixture models, one for each value
of `k` in `k_range`, and returns a table of fit indices to guide class
enumeration. The best model according to BIC is identified
automatically.

## Usage

``` r
compare_mixtures(
  X,
  k_range = 1:5,
  measurement = "binary",
  n_init = 10,
  n_steps = 1,
  ...
)
```

## Arguments

- X:

  A numeric matrix or data frame of indicator variables.

- k_range:

  Integer vector of class numbers to fit. All values must be \>= 1.
  Default is `1:5`.

- measurement:

  Character string specifying the measurement model type. See
  [`fit_mixture`](https://pdvalencia.github.io/mixtureEM/reference/fit_mixture.md)
  for accepted values. Default is `"binary"`.

- n_init:

  Positive integer. Number of random restarts per model. Default is
  `10`.

- n_steps:

  Integer. Estimation method: `1`, `2`, or `3`. Default is `1`.

- ...:

  Additional arguments passed to
  [`fit_mixture`](https://pdvalencia.github.io/mixtureEM/reference/fit_mixture.md).

## Value

A named list with three elements:

- `fit_table` Data frame with one row per K and columns `Classes`, `LL`,
  `Params`, `AIC`, `BIC`, `SABIC`, `Entropy`.

- `models` Named list of fitted `mixture_model` objects, one per K
  (names are `"K1"`, `"K2"`, etc.).

- `best_k` Integer. The value of K with the lowest BIC.

## Examples

``` r
set.seed(1)
X <- matrix(rbinom(500, 1, 0.5), nrow = 100)
result <- compare_mixtures(X, k_range = 1:4, measurement = "binary",
                           n_init = 5)
#> Running Model Selection across K = 1 to 4...
#> 
#> Fitting 1-class model...
#> Fitting 2-class model...
#> Fitting 3-class model...
#> Fitting 4-class model...
#> 
#> === Model Selection Summary ===
#>   Classes       LL Params     AIC     BIC   SABIC Entropy
#> 1       1 -342.102      5 694.203 707.229 691.438   1.000
#> 2       2 -340.085     11 702.169 730.826 696.085   0.303
#> 3       3 -337.018     17 708.037 752.325 698.634   0.458
#> 4       4 -334.542     23 715.084 775.003 702.364   0.590
#> 
#> -> Best model according to BIC: 1 classes
result$fit_table
#>   Classes        LL Params      AIC      BIC    SABIC   Entropy
#> 1       1 -342.1016      5 694.2032 707.2290 691.4378 1.0000000
#> 2       2 -340.0846     11 702.1693 730.8262 696.0854 0.3033087
#> 3       3 -337.0184     17 708.0368 752.3247 698.6345 0.4580078
#> 4       4 -334.5422     23 715.0844 775.0033 702.3636 0.5899555
result$best_k
#> [1] 1
```
