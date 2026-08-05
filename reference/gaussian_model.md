# Constructor for Gaussian models

Sets up the initial state and class structure for continuous emission
models (like Gaussian with diagonal or unit variance).

## Usage

``` r
gaussian_model(n_components, type = "gaussian_unit")
```

## Arguments

- n_components:

  Integer. The number of latent classes/components to estimate.

- type:

  Character. The specific variance structure, e.g., "gaussian_diag" or
  "gaussian_unit".

## Value

A list object containing the model state.
