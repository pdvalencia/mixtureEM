# Constructor for Gaussian models

Sets up the initial state and class structure for continuous emission
models (like Gaussian with diagonal or unit variance).

## Usage

``` r
gaussian_model(n_components, type = "gaussian_unit", variances_equal = FALSE)
```

## Arguments

- n_components:

  Integer. The number of latent classes/components to estimate.

- type:

  Character. The specific variance structure, e.g., "gaussian_diag" or
  "gaussian_unit".

- variances_equal:

  Logical. Hold each item's variance equal across the classes, so the
  classes differ in location only (the homoscedastic latent profile
  model, and the default parameterisation of several commercial
  programs). Ignored by the unit-variance types, which have no variances
  to estimate. The estimated variance is still stored once per class — a
  `K x J` matrix with identical rows — so every reader of the parameters
  (likelihood, plotting, alignment, boundary checks) is unchanged; the
  constraint lives in the M-step and in `n_parameters()`.

## Value

A list object containing the model state.
