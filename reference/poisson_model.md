# Constructor for Poisson (count) models

Sets up the initial state and class structure for count emission models,
in which each item is Poisson-distributed within class with its own
rate.

## Usage

``` r
poisson_model(n_components, type = "poisson", ...)
```

## Arguments

- n_components:

  Integer. The number of latent classes/components to estimate.

- type:

  Character. The specific distribution type, "poisson" or the
  missing-data variant "poisson_nan".

- ...:

  Additional arguments passed to the method.

## Value

A list object of class `c(type, "emission")` containing the model state.
