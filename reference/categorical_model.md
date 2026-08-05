# Constructor for Categorical models

Sets up the initial state and class structure for categorical emission
models (like Bernoulli or Multinoulli) before the EM algorithm runs.

## Usage

``` r
categorical_model(n_components, type = "bernoulli", max_val = NULL, ...)
```

## Arguments

- n_components:

  Integer. The number of latent classes/components to estimate.

- type:

  Character. The specific distribution type, usually "bernoulli".

- max_val:

  Integer or NULL. The maximum category value (used for multinoulli).

- ...:

  Additional arguments passed to the method.

## Value

A list object of class `c(type, "emission")` containing the model state.
