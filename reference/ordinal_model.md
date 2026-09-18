# Constructor for the ordinal emission model

Constructor for the ordinal emission model

## Usage

``` r
ordinal_model(n_components, n_features = NULL, cats, type = "ordinal", ...)
```

## Arguments

- n_components:

  Integer. Number of latent classes/components.

- n_features:

  Integer. Number of items (raw items, not expanded columns).

- cats:

  Integer vector of length `n_features`, each entry `>= 2`: the number
  of ordered categories for that item.

- type:

  Character. `"ordinal"` or `"ordinal_nan"`; both use the same methods
  (see the `_nan` aliasing below).

- ...:

  Unused; present for `build_emission()`'s uniform call shape.

## Value

A list object of class `c(type, "emission")`.
