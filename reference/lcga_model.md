# Constructor for latent class growth models

Sets up the emission state for a latent class growth model, in which
each class follows its own fixed-effect polynomial trajectory over the
occasions and there is no within-class random effect.

## Usage

``` r
lcga_model(n_components, design, family = "binomial", ...)
```

## Arguments

- n_components:

  Integer. The number of latent classes to estimate.

- design:

  Numeric matrix with one row per occasion and one column per growth
  coefficient, as built by the polynomial in the time scores.

- family:

  Character. `"binomial"`, `"gaussian"` or `"poisson"`.

- ...:

  Additional arguments, ignored.

## Value

A list object of class `c("lcga", "emission")`.
