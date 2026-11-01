# Constructor for latent class growth models

Sets up the emission state for a latent class growth model, in which
each class follows its own fixed-effect polynomial trajectory over the
occasions and there is no within-class random effect.

## Usage

``` r
lcga_model(
  n_components,
  design,
  family = "binomial",
  dispersion_structure = "occasion",
  inflation_design = NULL,
  ...
)
```

## Arguments

- n_components:

  Integer. The number of latent classes to estimate.

- design:

  Numeric matrix with one row per occasion and one column per growth
  coefficient, as built by the polynomial in the time scores.

- family:

  Character. `"binomial"`, `"gaussian"`, `"poisson"`, `"zip"` or
  `"negbin"`.

- dispersion_structure:

  Character. For `family = "negbin"`, whether the overdispersion is one
  value per occasion shared by the classes (`"occasion"`) or one per
  class constant over occasions (`"class"`). Ignored by the other
  families.

- inflation_design:

  For `family = "zip"`, the design of the zero-inflation logit, one row
  per occasion. `NULL`, the default, is a single constant column.
  Ignored by the other families.

- ...:

  Additional arguments, ignored.

## Value

A list object of class `c("lcga", "emission")`.
