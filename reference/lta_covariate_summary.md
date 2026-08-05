# Covariate Effects in a Latent Transition Model

Prints the multinomial-logit coefficients for the covariates predicting
the initial latent status and, where fitted, the transitions between
statuses. Coefficients are contrasts against the last latent status,
which is the reference category; exponentiating gives an odds ratio.

For transitions fitted with `transition_effects = "common"` the design
contains the origin-status dummies (labelled `from:k`) that supply the
row-specific intercepts, followed by the covariate slopes, which are
shared across origin statuses. With `"by_origin"` a separate table is
printed per origin status.

## Usage

``` r
lta_covariate_summary(object, digits = 3)
```

## Arguments

- object:

  A model fitted by
  [`fit_lta()`](https://pdvalencia.github.io/mixtureEM/reference/fit_lta.md)
  with covariates or a group.

- digits:

  Number of digits to print.

## Value

`object`, invisibly.
