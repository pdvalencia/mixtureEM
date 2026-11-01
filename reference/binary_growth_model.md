# Constructor for binary growth mixture models

Sets up the emission state for a growth mixture model of a binary
outcome: each class follows its own polynomial trajectory on the logit
scale, and cases vary about it through normally distributed random
growth factors that are integrated out by Gauss-Hermite quadrature.

## Usage

``` r
binary_growth_model(
  n_components,
  design,
  random_effects = "intercept_slope",
  psi = "equal",
  n_quadrature = 30,
  adaptive = FALSE,
  principal = TRUE,
  start_pool = NULL,
  ...
)
```

## Arguments

- n_components:

  Integer. The number of latent classes to estimate.

- design:

  Numeric matrix with one row per occasion and one column per growth
  factor, as built by the polynomial in the time scores.

- random_effects:

  Which growth factors vary within a class: `"intercept_slope"`,
  `"intercept"` or `"all"`.

- psi:

  `"equal"`, `"variances"` or `"free"`, as in
  [`fit_gmm()`](https://pdvalencia.github.io/mixtureEM/reference/fit_gmm.md).

- n_quadrature:

  The density of the quadrature rule: the number of nodes per random
  effect on the standard rule, and with two or more random effects on
  the principal-axis rule (`principal`), two thirds of this many nodes
  per unit of log-odds standard deviation along each axis, at least a
  tenth of it, and at most `n_quadrature^q` nodes in all.

- adaptive:

  Logical. Keep the standard normal's nodes for every response pattern
  (`FALSE`, the default) or centre them on each pattern's posterior.

- principal:

  Logical. With two or more random effects, lay the rule along each
  class's principal axes of the random effects (`TRUE`, the default,
  when `adaptive` is `FALSE`) or keep `n_quadrature` nodes on each
  standardised random effect.

- start_pool:

  Optional list of `n_components` x `ncol(design)` trajectory matrices
  that half the random starts are drawn around; `NULL` draws every start
  from the pooled trajectory and the observed response patterns.

- ...:

  Additional arguments, ignored.

## Value

A list object of class `c("binary_growth", "emission")`.
